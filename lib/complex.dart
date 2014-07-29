/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
library complex;

import 'dart:math' as math;
import './src/math.dart' as fastmath;

/**
 * Representation of a Complex number, i.e. a number which has both a
 * real and imaginary part.
 *
 * Implementations of arithmetic operations handle {@code NaN} and
 * infinite values according to the rules for {@link double}, i.e.
 * {@link #equals} is an equivalence relation for all instances that have
 * a {@code NaN} in either real or imaginary part, e.g. the following are
 * considered equal:
 * <ul>
 *  <li>{@code 1 + NaNi}</li>
 *  <li>{@code NaN + i}</li>
 *  <li>{@code NaN + NaNi}</li>
 * </ul>
 * Note that this is in contradiction with the IEEE-754 standard for floating
 * point numbers (according to which the test {@code x == x} must fail if
 * {@code x} is {@code NaN}).
 */
class Complex {

  /** The square root of -1. A number representing "0.0 + 1.0i" */
  static final I = new Complex(0.0, -1.0);
  /** A complex number representing "NaN + NaNi" */
  static final NAN = new Complex(double.NAN, double.NAN);
  /** A complex number representing "+INF + INFi" */
  static final INFINITY = new Complex(double.INFINITY, double.INFINITY);
  /** A complex number representing "1.0 + 0.0i" */
  static final ONE = new Complex(1.0);
  /** A complex number representing "0.0 + 0.0i" */
  static final ZERO = new Complex(0.0);

  double _imaginary, _real;

  /** The imaginary part. */
  double get imaginary => _imaginary;

  /** The real part. */
  double get real => _real;

  /**
   * Create a complex number given the real part and optionally the imaginary part.
   */
  Complex(this._real, [this._imaginary = 0.0]);

  /**
     * Get a hashCode for the complex number.
     * Any {@code Double.NaN} value in real or imaginary part produces
     * the same hash code {@code 7}.
     */
  int get hashCode {
    if (isNaN) {
      return 7;
    }
    return 37 * (17 * _imaginary.hashCode + _real.hashCode);
  }

  /**
   * True if the real and imaginary parts are finite; otherwise, false.
   */
  bool get isFinite => !isNaN && _real.isFinite && _imaginary.isFinite;

  /**
   * True if the real and imaginary parts are positive infinity or negative infinity; otherwise, false.
   */
  bool get isInfinite => !isNaN && (_real.isInfinite || _imaginary.isInfinite);

  /**
   * True if the real and imaginary parts are the double Not-a-Number value; otherwise, false.
   */
  bool get isNaN => _real.isNaN || _imaginary.isNaN;

  /**
   * Return the absolute value of this complex number.
   * Returns {@code NaN} if either real or imaginary part is {@code NaN}
   * and {@code double.INFINITY} if neither part is {@code NaN},
   * but at least one part is infinite.
   *
   * @return the absolute value.
   */
  double abs() {
    if (isNaN) {
      return double.NAN;
    }
    if (isInfinite) {
      return double.INFINITY;
    }
    if (_real.abs() < _imaginary.abs()) {
      if (_imaginary == 0.0) {
        return _real.abs();
      }
      double q = _real / _imaginary;
      return _imaginary.abs() * math.sqrt(1 + q * q);
    } else {
      if (_real == 0.0) {
        return _imaginary.abs();
      }
      double q = _imaginary / _real;
      return _real.abs() * math.sqrt(1 + q * q);
    }
  }

  /**
   * Returns a {@code Complex} whose value is
   * {@code (this + addend)}.
   * Uses the definitional formula
   * <pre>
   *  <code>
   *   (a + bi) + (c + di) = (a+c) + (b+d)i
   *  </code>
   * </pre>
   * <br/>
   * If either {@code this} or {@code addend} has a {@code NaN} value in
   * either part, {@link #NaN} is returned; otherwise {@code Infinite}
   * and {@code NaN} values are returned in the parts of the result
   * according to the rules for {@link double} arithmetic.
   */
  Complex operator +(Object addend) {
    if (addend is Complex) {
      if (isNaN || addend.isNaN) {
        return Complex.NAN;
      }
      return complex(_real + addend.real, _imaginary + addend.imaginary);
    } else if (addend is num) {
      if (isNaN || addend.isNaN) {
        return Complex.NAN;
      }
      return complex(_real + addend, _imaginary);
    } else {
      throw new ArgumentError('factor must be a num or a Complex');
    }
  }


  /**
   * Return the conjugate of this complex number.
   * The conjugate of {@code a + bi} is {@code a - bi}.
   * <br/>
   * {@link #NaN} is returned if either the real or imaginary
   * part of this Complex number equals {@code double.NAN}.
   * <br/>
   * If the imaginary part is infinite, and the real part is not
   * {@code NaN}, the returned value has infinite imaginary part
   * of the opposite sign, e.g. the conjugate of
   * {@code 1 + INFINITY i} is {@code 1 - NEGATIVE_INFINITY i}.
   *
   * @return the conjugate of this Complex object.
   */
  Complex conjugate() {
    if (isNaN) {
      return Complex.NAN;
    }
    return complex(real, -imaginary);
  }

  /**
   * Returns a {@code Complex} whose value is
   * {@code (this / divisor)}.
   * Implements the definitional formula
   * <pre>
   *  <code>
   *    a + bi          ac + bd + (bc - ad)i
   *    ----------- = -------------------------
   *    c + di         c<sup>2</sup> + d<sup>2</sup>
   *  </code>
   * </pre>
   * but uses
   * <a href="http://doi.acm.org/10.1145/1039813.1039814">
   * prescaling of operands</a> to limit the effects of overflows and
   * underflows in the computation.
   * <br/>
   * {@code Infinite} and {@code NaN} values are handled according to the
   * following rules, applied in the order presented:
   * <ul>
   *  <li>If either {@code this} or {@code divisor} has a {@code NaN} value
   *   in either part, {@link #NaN} is returned.
   *  </li>
   *  <li>If {@code divisor} equals {@link #ZERO}, {@link #NaN} is returned.
   *  </li>
   *  <li>If {@code this} and {@code divisor} are both infinite,
   *   {@link #NaN} is returned.
   *  </li>
   *  <li>If {@code this} is finite (i.e., has no {@code Infinite} or
   *   {@code NaN} parts) and {@code divisor} is infinite (one or both parts
   *   infinite), {@link #ZERO} is returned.
   *  </li>
   *  <li>If {@code this} is infinite and {@code divisor} is finite,
   *   {@code NaN} values are returned in the parts of the result if the
   *   {@link java.lang.Double} rules applied to the definitional formula
   *   force {@code NaN} results.
   *  </li>
   * </ul>
   */
  Complex operator /(Object divisor) {
    if (divisor is Complex) {
      if (isNaN || divisor.isNaN) {
        return NAN;
      }

      final double c = divisor.real;
      final double d = divisor.imaginary;
      if (c == 0.0 && d == 0.0) {
        return NAN;
      }

      if (divisor.isInfinite && !isInfinite) {
        return ZERO;
      }

      if (c.abs() < d.abs()) {
        double q = c / d;
        double denominator = c * q + d;
        return complex((_real * q + _imaginary) / denominator, (_imaginary * q - _real) / denominator);
      } else {
        double q = d / c;
        double denominator = d * q + c;
        return complex((_imaginary * q + _real) / denominator, (_imaginary - _real * q) / denominator);
      }
    } else if (divisor is num) {

      if (isNaN || divisor.isNaN) {
        return NAN;
      }
      if (divisor == 0) {
        return NAN;
      }
      if (divisor.isInfinite) {
        return !isInfinite ? ZERO : NAN;
      }
      return complex(_real / divisor, _imaginary / divisor);
    } else {
      throw new ArgumentError('factor must be a num or a Complex');
    }
  }

  /**
   * Returns the multiplicative inverse of {@code this} element.
   */
  Complex reciprocal() {
    if (isNaN) {
      return NAN;
    }

    if (_real == 0.0 && _imaginary == 0.0) {
      return INFINITY;
    }

    if (isInfinite) {
      return ZERO;
    }

    if (_real.abs() < _imaginary.abs()) {
      double q = _real / _imaginary;
      double scale = 1.0 / (_real * q + _imaginary);
      return complex(scale * q, -scale);
    } else {
      double q = _imaginary / _real;
      double scale = 1.0 / (_imaginary * q + _real);
      return complex(scale, -scale * q);
    }
  }

  /**
   * Test for equality with another object.
   * If both the real and imaginary parts of two complex numbers
   * are exactly the same, and neither is {@code Double.NaN}, the two
   * Complex objects are considered to be equal.
   * The behavior is the same as for JDK's {@link Double#equals(Object)
   * Double}:
   * <ul>
   *  <li>All {@code NaN} values are considered to be equal,
   *   i.e, if either (or both) real and imaginary parts of the complex
   *   number are equal to {@code Double.NaN}, the complex number is equal
   *   to {@code NaN}.
   *  </li>
   *  <li>
   *   Instances constructed with different representations of zero (i.e.
   *   either "0" or "-0") are <em>not</em> considered to be equal.
   *  </li>
   * </ul>
   */
  bool operator ==(Object other) {
    if (this == other) {
      return true;
    }
    if (other is Complex) {
      Complex c = other;
      if (c.isNaN) {
        return isNaN;
      } else {
        return _real == c.real && _imaginary == c.imaginary;
      }
    }
    return false;
  }


  /**
     * Returns a {@code Complex} whose value is {@code this * factor}.
     * Implements preliminary checks for {@code NaN} and infinity followed by
     * the definitional formula:
     * <pre>
     *  <code>
     *   (a + bi)(c + di) = (ac - bd) + (ad + bc)i
     *  </code>
     * </pre>
     * Returns {@link #NaN} if either {@code this} or {@code factor} has one or
     * more {@code NaN} parts.
     * <br/>
     * Returns {@link #INF} if neither {@code this} nor {@code factor} has one
     * or more {@code NaN} parts and if either {@code this} or {@code factor}
     * has one or more infinite parts (same result is returned regardless of
     * the sign of the components).
     * <br/>
     * Returns finite values in components of the result per the definitional
     * formula in all remaining cases.
     *
     * @param  factor value to be multiplied by this {@code Complex}.
     * @return {@code this * factor}.
     * @throws NullArgumentException if {@code factor} is {@code null}.
     */
  Complex operator *(Object factor) {
    if (factor is Complex) {
      if (isNaN || factor.isNaN) {
        return NAN;
      }
      if (_real.isInfinite || _imaginary.isInfinite || factor._real.isInfinite || factor._imaginary.isInfinite) {
        // we don't use isInfinite to avoid testing for NaN again
        return INFINITY;
      }
      return complex(_real * factor._real - _imaginary * factor._imaginary, _real * factor._imaginary + _imaginary * factor._real);
    } else if (factor is num) {
      if (isNaN || factor.isNaN) {
        return NAN;
      }
      if (_real.isInfinite || _imaginary.isInfinite || factor.isInfinite) {
        // we don't use isInfinite to avoid testing for NaN again
        return INFINITY;
      }
      return complex(_real * factor, _imaginary * factor);
    } else {
      throw new ArgumentError('factor must be a num or a Complex');
    }
  }

    /**
     * Returns a {@code Complex} whose value is {@code (-this)}.
     * Returns {@code NAN} if either real or imaginary
     * part of this Complex number equals {@code double.NAN}.
     */
    //Complex operator unary-() {
    Complex negate() {
        if (isNaN) {
            return NAN;
        }

        return complex(-_real, -_imaginary);
    }

    /**
     * Returns a {@code Complex} whose value is
     * {@code (this - subtrahend)}.
     * Uses the definitional formula
     * <pre>
     *  <code>
     *   (a + bi) - (c + di) = (a-c) + (b-d)i
     *  </code>
     * </pre>
     * If either {@code this} or {@code subtrahend} has a {@code NaN]} value in either part,
     * {@link #NaN} is returned; otherwise infinite and {@code NaN} values are
     * returned in the parts of the result according to the rules for
     * {@link java.lang.Double} arithmetic.
     *
     * @param  subtrahend value to be subtracted from this {@code Complex}.
     * @return {@code this - subtrahend}.
     * @throws NullArgumentException if {@code subtrahend} is {@code null}.
     */
    Complex operator -(Object subtrahend) {
      if (subtrahend is Complex) {
        if (isNaN || subtrahend.isNaN) {
            return NAN;
        }

        return complex(real - subtrahend._real,
                             imaginary - subtrahend._imaginary);
      } else if (subtrahend is num) {
        if (isNaN || subtrahend.isNaN) {
            return NAN;
        }
        return complex(real - subtrahend, imaginary);
      } else {
        throw new ArgumentError('factor must be a num or a Complex');
      }
    }


    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/InverseCosine.html" TARGET="_top">
     * inverse cosine</a> of this complex number.
     * Implements the formula:
     * <pre>
     *  <code>
     *   acos(z) = -i (log(z + i (sqrt(1 - z<sup>2</sup>))))
     *  </code>
     * </pre>
     * Returns {@link Complex#NaN} if either real or imaginary part of the
     * input argument is {@code NaN} or infinite.
     *
     * @return the inverse cosine of this complex number.
     * @since 1.2
     */
    Complex acos() {
        if (isNaN) {
            return NAN;
        }
        return this + (this.sqrt1z() * I).log() * (I.negate());
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/InverseSine.html" TARGET="_top">
     * inverse sine</a> of this complex number.
     * Implements the formula:
     * <pre>
     *  <code>
     *   asin(z) = -i (log(sqrt(1 - z<sup>2</sup>) + iz))
     *  </code>
     * </pre>
     * Returns {@link Complex#NaN} if either real or imaginary part of the
     * input argument is {@code NaN} or infinite.
     *
     * @return the inverse sine of this complex number.
     * @since 1.2
     */
    Complex asin() {
        if (isNaN) {
            return NAN;
        }

        return sqrt1z() + (this * I).log() * I.negate();
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/InverseTangent.html" TARGET="_top">
     * inverse tangent</a> of this complex number.
     * Implements the formula:
     * <pre>
     *  <code>
     *   atan(z) = (i/2) log((i + z)/(i - z))
     *  </code>
     * </pre>
     * Returns {@link Complex#NaN} if either real or imaginary part of the
     * input argument is {@code NaN} or infinite.
     *
     * @return the inverse tangent of this complex number
     * @since 1.2
     */
    Complex atan() {
        if (isNaN) {
            return NAN;
        }

        return (this + I) / (I - this).log() * (I / complex(2.0, 0.0));
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/Cosine.html" TARGET="_top">
     * cosine</a>
     * of this complex number.
     * Implements the formula:
     * <pre>
     *  <code>
     *   cos(a + bi) = cos(a)cosh(b) - sin(a)sinh(b)i
     *  </code>
     * </pre>
     * where the (real) functions on the right-hand side are
     * {@link FastMath#sin}, {@link FastMath#cos},
     * {@link FastMath#cosh} and {@link FastMath#sinh}.
     * <br/>
     * Returns {@link Complex#NaN} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * <br/>
     * Infinite values in real or imaginary parts of the input may result in
     * infinite or NaN values returned in parts of the result.
     * <pre>
     *  Examples:
     *  <code>
     *   cos(1 &plusmn; INFINITY i) = 1 &#x2213; INFINITY i
     *   cos(&plusmn;INFINITY + i) = NaN + NaN i
     *   cos(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
     *  </code>
     * </pre>
     *
     * @return the cosine of this complex number.
     * @since 1.2
     */
    Complex cos() {
        if (isNaN) {
            return NAN;
        }

        return complex(math.cos(real) * fastmath.cosh(imaginary),
                             -math.sin(real) * fastmath.sinh(imaginary));
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/HyperbolicCosine.html" TARGET="_top">
     * hyperbolic cosine</a> of this complex number.
     * Implements the formula:
     * <pre>
     *  <code>
     *   cosh(a + bi) = cosh(a)cos(b) + sinh(a)sin(b)i}
     *  </code>
     * </pre>
     * where the (real) functions on the right-hand side are
     * {@link FastMath#sin}, {@link FastMath#cos},
     * {@link FastMath#cosh} and {@link FastMath#sinh}.
     * <br/>
     * Returns {@link Complex#NaN} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * <br/>
     * Infinite values in real or imaginary parts of the input may result in
     * infinite or NaN values returned in parts of the result.
     * <pre>
     *  Examples:
     *  <code>
     *   cosh(1 &plusmn; INFINITY i) = NaN + NaN i
     *   cosh(&plusmn;INFINITY + i) = INFINITY &plusmn; INFINITY i
     *   cosh(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
     *  </code>
     * </pre>
     *
     * @return the hyperbolic cosine of this complex number.
     * @since 1.2
     */
    Complex cosh() {
        if (isNaN) {
            return NAN;
        }

        return complex(fastmath.cosh(real) * math.cos(imaginary),
                             fastmath.sinh(real) * math.sin(imaginary));
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/ExponentialFunction.html" TARGET="_top">
     * exponential function</a> of this complex number.
     * Implements the formula:
     * <pre>
     *  <code>
     *   exp(a + bi) = exp(a)cos(b) + exp(a)sin(b)i
     *  </code>
     * </pre>
     * where the (real) functions on the right-hand side are
     * {@link FastMath#exp}, {@link FastMath#cos}, and
     * {@link FastMath#sin}.
     * <br/>
     * Returns {@link Complex#NaN} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * <br/>
     * Infinite values in real or imaginary parts of the input may result in
     * infinite or NaN values returned in parts of the result.
     * <pre>
     *  Examples:
     *  <code>
     *   exp(1 &plusmn; INFINITY i) = NaN + NaN i
     *   exp(INFINITY + i) = INFINITY + INFINITY i
     *   exp(-INFINITY + i) = 0 + 0i
     *   exp(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
     *  </code>
     * </pre>
     *
     * @return <code><i>e</i><sup>this</sup></code>.
     * @since 1.2
     */
    Complex exp() {
        if (isNaN) {
            return NAN;
        }

        double expReal = math.exp(real);
        return complex(expReal *  math.cos(imaginary),
                             expReal * math.sin(imaginary));
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/NaturalLogarithm.html" TARGET="_top">
     * natural logarithm</a> of this complex number.
     * Implements the formula:
     * <pre>
     *  <code>
     *   log(a + bi) = ln(|a + bi|) + arg(a + bi)i
     *  </code>
     * </pre>
     * where ln on the right hand side is {@link FastMath#log},
     * {@code |a + bi|} is the modulus, {@link Complex#abs},  and
     * {@code arg(a + bi) = }{@link FastMath#atan2}(b, a).
     * <br/>
     * Returns {@link Complex#NaN} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * <br/>
     * Infinite (or critical) values in real or imaginary parts of the input may
     * result in infinite or NaN values returned in parts of the result.
     * <pre>
     *  Examples:
     *  <code>
     *   log(1 &plusmn; INFINITY i) = INFINITY &plusmn; (&pi;/2)i
     *   log(INFINITY + i) = INFINITY + 0i
     *   log(-INFINITY + i) = INFINITY + &pi;i
     *   log(INFINITY &plusmn; INFINITY i) = INFINITY &plusmn; (&pi;/4)i
     *   log(-INFINITY &plusmn; INFINITY i) = INFINITY &plusmn; (3&pi;/4)i
     *   log(0 + 0i) = -INFINITY + 0i
     *  </code>
     * </pre>
     *
     * @return the value <code>ln &nbsp; this</code>, the natural logarithm
     * of {@code this}.
     * @since 1.2
     */
    Complex log() {
        if (isNaN) {
            return NAN;
        }

        return complex(math.log(abs()),
                             fastmath.atan2(imaginary, real));
    }

    /**
     * Returns of value of this complex number raised to the power of {@code x}.
     * Implements the formula:
     * <pre>
     *  <code>
     *   y<sup>x</sup> = exp(x&middot;log(y))
     *  </code>
     * </pre>
     * where {@code exp} and {@code log} are {@link #exp} and
     * {@link #log}, respectively.
     * <br/>
     * Returns {@link Complex#NaN} if either real or imaginary part of the
     * input argument is {@code NaN} or infinite, or if {@code y}
     * equals {@link Complex#ZERO}.
     *
     * @param  x exponent to which this {@code Complex} is to be raised.
     * @return <code> this<sup>{@code x}</sup></code>.
     * @throws NullArgumentException if x is {@code null}.
     * @since 1.2
     */
    Complex power(Complex x) {
        return (this.log() * x).exp();
    }

    /**
     * Returns of value of this complex number raised to the power of {@code x}.
     *
     * @param  x exponent to which this {@code Complex} is to be raised.
     * @return <code>this<sup>x</sup></code>.
     * @see #pow(Complex)
     */
     Complex pow(double x) {
        return (this.log() * x).exp();
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/Sine.html" TARGET="_top">
     * sine</a>
     * of this complex number.
     * Implements the formula:
     * <pre>
     *  <code>
     *   sin(a + bi) = sin(a)cosh(b) - cos(a)sinh(b)i
     *  </code>
     * </pre>
     * where the (real) functions on the right-hand side are
     * {@link FastMath#sin}, {@link FastMath#cos},
     * {@link FastMath#cosh} and {@link FastMath#sinh}.
     * <br/>
     * Returns {@link Complex#NaN} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * <br/>
     * Infinite values in real or imaginary parts of the input may result in
     * infinite or {@code NaN} values returned in parts of the result.
     * <pre>
     *  Examples:
     *  <code>
     *   sin(1 &plusmn; INFINITY i) = 1 &plusmn; INFINITY i
     *   sin(&plusmn;INFINITY + i) = NaN + NaN i
     *   sin(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
     *  </code>
     * </pre>
     *
     * @return the sine of this complex number.
     * @since 1.2
     */
    Complex sin() {
        if (isNaN) {
            return NAN;
        }

        return complex(math.sin(real) * fastmath.cosh(imaginary),
                             math.cos(real) * fastmath.sinh(imaginary));
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/HyperbolicSine.html" TARGET="_top">
     * hyperbolic sine</a> of this complex number.
     * Implements the formula:
     * <pre>
     *  <code>
     *   sinh(a + bi) = sinh(a)cos(b)) + cosh(a)sin(b)i
     *  </code>
     * </pre>
     * where the (real) functions on the right-hand side are
     * {@link FastMath#sin}, {@link FastMath#cos},
     * {@link FastMath#cosh} and {@link FastMath#sinh}.
     * <br/>
     * Returns {@link Complex#NaN} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * <br/>
     * Infinite values in real or imaginary parts of the input may result in
     * infinite or NaN values returned in parts of the result.
     * <pre>
     *  Examples:
     *  <code>
     *   sinh(1 &plusmn; INFINITY i) = NaN + NaN i
     *   sinh(&plusmn;INFINITY + i) = &plusmn; INFINITY + INFINITY i
     *   sinh(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
     *  </code>
     * </pre>
     *
     * @return the hyperbolic sine of {@code this}.
     * @since 1.2
     */
    Complex sinh() {
        if (isNaN) {
            return NAN;
        }

        return complex(fastmath.sinh(real) * math.cos(imaginary),
            fastmath.cosh(real) * math.sin(imaginary));
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/SquareRoot.html" TARGET="_top">
     * square root</a> of this complex number.
     * Implements the following algorithm to compute {@code sqrt(a + bi)}:
     * <ol><li>Let {@code t = sqrt((|a| + |a + bi|) / 2)}</li>
     * <li><pre>if {@code  a &#8805; 0} return {@code t + (b/2t)i}
     *  else return {@code |b|/2t + sign(b)t i }</pre></li>
     * </ol>
     * where <ul>
     * <li>{@code |a| = }{@link FastMath#abs}(a)</li>
     * <li>{@code |a + bi| = }{@link Complex#abs}(a + bi)</li>
     * <li>{@code sign(b) =  }{@link FastMath#copySign(double,double) copySign(1d, b)}
     * </ul>
     * <br/>
     * Returns {@link Complex#NaN} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * <br/>
     * Infinite values in real or imaginary parts of the input may result in
     * infinite or NaN values returned in parts of the result.
     * <pre>
     *  Examples:
     *  <code>
     *   sqrt(1 &plusmn; INFINITY i) = INFINITY + NaN i
     *   sqrt(INFINITY + i) = INFINITY + 0i
     *   sqrt(-INFINITY + i) = 0 + INFINITY i
     *   sqrt(INFINITY &plusmn; INFINITY i) = INFINITY + NaN i
     *   sqrt(-INFINITY &plusmn; INFINITY i) = NaN &plusmn; INFINITY i
     *  </code>
     * </pre>
     *
     * @return the square root of {@code this}.
     * @since 1.2
     */
    Complex sqrt() {
        if (isNaN) {
            return NAN;
        }

        if (real == 0.0 && imaginary == 0.0) {
            return complex(0.0, 0.0);
        }

        double t = math.sqrt((real.abs() + abs()) / 2.0);
        if (real >= 0.0) {
            return complex(t, imaginary / (2.0 * t));
        } else {
            return complex(imaginary.abs() / (2.0 * t),
                                 fastmath.copySign(1.0, imaginary) * t);
        }
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/SquareRoot.html" TARGET="_top">
     * square root</a> of <code>1 - this<sup>2</sup></code> for this complex
     * number.
     * Computes the result directly as
     * {@code sqrt(ONE.subtract(z.multiply(z)))}.
     * <br/>
     * Returns {@link Complex#NaN} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * <br/>
     * Infinite values in real or imaginary parts of the input may result in
     * infinite or NaN values returned in parts of the result.
     *
     * @return the square root of <code>1 - this<sup>2</sup></code>.
     * @since 1.2
     */
    Complex sqrt1z() {
        return (complex(1.0, 0.0) - (this * this)).sqrt();
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/Tangent.html" TARGET="_top">
     * tangent</a> of this complex number.
     * Implements the formula:
     * <pre>
     *  <code>
     *   tan(a + bi) = sin(2a)/(cos(2a)+cosh(2b)) + [sinh(2b)/(cos(2a)+cosh(2b))]i
     *  </code>
     * </pre>
     * where the (real) functions on the right-hand side are
     * {@link FastMath#sin}, {@link FastMath#cos}, {@link FastMath#cosh} and
     * {@link FastMath#sinh}.
     * <br/>
     * Returns {@link Complex#NaN} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * <br/>
     * Infinite (or critical) values in real or imaginary parts of the input may
     * result in infinite or NaN values returned in parts of the result.
     * <pre>
     *  Examples:
     *  <code>
     *   tan(a &plusmn; INFINITY i) = 0 &plusmn; i
     *   tan(&plusmn;INFINITY + bi) = NaN + NaN i
     *   tan(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
     *   tan(&plusmn;&pi;/2 + 0 i) = &plusmn;INFINITY + NaN i
     *  </code>
     * </pre>
     *
     * @return the tangent of {@code this}.
     * @since 1.2
     */
    Complex tan() {
        if (isNaN || real.isInfinite) {
            return NAN;
        }
        if (imaginary > 20.0) {
            return complex(0.0, 1.0);
        }
        if (imaginary < -20.0) {
            return complex(0.0, -1.0);
        }

        double real2 = 2.0 * real;
        double imaginary2 = 2.0 * imaginary;
        double d = math.cos(real2) + fastmath.cosh(imaginary2);

        return complex(math.sin(real2) / d,
                             fastmath.sinh(imaginary2) / d);
    }

    /**
     * Compute the
     * <a href="http://mathworld.wolfram.com/HyperbolicTangent.html" TARGET="_top">
     * hyperbolic tangent</a> of this complex number.
     * Implements the formula:
     * <pre>
     *  <code>
     *   tan(a + bi) = sinh(2a)/(cosh(2a)+cos(2b)) + [sin(2b)/(cosh(2a)+cos(2b))]i
     *  </code>
     * </pre>
     * where the (real) functions on the right-hand side are
     * {@link FastMath#sin}, {@link FastMath#cos}, {@link FastMath#cosh} and
     * {@link FastMath#sinh}.
     * <br/>
     * Returns {@link Complex#NaN} if either real or imaginary part of the
     * input argument is {@code NaN}.
     * <br/>
     * Infinite values in real or imaginary parts of the input may result in
     * infinite or NaN values returned in parts of the result.
     * <pre>
     *  Examples:
     *  <code>
     *   tanh(a &plusmn; INFINITY i) = NaN + NaN i
     *   tanh(&plusmn;INFINITY + bi) = &plusmn;1 + 0 i
     *   tanh(&plusmn;INFINITY &plusmn; INFINITY i) = NaN + NaN i
     *   tanh(0 + (&pi;/2)i) = NaN + INFINITY i
     *  </code>
     * </pre>
     *
     * @return the hyperbolic tangent of {@code this}.
     * @since 1.2
     */
    Complex tanh() {
        if (isNaN || imaginary.isInfinite) {
            return NAN;
        }
        if (real > 20.0) {
            return complex(1.0, 0.0);
        }
        if (real < -20.0) {
            return complex(-1.0, 0.0);
        }
        double real2 = 2.0 * real;
        double imaginary2 = 2.0 * imaginary;
        double d = fastmath.cosh(real2) + math.cos(imaginary2);

        return complex(fastmath.sinh(real2) / d,
                             math.sin(imaginary2) / d);
    }



    /**
     * Compute the argument of this complex number.
     * The argument is the angle phi between the positive real axis and
     * the point representing this number in the complex plane.
     * The value returned is between -PI (not inclusive)
     * and PI (inclusive), with negative values returned for numbers with
     * negative imaginary parts.
     * <br/>
     * If either real or imaginary part (or both) is NaN, NaN is returned.
     * Infinite parts are handled as {@code Math.atan2} handles them,
     * essentially treating finite parts as zero in the presence of an
     * infinite coordinate and returning a multiple of pi/4 depending on
     * the signs of the infinite parts.
     * See the javadoc for {@code Math.atan2} for full details.
     *
     * @return the argument of {@code this}.
     */
    double argument() {
        return fastmath.atan2(imaginary, real);
    }

    /**
     * Computes the n-th roots of this complex number.
     * The nth roots are defined by the formula:
     * <pre>
     *  <code>
     *   z<sub>k</sub> = abs<sup>1/n</sup> (cos(phi + 2&pi;k/n) + i (sin(phi + 2&pi;k/n))
     *  </code>
     * </pre>
     * for <i>{@code k=0, 1, ..., n-1}</i>, where {@code abs} and {@code phi}
     * are respectively the {@link #abs() modulus} and
     * {@link #getArgument() argument} of this complex number.
     * <br/>
     * If one or both parts of this complex number is NaN, a list with just
     * one element, {@link #NaN} is returned.
     * if neither part is NaN, but at least one part is infinite, the result
     * is a one-element list containing {@link #INF}.
     *
     * @param n Degree of root.
     * @return a List<Complex> of all {@code n}-th roots of {@code this}.
     * @throws NotPositiveException if {@code n <= 0}.
     * @since 2.0
     */
    List<Complex> nthRoot(int n) {
        if (n <= 0) {
            throw new ArgumentError("Can't compute nth root for negative n");
        }

        final result = new List<Complex>();

        if (isNaN) {
            result.add(NAN);
            return result;
        }
        if (isInfinite) {
            result.add(INFINITY);
            return result;
        }

        // nth root of abs -- faster / more accurate to use a solver here?
        final double nthRootOfAbs = math.pow(abs(), 1.0 / n);

        // Compute nth roots of complex number with k = 0, 1, ... n-1
        final double nthPhi = argument() / n;
        final double slice = 2 * math.PI / n;
        double innerPart = nthPhi;
        for (int k = 0; k < n ; k++) {
            // inner part
            final double realPart = nthRootOfAbs * math.cos(innerPart);
            final double imaginaryPart = nthRootOfAbs * math.sin(innerPart);
            result.add(complex(realPart, imaginaryPart));
            innerPart += slice;
        }

        return result;
    }

    String toString() {
        return "($real, $imaginary)";
    }

}

/**
 * Create a complex number given the real part and optionally the imaginary part.
 */
Complex complex(double real, [double imaginary = 0.0]) {
  return new Complex(real, imaginary);
}
