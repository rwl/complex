part of '../complex.dart';

/// Representation of a complex number, i.e. a number which has both a
/// real and imaginary part.
///
/// Implementations of arithmetic operations handle `NaN` and
/// infinite values according to the rules for [double], i.e.
/// [==] is an equivalence relation for all instances that have
/// a `NaN` in either real or imaginary part, e.g. the following are
/// considered equal:
///
/// * `1 + NaNi`
/// * `NaN + i`
/// * `NaN + NaNi`
///
/// Note that this is in contradiction with the IEEE-754 standard for floating
/// point numbers (according to which the test `x == x` must fail if
/// `x` is `NaN`).
abstract class Complex {
  /// The square root of -1. A number representing "0.0 + 1.0i"
  static const i = Cartesian(0.0, 1.0);

  /// A complex number representing "NaN + NaNi"
  static const nan = Cartesian(double.nan, double.nan);

  /// A complex number representing "+INF + INFi"
  static const infinity = Cartesian(double.infinity, double.infinity);

  /// A complex number representing "1.0 + 0.0i"
  static const one = Cartesian(1.0);

  /// A complex number representing "0.0 + 0.0i"
  static const zero = Cartesian(0.0);

  /// [Complex] as a [Cartesian] with optional `imaginary` part.
  ///
  /// Example:
  /// ```dart
  /// const realComplex = Complex(3.0);
  /// const complexWithImaginary = Complex(5.3, 4.5);
  /// ```
  const factory Complex(double real, [double imaginary]) = Cartesian;

  factory Complex.imag(double imaginary) => Cartesian(0.0, imaginary);

  factory Complex.polar(double r, double phase) {
    r ??= 0.0;
    phase ??= 0.0;
    return Cartesian(r * math.cos(phase), r * math.sin(phase));
  }

  /// The real part.
  double get real;

  /// The imaginary part.
  double get imaginary;

  /// True if the real and imaginary parts are finite; otherwise, false.
  bool get isFinite;

  /// True if the real and imaginary parts are positive infinity or negative
  /// infinity; otherwise, false.
  bool get isInfinite;

  /// True if the real and imaginary parts are the double Not-a-Number value;
  /// otherwise, false.
  bool get isNaN;

  /// Returns a `Complex` whose value is `this * factor`.
  /// Implements preliminary checks for `NaN` and infinity followed by
  /// the definitional formula:
  ///
  ///     (a + bi)(c + di) = (ac - bd) + (ad + bc)i
  ///
  /// Returns [nan] if either `this` or `factor` has one or
  /// more `NaN` parts.
  ///
  /// Returns [infinity] if neither `this` nor `factor` has one
  /// or more `NaN` parts and if either `this` or `factor`
  /// has one or more infinite parts (same result is returned regardless of
  /// the sign of the components).
  ///
  /// Returns finite values in components of the result per the definitional
  /// formula in all remaining cases.
  Complex operator *(Object factor);

  /// Returns a `Complex` whose value is
  /// `(this + addend)`.
  /// Uses the definitional formula
  ///
  ///     (a + bi) + (c + di) = (a+c) + (b+d)i
  ///
  /// If either `this` or [addend] has a `NaN` value in
  /// either part, [NaN] is returned; otherwise `Infinite`
  /// and `NaN` values are returned in the parts of the result
  /// according to the rules for [double] arithmetic.
  Complex operator +(Object addend);

  /// Returns a `Complex` whose value is
  /// `this - subtrahend`.
  /// Uses the definitional formula
  ///
  ///     (a + bi) - (c + di) = (a-c) + (b-d)i
  ///
  /// If either `this` or `subtrahend` has a `NaN` value in either part,
  /// [nan] is returned; otherwise infinite and `NaN` values are
  /// returned in the parts of the result according to the rules for
  /// [double] arithmetic.
  Complex operator -(Object subtrahend);

  /// Negate operator. Returns a `Complex` whose value is `-this`.
  /// Returns `NAN` if either real or imaginary
  /// part of this complex number equals `double.nan`.
  Complex operator -();

  /// Returns a `Complex` whose value is
  /// `(this / divisor)`.
  /// Implements the definitional formula
  ///
  ///       a + bi       ac + bd + (bc - ad)i
  ///     ---------- = ------------------------
  ///       c + di           c^2 + d^2
  ///
  /// but uses [prescaling of operands](http://doi.acm.org/10.1145/1039813.1039814)
  /// to limit the effects of overflows and underflows in the computation.
  ///
  /// `Infinite` and `NaN` values are handled according to the
  /// following rules, applied in the order presented:
  ///
  /// * If either `this` or `divisor` has a `NaN` value
  ///   in either part, [nan] is returned.
  /// * If `divisor` equals [zero], [nan] is returned.
  /// * If `this` and `divisor` are both infinite,
  ///   [nan] is returned.
  /// * If `this` is finite (i.e., has no `Infinite` or
  ///   `NAN` parts) and `divisor` is infinite (one or both parts
  ///   infinite), [zero] is returned.
  /// * If `this` is infinite and `divisor` is finite,
  /// `NAN` values are returned in the parts of the result if the
  /// [double] rules applied to the definitional formula
  /// force `NaN` results.
  Complex operator /(Object divisor);

  /// Test for equality with another object.
  ///
  /// If both the real and imaginary parts of two complex numbers
  /// are exactly the same, and neither is `double.nan`, the two
  /// Complex objects are considered to be equal.
  ///
  /// * All `NaN` values are considered to be equal,
  ///   i.e, if either (or both) real and imaginary parts of the complex
  ///   number are equal to `double.nan`, the complex number is equal
  ///   to `NaN`.
  /// * Instances constructed with different representations of zero (i.e.
  ///   either "0" or "-0") are *not* considered to be equal.
  bool operator ==(Object other);

  /// Return the absolute value of this complex number.
  /// Returns `NaN` if either real or imaginary part is `NaN`
  /// and `double.INFINITY` if neither part is `NaN`,
  /// but at least one part is infinite.
  double abs();

  /// Compute the argument of this complex number.
  ///
  /// The argument is the angle phi between the positive real axis and
  /// the point representing this number in the complex plane.
  /// The value returned is between -PI (not inclusive)
  /// and PI (inclusive), with negative values returned for numbers with
  /// negative imaginary parts.
  ///
  /// If either real or imaginary part (or both) is NaN, NaN is returned.
  /// Infinite parts are handled as `Math.atan2} handles them,
  /// essentially treating finite parts as zero in the presence of an
  /// infinite coordinate and returning a multiple of pi/4 depending on
  /// the signs of the infinite parts.
  double argument();

  /// Return the conjugate of this complex number.
  /// The conjugate of `a + bi` is `a - bi`.
  ///
  /// [NaN] is returned if either the real or imaginary
  /// part of this Complex number equals `double.nan`.
  ///
  /// If the imaginary part is infinite, and the real part is not
  /// `NaN`, the returned value has infinite imaginary part
  /// of the opposite sign, e.g. the conjugate of
  /// `1 + INFINITY i is `1 - NEGATIVE_INFINITY i`.
  Complex conjugate();

  /// Compute the [exponential function](http://mathworld.wolfram.com/ExponentialFunction.html)
  /// of this complex number.
  ///
  /// Implements the formula:
  ///
  ///     exp(a + bi) = exp(a)cos(b) + exp(a)sin(b)i
  ///
  /// where the (real) functions on the right-hand side are
  /// [math.exp], [math.cos], and [math.sin].
  ///
  /// Returns [nan] if either real or imaginary part of the
  /// input argument is `NaN`.
  ///
  /// Infinite values in real or imaginary parts of the input may result in
  /// infinite or `NaN` values returned in parts of the result.
  ///
  /// Examples:
  ///
  ///     exp(1 ± INFINITY i) = NaN + NaN i
  ///     exp(INFINITY + i) = INFINITY + INFINITY i
  ///     exp(-INFINITY + i) = 0 + 0i
  ///     exp(±INFINITY ± INFINITY i) = NaN + NaN i
  Complex exp();

  /// Compute the [natural logarithm](http://mathworld.wolfram.com/NaturalLogarithm.html)
  /// of this complex number.
  ///
  /// Implements the formula:
  ///
  ///     log(a + bi) = ln(|a + bi|) + arg(a + bi)i
  ///
  /// where ln on the right hand side is [math.log],
  /// `|a + bi|` is the modulus, [abs],  and
  /// `arg(a + bi) = atan2(b, a).
  ///
  /// Returns [nan] if either real or imaginary part of the
  /// input argument is `NaN`.
  ///
  /// Infinite (or critical) values in real or imaginary parts of the input may
  /// result in infinite or NaN values returned in parts of the result.
  ///
  /// Examples:
  ///
  ///     log(1 ± INFINITY i) = INFINITY ± (π/2)i
  ///     log(INFINITY + i) = INFINITY + 0i
  ///     log(-INFINITY + i) = INFINITY + &pi;i
  ///     log(INFINITY ± INFINITY i) = INFINITY ± (π/4)i
  ///     log(-INFINITY ± INFINITY i) = INFINITY ± (3π/4)i
  ///     log(0 + 0i) = -INFINITY + 0i
  Complex log();

  /// Returns of value of this complex number raised to the power of `x`.
  Complex pow(num x);

  /// Returns of value of this complex number raised to the power of `x`.
  ///
  /// Implements the formula:
  ///
  ///     y^x = exp(x·log(y))
  ///
  /// where `exp` and `log` are [exp] and
  /// [log], respectively.
  ///
  /// Returns [nan] if either real or imaginary part of the
  /// input argument is `NaN` or infinite, or if `y`
  /// equals [zero].
  Complex power<T extends Complex>(T x);

  /// Returns the multiplicative inverse of `this` element.
  Complex reciprocal();

  /// Compute the [square root](http://mathworld.wolfram.com/SquareRoot.html)
  /// of this complex number.
  ///
  /// Implements the following algorithm to compute `sqrt(a + bi)}:
  ///
  /// 1. Let `t = sqrt((|a| + |a + bi|) / 2)`
  /// 2. if `a ≥ 0` return `t + (b/2t)i`
  ///    else return `|b|/2t + sign(b)t i`
  ///
  /// where:
  ///
  /// * `|a| = abs(a)`
  /// * `|a + bi| = abs(a + bi)`
  /// * `sign(b) = copySign(double, double) copySign(1d, b)`
  ///
  /// Returns [nan] if either real or imaginary part of the
  /// input argument is `NaN`.
  ///
  /// Infinite values in real or imaginary parts of the input may result in
  /// infinite or NaN values returned in parts of the result.
  ///
  /// Examples:
  ///
  ///     sqrt(1 ± INFINITY i) = INFINITY + NaN i
  ///     sqrt(INFINITY + i) = INFINITY + 0i
  ///     sqrt(-INFINITY + i) = 0 + INFINITY i
  ///     sqrt(INFINITY ± INFINITY i) = INFINITY + NaN i
  ///     sqrt(-INFINITY ± INFINITY i) = NaN ± INFINITY i
  Complex sqrt();

  /// Compute the [square root](http://mathworld.wolfram.com/SquareRoot.html)
  /// of `1 - this^2` for this complex number.
  ///
  /// Computes the result directly as `sqrt(ONE - (z * z))`.
  ///
  /// Returns [nan] if either real or imaginary part of the
  /// input argument is `NaN`.
  ///
  /// Infinite values in real or imaginary parts of the input may result in
  /// infinite or NaN values returned in parts of the result.
  Complex sqrt1z();
}

extension ComplexerX on num {
  Complex get real => Cartesian(this.toDouble());
  Complex get imag => Cartesian(0, this.toDouble());
}
