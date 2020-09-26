part of '../complex.dart';

extension ComplexTrigonometricX<T extends Complex> on T {
  /// Compute the [sine](http://mathworld.wolfram.com/Sine.html)
  /// of this complex number.
  ///
  /// Implements the formula:
  ///
  ///     sin(a + bi) = sin(a)cosh(b) - cos(a)sinh(b)i
  ///
  /// where the (real) functions on the right-hand side are
  /// [math.sin], [math.cos], [fastmath.cosh] and [fastmath.sinh].
  ///
  /// Returns [nan] if either real or imaginary part of the
  /// input argument is `NaN`.
  ///
  /// Infinite values in real or imaginary parts of the input may result in
  /// infinite or `NaN` values returned in parts of the result.
  ///
  /// Examples:
  ///
  ///     sin(1 ± INFINITY i) = 1 ± INFINITY i
  ///     sin(±INFINITY + i) = NaN + NaN i
  ///     sin(±INFINITY ± INFINITY i) = NaN + NaN i
  Complex sin() {
    if (isNaN) {
      return Complex.nan;
    }
    final _real = math.sin(real) * fastmath.cosh(imaginary);
    final _imag = math.cos(real) * fastmath.sinh(imaginary);
    return Cartesian(_real, _imag);
  }

  /// Compute the [cosine](http://mathworld.wolfram.com/Cosine.html)
  /// of this complex number.
  ///
  /// Implements the formula:
  ///
  ///     cos(a + bi) = cos(a)cosh(b) - sin(a)sinh(b)i
  ///
  /// where the (real) functions on the right-hand side are
  /// [math.sin], [math.cos], [fastmath.cosh] and [fastmath.sinh].
  ///
  /// Returns [nan] if either real or imaginary part of the
  /// input argument is `NaN`.
  ///
  /// Infinite values in real or imaginary parts of the input may result in
  /// infinite or `NaN` values returned in parts of the result.
  ///
  /// Examples:
  ///
  ///     cos(1 ± INFINITY i) = 1 ∓ INFINITY i
  ///     cos(±INFINITY + i) = NaN + NaN i
  ///     cos(±INFINITY ± INFINITY i) = NaN + NaN i
  Complex cos() {
    if (isNaN) return Complex.nan;

    final _real = math.cos(real) * fastmath.cosh(imaginary);
    final _imag = -math.sin(real) * fastmath.sinh(imaginary);
    return Cartesian(_real, _imag);
  }

  /// Compute the [tangent](http://mathworld.wolfram.com/Tangent.html)
  /// of this complex number.
  ///
  /// Implements the formula:
  ///
  ///     tan(a + bi) = sin(2a)/(cos(2a)+cosh(2b)) + [sinh(2b)/(cos(2a)+cosh(2b))]i
  ///
  /// where the (real) functions on the right-hand side are
  /// [math.sin], [math.cos], [fastmath.cosh] and [fastmath.sinh].
  ///
  /// Returns [nan] if either real or imaginary part of the
  /// input argument is `NaN`.
  ///
  /// Infinite (or critical) values in real or imaginary parts of the input may
  /// result in infinite or `NaN` values returned in parts of the result.
  ///
  /// Examples:
  ///
  ///     tan(a ± INFINITY i) = 0 ± i
  ///     tan(±INFINITY + bi) = NaN + NaN i
  ///     tan(±INFINITY ± INFINITY i) = NaN + NaN i
  ///     tan(±π/2 + 0 i) = ±INFINITY + NaN i
  Complex tan() {
    if (isNaN || real.isInfinite) {
      return Complex.nan;
    }
    if (imaginary > 20.0) {
      return Cartesian(0.0, 1.0);
    }
    if (imaginary < -20.0) {
      return Cartesian(0.0, -1.0);
    }

    double real2 = 2.0 * real;
    double imaginary2 = 2.0 * imaginary;
    double d = math.cos(real2) + fastmath.cosh(imaginary2);

    return Cartesian(math.sin(real2) / d, fastmath.sinh(imaginary2) / d);
  }

  /// Compute the [inverse cosine](http://mathworld.wolfram.com/InverseCosine.html)
  /// of this complex number.
  ///
  /// Implements the formula:
  ///
  ///     acos(z) = -i (log(z + i (sqrt(1 - z^2))))
  ///
  /// Returns [nan] if either real or imaginary part of the
  /// input argument is `NaN` or infinite.
  Complex acos() {
    if (isNaN) return Complex.nan;
    return (this + (this.sqrt1z() * Complex.i)).log() * (-Complex.i);
  }

  /// Compute the [inverse sine](http://mathworld.wolfram.com/InverseSine.html)
  /// of this complex number.
  ///
  /// Implements the formula:
  ///
  ///     asin(z) = -i (log(sqrt(1 - z^2) + iz))
  ///
  /// Returns [nan] if either real or imaginary part of the
  /// input argument is `NaN` or infinite.
  Complex asin() {
    if (isNaN) return Complex.nan;
    return (sqrt1z() + (this * Complex.i)).log() * -Complex.i;
  }

  /// Compute the [inverse tangent](http://mathworld.wolfram.com/InverseTangent.html)
  /// of this complex number.
  ///
  /// Implements the formula:
  ///
  ///     atan(z) = (i/2) log((i + z)/(i - z))
  ///
  /// Returns [nan] if either real or imaginary part of the
  /// input argument is `NaN` or infinite.
  Complex atan() {
    if (isNaN) return Complex.nan;
    // TODO(kranfix): Optimize with Imaginary
    return ((this + Complex.i) / (Complex.i - this)).log() *
        (Complex.i / Cartesian(2.0, 0.0));
  }
}
