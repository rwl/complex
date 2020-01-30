part of '../complex.dart';

extension ComplexHyperbolicX<T extends Complex> on T {
  /// Compute the [hyperbolic sine](http://mathworld.wolfram.com/HyperbolicSine.html)
  /// of this complex number.
  ///
  /// Implements the formula:
  ///
  ///     sinh(a + bi) = sinh(a)cos(b)) + cosh(a)sin(b)i
  ///
  /// where the (real) functions on the right-hand side are
  /// [math.sin], [math.cos], [fastmath.cosh] and [fastmath.sinh].
  ///
  /// Returns [nan] if either real or imaginary part of the
  /// input argument is `NaN`.
  ///
  /// Infinite values in real or imaginary parts of the input may result in
  /// infinite or NaN values returned in parts of the result.
  ///
  /// Examples:
  ///
  ///     sinh(1 ± INFINITY i) = NaN + NaN i
  ///     sinh(±INFINITY + i) = ± INFINITY + INFINITY i
  ///     sinh(±INFINITY ± INFINITY i) = NaN + NaN i
  Complex sinh() {
    if (isNaN) return Complex.nan;

    final _real = fastmath.sinh(real) * math.cos(imaginary);
    final _imag = fastmath.cosh(real) * math.sin(imaginary);
    return Cartesian(_real, _imag);
  }

  /// Compute the [hyperbolic cosine](http://mathworld.wolfram.com/HyperbolicCosine.html)
  /// of this complex number.
  ///
  /// Implements the formula:
  ///
  ///     cosh(a + bi) = cosh(a)cos(b) + sinh(a)sin(b)i
  ///
  /// where the (real) functions on the right-hand side are
  /// [math.sin], [math.cos], [fastmath.cosh] and [fastmath.sinh].
  ///
  /// Returns [nan] if either real or imaginary part of the
  /// input argument is `NaN`.
  ///
  /// Infinite values in real or imaginary parts of the input may result in
  /// infinite or NaN values returned in parts of the result.
  ///
  /// Examples:
  ///
  ///     cosh(1 ± INFINITY i) = NaN + NaN i
  ///     cosh(±INFINITY + i) = INFINITY ± INFINITY i
  ///     cosh±INFINITY &plusmn; INFINITY i) = NaN + NaN i
  Complex cosh() {
    if (isNaN) {
      return Complex.nan;
    }
    final _real = fastmath.cosh(real) * math.cos(imaginary);
    final _imag = fastmath.sinh(real) * math.sin(imaginary);
    return Cartesian(_real, _imag);
  }

  /// Compute the [hyperbolic tangent](http://mathworld.wolfram.com/HyperbolicTangent.html)
  /// of this complex number.
  ///
  /// Implements the formula:
  ///
  ///     tan(a + bi) = sinh(2a)/(cosh(2a)+cos(2b)) + [sin(2b)/(cosh(2a)+cos(2b))]i
  ///
  /// where the (real) functions on the right-hand side are
  /// [math.sin], [math.cos], [fastmath.cosh] and [fastmath.sinh].
  ///
  /// Returns [nan] if either real or imaginary part of the
  /// input argument is `NaN`.
  ///
  /// Infinite values in real or imaginary parts of the input may result in
  /// infinite or NaN values returned in parts of the result.
  ///
  /// Examples:
  ///
  ///     tanh(a ± INFINITY i) = NaN + NaN i
  ///     tanh(±INFINITY + bi) = ±1 + 0 i
  ///     tanh(±INFINITY ± INFINITY i) = NaN + NaN i
  ///     tanh(0 + (π/2)i) = NaN + INFINITY i
  Complex tanh() {
    if (isNaN || imaginary.isInfinite) {
      return Complex.nan;
    }
    if (real > 20.0) {
      return Cartesian(1.0, 0.0);
    }
    if (real < -20.0) {
      return Cartesian(-1.0, 0.0);
    }
    double real2 = 2.0 * real;
    double imaginary2 = 2.0 * imaginary;
    double d = fastmath.cosh(real2) + math.cos(imaginary2);

    return Cartesian(fastmath.sinh(real2) / d, math.sin(imaginary2) / d);
  }
}
