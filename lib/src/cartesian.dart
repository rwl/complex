part of '../complex.dart';

class Cartesian implements Complex {
  final double imaginary, real;

  /// Create a complex number given the real part and optionally the imaginary
  /// part.
  const Cartesian(this.real, [this.imaginary = 0]);

  @override
  bool get isFinite => !isNaN && real.isFinite && imaginary.isFinite;

  @override
  bool get isInfinite => !isNaN && (real.isInfinite || imaginary.isInfinite);

  @override
  bool get isNaN => real.isNaN || imaginary.isNaN;

  @override
  double abs() {
    if (isNaN) return double.nan;
    if (isInfinite) return double.infinity;

    double x = real.abs();
    double y = imaginary.abs();

    if (x > y) {
      final z = x;
      x = y;
      y = z;
    }
    if (x == 0.0) return y;
    double q = x / y;
    return y * math.sqrt(1 + q * q);
  }

  @override
  Cartesian operator +(Object addend) {
    if (isNaN) return Complex.nan;
    if (addend is Complex) {
      if (addend.isNaN) return Complex.nan;
      return Cartesian(real + addend.real, imaginary + addend.imaginary);
    } else if (addend is num) {
      if (addend.isNaN) return Complex.nan;
      return Cartesian(real + addend, imaginary);
    } else {
      throw ArgumentError('factor must be a num or a Complex');
    }
  }

  @override
  Cartesian conjugate() {
    if (isNaN) return Complex.nan;
    return Cartesian(real, -imaginary);
  }

  @override
  Cartesian operator /(Object divisor) {
    if (divisor is Complex) {
      if (isNaN || divisor.isNaN) {
        return Complex.nan;
      }

      final double c = divisor.real;
      final double d = divisor.imaginary;
      if (c == 0.0 && d == 0.0) {
        return Complex.nan;
      }

      if (divisor.isInfinite && !isInfinite) {
        return Complex.zero;
      }

      if (c.abs() < d.abs()) {
        double q = c / d;
        double denominator = c * q + d;
        return Cartesian(
          (real * q + imaginary) / denominator,
          (imaginary * q - real) / denominator,
        );
      } else {
        double q = d / c;
        double denominator = d * q + c;
        return Cartesian(
          (imaginary * q + real) / denominator,
          (imaginary - real * q) / denominator,
        );
      }
    } else if (divisor is num) {
      if (isNaN || divisor.isNaN) {
        return Complex.nan;
      }
      if (divisor == 0) {
        return Complex.nan;
      }
      if (divisor.isInfinite) {
        return isFinite ? Complex.zero : Complex.nan;
      }
      return Cartesian(real / divisor, imaginary / divisor);
    } else {
      throw ArgumentError('factor must be a num or a Complex');
    }
  }

  @override
  Cartesian reciprocal() {
    if (isNaN) {
      return Complex.nan;
    }

    if (real == 0.0 && imaginary == 0.0) {
      return Complex.infinity;
    }

    if (isInfinite) {
      return Complex.zero;
    }

    if (real.abs() < imaginary.abs()) {
      double q = real / imaginary;
      double scale = 1.0 / (real * q + imaginary);
      return Cartesian(scale * q, -scale);
    } else {
      double q = imaginary / real;
      double scale = 1.0 / (imaginary * q + real);
      return Cartesian(scale, -scale * q);
    }
  }

  @override
  bool operator ==(Object other) {
    if (identical(this, other)) {
      return true;
    }
    if (other is Cartesian) {
      Cartesian c = other;
      if (c.isNaN) {
        return isNaN;
      } else {
        return real == c.real && imaginary == c.imaginary;
      }
    }
    return false;
  }

  @override
  Cartesian operator *(Object factor) {
    if (factor is Cartesian) {
      if (isNaN || factor.isNaN) {
        return Complex.nan;
      }
      if (real.isInfinite ||
          imaginary.isInfinite ||
          factor.real.isInfinite ||
          factor.imaginary.isInfinite) {
        // we don't use isInfinite to avoid testing for NaN again
        return Complex.infinity;
      }
      return Cartesian(real * factor.real - imaginary * factor.imaginary,
          real * factor.imaginary + imaginary * factor.real);
    } else if (factor is num) {
      if (isNaN || factor.isNaN) {
        return Complex.nan;
      }
      if (real.isInfinite || imaginary.isInfinite || factor.isInfinite) {
        // we don't use isInfinite to avoid testing for NaN again
        return Complex.infinity;
      }
      return Cartesian(real * factor, imaginary * factor);
    } else {
      throw ArgumentError('factor must be a num or a Complex');
    }
  }

  @override
  Cartesian operator -() {
    if (isNaN) {
      return Complex.nan;
    }

    return Cartesian(-real, -imaginary);
  }

  @override
  Cartesian operator -(Object subtrahend) {
    if (subtrahend is Cartesian) {
      if (isNaN || subtrahend.isNaN) {
        return Complex.nan;
      }

      return Cartesian(
          real - subtrahend.real, imaginary - subtrahend.imaginary);
    } else if (subtrahend is num) {
      if (isNaN || subtrahend.isNaN) {
        return Complex.nan;
      }
      return Cartesian(real - subtrahend, imaginary);
    } else {
      throw ArgumentError('factor must be a num or a Complex');
    }
  }

  @override
  Cartesian exp() {
    if (isNaN) return Complex.nan;

    double expReal = math.exp(real);
    return Cartesian(
        expReal * math.cos(imaginary), expReal * math.sin(imaginary));
  }

  @override
  Cartesian log() {
    if (isNaN) Complex.nan;
    return Cartesian(math.log(abs()), fastmath.atan2(imaginary, real));
  }

  @override
  Cartesian power<T extends Complex>(T z) => (this.log() * z).exp();

  @override
  Cartesian pow(num x) => (this.log() * x).exp();

  @override
  Cartesian sqrt() {
    if (isNaN) Complex.nan;

    if (real == 0.0 && imaginary == 0.0) {
      return Complex.zero;
    }
    double t = math.sqrt((real.abs() + abs()) / 2.0);
    if (real >= 0.0) {
      return Cartesian(t, imaginary / (2.0 * t));
    } else {
      return Cartesian(
          imaginary.abs() / (2.0 * t), fastmath.copySign(1.0, imaginary) * t);
    }
  }

  @override
  Cartesian sqrt1z() {
    final a = (Complex.one - (this * this)).sqrt();
    return a;
  }

  @override
  double argument() => fastmath.atan2(imaginary, real);

  /// Get a hashCode for the complex number.
  ///
  /// Any [double.nan] value in real or imaginary part produces
  /// the same hash code `7`.
  @override
  int get hashCode {
    if (isNaN) return 7;
    return 37 * (17 * imaginary.hashCode + real.hashCode);
  }

  @override
  String toString() => "($real, $imaginary)";
}
