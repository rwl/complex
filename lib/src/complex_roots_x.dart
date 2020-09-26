part of '../complex.dart';

extension ComplexRootsX<T extends Complex> on T {
  /// Computes the n-th roots of this complex number.
  ///
  /// The nth roots are defined by the formula:
  ///
  ///     z_k = abs^(1/n) (cos(phi + 2πk/n) + i (sin(phi + 2πk/n))
  ///
  /// for `k=0, 1, ..., n-1`, where `abs` and `phi`
  /// are respectively the modulus and argument of this complex number.
  ///
  /// If one or both parts of this complex number is NaN, a list with just
  /// one element, [nan] is returned.
  /// If neither part is `NaN`, but at least one part is infinite, the result
  /// is a one-element list containing [infinity].
  List<Cartesian> nthRoot(int n) {
    if (n <= 0) {
      throw ArgumentError("Can't compute nth root for negative n");
    }

    final result = List<Cartesian>();

    if (isNaN) {
      result.add(Complex.nan);
      return result;
    }
    if (isInfinite) {
      result.add(Complex.infinity);
      return result;
    }

    // nth root of abs -- faster / more accurate to use a solver here?
    final double nthRootOfAbs = math.pow(abs(), 1.0 / n);

    // Compute nth roots of complex number with k = 0, 1, ... n-1
    final double nthPhi = argument() / n;
    final double slice = 2 * math.pi / n;
    double innerPart = nthPhi;
    for (int k = 0; k < n; k++) {
      // inner part
      final double realPart = nthRootOfAbs * math.cos(innerPart);
      final double imaginaryPart = nthRootOfAbs * math.sin(innerPart);
      result.add(Cartesian(realPart, imaginaryPart));
      innerPart += slice;
    }

    return result;
  }
}
