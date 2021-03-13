part of '../complex.dart';

/// # ComplexRootsX
///
/// {@template ComplexRootsX_nthRoot}
/// ## ComplexRootsX
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
/// {@endtemplate}
extension ComplexRootsX<T extends Complex> on T {
  /// {@macro ComplexRootsX_nthRoot}
  List<Complex> nthRoot(int n) {
    if (n <= 0) {
      throw ArgumentError("Can't compute nth root for negative n");
    }

    if (isNaN) {
      return [Complex.nan];
    } else if (isInfinite) {
      return [Complex.infinity];
    }

    // nth root of abs -- faster / more accurate to use a solver here?
    final nthRootOfAbs = math.pow(abs(), 1.0 / n).toDouble();

    // Compute nth roots of complex number with k = 0, 1, ... n-1
    final nthPhi = argument() / n;
    final slice = 2 * math.pi / n;
    var innerPart = nthPhi;
    return List.generate(n, (_) {
      final c = Complex.polar(nthRootOfAbs, innerPart);
      innerPart += slice;
      return c;
    });
  }
}
