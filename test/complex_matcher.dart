import 'package:test/test.dart';
import 'package:complex/complex.dart';

/// Class to test extending Complex
class TestComplex extends Cartesian {
  TestComplex(double real, double imaginary) : super(real, imaginary);

  factory TestComplex.from(Complex other) {
    return TestComplex(other.real, other.imaginary);
  }

  @override
  String toString() {
    return "$real ${imaginary}j";
  }
}

/// Returns a matcher which matches if the match argument is within [delta]
/// of some [value]; i.e. if the match argument is greater than
/// than or equal [value]-[delta] and less than or equal to [value]+[delta].
Matcher closeToZ(Complex value, num delta) => _IsCloseToZ(value, delta);

class _IsCloseToZ extends Matcher {
  final Complex _value;
  final num _delta;

  const _IsCloseToZ(this._value, this._delta);

  bool matches(item, Map matchState) {
    if (item is! Complex) {
      return false;
    }
    var re_diff = item.real - _value.real;
    if (re_diff < 0) re_diff = -re_diff;
    if (re_diff > _delta) {
      return false;
    }
    var im_diff = item.imaginary - _value.imaginary;
    if (im_diff < 0) im_diff = -im_diff;
    if (im_diff > _delta) {
      return false;
    }
    return true;
  }

  Description describe(Description description) => description
      .add('a complex value within ')
      .addDescriptionOf(_delta)
      .add(' of ')
      .addDescriptionOf(_value);

  Description describeMismatch(
      item, Description mismatchDescription, Map matchState, bool verbose) {
    if (item is! Complex) {
      return mismatchDescription.add(' not complex');
    } else {
      var diff = item.abs() - _value.abs();
      if (diff < 0) diff = -diff;
      return mismatchDescription.add(' differs by ').addDescriptionOf(diff);
    }
  }
}
