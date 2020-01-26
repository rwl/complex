// Licensed to the Apache Software Foundation (ASF) under one or more
// contributor license agreements.  See the NOTICE file distributed with
// this work for additional information regarding copyright ownership.
// The ASF licenses this file to You under the Apache License, Version 2.0
// (the "License"); you may not use this file except in compliance with
// the License.  You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

import 'dart:math' as math;
import 'package:test/test.dart';
import 'package:complex/complex.dart';

main() {
  double inf = double.infinity;
  double neginf = double.negativeInfinity;
  double nan = double.nan;
  double pi = math.pi;
  Complex oneInf = Complex(1, inf);
  Complex oneNegInf = Complex(1, neginf);
  Complex infOne = Complex(inf, 1);
  Complex infZero = Complex(inf, 0);
  Complex infNaN = Complex(inf, nan);
  Complex infNegInf = Complex(inf, neginf);
  Complex infInf = Complex(inf, inf);
  Complex negInfInf = Complex(neginf, inf);
  Complex negInfZero = Complex(neginf, 0);
  Complex negInfOne = Complex(neginf, 1);
  Complex negInfNaN = Complex(neginf, nan);
  Complex negInfNegInf = Complex(neginf, neginf);
  Complex oneNaN = Complex(1, nan);
  Complex zeroInf = Complex(0, inf);
  Complex zeroNaN = Complex(0, nan);
  Complex nanInf = Complex(nan, inf);
  Complex nanNegInf = Complex(nan, neginf);
  Complex nanZero = Complex(nan, 0);

  test('Constructor', () {
    Complex z = Complex(3.0, 4.0);
    expect(3.0, closeTo(z.real, 1.0e-5));
    expect(4.0, closeTo(z.imaginary, 1.0e-5));
  });

  test('ConstructorNaN', () {
    Complex z = Complex(3.0, double.nan);
    expect(z.isNaN, isTrue);

    z = Complex(nan, 4.0);
    expect(z.isNaN, isTrue);

    z = Complex(3.0, 4.0);
    expect(z.isNaN, isFalse);
  });

  test('Abs', () {
    Complex z = Complex(3.0, 4.0);
    expect(5.0, closeTo(z.abs(), 1.0e-5));
  });

  test('AbsNaN', () {
    expect(Complex.nan.abs().isNaN, isTrue);
    Complex z = Complex(inf, nan);
    expect(z.abs().isNaN, isTrue);
  });

  test('AbsInfinite', () {
    Complex z = Complex(inf, 0);
    expect(inf, equals(z.abs()));
    z = Complex(0, neginf);
    expect(inf, equals(z.abs()));
    z = Complex(inf, neginf);
    expect(inf, equals(z.abs()));
  });

  test('Add', () {
    Complex x = Complex(3.0, 4.0);
    Complex y = Complex(5.0, 6.0);
    Complex z = x + y;
    expect(8.0, closeTo(z.real, 1.0e-5));
    expect(10.0, closeTo(z.imaginary, 1.0e-5));
  });

  test('AddNaN', () {
    Complex x = Complex(3.0, 4.0);
    Complex z = x + Complex.nan;
    expect(Complex.nan, equals(z));
    z = Complex(1, nan);
    Complex w = x + z;
    expect(Complex.nan, equals(w));
  });

  test('AddInf', () {
    Complex x = Complex(1, 1);
    Complex z = Complex(inf, 0);
    Complex w = x + z;
    expect(w.imaginary, equals(1));
    expect(inf, equals(w.real));

    x = Complex(neginf, 0);
    expect((x + z).real.isNaN, isTrue);
  });

  test('ScalarAdd', () {
    Complex x = Complex(3.0, 4.0);
    double yDouble = 2.0;
    Complex yComplex = Complex(yDouble);
    expect(x + yComplex, equals(x + yDouble));
  });

  test('ScalarAddNaN', () {
    Complex x = Complex(3.0, 4.0);
    double yDouble = double.nan;
    Complex yComplex = Complex(yDouble);
    expect(x + yComplex, equals(x + yDouble));
  });

  test('ScalarAddInf', () {
    Complex x = Complex(1, 1);
    double yDouble = double.infinity;

    Complex yComplex = Complex(yDouble);
    expect(x + yComplex, x + yDouble);

    x = Complex(neginf, 0);
    expect(x + yComplex, x + yDouble);
  });

  test('Conjugate', () {
    Complex x = Complex(3.0, 4.0);
    Complex z = x.conjugate();
    expect(3.0, closeTo(z.real, 1.0e-5));
    expect(-4.0, closeTo(z.imaginary, 1.0e-5));
  });

  test('ConjugateNaN', () {
    Complex z = Complex.nan.conjugate();
    expect(z.isNaN, isTrue);
  });

  test('ConjugateInfinite', () {
    Complex z = Complex(0, inf);
    expect(neginf, equals(z.conjugate().imaginary));
    z = Complex(0, neginf);
    expect(inf, equals(z.conjugate().imaginary));
  });

  test('Divide', () {
    Complex x = Complex(3.0, 4.0);
    Complex y = Complex(5.0, 6.0);
    Complex z = x / y;
    expect(39.0 / 61.0, closeTo(z.real, 1.0e-5));
    expect(2.0 / 61.0, closeTo(z.imaginary, 1.0e-5));
  });

  test('DivideReal', () {
    Complex x = Complex(2, 3);
    Complex y = Complex(2, 0);
    expect(Complex(1, 1.5), equals(x / y));
  });

  test('DivideImaginary', () {
    Complex x = Complex(2, 3);
    Complex y = Complex(0, 2);
    expect(Complex(1.5, -1), equals(x / y));
  });

  test('DivideInf', () {
    Complex x = Complex(3, 4);
    Complex w = Complex(neginf, inf);
    expect(x / w == Complex.zero, isTrue);

    Complex z = w / x;
    expect(z.real.isNaN, isTrue);
    expect(inf, equals(z.imaginary));

    w = Complex(inf, inf);
    z = w / x;
    expect(z.imaginary.isNaN, isTrue);
    expect(inf, equals(z.real));

    w = Complex(1, inf);
    z = w / w;
    expect(z.real.isNaN, isTrue);
    expect(z.imaginary.isNaN, isTrue);
  });

  test('DivideZero', () {
    Complex x = Complex(3.0, 4.0);
    Complex z = x / Complex.zero;
    // expect(z, Complex.INF); // See MATH-657
    expect(z, equals(Complex.nan));
  });

  test('DivideZeroZero', () {
    Complex x = Complex(0.0, 0.0);
    Complex z = x / Complex.zero;
    expect(z, equals(Complex.nan));
  });

  test('DivideNaN', () {
    Complex x = Complex(3.0, 4.0);
    Complex z = x / Complex.nan;
    expect(z.isNaN, isTrue);
  });

  test('DivideNaNInf', () {
    Complex z = oneInf / Complex.one;
    expect(z.real.isNaN, isTrue);
    expect(inf, equals(z.imaginary));

    z = negInfNegInf / oneNaN;
    expect(z.real.isNaN, isTrue);
    expect(z.imaginary.isNaN, isTrue);

    z = negInfInf / Complex.one;
    expect(z.real.isNaN, isTrue);
    expect(z.imaginary.isNaN, isTrue);
  });

  test('ScalarDivide', () {
    Complex x = Complex(3.0, 4.0);
    double yDouble = 2.0;
    Complex yComplex = Complex(yDouble);
    expect(x / yComplex, equals(x / yDouble));
  });

  test('ScalarDivideNaN', () {
    Complex x = Complex(3.0, 4.0);
    double yDouble = double.nan;
    Complex yComplex = Complex(yDouble);
    expect(x / yComplex, equals(x / yDouble));
  });

  test('ScalarDivideInf', () {
    Complex x = Complex(1, 1);
    double yDouble = double.infinity;
    Complex yComplex = Complex(yDouble);
    expect(x / yComplex, equals(x / yDouble));

    yDouble = double.negativeInfinity;
    yComplex = Complex(yDouble);
    expect(x / yComplex, equals(x / yDouble));

    x = Complex(1, double.negativeInfinity);
    expect(x / yComplex, equals(x / yDouble));
  });

  test('ScalarDivideZero', () {
    Complex x = Complex(1, 1);
    expect(x / Complex.zero, equals(x / 0));
  });

  test('Reciprocal', () {
    Complex z = Complex(5.0, 6.0);
    Complex act = z.reciprocal();
    double expRe = 5.0 / 61.0;
    double expIm = -6.0 / 61.0;
    expect(expRe, closeTo(act.real, /*FastMath.ulp(expRe)*/ 1.0e-12));
    expect(expIm, closeTo(act.imaginary, /*FastMath.ulp(expIm)*/ 1.0e-12));
  });

  test('ReciprocalReal', () {
    Complex z = Complex(-2.0, 0.0);
    //expect(Complex.equals(Complex(-0.5, 0.0), z.reciprocal()), isTrue);
    expect(Complex(-0.5, 0.0), equals(z.reciprocal()));
  });

  test('ReciprocalImaginary', () {
    Complex z = Complex(0.0, -2.0);
    expect(Complex(0.0, 0.5), equals(z.reciprocal()));
  });

  test('ReciprocalInf', () {
    Complex z = Complex(neginf, inf);
    expect(z.reciprocal() == Complex.zero, isTrue);

    z = Complex(1, inf).reciprocal();
    expect(z, equals(Complex.zero));
  });

  test('ReciprocalZero', () {
    expect(Complex.zero.reciprocal(), equals(Complex.infinity));
  });

  test('ReciprocalNaN', () {
    expect(Complex.nan.reciprocal().isNaN, isTrue);
  });

  test('Multiply', () {
    Complex x = Complex(3.0, 4.0);
    Complex y = Complex(5.0, 6.0);
    Complex z = x * y;
    expect(-9.0, closeTo(z.real, 1.0e-5));
    expect(38.0, closeTo(z.imaginary, 1.0e-5));
  });

  test('MultiplyNaN', () {
    Complex x = Complex(3.0, 4.0);
    Complex z = x * Complex.nan;
    expect(Complex.nan, equals(z));
    z = Complex.nan * 5;
    expect(Complex.nan, equals(z));
  });

  test('MultiplyInfInf', () {
    // expect(infInf.multiply(infInf).isNaN()); // MATH-620
    expect((infInf * infInf).isInfinite, isTrue);
  });

  test('MultiplyNaNInf', () {
    Complex z = Complex(1, 1);
    Complex w = z * infOne;
    expect(w.real, equals(inf));
    expect(w.imaginary, equals(inf));

    // [MATH-164]
    expect(Complex(1, 0) * infInf == Complex.infinity, isTrue);
    expect(Complex(-1, 0) * infInf == Complex.infinity, isTrue);
    expect(Complex(1, 0) * negInfZero == Complex.infinity, isTrue);

    w = oneInf * oneNegInf;
    expect(w.real, equals(inf));
    expect(w.imaginary, equals(inf));

    w = negInfNegInf * oneNaN;
    expect(w.real.isNaN, isTrue);
    expect(w.imaginary.isNaN, isTrue);

    z = Complex(1, neginf);
    expect(Complex.infinity, equals(z * z));
  });

  test('ScalarMultiply', () {
    Complex x = Complex(3.0, 4.0);
    double yDouble = 2.0;
    Complex yComplex = Complex(yDouble);
    expect(x * yComplex, equals(x * yDouble));
    int zInt = -5;
    Complex zComplex = Complex(zInt);
    expect(x * zComplex, equals(x * zInt));
  });

  test('ScalarMultiplyNaN', () {
    Complex x = Complex(3.0, 4.0);
    double yDouble = double.nan;
    Complex yComplex = Complex(yDouble);
    expect(x * yComplex, equals(x * yDouble));
  });

  test('ScalarMultiplyInf', () {
    Complex x = Complex(1, 1);
    double yDouble = double.infinity;
    Complex yComplex = Complex(yDouble);
    expect(x * yComplex, equals(x * yDouble));

    yDouble = double.negativeInfinity;
    yComplex = Complex(yDouble);
    expect(x * yComplex, equals(x * yDouble));
  });

  test('Negate', () {
    Complex x = Complex(3.0, 4.0);
    Complex z = -x;
    expect(-3.0, closeTo(z.real, 1.0e-5));
    expect(-4.0, closeTo(z.imaginary, 1.0e-5));
  });

  test('NegateNaN', () {
    Complex z = -Complex.nan;
    expect(z.isNaN, isTrue);
  });

  test('Subtract', () {
    Complex x = Complex(3.0, 4.0);
    Complex y = Complex(5.0, 6.0);
    Complex z = x - y;
    expect(-2.0, closeTo(z.real, 1.0e-5));
    expect(-2.0, closeTo(z.imaginary, 1.0e-5));
  });

  test('SubtractNaN', () {
    Complex x = Complex(3.0, 4.0);
    Complex z = x - Complex.nan;
    expect(Complex.nan, equals(z));
    z = Complex(1, nan);
    Complex w = x - z;
    expect(Complex.nan, equals(w));
  });

  test('SubtractInf', () {
    Complex x = Complex(1, 1);
    Complex z = Complex(neginf, 0);
    Complex w = x - z;
    expect(w.imaginary, equals(1));
    expect(inf, equals(w.real));

    x = Complex(neginf, 0);
    expect((x - z).real.isNaN, isTrue);
  });

  test('ScalarSubtract', () {
    Complex x = Complex(3.0, 4.0);
    double yDouble = 2.0;
    Complex yComplex = Complex(yDouble);
    expect(x - yComplex, equals(x - yDouble));
  });

  test('ScalarSubtractNaN', () {
    Complex x = Complex(3.0, 4.0);
    double yDouble = double.nan;
    Complex yComplex = Complex(yDouble);
    expect(x - yComplex, equals(x - yDouble));
  });

  test('ScalarSubtractInf', () {
    Complex x = Complex(1, 1);
    double yDouble = double.infinity;
    Complex yComplex = Complex(yDouble);
    expect(x - yComplex, equals(x - yDouble));

    x = Complex(neginf, 0);
    expect(x - yComplex, equals(x - yDouble));
  });

  test('EqualsNull', () {
    Complex x = Complex(3.0, 4.0);
    expect(x == null, isFalse);
  });

  /*test('FloatingPointEqualsPrecondition1', () {
    expect(() {
      Complex.equals(Complex(3.0, 4.0), null, 3);
    }, throwsA(NullPointerException));
  });

  test('FloatingPointEqualsPrecondition2', () {
    expect(() {
      Complex.equals(null, Complex(3.0, 4.0), 3);
    }, throwsA(NullPointerException));
  });*/

  test('EqualsClass', () {
    Complex x = Complex(3.0, 4.0);
    expect(x == TestComplex, isFalse);
  });

  test('EqualsSame', () {
    Complex x = Complex(3.0, 4.0);
    expect(x == x, isTrue);
  });

  test('FloatingPointEquals', () {
    double re = -3.21;
    double im = 456789e10;

    final Complex x = Complex(re, im);
    Complex y = Complex(re, im);

    expect(x == y, isTrue);
    /*expect(Complex.equals(x, y), isTrue);

    final int maxUlps = 5;
    for (int i = 0; i < maxUlps; i++) {
      re = FastMath.nextUp(re);
      im = FastMath.nextUp(im);
    }
    y = Complex(re, im);
    expect(Complex.equals(x, y, maxUlps), isTrue);

    re = FastMath.nextUp(re);
    im = FastMath.nextUp(im);
    y = Complex(re, im);
    expect(Complex.equals(x, y, maxUlps), isFalse);*/
  });

  /*test('FloatingPointEqualsNaN', () {
    Complex c = Complex(double.nan, 1);
    expect(Complex.equals(c, c), isFalse);

    c = Complex(1, double.nan);
    expect(Complex.equals(c, c), isFalse);
  });

  test('FloatingPointEqualsWithAllowedDelta', () {
    final double re = 153.0000;
    final double im = 152.9375;
    final double tol1 = 0.0625;
    final Complex x = Complex(re, im);
    final Complex y = Complex(re + tol1, im + tol1);
    expect(Complex.equals(x, y, tol1), isTrue);

    final double tol2 = 0.0624;
    expect(Complex.equals(x, y, tol2), isFalse);
  });

  test('FloatingPointEqualsWithRelativeTolerance', () {
    final double tol = 1e-4;
    final double re = 1;
    final double im = 1e10;

    final double f = 1 + tol;
    final Complex x = Complex(re, im);
    final Complex y = Complex(re * f, im * f);
    expect(Complex.equalsWithRelativeTolerance(x, y, tol), isTrue);
  });*/

  test('EqualsTrue', () {
    Complex x = Complex(3.0, 4.0);
    Complex y = Complex(3.0, 4.0);
    expect(x == y, isTrue);
  });

  test('EqualsRealDifference', () {
    Complex x = Complex(0.0, 0.0);
    Complex y = Complex(0.0 + double.minPositive, 0.0);
    expect(x == y, isFalse);
  });

  test('EqualsImaginaryDifference', () {
    Complex x = Complex(0.0, 0.0);
    Complex y = Complex(0.0, 0.0 + double.minPositive);
    expect(x == y, isFalse);
  });

  test('EqualsNaN', () {
    Complex realNaN = Complex(double.nan, 0.0);
    Complex imaginaryNaN = Complex(0.0, double.nan);
    Complex complexNaN = Complex.nan;
    expect(realNaN == imaginaryNaN, isTrue);
    expect(imaginaryNaN == complexNaN, isTrue);
    expect(realNaN == complexNaN, isTrue);
  });

  test('HashCode', () {
    Complex x = Complex(0.0, 0.0);
    Complex y = Complex(0.0, 0.0 + double.minPositive);
    // TODO expect(x.hashCode == y.hashCode, isFalse);
    y = Complex(0.0 + double.minPositive, 0.0);
    // TODO expect(x.hashCode == y.hashCode, isFalse);
    Complex realNaN = Complex(double.nan, 0.0);
    Complex imaginaryNaN = Complex(0.0, double.nan);
    expect(realNaN.hashCode, equals(imaginaryNaN.hashCode));
    expect(imaginaryNaN.hashCode, equals(Complex.nan.hashCode));

    // MATH-1118
    // "equals" and "hashCode" must be compatible: if two objects have
    // different hash codes, "equals" must return false.
    final String msg = "'==' not compatible with 'hashCode'";

    x = Complex(0.0, 0.0);
    y = Complex(0.0, -0.0);
    // TODO expect(x.hashCode != y.hashCode, isTrue);
    // TODO expect(x == y, isFalse, reason:msg);

    x = Complex(0.0, 0.0);
    y = Complex(-0.0, 0.0);
    // TODO expect(x.hashCode != y.hashCode, isTrue);
    // TODO expect(x == y, isFalse, reason:msg);
  });

  test('Acos', () {
    Complex z = Complex(3, 4);
    Complex expected = Complex(0.936812, -2.30551);
    expect(expected, closeToZ(z.acos(), 1.0e-5));
    expect(Complex(math.acos(0), 0), closeToZ(Complex.zero.acos(), 1.0e-12));
  });

  test('AcosInf', () {
    expect(Complex.nan, equals(oneInf.acos()));
    expect(Complex.nan, equals(oneNegInf.acos()));
    expect(Complex.nan, equals(infOne.acos()));
    expect(Complex.nan, equals(negInfOne.acos()));
    expect(Complex.nan, equals(infInf.acos()));
    expect(Complex.nan, equals(infNegInf.acos()));
    expect(Complex.nan, equals(negInfInf.acos()));
    expect(Complex.nan, equals(negInfNegInf.acos()));
  });

  test('AcosNaN', () {
    expect(Complex.nan.acos().isNaN, isTrue);
  });

  test('Asin', () {
    Complex z = Complex(3, 4);
    Complex expected = Complex(0.633984, 2.30551);
    expect(expected, closeToZ(z.asin(), 1.0e-5));
  });

  test('AsinNaN', () {
    expect(Complex.nan.asin().isNaN, isTrue);
  });

  test('AsinInf', () {
    expect(Complex.nan, equals(oneInf.asin()));
    expect(Complex.nan, equals(oneNegInf.asin()));
    expect(Complex.nan, equals(infOne.asin()));
    expect(Complex.nan, equals(negInfOne.asin()));
    expect(Complex.nan, equals(infInf.asin()));
    expect(Complex.nan, equals(infNegInf.asin()));
    expect(Complex.nan, equals(negInfInf.asin()));
    expect(Complex.nan, equals(negInfNegInf.asin()));
  });

  test('Atan', () {
    Complex z = Complex(3, 4);
    Complex expected = Complex(1.44831, 0.158997);
    expect(expected, closeToZ(z.atan(), 1.0e-5));
  });

  test('AtanInf', () {
    expect(Complex.nan, equals(oneInf.atan()));
    expect(Complex.nan, equals(oneNegInf.atan()));
    expect(Complex.nan, equals(infOne.atan()));
    expect(Complex.nan, equals(negInfOne.atan()));
    expect(Complex.nan, equals(infInf.atan()));
    expect(Complex.nan, equals(infNegInf.atan()));
    expect(Complex.nan, equals(negInfInf.atan()));
    expect(Complex.nan, equals(negInfNegInf.atan()));
  });

  test('AtanI', () {
    expect(Complex.i.atan().isNaN, isTrue);
  });

  test('AtanNaN', () {
    expect(Complex.nan.atan().isNaN, isTrue);
  });

  test('Cos', () {
    Complex z = Complex(3, 4);
    Complex expected = Complex(-27.03495, -3.851153);
    expect(expected, closeToZ(z.cos(), 1.0e-5));
  });

  test('CosNaN', () {
    expect(Complex.nan.cos().isNaN, isTrue);
  });

  test('CosInf', () {
    expect(infNegInf, equals(oneInf.cos()));
    expect(infInf, equals(oneNegInf.cos()));
    expect(Complex.nan, equals(infOne.cos()));
    expect(Complex.nan, equals(negInfOne.cos()));
    expect(Complex.nan, equals(infInf.cos()));
    expect(Complex.nan, equals(infNegInf.cos()));
    expect(Complex.nan, equals(negInfInf.cos()));
    expect(Complex.nan, equals(negInfNegInf.cos()));
  });

  test('Cosh', () {
    Complex z = Complex(3, 4);
    Complex expected = Complex(-6.58066, -7.58155);
    expect(expected, closeToZ(z.cosh(), 1.0e-5));
  });

  test('CoshNaN', () {
    expect(Complex.nan.cosh().isNaN, isTrue);
  });

  test('CoshInf', () {
    expect(Complex.nan, equals(oneInf.cosh()));
    expect(Complex.nan, equals(oneNegInf.cosh()));
    expect(infInf, equals(infOne.cosh()));
    expect(infNegInf, equals(negInfOne.cosh()));
    expect(Complex.nan, equals(infInf.cosh()));
    expect(Complex.nan, equals(infNegInf.cosh()));
    expect(Complex.nan, equals(negInfInf.cosh()));
    expect(Complex.nan, equals(negInfNegInf.cosh()));
  });

  test('Exp', () {
    Complex z = Complex(3, 4);
    Complex expected = Complex(-13.12878, -15.20078);
    expect(expected, closeToZ(z.exp(), 1.0e-5));
    expect(Complex.one, closeToZ(Complex.zero.exp(), 10e-12));
    Complex iPi = Complex.i * Complex(pi, 0);
    expect(-Complex.one, closeToZ(iPi.exp(), 10e-12));
  });

  test('ExpNaN', () {
    expect(Complex.nan.exp().isNaN, isTrue);
  });

  test('ExpInf', () {
    expect(Complex.nan, equals(oneInf.exp()));
    expect(Complex.nan, equals(oneNegInf.exp()));
    expect(infInf, equals(infOne.exp()));
    expect(Complex.zero, equals(negInfOne.exp()));
    expect(Complex.nan, equals(infInf.exp()));
    expect(Complex.nan, equals(infNegInf.exp()));
    expect(Complex.nan, equals(negInfInf.exp()));
    expect(Complex.nan, equals(negInfNegInf.exp()));
  });

  test('Log', () {
    Complex z = Complex(3, 4);
    Complex expected = Complex(1.60944, 0.927295);
    expect(expected, closeToZ(z.log(), 1.0e-5));
  });

  test('LogNaN', () {
    expect(Complex.nan.log().isNaN, isTrue);
  });

  test('LogInf', () {
    expect(Complex(inf, pi / 2), closeToZ(oneInf.log(), 10e-12));
    expect(Complex(inf, -pi / 2), closeToZ(oneNegInf.log(), 10e-12));
    expect(infZero, closeToZ(infOne.log(), 10e-12));
    expect(Complex(inf, pi), closeToZ(negInfOne.log(), 10e-12));
    expect(Complex(inf, pi / 4), closeToZ(infInf.log(), 10e-12));
    expect(Complex(inf, -pi / 4), closeToZ(infNegInf.log(), 10e-12));
    expect(Complex(inf, 3 * pi / 4), closeToZ(negInfInf.log(), 10e-12));
    expect(Complex(inf, -3 * pi / 4), closeToZ(negInfNegInf.log(), 10e-12));
  });

  test('LogZero', () {
    expect(negInfZero, equals(Complex.zero.log()));
  });

  test('Pow', () {
    Complex x = Complex(3, 4);
    Complex y = Complex(5, 6);
    Complex expected = Complex(-1.860893, 11.83677);
    expect(expected, closeToZ(x.power(y), 1.0e-5));
  });

  test('PowNaNBase', () {
    Complex x = Complex(3, 4);
    expect(Complex.nan.power(x).isNaN, isTrue);
  });

  test('PowNaNExponent', () {
    Complex x = Complex(3, 4);
    expect(x.power(Complex.nan).isNaN, isTrue);
  });

  test('PowInf', () {
    expect(Complex.nan, equals(Complex.one.power(oneInf)));
    expect(Complex.nan, equals(Complex.one.power(oneNegInf)));
    expect(Complex.nan, equals(Complex.one.power(infOne)));
    expect(Complex.nan, equals(Complex.one.power(infInf)));
    expect(Complex.nan, equals(Complex.one.power(infNegInf)));
    expect(Complex.nan, equals(Complex.one.power(negInfInf)));
    expect(Complex.nan, equals(Complex.one.power(negInfNegInf)));
    expect(Complex.nan, equals(infOne.power(Complex.one)));
    expect(Complex.nan, equals(negInfOne.power(Complex.one)));
    expect(Complex.nan, equals(infInf.power(Complex.one)));
    expect(Complex.nan, equals(infNegInf.power(Complex.one)));
    expect(Complex.nan, equals(negInfInf.power(Complex.one)));
    expect(Complex.nan, equals(negInfNegInf.power(Complex.one)));
    expect(Complex.nan, equals(negInfNegInf.power(infNegInf)));
    expect(Complex.nan, equals(negInfNegInf.power(negInfNegInf)));
    expect(Complex.nan, equals(negInfNegInf.power(infInf)));
    expect(Complex.nan, equals(infInf.power(infNegInf)));
    expect(Complex.nan, equals(infInf.power(negInfNegInf)));
    expect(Complex.nan, equals(infInf.power(infInf)));
    expect(Complex.nan, equals(infNegInf.power(infNegInf)));
    expect(Complex.nan, equals(infNegInf.power(negInfNegInf)));
    expect(Complex.nan, equals(infNegInf.power(infInf)));
  });

  test('PowZero', () {
    expect(Complex.nan, equals(Complex.zero.power(Complex.one)));
    expect(Complex.nan, equals(Complex.zero.power(Complex.zero)));
    expect(Complex.nan, equals(Complex.zero.power(Complex.i)));
    expect(Complex.one, closeToZ(Complex.one.power(Complex.zero), 10e-12));
    expect(Complex.one, closeToZ(Complex.i.power(Complex.zero), 10e-12));
    expect(Complex.one, closeToZ(Complex(-1, 3).power(Complex.zero), 10e-12));
  });

  test('ScalarPow', () {
    Complex x = Complex(3, 4);
    double yDouble = 5.0;
    Complex yComplex = Complex(yDouble);
    expect(x.power(yComplex), equals(x.pow(yDouble)));
  });

  test('ScalarPowNaNBase', () {
    Complex x = Complex.nan;
    double yDouble = 5.0;
    Complex yComplex = Complex(yDouble);
    expect(x.power(yComplex), equals(x.pow(yDouble)));
  });

  test('ScalarPowNaNExponent', () {
    Complex x = Complex(3, 4);
    double yDouble = double.nan;
    Complex yComplex = Complex(yDouble);
    expect(x.power(yComplex), x.pow(yDouble));
  });

  test('ScalarPowInf', () {
    expect(Complex.nan, equals(Complex.one.pow(double.infinity)));
    expect(Complex.nan, equals(Complex.one.pow(double.negativeInfinity)));
    expect(Complex.nan, equals(infOne.pow(1.0)));
    expect(Complex.nan, equals(negInfOne.pow(1.0)));
    expect(Complex.nan, equals(infInf.pow(1.0)));
    expect(Complex.nan, equals(infNegInf.pow(1.0)));
    expect(Complex.nan, equals(negInfInf.pow(10)));
    expect(Complex.nan, equals(negInfNegInf.pow(1.0)));
    expect(Complex.nan, equals(negInfNegInf.pow(double.infinity)));
    expect(Complex.nan, equals(negInfNegInf.pow(double.infinity)));
    expect(Complex.nan, equals(infInf.pow(double.infinity)));
    expect(Complex.nan, equals(infInf.pow(double.negativeInfinity)));
    expect(Complex.nan, equals(infNegInf.pow(double.negativeInfinity)));
    expect(Complex.nan, equals(infNegInf.pow(double.infinity)));
  });

  test('ScalarPowZero', () {
    expect(Complex.nan, equals(Complex.zero.pow(1.0)));
    expect(Complex.nan, equals(Complex.zero.pow(0.0)));
    expect(Complex.one, closeToZ(Complex.one.pow(0.0), 10e-12));
    expect(Complex.one, closeToZ(Complex.i.pow(0.0), 10e-12));
    expect(Complex.one, closeToZ(Complex(-1, 3).pow(0.0), 10e-12));
  });

  test('powNull', () {
    expect(() => Complex.one.pow(null), throwsArgumentError);
  });

  test('Sin', () {
    Complex z = Complex(3, 4);
    Complex expected = Complex(3.853738, -27.01681);
    expect(expected, closeToZ(z.sin(), 1.0e-5));
  });

  test('SinInf', () {
    expect(infInf, equals(oneInf.sin()));
    expect(infNegInf, equals(oneNegInf.sin()));
    expect(Complex.nan, equals(infOne.sin()));
    expect(Complex.nan, equals(negInfOne.sin()));
    expect(Complex.nan, equals(infInf.sin()));
    expect(Complex.nan, equals(infNegInf.sin()));
    expect(Complex.nan, equals(negInfInf.sin()));
    expect(Complex.nan, equals(negInfNegInf.sin()));
  });

  test('SinNaN', () {
    expect(Complex.nan.sin().isNaN, isTrue);
  });

  test('Sinh', () {
    Complex z = Complex(3, 4);
    Complex expected = Complex(-6.54812, -7.61923);
    expect(expected, closeToZ(z.sinh(), 1.0e-5));
  });

  test('SinhNaN', () {
    expect(Complex.nan.sinh().isNaN, isTrue);
  });

  test('SinhInf', () {
    expect(Complex.nan, equals(oneInf.sinh()));
    expect(Complex.nan, equals(oneNegInf.sinh()));
    expect(infInf, equals(infOne.sinh()));
    expect(negInfInf, equals(negInfOne.sinh()));
    expect(Complex.nan, equals(infInf.sinh()));
    expect(Complex.nan, equals(infNegInf.sinh()));
    expect(Complex.nan, equals(negInfInf.sinh()));
    expect(Complex.nan, equals(negInfNegInf.sinh()));
  });

  test('SqrtRealPositive', () {
    Complex z = Complex(3, 4);
    Complex expected = Complex(2, 1);
    expect(expected, closeToZ(z.sqrt(), 1.0e-5));
  });

  test('SqrtRealZero', () {
    Complex z = Complex(0.0, 4);
    Complex expected = Complex(1.41421, 1.41421);
    expect(expected, closeToZ(z.sqrt(), 1.0e-5));
  });

  test('SqrtRealNegative', () {
    Complex z = Complex(-3.0, 4);
    Complex expected = Complex(1, 2);
    expect(expected, closeToZ(z.sqrt(), 1.0e-5));
  });

  test('SqrtImaginaryZero', () {
    Complex z = Complex(-3.0, 0.0);
    Complex expected = Complex(0.0, 1.73205);
  });

  test('SqrtImaginaryNegative', () {
    Complex z = Complex(-3.0, -4.0);
    Complex expected = Complex(1.0, -2.0);
    expect(expected, closeToZ(z.sqrt(), 1.0e-5));
  });

  test('SqrtPolar', () {
    double r = 1.0;
    for (int i = 0; i < 5; i++) {
      r += i;
      double theta = 0.0;
      for (int j = 0; j < 11; j++) {
        theta += pi / 12;
        Complex z = Complex.polar(r, theta);
        Complex sqrtz = Complex.polar(math.sqrt(r), theta / 2);
        expect(sqrtz, closeToZ(z.sqrt(), 10e-12));
      }
    }
  });

  test('SqrtNaN', () {
    expect(Complex.nan.sqrt().isNaN, isTrue);
  });

  test('SqrtInf', () {
    expect(infNaN, equals(oneInf.sqrt()));
    expect(infNaN, equals(oneNegInf.sqrt()));
    expect(infZero, equals(infOne.sqrt()));
    expect(zeroInf, equals(negInfOne.sqrt()));
    expect(infNaN, equals(infInf.sqrt()));
    expect(infNaN, equals(infNegInf.sqrt()));
    expect(nanInf, equals(negInfInf.sqrt()));
    expect(nanNegInf, equals(negInfNegInf.sqrt()));
  });

  test('Sqrt1z', () {
    Complex z = Complex(3, 4);
    Complex expected = Complex(4.08033, -2.94094);
    expect(expected, closeToZ(z.sqrt1z(), 1.0e-5));
  });

  test('Sqrt1zNaN', () {
    expect(Complex.nan.sqrt1z().isNaN, isTrue);
  });

  test('Tan', () {
    Complex z = Complex(3, 4);
    Complex expected = Complex(-0.000187346, 0.999356);
    expect(expected, closeToZ(z.tan(), 1.0e-5));
    /* Check that no overflow occurs (MATH-722) */
    Complex actual = Complex(3.0, 1E10).tan();
    expected = Complex(0, 1);
    expect(expected, closeToZ(actual, 1.0e-5));
    actual = Complex(3.0, -1E10).tan();
    expected = Complex(0, -1);
    expect(expected, closeToZ(actual, 1.0e-5));
  });

  test('TanNaN', () {
    expect(Complex.nan.tan().isNaN, isTrue);
  });

  test('TanInf', () {
    expect(Complex(0.0, 1.0), equals(oneInf.tan()));
    expect(Complex(0.0, -1.0), equals(oneNegInf.tan()));
    expect(Complex.nan, equals(infOne.tan()));
    expect(Complex.nan, equals(negInfOne.tan()));
    expect(Complex.nan, equals(infInf.tan()));
    expect(Complex.nan, equals(infNegInf.tan()));
    expect(Complex.nan, equals(negInfInf.tan()));
    expect(Complex.nan, equals(negInfNegInf.tan()));
  });

  test('TanCritical', () {
    expect(infNaN, equals(Complex(pi / 2, 0).tan()));
    expect(negInfNaN, equals(Complex(-pi / 2, 0).tan()));
  });

  test('Tanh', () {
    Complex z = Complex(3, 4);
    Complex expected = Complex(1.00071, 0.00490826);
    expect(expected, closeToZ(z.tanh(), 1.0e-5));
    /* Check that no overflow occurs (MATH-722) */
    Complex actual = Complex(1E10, 3.0).tanh();
    expected = Complex(1, 0);
    expect(expected, closeToZ(actual, 1.0e-5));
    actual = Complex(-1E10, 3.0).tanh();
    expected = Complex(-1, 0);
    expect(expected, closeToZ(actual, 1.0e-5));
  });

  test('TanhNaN', () {
    expect(Complex.nan.tanh().isNaN, isTrue);
  });

  test('TanhInf', () {
    expect(Complex.nan, equals(oneInf.tanh()));
    expect(Complex.nan, equals(oneNegInf.tanh()));
    expect(Complex(1.0, 0.0), equals(infOne.tanh()));
    expect(Complex(-1.0, 0.0), equals(negInfOne.tanh()));
    expect(Complex.nan, equals(infInf.tanh()));
    expect(Complex.nan, equals(infNegInf.tanh()));
    expect(Complex.nan, equals(negInfInf.tanh()));
    expect(Complex.nan, equals(negInfNegInf.tanh()));
  });

  test('TanhCritical', () {
    expect(nanInf, equals(Complex(0, pi / 2).tanh()));
  });

  /// test issue MATH-221

  test('Math221', () {
    //expect(Complex.equals(Complex(0, -1), Complex(0, 1) * Complex(-1, 0)), isTrue);
    expect(Complex(0, -1) == Complex(0, 1) * Complex(-1, 0), isTrue);
  });

  /// Test: computing <b>third roots</b> of z.
  ///
  ///     z = -2 + 2 * i
  ///      => z_0 =  1      +          i
  ///      => z_1 = -1.3660 + 0.3660 * i
  ///      => z_2 =  0.3660 - 1.3660 * i

  test('NthRoot_normal_thirdRoot', () {
    // The complex number we want to compute all third-roots for.
    Complex z = Complex(-2, 2);
    // The List holding all third roots
    List<Complex> thirdRootsOfZ = z.nthRoot(3); //.toArray(Complex[0]);
    // Returned Collection must not be empty!
    expect(3, equals(thirdRootsOfZ.length));
    // test z_0
    expect(1.0, closeTo(thirdRootsOfZ[0].real, 1.0e-5));
    expect(1.0, closeTo(thirdRootsOfZ[0].imaginary, 1.0e-5));
    // test z_1
    expect(-1.3660254037844386, closeTo(thirdRootsOfZ[1].real, 1.0e-5));
    expect(0.36602540378443843, closeTo(thirdRootsOfZ[1].imaginary, 1.0e-5));
    // test z_2
    expect(0.366025403784439, closeTo(thirdRootsOfZ[2].real, 1.0e-5));
    expect(-1.3660254037844384, closeTo(thirdRootsOfZ[2].imaginary, 1.0e-5));
  });

  /// Test: computing <b>fourth roots</b> of z.
  ///
  ///     z = 5 - 2 * i
  ///      => z_0 =  1.5164 - 0.1446 * i
  ///      => z_1 =  0.1446 + 1.5164 * i
  ///      => z_2 = -1.5164 + 0.1446 * i
  ///      => z_3 = -1.5164 - 0.1446 * i

  test('NthRoot_normal_fourthRoot', () {
    // The complex number we want to compute all third-roots for.
    Complex z = Complex(5, -2);
    // The List holding all fourth roots
    List<Complex> fourthRootsOfZ = z.nthRoot(4); //.toArray(Complex[0]);
    // Returned Collection must not be empty!
    expect(4, equals(fourthRootsOfZ.length));
    // test z_0
    expect(1.5164629308487783, closeTo(fourthRootsOfZ[0].real, 1.0e-5));
    expect(-0.14469266210702247, closeTo(fourthRootsOfZ[0].imaginary, 1.0e-5));
    // test z_1
    expect(0.14469266210702256, closeTo(fourthRootsOfZ[1].real, 1.0e-5));
    expect(1.5164629308487783, closeTo(fourthRootsOfZ[1].imaginary, 1.0e-5));
    // test z_2
    expect(-1.5164629308487783, closeTo(fourthRootsOfZ[2].real, 1.0e-5));
    expect(0.14469266210702267, closeTo(fourthRootsOfZ[2].imaginary, 1.0e-5));
    // test z_3
    expect(-0.14469266210702275, closeTo(fourthRootsOfZ[3].real, 1.0e-5));
    expect(-1.5164629308487783, closeTo(fourthRootsOfZ[3].imaginary, 1.0e-5));
  });

  /// Test: computing <b>third roots</b> of z.
  ///
  ///     z = 8
  ///      => z_0 =  2
  ///      => z_1 = -1 + 1.73205 * i
  ///      => z_2 = -1 - 1.73205 * i

  test('NthRoot_cornercase_thirdRoot_imaginaryPartEmpty', () {
    // The number 8 has three third roots. One we all already know is the number 2.
    // But there are two more complex roots.
    Complex z = Complex(8, 0);
    // The List holding all third roots
    List<Complex> thirdRootsOfZ = z.nthRoot(3); //.toArray(Complex[0]);
    // Returned Collection must not be empty!
    expect(3, equals(thirdRootsOfZ.length));
    // test z_0
    expect(2.0, closeTo(thirdRootsOfZ[0].real, 1.0e-5));
    expect(0.0, closeTo(thirdRootsOfZ[0].imaginary, 1.0e-5));
    // test z_1
    expect(-1.0, closeTo(thirdRootsOfZ[1].real, 1.0e-5));
    expect(1.7320508075688774, closeTo(thirdRootsOfZ[1].imaginary, 1.0e-5));
    // test z_2
    expect(-1.0, closeTo(thirdRootsOfZ[2].real, 1.0e-5));
    expect(-1.732050807568877, closeTo(thirdRootsOfZ[2].imaginary, 1.0e-5));
  });

  /// Test: computing <b>third roots</b> of z with real part 0.
  ///
  ///     z = 2 * i
  ///      => z_0 =  1.0911 + 0.6299 * i
  ///      => z_1 = -1.0911 + 0.6299 * i
  ///      => z_2 = -2.3144 - 1.2599 * i

  test('NthRoot_cornercase_thirdRoot_realPartZero', () {
    // complex number with only imaginary part
    Complex z = Complex(0, 2);
    // The List holding all third roots
    List<Complex> thirdRootsOfZ = z.nthRoot(3); //.toArray(Complex[0]);
    // Returned Collection must not be empty!
    expect(3, equals(thirdRootsOfZ.length));
    // test z_0
    expect(1.0911236359717216, closeTo(thirdRootsOfZ[0].real, 1.0e-5));
    expect(0.6299605249474365, closeTo(thirdRootsOfZ[0].imaginary, 1.0e-5));
    // test z_1
    expect(-1.0911236359717216, closeTo(thirdRootsOfZ[1].real, 1.0e-5));
    expect(0.6299605249474365, closeTo(thirdRootsOfZ[1].imaginary, 1.0e-5));
    // test z_2
    expect(-2.3144374213981936E-16, closeTo(thirdRootsOfZ[2].real, 1.0e-5));
    expect(-1.2599210498948732, closeTo(thirdRootsOfZ[2].imaginary, 1.0e-5));
  });

  /// Test cornercases with NaN and Infinity.

  test('NthRoot_cornercase_NAN_Inf', () {
    // NaN + finite -> NaN
    List<Complex> roots = oneNaN.nthRoot(3);
    expect(1, equals(roots.length));
    expect(Complex.nan, equals(roots[0]));

    roots = nanZero.nthRoot(3);
    expect(1, equals(roots.length));
    expect(Complex.nan, equals(roots[0]));

    // NaN + infinite -> NaN
    roots = nanInf.nthRoot(3);
    expect(1, equals(roots.length));
    expect(Complex.nan, equals(roots[0]));

    // finite + infinite -> Inf
    roots = oneInf.nthRoot(3);
    expect(1, equals(roots.length));
    expect(Complex.infinity, equals(roots[0]));

    // infinite + infinite -> Inf
    roots = negInfInf.nthRoot(3);
    expect(1, equals(roots.length));
    expect(Complex.infinity, equals(roots[0]));
  });

  /// Test standard values

  test('argument', () {
    Complex z = Complex(1, 0);
    expect(0.0, closeTo(z.argument(), 1.0e-12));

    z = Complex(1, 1);
    expect(math.pi / 4, closeTo(z.argument(), 1.0e-12));

    z = Complex(0, 1);
    expect(math.pi / 2, closeTo(z.argument(), 1.0e-12));

    z = Complex(-1, 1);
    expect(3 * math.pi / 4, closeTo(z.argument(), 1.0e-12));

    z = Complex(-1, 0);
    expect(math.pi, closeTo(z.argument(), 1.0e-12));

    z = Complex(-1, -1);
    expect(-3 * math.pi / 4, closeTo(z.argument(), 1.0e-12));

    z = Complex(0, -1);
    expect(-math.pi / 2, closeTo(z.argument(), 1.0e-12));

    z = Complex(1, -1);
    expect(-math.pi / 4, closeTo(z.argument(), 1.0e-12));
  });

  /// Verify atan2-style handling of infinite parts

  test('argumentInf', () {
    expect(math.pi / 4, closeTo(infInf.argument(), 1.0e-12));
    expect(math.pi / 2, closeTo(oneInf.argument(), 1.0e-12));
    expect(0.0, closeTo(infOne.argument(), 1.0e-12));
    expect(math.pi / 2, closeTo(zeroInf.argument(), 1.0e-12));
    expect(0.0, closeTo(infZero.argument(), 1.0e-12));
    expect(math.pi, closeTo(negInfOne.argument(), 1.0e-12));
    expect(-3.0 * math.pi / 4, closeTo(negInfNegInf.argument(), 1.0e-12));
    expect(-math.pi / 2, closeTo(oneNegInf.argument(), 1.0e-12));
  });

  /// Verify that either part NaN results in NaN

  test('GetArgumentNaN', () {
    expect(nanZero.argument().isNaN, isTrue);
    expect(zeroNaN.argument().isNaN, isTrue);
    expect(Complex.nan.argument().isNaN, isTrue);
  });

  /*test('Serial', () {
    Complex z = Complex(3.0, 4.0);
    expect(z, equals(TestUtils.serializeAndRecover(z)));
    Complex ncmplx = TestUtils.serializeAndRecover(oneNaN) as Complex;
    expect(nanZero, equals(ncmplx));
    expect(ncmplx.isNaN, isTrue);
    Complex infcmplx = TestUtils.serializeAndRecover(infInf) as Complex;
    expect(infInf, equals(infcmplx));
    expect(infcmplx.isInfinite, isTrue);
    TestComplex tz = TestComplex(3.0, 4.0);
    expect(tz, TestUtils.serializeAndRecover(tz));
    TestComplex ntcmplx = TestUtils.serializeAndRecover(TestComplex.from(oneNaN)) as TestComplex;
    expect(nanZero, equals(ntcmplx));
    expect(ntcmplx.isNaN, isTrue);
    TestComplex inftcmplx = TestUtils.serializeAndRecover(TestComplex.from(infInf)) as TestComplex;
    expect(infInf, inftcmplx);
    expect(inftcmplx.isInfinite, isTrue);
  });*/
}

/// Class to test extending Complex
class TestComplex extends Complex {
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
