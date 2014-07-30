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
import 'dart:math' as math;
import 'package:unittest/unittest.dart';
import 'package:complex/complex.dart';

main() {
  double inf = double.INFINITY;
  double neginf = double.NEGATIVE_INFINITY;
  double nan = double.NAN;
  double pi = math.PI;
  Complex oneInf = new Complex(1, inf);
  Complex oneNegInf = new Complex(1, neginf);
  Complex infOne = new Complex(inf, 1);
  Complex infZero = new Complex(inf, 0);
  Complex infNaN = new Complex(inf, nan);
  Complex infNegInf = new Complex(inf, neginf);
  Complex infInf = new Complex(inf, inf);
  Complex negInfInf = new Complex(neginf, inf);
  Complex negInfZero = new Complex(neginf, 0);
  Complex negInfOne = new Complex(neginf, 1);
  Complex negInfNaN = new Complex(neginf, nan);
  Complex negInfNegInf = new Complex(neginf, neginf);
  Complex oneNaN = new Complex(1, nan);
  Complex zeroInf = new Complex(0, inf);
  Complex zeroNaN = new Complex(0, nan);
  Complex nanInf = new Complex(nan, inf);
  Complex nanNegInf = new Complex(nan, neginf);
  Complex nanZero = new Complex(nan, 0);

  test('Constructor', () {
    Complex z = new Complex(3.0, 4.0);
    Assert.assertEquals(3.0, closeTo(z.real, 1.0e-5));
    Assert.assertEquals(4.0, closeTo(z.imaginary, 1.0e-5));
  });

  test('ConstructorNaN', () {
    Complex z = new Complex(3.0, double.NAN);
    Assert.assertTrue(z.isNaN, isTrue);

    z = new Complex(nan, 4.0);
    Assert.assertTrue(z.isNaN, isTrue);

    z = new Complex(3.0, 4.0);
    Assert.assertFalse(z.isNaN, isFalse);
  });

  test('Abs', () {
    Complex z = new Complex(3.0, 4.0);
    Assert.assertEquals(5.0, closeTo(z.abs(), 1.0e-5));
  });

  test('AbsNaN', () {
    Assert.assertTrue(Complex.NAN.abs().isNaN, isTrue);
    Complex z = new Complex(inf, nan);
    Assert.assertTrue(z.abs().isNaN, isTrue);
  });

  test('AbsInfinite', () {
    Complex z = new Complex(inf, 0);
    Assert.assertEquals(inf, equals(z.abs()));
    z = new Complex(0, neginf);
    Assert.assertEquals(inf, equals(z.abs()));
    z = new Complex(inf, neginf);
    Assert.assertEquals(inf, equals(z.abs()));
  });

  test('Add', () {
    Complex x = new Complex(3.0, 4.0);
    Complex y = new Complex(5.0, 6.0);
    Complex z = x + y;
    Assert.assertEquals(8.0, closeTo(z.real, 1.0e-5));
    Assert.assertEquals(10.0, closeTo(z.imaginary, 1.0e-5));
  });

  test('AddNaN', () {
    Complex x = new Complex(3.0, 4.0);
    Complex z = x + Complex.NAN;
    Assert.assertSame(Complex.NAN, same(z));
    z = new Complex(1, nan);
    Complex w = x + z;
    Assert.assertSame(Complex.NAN, same(w));
  });

  test('AddInf', () {
    Complex x = new Complex(1, 1);
    Complex z = new Complex(inf, 0);
    Complex w = x + z;
    Assert.assertEquals(w.imaginary, equals(1));
    Assert.assertEquals(inf, equals(w.real));

    x = new Complex(neginf, 0);
    Assert.assertTrue((x + z).real.isNaN, isTrue);
  });

  test('ScalarAdd', () {
    Complex x = new Complex(3.0, 4.0);
    double yDouble = 2.0;
    Complex yComplex = new Complex(yDouble);
    Assert.assertEquals(x + yComplex, equals(x + yDouble));
  });

  test('ScalarAddNaN', () {
    Complex x = new Complex(3.0, 4.0);
    double yDouble = double.NAN;
    Complex yComplex = new Complex(yDouble);
    Assert.assertEquals(x + yComplex, equals(x + yDouble));
  });

  test('ScalarAddInf', () {
    Complex x = new Complex(1, 1);
    double yDouble = double.INFINITY;

    Complex yComplex = new Complex(yDouble);
    Assert.assertEquals(x + yComplex, x + yDouble);

    x = new Complex(neginf, 0);
    Assert.assertEquals(x + yComplex, x + yDouble);
  });

  test('Conjugate', () {
    Complex x = new Complex(3.0, 4.0);
    Complex z = x.conjugate();
    Assert.assertEquals(3.0, closeTo(z.real, 1.0e-5));
    Assert.assertEquals(-4.0, closeTo(z.imaginary, 1.0e-5));
  });

  test('ConjugateNaN', () {
    Complex z = Complex.NaN.conjugate();
    Assert.assertTrue(z.isNaN);
  });

  test('ConjugateInfiinite', () {
    Complex z = new Complex(0, inf);
    Assert.assertEquals(neginf, equals(z.conjugate().imaginary));
    z = new Complex(0, neginf);
    Assert.assertEquals(inf, equals(z.conjugate().imaginary));
  });

  test('Divide', () {
    Complex x = new Complex(3.0, 4.0);
    Complex y = new Complex(5.0, 6.0);
    Complex z = x / y;
    Assert.assertEquals(39.0 / 61.0, closeTo(z.real, 1.0e-5));
    Assert.assertEquals(2.0 / 61.0, closeTo(z.imaginary, 1.0e-5));
  });

  test('DivideReal', () {
    Complex x = new Complex(2, 3);
    Complex y = new Complex(2, 0);
    Assert.assertEquals(new Complex(1, 1.5), equals(x / y));
  });

  test('DivideImaginary', () {
    Complex x = new Complex(2, 3);
    Complex y = new Complex(0, 2);
    Assert.assertEquals(new Complex(1.5, -1), equals(x / y));
  });

  test('DivideInf', () {
    Complex x = new Complex(3, 4);
    Complex w = new Complex(neginf, inf);
    Assert.assertTrue(x / w == Complex.ZERO, isTrue);

    Complex z = w / x;
    Assert.assertTrue(z.real.isNaN, isTrue);
    Assert.assertEquals(inf, equals(z.imaginary));

    w = new Complex(inf, inf);
    z = w / x;
    Assert.assertTrue(z.imaginary.isNaN, isTrue);
    Assert.assertEquals(inf, equals(z.real));

    w = new Complex(1, inf);
    z = w / w;
    Assert.assertTrue(z.real.isNaN);
    Assert.assertTrue(z.imaginary.isNaN);
  });

  test('DivideZero', () {
    Complex x = new Complex(3.0, 4.0);
    Complex z = x / Complex.ZERO;
    // Assert.assertEquals(z, Complex.INF); // See MATH-657
    Assert.assertEquals(z, equals(Complex.NAN));
  });

  test('DivideZeroZero', () {
    Complex x = new Complex(0.0, 0.0);
    Complex z = x / Complex.ZERO;
    Assert.assertEquals(z, equals(Complex.NAN));
  });

  test('DivideNaN', () {
    Complex x = new Complex(3.0, 4.0);
    Complex z = x / Complex.NAN;
    Assert.assertTrue(z.isNaN, isTrue);
  });

  test('DivideNaNInf', () {
    Complex z = oneInf / Complex.ONE;
    Assert.assertTrue(z.real.isNaN, isTrue);
    Assert.assertEquals(inf, equals(z.imaginary));

    z = negInfNegInf / oneNaN;
    Assert.assertTrue(z.real.isNaN, isTrue);
    Assert.assertTrue(z.imaginary.isNaN, isTrue);

    z = negInfInf / Complex.ONE;
    Assert.assertTrue(z.real.isNaN, isTrue);
    Assert.assertTrue(z.imaginary.isNaN, isTrue);
  });

  test('ScalarDivide', () {
    Complex x = new Complex(3.0, 4.0);
    double yDouble = 2.0;
    Complex yComplex = new Complex(yDouble);
    Assert.assertEquals(x / yComplex, equals(x / yDouble));
  });

  test('ScalarDivideNaN', () {
    Complex x = new Complex(3.0, 4.0);
    double yDouble = double.NAN;
    Complex yComplex = new Complex(yDouble);
    Assert.assertEquals(x / yComplex, equals(x / yDouble));
  });

  test('ScalarDivideInf', () {
    Complex x = new Complex(1, 1);
    double yDouble = double.INFINITY;
    Complex yComplex = new Complex(yDouble);
    TestUtils.assertEquals(x / yComplex, equals(x / yDouble));

    yDouble = double.NEGATIVE_INFINITY;
    yComplex = new Complex(yDouble);
    TestUtils.assertEquals(x / yComplex, equals(x / yDouble));

    x = new Complex(1, double.NEGATIVE_INFINITY);
    TestUtils.assertEquals(x / yComplex, equals(x / yDouble));
  });

  test('ScalarDivideZero', () {
    Complex x = new Complex(1, 1);
    TestUtils.assertEquals(x / Complex.ZERO, equals(x / 0));
  });

  test('Reciprocal', () {
    Complex z = new Complex(5.0, 6.0);
    Complex act = z.reciprocal();
    double expRe = 5.0 / 61.0;
    double expIm = -6.0 / 61.0;
    Assert.assertEquals(expRe, closeTo(act.real, FastMath.ulp(expRe)));
    Assert.assertEquals(expIm, closeTo(act.imaginary, FastMath.ulp(expIm)));
  });

  test('ReciprocalReal', () {
    Complex z = new Complex(-2.0, 0.0);
    Assert.assertTrue(Complex.equals(new Complex(-0.5, 0.0), z.reciprocal()), isTrue);
  });


  test('ReciprocalImaginary', () {
    Complex z = new Complex(0.0, -2.0);
    Assert.assertEquals(new Complex(0.0, 0.5), equals(z.reciprocal()));
  });


  test('ReciprocalInf', () {
    Complex z = new Complex(neginf, inf);
    Assert.assertTrue(z.reciprocal() == Complex.ZERO, isTrue);

    z = new Complex(1, inf).reciprocal();
    Assert.assertEquals(z, equals(Complex.ZERO));
  });


  test('ReciprocalZero', () {
    Assert.assertEquals(Complex.ZERO.reciprocal(), equals(Complex.INFINITY));
  });


  test('ReciprocalNaN', () {
    Assert.assertTrue(Complex.NAN.reciprocal().isNaN, isTrue);
  });


  test('Multiply', () {
    Complex x = new Complex(3.0, 4.0);
    Complex y = new Complex(5.0, 6.0);
    Complex z = x * y;
    Assert.assertEquals(-9.0, closeTo(z.real, 1.0e-5));
    Assert.assertEquals(38.0, closeTo(z.imaginary, 1.0e-5));
  });


  test('MultiplyNaN', () {
    Complex x = new Complex(3.0, 4.0);
    Complex z = x * Complex.NAN;
    Assert.assertSame(Complex.NaN, same(z));
    z = Complex.NaN * 5;
    Assert.assertSame(Complex.NaN, same(z));
  });


  test('MultiplyInfInf', () {
    // Assert.assertTrue(infInf.multiply(infInf).isNaN()); // MATH-620
    Assert.assertTrue((infInf * infInf).isInfinite, isTrue);
  });


  test('MultiplyNaNInf', () {
    Complex z = new Complex(1, 1);
    Complex w = z * infOne;
    Assert.assertEquals(w.real, equals(inf));
    Assert.assertEquals(w.imaginary, equals(inf));

    // [MATH-164]
    Assert.assertTrue(new Complex(1, 0) * infInf == Complex.INFINITY, isTrue);
    Assert.assertTrue(new Complex(-1, 0) * infInf == Complex.INFINITY, isTrue);
    Assert.assertTrue(new Complex(1, 0) * negInfZero == Complex.INFINITY, isTrue);

    w = oneInf * oneNegInf;
    Assert.assertEquals(w.real, equals(inf));
    Assert.assertEquals(w.imaginary, equals(inf));

    w = negInfNegInf * oneNaN;
    Assert.assertTrue(w.real.isNaN, isTrue);
    Assert.assertTrue(w.imaginary.isNaN, isTrue);

    z = new Complex(1, neginf);
    Assert.assertSame(Complex.INFINITY, same(z * z));
  });


  test('ScalarMultiply', () {
    Complex x = new Complex(3.0, 4.0);
    double yDouble = 2.0;
    Complex yComplex = new Complex(yDouble);
    Assert.assertEquals(x * yComplex, equals(x * yDouble));
    int zInt = -5;
    Complex zComplex = new Complex(zInt);
    Assert.assertEquals(x * zComplex, equals(x * zInt));
  });


  test('ScalarMultiplyNaN', () {
    Complex x = new Complex(3.0, 4.0);
    double yDouble = double.NAN;
    Complex yComplex = new Complex(yDouble);
    Assert.assertEquals(x * yComplex, equals(x * yDouble));
  });


  test('ScalarMultiplyInf', () {
    Complex x = new Complex(1, 1);
    double yDouble = double.INFINITY;
    Complex yComplex = new Complex(yDouble);
    Assert.assertEquals(x * yComplex, equals(x * yDouble));

    yDouble = double.NEGATIVE_INFINITY;
    yComplex = new Complex(yDouble);
    Assert.assertEquals(x * yComplex, equals(x * yDouble));
  });


  test('Negate', () {
    Complex x = new Complex(3.0, 4.0);
    Complex z = x.negate();
    Assert.assertEquals(-3.0, closeTo(z.real, 1.0e-5));
    Assert.assertEquals(-4.0, closeTo(z.imaginary, 1.0e-5));
  });


  test('NegateNaN', () {
    Complex z = Complex.NaN.negate();
    Assert.assertTrue(z.isNaN, isTrue);
  });


  test('Subtract', () {
    Complex x = new Complex(3.0, 4.0);
    Complex y = new Complex(5.0, 6.0);
    Complex z = x - y;
    Assert.assertEquals(-2.0, closeTo(z.real, 1.0e-5));
    Assert.assertEquals(-2.0, closeTo(z.imaginary, 1.0e-5));
  });


  test('SubtractNaN', () {
    Complex x = new Complex(3.0, 4.0);
    Complex z = x - Complex.NAN;
    Assert.assertSame(Complex.NAN, same(z));
    z = new Complex(1, nan);
    Complex w = x - z;
    Assert.assertSame(Complex.NAN, same(w));
  });


  test('SubtractInf', () {
    Complex x = new Complex(1, 1);
    Complex z = new Complex(neginf, 0);
    Complex w = x - z;
    Assert.assertEquals(w.imaginary, equals(1));
    Assert.assertEquals(inf, equals(w.real));

    x = new Complex(neginf, 0);
    Assert.assertTrue((x - z).real.isNaN, isTrue);
  });


  test('ScalarSubtract', () {
    Complex x = new Complex(3.0, 4.0);
    double yDouble = 2.0;
    Complex yComplex = new Complex(yDouble);
    Assert.assertEquals(x - yComplex, equals(x - yDouble));
  });


  test('ScalarSubtractNaN', () {
    Complex x = new Complex(3.0, 4.0);
    double yDouble = double.NAN;
    Complex yComplex = new Complex(yDouble);
    Assert.assertEquals(x - yComplex, equals(x - yDouble));
  });


  test('ScalarSubtractInf', () {
    Complex x = new Complex(1, 1);
    double yDouble = double.INFINITY;
    Complex yComplex = new Complex(yDouble);
    Assert.assertEquals(x - yComplex, equals(x - yDouble));

    x = new Complex(neginf, 0);
    Assert.assertEquals(x - yComplex, equals(x - yDouble));
  });



  test('EqualsNull', () {
    Complex x = new Complex(3.0, 4.0);
    Assert.assertFalse(x == null, isFalse);
  });

  test('FloatingPointEqualsPrecondition1', () {
    expect(() {
      Complex.equals(new Complex(3.0, 4.0), null, 3);
    }, throwsA(NullPointerException));
  });

  test('FloatingPointEqualsPrecondition2', () {
    expect(() {
      Complex.equals(null, new Complex(3.0, 4.0), 3);
    }, throwsA(NullPointerException));
  });


  test('EqualsClass', () {
    Complex x = new Complex(3.0, 4.0);
    Assert.assertFalse(x == TextComplex, isFalse);
  });


  test('EqualsSame', () {
    Complex x = new Complex(3.0, 4.0);
    Assert.assertTrue(x == x, isTrue);
  });


  test('FloatingPointEquals', () {
    double re = -3.21;
    double im = 456789e10;

    final Complex x = new Complex(re, im);
    Complex y = new Complex(re, im);

    Assert.assertTrue(x == y, isTrue);
    Assert.assertTrue(Complex.equals(x, y), isTrue);

    final int maxUlps = 5;
    for (int i = 0; i < maxUlps; i++) {
      re = FastMath.nextUp(re);
      im = FastMath.nextUp(im);
    }
    y = new Complex(re, im);
    Assert.assertTrue(Complex.equals(x, y, maxUlps), isTrue);

    re = FastMath.nextUp(re);
    im = FastMath.nextUp(im);
    y = new Complex(re, im);
    Assert.assertFalse(Complex.equals(x, y, maxUlps), isFalse);
  });


  test('FloatingPointEqualsNaN', () {
    Complex c = new Complex(double.NAN, 1);
    Assert.assertFalse(Complex.equals(c, c), isFalse);

    c = new Complex(1, double.NAN);
    Assert.assertFalse(Complex.equals(c, c), isFalse);
  });


  test('FloatingPointEqualsWithAllowedDelta', () {
    final double re = 153.0000;
    final double im = 152.9375;
    final double tol1 = 0.0625;
    final Complex x = new Complex(re, im);
    final Complex y = new Complex(re + tol1, im + tol1);
    Assert.assertTrue(Complex.equals(x, y, tol1), isTrue);

    final double tol2 = 0.0624;
    Assert.assertFalse(Complex.equals(x, y, tol2), isFalse);
  });


  test('FloatingPointEqualsWithRelativeTolerance', () {
    final double tol = 1e-4;
    final double re = 1;
    final double im = 1e10;

    final double f = 1 + tol;
    final Complex x = new Complex(re, im);
    final Complex y = new Complex(re * f, im * f);
    Assert.assertTrue(Complex.equalsWithRelativeTolerance(x, y, tol), isTrue);
  });


  test('EqualsTrue', () {
    Complex x = new Complex(3.0, 4.0);
    Complex y = new Complex(3.0, 4.0);
    Assert.assertTrue(x == y, isTrue);
  });


  test('EqualsRealDifference', () {
    Complex x = new Complex(0.0, 0.0);
    Complex y = new Complex(0.0 + double.MIN_POSITIVE, 0.0);
    Assert.assertFalse(x == y, isFalse);
  });


  test('EqualsImaginaryDifference', () {
    Complex x = new Complex(0.0, 0.0);
    Complex y = new Complex(0.0, 0.0 + double.MIN_POSITIVE);
    Assert.assertFalse(x == y, isFalse);
  });


  test('EqualsNaN', () {
    Complex realNaN = new Complex(double.NAN, 0.0);
    Complex imaginaryNaN = new Complex(0.0, double.NAN);
    Complex complexNaN = Complex.NAN;
    Assert.assertTrue(realNaN == imaginaryNaN, isTrue);
    Assert.assertTrue(imaginaryNaN == complexNaN, isTrue);
    Assert.assertTrue(realNaN == complexNaN, isTrue);
  });


  test('HashCode', () {
    Complex x = new Complex(0.0, 0.0);
    Complex y = new Complex(0.0, 0.0 + double.MIN_POSITIVE);
    Assert.assertFalse(x.hashCode == y.hashCode, isFalse);
    y = new Complex(0.0 + double.MIN_POSITIVE, 0.0);
    Assert.assertFalse(x.hashCode == y.hashCode, isFalse);
    Complex realNaN = new Complex(double.NAN, 0.0);
    Complex imaginaryNaN = new Complex(0.0, double.NAN);
    Assert.assertEquals(realNaN.hashCode, equals(imaginaryNaN.hashCode));
    Assert.assertEquals(imaginaryNaN.hashCode, equals(Complex.NAN.hashCode));

    // MATH-1118
    // "equals" and "hashCode" must be compatible: if two objects have
    // different hash codes, "equals" must return false.
    final String msg = "'equals' not compatible with 'hashCode'";

    x = new Complex(0.0, 0.0);
    y = new Complex(0.0, -0.0);
    Assert.assertTrue(x.hashCode != y.hashCode, isTrue);
    Assert.assertFalse(x == y, isFalse, msg);

    x = new Complex(0.0, 0.0);
    y = new Complex(-0.0, 0.0);
    Assert.assertTrue(x.hashCode != y.hashCode, isTrue);
    Assert.assertFalse(x == y, isFalse, msg);
  });


  test('Acos', () {
    Complex z = new Complex(3, 4);
    Complex expected = new Complex(0.936812, -2.30551);
    TestUtils.assertEquals(expected, closeTo(z.acos(), 1.0e-5));
    TestUtils.assertEquals(new Complex(fastmath.acos(0), 0), closeTo(Complex.ZERO.acos(), 1.0e-12));
  });


  test('AcosInf', () {
    TestUtils.assertSame(Complex.NaN, oneInf.acos());
    TestUtils.assertSame(Complex.NaN, oneNegInf.acos());
    TestUtils.assertSame(Complex.NaN, infOne.acos());
    TestUtils.assertSame(Complex.NaN, negInfOne.acos());
    TestUtils.assertSame(Complex.NaN, infInf.acos());
    TestUtils.assertSame(Complex.NaN, infNegInf.acos());
    TestUtils.assertSame(Complex.NaN, negInfInf.acos());
    TestUtils.assertSame(Complex.NaN, negInfNegInf.acos());
  });


  test('AcosNaN', () {
    Assert.assertTrue(Complex.NaN.acos().isNaN());
  });


  test('Asin', () {
    Complex z = new Complex(3, 4);
    Complex expected = new Complex(0.633984, 2.30551);
    TestUtils.assertEquals(expected, closeTo(z.asin(), 1.0e-5));
  });


  test('AsinNaN', () {
    Assert.assertTrue(Complex.NaN.asin().isNaN());
  });


  test('AsinInf', () {
    TestUtils.assertSame(Complex.NaN, oneInf.asin());
    TestUtils.assertSame(Complex.NaN, oneNegInf.asin());
    TestUtils.assertSame(Complex.NaN, infOne.asin());
    TestUtils.assertSame(Complex.NaN, negInfOne.asin());
    TestUtils.assertSame(Complex.NaN, infInf.asin());
    TestUtils.assertSame(Complex.NaN, infNegInf.asin());
    TestUtils.assertSame(Complex.NaN, negInfInf.asin());
    TestUtils.assertSame(Complex.NaN, negInfNegInf.asin());
  });



  test('Atan', () {
    Complex z = new Complex(3, 4);
    Complex expected = new Complex(1.44831, 0.158997);
    TestUtils.assertEquals(expected, z.atan(), 1.0e-5);
  });


  test('AtanInf', () {
    TestUtils.assertSame(Complex.NaN, oneInf.atan());
    TestUtils.assertSame(Complex.NaN, oneNegInf.atan());
    TestUtils.assertSame(Complex.NaN, infOne.atan());
    TestUtils.assertSame(Complex.NaN, negInfOne.atan());
    TestUtils.assertSame(Complex.NaN, infInf.atan());
    TestUtils.assertSame(Complex.NaN, infNegInf.atan());
    TestUtils.assertSame(Complex.NaN, negInfInf.atan());
    TestUtils.assertSame(Complex.NaN, negInfNegInf.atan());
  });


  test('AtanI', () {
    Assert.assertTrue(Complex.I.atan().isNaN());
  });


  test('AtanNaN', () {
    Assert.assertTrue(Complex.NaN.atan().isNaN());
  });


  test('Cos', () {
    Complex z = new Complex(3, 4);
    Complex expected = new Complex(-27.03495, -3.851153);
    TestUtils.assertEquals(expected, closeTo(z.cos(), 1.0e-5));
  });


  test('CosNaN', () {
    Assert.assertTrue(Complex.NaN.cos().isNaN());
  });


  test('CosInf', () {
    TestUtils.assertSame(infNegInf, oneInf.cos());
    TestUtils.assertSame(infInf, oneNegInf.cos());
    TestUtils.assertSame(Complex.NaN, infOne.cos());
    TestUtils.assertSame(Complex.NaN, negInfOne.cos());
    TestUtils.assertSame(Complex.NaN, infInf.cos());
    TestUtils.assertSame(Complex.NaN, infNegInf.cos());
    TestUtils.assertSame(Complex.NaN, negInfInf.cos());
    TestUtils.assertSame(Complex.NaN, negInfNegInf.cos());
  });


  test('Cosh', () {
    Complex z = new Complex(3, 4);
    Complex expected = new Complex(-6.58066, -7.58155);
    TestUtils.assertEquals(expected, closeTo(z.cosh(), 1.0e-5));
  });


  test('CoshNaN', () {
    Assert.assertTrue(Complex.NaN.cosh().isNaN());
  });


  test('CoshInf', () {
    TestUtils.assertSame(Complex.NaN, oneInf.cosh());
    TestUtils.assertSame(Complex.NaN, oneNegInf.cosh());
    TestUtils.assertSame(infInf, infOne.cosh());
    TestUtils.assertSame(infNegInf, negInfOne.cosh());
    TestUtils.assertSame(Complex.NaN, infInf.cosh());
    TestUtils.assertSame(Complex.NaN, infNegInf.cosh());
    TestUtils.assertSame(Complex.NaN, negInfInf.cosh());
    TestUtils.assertSame(Complex.NaN, negInfNegInf.cosh());
  });


  test('Exp', () {
    Complex z = new Complex(3, 4);
    Complex expected = new Complex(-13.12878, -15.20078);
    TestUtils.assertEquals(expected, closeTo(z.exp(), 1.0e-5));
    TestUtils.assertEquals(Complex.ONE, closeTo(Complex.ZERO.exp(), 10e-12));
    Complex iPi = Complex.I.multiply(new Complex(pi, 0));
    TestUtils.assertEquals(Complex.ONE.negate(), closeTo(iPi.exp(), 10e-12));
  });


  test('ExpNaN', () {
    Assert.assertTrue(Complex.NaN.exp().isNaN());
  });


  test('ExpInf', () {
    TestUtils.assertSame(Complex.NaN, oneInf.exp());
    TestUtils.assertSame(Complex.NaN, oneNegInf.exp());
    TestUtils.assertSame(infInf, infOne.exp());
    TestUtils.assertSame(Complex.ZERO, negInfOne.exp());
    TestUtils.assertSame(Complex.NaN, infInf.exp());
    TestUtils.assertSame(Complex.NaN, infNegInf.exp());
    TestUtils.assertSame(Complex.NaN, negInfInf.exp());
    TestUtils.assertSame(Complex.NaN, negInfNegInf.exp());
  });


  test('Log', () {
    Complex z = new Complex(3, 4);
    Complex expected = new Complex(1.60944, 0.927295);
    TestUtils.assertEquals(expected, closeTo(z.log(), 1.0e-5));
  });


  test('LogNaN', () {
    Assert.assertTrue(Complex.NaN.log().isNaN());
  });


  test('LogInf', () {
    TestUtils.assertEquals(new Complex(inf, pi / 2), closeTo(oneInf.log(), 10e-12));
    TestUtils.assertEquals(new Complex(inf, -pi / 2), closeTo(oneNegInf.log(), 10e-12));
    TestUtils.assertEquals(infZero, closeTo(infOne.log(), 10e-12));
    TestUtils.assertEquals(new Complex(inf, pi), closeTo(negInfOne.log(), 10e-12));
    TestUtils.assertEquals(new Complex(inf, pi / 4), closeTo(infInf.log(), 10e-12));
    TestUtils.assertEquals(new Complex(inf, -pi / 4), closeTo(infNegInf.log(), 10e-12));
    TestUtils.assertEquals(new Complex(inf, 3 * pi / 4), closeTo(negInfInf.log(), 10e-12));
    TestUtils.assertEquals(new Complex(inf, -3 * pi / 4), closeTo(negInfNegInf.log(), 10e-12));
  });


  test('LogZero', () {
    TestUtils.assertSame(negInfZero, Complex.ZERO.log());
  });


  test('Pow', () {
    Complex x = new Complex(3, 4);
    Complex y = new Complex(5, 6);
    Complex expected = new Complex(-1.860893, 11.83677);
    TestUtils.assertEquals(expected, closeTo(x.pow(y), 1.0e-5));
  });


  test('PowNaNBase', () {
    Complex x = new Complex(3, 4);
    Assert.assertTrue(Complex.NaN.pow(x).isNaN());
  });


  test('PowNaNExponent', () {
    Complex x = new Complex(3, 4);
    Assert.assertTrue(x.pow(Complex.NaN).isNaN());
  });


  test('PowInf', () {
    TestUtils.assertSame(Complex.NaN, Complex.ONE.pow(oneInf));
    TestUtils.assertSame(Complex.NaN, Complex.ONE.pow(oneNegInf));
    TestUtils.assertSame(Complex.NaN, Complex.ONE.pow(infOne));
    TestUtils.assertSame(Complex.NaN, Complex.ONE.pow(infInf));
    TestUtils.assertSame(Complex.NaN, Complex.ONE.pow(infNegInf));
    TestUtils.assertSame(Complex.NaN, Complex.ONE.pow(negInfInf));
    TestUtils.assertSame(Complex.NaN, Complex.ONE.pow(negInfNegInf));
    TestUtils.assertSame(Complex.NaN, infOne.pow(Complex.ONE));
    TestUtils.assertSame(Complex.NaN, negInfOne.pow(Complex.ONE));
    TestUtils.assertSame(Complex.NaN, infInf.pow(Complex.ONE));
    TestUtils.assertSame(Complex.NaN, infNegInf.pow(Complex.ONE));
    TestUtils.assertSame(Complex.NaN, negInfInf.pow(Complex.ONE));
    TestUtils.assertSame(Complex.NaN, negInfNegInf.pow(Complex.ONE));
    TestUtils.assertSame(Complex.NaN, negInfNegInf.pow(infNegInf));
    TestUtils.assertSame(Complex.NaN, negInfNegInf.pow(negInfNegInf));
    TestUtils.assertSame(Complex.NaN, negInfNegInf.pow(infInf));
    TestUtils.assertSame(Complex.NaN, infInf.pow(infNegInf));
    TestUtils.assertSame(Complex.NaN, infInf.pow(negInfNegInf));
    TestUtils.assertSame(Complex.NaN, infInf.pow(infInf));
    TestUtils.assertSame(Complex.NaN, infNegInf.pow(infNegInf));
    TestUtils.assertSame(Complex.NaN, infNegInf.pow(negInfNegInf));
    TestUtils.assertSame(Complex.NaN, infNegInf.pow(infInf));
  });


  test('PowZero', () {
    TestUtils.assertSame(Complex.NaN, Complex.ZERO.pow(Complex.ONE));
    TestUtils.assertSame(Complex.NaN, Complex.ZERO.pow(Complex.ZERO));
    TestUtils.assertSame(Complex.NaN, Complex.ZERO.pow(Complex.I));
    TestUtils.assertEquals(Complex.ONE, closeTo(Complex.ONE.pow(Complex.ZERO), 10e-12));
    TestUtils.assertEquals(Complex.ONE, closeTo(Complex.I.pow(Complex.ZERO), 10e-12));
    TestUtils.assertEquals(Complex.ONE, closeTo(new Complex(-1, 3).pow(Complex.ZERO), 10e-12));
  });


  test('ScalarPow', () {
    Complex x = new Complex(3, 4);
    double yDouble = 5.0;
    Complex yComplex = new Complex(yDouble);
    Assert.assertEquals(x.pow(yComplex), x.pow(yDouble));
  });


  test('ScalarPowNaNBase', () {
    Complex x = Complex.NaN;
    double yDouble = 5.0;
    Complex yComplex = new Complex(yDouble);
    Assert.assertEquals(x.pow(yComplex), x.pow(yDouble));
  });


  test('ScalarPowNaNExponent', () {
    Complex x = new Complex(3, 4);
    double yDouble = Double.NaN;
    Complex yComplex = new Complex(yDouble);
    Assert.assertEquals(x.pow(yComplex), x.pow(yDouble));
  });


  test('ScalarPowInf', () {
    TestUtils.assertSame(Complex.NaN, Complex.ONE.pow(Double.POSITIVE_INFINITY));
    TestUtils.assertSame(Complex.NaN, Complex.ONE.pow(Double.NEGATIVE_INFINITY));
    TestUtils.assertSame(Complex.NaN, infOne.pow(1.0));
    TestUtils.assertSame(Complex.NaN, negInfOne.pow(1.0));
    TestUtils.assertSame(Complex.NaN, infInf.pow(1.0));
    TestUtils.assertSame(Complex.NaN, infNegInf.pow(1.0));
    TestUtils.assertSame(Complex.NaN, negInfInf.pow(10));
    TestUtils.assertSame(Complex.NaN, negInfNegInf.pow(1.0));
    TestUtils.assertSame(Complex.NaN, negInfNegInf.pow(Double.POSITIVE_INFINITY));
    TestUtils.assertSame(Complex.NaN, negInfNegInf.pow(Double.POSITIVE_INFINITY));
    TestUtils.assertSame(Complex.NaN, infInf.pow(Double.POSITIVE_INFINITY));
    TestUtils.assertSame(Complex.NaN, infInf.pow(Double.NEGATIVE_INFINITY));
    TestUtils.assertSame(Complex.NaN, infNegInf.pow(Double.NEGATIVE_INFINITY));
    TestUtils.assertSame(Complex.NaN, infNegInf.pow(Double.POSITIVE_INFINITY));
  });


  test('ScalarPowZero', () {
    TestUtils.assertSame(Complex.NaN, Complex.ZERO.pow(1.0));
    TestUtils.assertSame(Complex.NaN, Complex.ZERO.pow(0.0));
    TestUtils.assertEquals(Complex.ONE, closeTo(Complex.ONE.pow(0.0), 10e-12));
    TestUtils.assertEquals(Complex.ONE, closeTo(Complex.I.pow(0.0), 10e-12));
    TestUtils.assertEquals(Complex.ONE, closeTo(new Complex(-1, 3).pow(0.0), 10e-12));
  });

  //(expected=NullArgumentException.class)
  test('powNull', () {
    Complex.ONE.pow(null);
  });


  test('Sin', () {
    Complex z = new Complex(3, 4);
    Complex expected = new Complex(3.853738, -27.01681);
    TestUtils.assertEquals(expected, closeTo(z.sin(), 1.0e-5));
  });


  test('SinInf', () {
    TestUtils.assertSame(infInf, oneInf.sin());
    TestUtils.assertSame(infNegInf, oneNegInf.sin());
    TestUtils.assertSame(Complex.NaN, infOne.sin());
    TestUtils.assertSame(Complex.NaN, negInfOne.sin());
    TestUtils.assertSame(Complex.NaN, infInf.sin());
    TestUtils.assertSame(Complex.NaN, infNegInf.sin());
    TestUtils.assertSame(Complex.NaN, negInfInf.sin());
    TestUtils.assertSame(Complex.NaN, negInfNegInf.sin());
  });


  test('SinNaN', () {
    Assert.assertTrue(Complex.NaN.sin().isNaN());
  });


  test('Sinh', () {
    Complex z = new Complex(3, 4);
    Complex expected = new Complex(-6.54812, -7.61923);
    TestUtils.assertEquals(expected, closeTo(z.sinh(), 1.0e-5));
  });


  test('SinhNaN', () {
    Assert.assertTrue(Complex.NaN.sinh().isNaN());
  });


  test('SinhInf', () {
    TestUtils.assertSame(Complex.NaN, oneInf.sinh());
    TestUtils.assertSame(Complex.NaN, oneNegInf.sinh());
    TestUtils.assertSame(infInf, infOne.sinh());
    TestUtils.assertSame(negInfInf, negInfOne.sinh());
    TestUtils.assertSame(Complex.NaN, infInf.sinh());
    TestUtils.assertSame(Complex.NaN, infNegInf.sinh());
    TestUtils.assertSame(Complex.NaN, negInfInf.sinh());
    TestUtils.assertSame(Complex.NaN, negInfNegInf.sinh());
  });


  test('SqrtRealPositive', () {
    Complex z = new Complex(3, 4);
    Complex expected = new Complex(2, 1);
    TestUtils.assertEquals(expected, closeTo(z.sqrt(), 1.0e-5));
  });


  test('SqrtRealZero', () {
    Complex z = new Complex(0.0, 4);
    Complex expected = new Complex(1.41421, 1.41421);
    TestUtils.assertEquals(expected, closeTo(z.sqrt(), 1.0e-5));
  });


  test('SqrtRealNegative', () {
    Complex z = new Complex(-3.0, 4);
    Complex expected = new Complex(1, 2);
    TestUtils.assertEquals(expected, closeTo(z.sqrt(), 1.0e-5));
  });


  test('SqrtImaginaryZero', () {
    Complex z = new Complex(-3.0, 0.0);
    Complex expected = new Complex(0.0, 1.73205);
    TestUtils.assertEquals(expected, closeTo(z.sqrt(), 1.0e-5));
  });


  test('SqrtImaginaryNegative', () {
    Complex z = new Complex(-3.0, -4.0);
    Complex expected = new Complex(1.0, -2.0);
    TestUtils.assertEquals(expected, closeTo(z.sqrt(), 1.0e-5));
  });


  test('SqrtPolar', () {
    double r = 1;
    for (int i = 0; i < 5; i++) {
      r += i;
      double theta = 0;
      for (int j = 0; j < 11; j++) {
        theta += pi / 12;
        Complex z = ComplexUtils.polar2Complex(r, theta);
        Complex sqrtz = ComplexUtils.polar2Complex(FastMath.sqrt(r), theta / 2);
        TestUtils.assertEquals(sqrtz, closeTo(z.sqrt(), 10e-12));
      }
    }
  });


  test('SqrtNaN', () {
    Assert.assertTrue(Complex.NaN.sqrt().isNaN());
  });


  test('SqrtInf', () {
    TestUtils.assertSame(infNaN, oneInf.sqrt());
    TestUtils.assertSame(infNaN, oneNegInf.sqrt());
    TestUtils.assertSame(infZero, infOne.sqrt());
    TestUtils.assertSame(zeroInf, negInfOne.sqrt());
    TestUtils.assertSame(infNaN, infInf.sqrt());
    TestUtils.assertSame(infNaN, infNegInf.sqrt());
    TestUtils.assertSame(nanInf, negInfInf.sqrt());
    TestUtils.assertSame(nanNegInf, negInfNegInf.sqrt());
  });


  test('Sqrt1z', () {
    Complex z = new Complex(3, 4);
    Complex expected = new Complex(4.08033, -2.94094);
    TestUtils.assertEquals(expected, closeTo(z.sqrt1z(), 1.0e-5));
  });


  test('Sqrt1zNaN', () {
    Assert.assertTrue(Complex.NaN.sqrt1z().isNaN());
  });


  test('Tan', () {
    Complex z = new Complex(3, 4);
    Complex expected = new Complex(-0.000187346, 0.999356);
    TestUtils.assertEquals(expected, closeTo(z.tan(), 1.0e-5));
    /* Check that no overflow occurs (MATH-722) */
    Complex actual = new Complex(3.0, 1E10).tan();
    expected = new Complex(0, 1);
    TestUtils.assertEquals(expected, closeTo(actual, 1.0e-5));
    actual = new Complex(3.0, -1E10).tan();
    expected = new Complex(0, -1);
    TestUtils.assertEquals(expected, closeTo(actual, 1.0e-5));
  });


  test('TanNaN', () {
    Assert.assertTrue(Complex.NaN.tan().isNaN());
  });


  test('TanInf', () {
    TestUtils.assertSame(Complex.valueOf(0.0, 1.0), oneInf.tan());
    TestUtils.assertSame(Complex.valueOf(0.0, -1.0), oneNegInf.tan());
    TestUtils.assertSame(Complex.NaN, infOne.tan());
    TestUtils.assertSame(Complex.NaN, negInfOne.tan());
    TestUtils.assertSame(Complex.NaN, infInf.tan());
    TestUtils.assertSame(Complex.NaN, infNegInf.tan());
    TestUtils.assertSame(Complex.NaN, negInfInf.tan());
    TestUtils.assertSame(Complex.NaN, negInfNegInf.tan());
  });


  test('TanCritical', () {
    TestUtils.assertSame(infNaN, new Complex(pi / 2, 0).tan());
    TestUtils.assertSame(negInfNaN, new Complex(-pi / 2, 0).tan());
  });


  test('Tanh', () {
    Complex z = new Complex(3, 4);
    Complex expected = new Complex(1.00071, 0.00490826);
    TestUtils.assertEquals(expected, closeTo(z.tanh(), 1.0e-5));
    /* Check that no overflow occurs (MATH-722) */
    Complex actual = new Complex(1E10, 3.0).tanh();
    expected = new Complex(1, 0);
    TestUtils.assertEquals(expected, closeTo(actual, 1.0e-5));
    actual = new Complex(-1E10, 3.0).tanh();
    expected = new Complex(-1, 0);
    TestUtils.assertEquals(expected, closeTo(actual, 1.0e-5));
  });


  test('TanhNaN', () {
    Assert.assertTrue(Complex.NaN.tanh().isNaN());
  });


  test('TanhInf', () {
    TestUtils.assertSame(Complex.NaN, oneInf.tanh());
    TestUtils.assertSame(Complex.NaN, oneNegInf.tanh());
    TestUtils.assertSame(Complex.valueOf(1.0, 0.0), infOne.tanh());
    TestUtils.assertSame(Complex.valueOf(-1.0, 0.0), negInfOne.tanh());
    TestUtils.assertSame(Complex.NaN, infInf.tanh());
    TestUtils.assertSame(Complex.NaN, infNegInf.tanh());
    TestUtils.assertSame(Complex.NaN, negInfInf.tanh());
    TestUtils.assertSame(Complex.NaN, negInfNegInf.tanh());
  });


  test('TanhCritical', () {
    TestUtils.assertSame(nanInf, new Complex(0, pi / 2).tanh());
  });

  /** test issue MATH-221 */

  test('Math221', () {
    Assert.assertTrue(Complex.equals(new Complex(0, -1), new Complex(0, 1).multiply(new Complex(-1, 0))));
  });

  /**
     * Test: computing <b>third roots</b> of z.
     * <pre>
     * <code>
     * <b>z = -2 + 2 * i</b>
     *   => z_0 =  1      +          i
     *   => z_1 = -1.3660 + 0.3660 * i
     *   => z_2 =  0.3660 - 1.3660 * i
     * </code>
     * </pre>
     */

  test('NthRoot_normal_thirdRoot', () {
    // The complex number we want to compute all third-roots for.
    Complex z = new Complex(-2, 2);
    // The List holding all third roots
    List<Complex> thirdRootsOfZ = z.nthRoot(3);//.toArray(new Complex[0]);
    // Returned Collection must not be empty!
    Assert.assertEquals(3, thirdRootsOfZ.length);
    // test z_0
    Assert.assertEquals(1.0, closeTo(thirdRootsOfZ[0].real, 1.0e-5));
    Assert.assertEquals(1.0, closeTo(thirdRootsOfZ[0].imaginary, 1.0e-5));
    // test z_1
    Assert.assertEquals(-1.3660254037844386, closeTo(thirdRootsOfZ[1].real, 1.0e-5));
    Assert.assertEquals(0.36602540378443843, closeTo(thirdRootsOfZ[1].imaginary, 1.0e-5));
    // test z_2
    Assert.assertEquals(0.366025403784439, closeTo(thirdRootsOfZ[2].real, 1.0e-5));
    Assert.assertEquals(-1.3660254037844384, closeTo(thirdRootsOfZ[2].imaginary, 1.0e-5));
  });


  /**
     * Test: computing <b>fourth roots</b> of z.
     * <pre>
     * <code>
     * <b>z = 5 - 2 * i</b>
     *   => z_0 =  1.5164 - 0.1446 * i
     *   => z_1 =  0.1446 + 1.5164 * i
     *   => z_2 = -1.5164 + 0.1446 * i
     *   => z_3 = -1.5164 - 0.1446 * i
     * </code>
     * </pre>
     */

  test('NthRoot_normal_fourthRoot', () {
    // The complex number we want to compute all third-roots for.
    Complex z = new Complex(5, -2);
    // The List holding all fourth roots
    List<Complex> fourthRootsOfZ = z.nthRoot(4);//.toArray(new Complex[0]);
    // Returned Collection must not be empty!
    Assert.assertEquals(4, fourthRootsOfZ.length);
    // test z_0
    Assert.assertEquals(1.5164629308487783, closeTo(fourthRootsOfZ[0].real, 1.0e-5));
    Assert.assertEquals(-0.14469266210702247, closeTo(fourthRootsOfZ[0].imaginary, 1.0e-5));
    // test z_1
    Assert.assertEquals(0.14469266210702256, closeTo(fourthRootsOfZ[1].real, 1.0e-5));
    Assert.assertEquals(1.5164629308487783, closeTo(fourthRootsOfZ[1].imaginary, 1.0e-5));
    // test z_2
    Assert.assertEquals(-1.5164629308487783, closeTo(fourthRootsOfZ[2].real, 1.0e-5));
    Assert.assertEquals(0.14469266210702267, closeTo(fourthRootsOfZ[2].imaginary, 1.0e-5));
    // test z_3
    Assert.assertEquals(-0.14469266210702275, closeTo(fourthRootsOfZ[3].real, 1.0e-5));
    Assert.assertEquals(-1.5164629308487783, closeTo(fourthRootsOfZ[3].imaginary, 1.0e-5));
  });

  /**
     * Test: computing <b>third roots</b> of z.
     * <pre>
     * <code>
     * <b>z = 8</b>
     *   => z_0 =  2
     *   => z_1 = -1 + 1.73205 * i
     *   => z_2 = -1 - 1.73205 * i
     * </code>
     * </pre>
     */

  test('NthRoot_cornercase_thirdRoot_imaginaryPartEmpty', () {
    // The number 8 has three third roots. One we all already know is the number 2.
    // But there are two more complex roots.
    Complex z = new Complex(8, 0);
    // The List holding all third roots
    List<Complex> thirdRootsOfZ = z.nthRoot(3);//.toArray(new Complex[0]);
    // Returned Collection must not be empty!
    Assert.assertEquals(3, thirdRootsOfZ.length);
    // test z_0
    Assert.assertEquals(2.0, closeTo(thirdRootsOfZ[0].real, 1.0e-5));
    Assert.assertEquals(0.0, closeTo(thirdRootsOfZ[0].imaginary, 1.0e-5));
    // test z_1
    Assert.assertEquals(-1.0, closeTo(thirdRootsOfZ[1].real, 1.0e-5));
    Assert.assertEquals(1.7320508075688774, closeTo(thirdRootsOfZ[1].imaginary, 1.0e-5));
    // test z_2
    Assert.assertEquals(-1.0, closeTo(thirdRootsOfZ[2].real, 1.0e-5));
    Assert.assertEquals(-1.732050807568877, closeTo(thirdRootsOfZ[2].imaginary, 1.0e-5));
  });


  /**
     * Test: computing <b>third roots</b> of z with real part 0.
     * <pre>
     * <code>
     * <b>z = 2 * i</b>
     *   => z_0 =  1.0911 + 0.6299 * i
     *   => z_1 = -1.0911 + 0.6299 * i
     *   => z_2 = -2.3144 - 1.2599 * i
     * </code>
     * </pre>
     */

  test('NthRoot_cornercase_thirdRoot_realPartZero', () {
    // complex number with only imaginary part
    Complex z = new Complex(0, 2);
    // The List holding all third roots
    List<Complex> thirdRootsOfZ = z.nthRoot(3);//.toArray(new Complex[0]);
    // Returned Collection must not be empty!
    Assert.assertEquals(3, thirdRootsOfZ.length);
    // test z_0
    Assert.assertEquals(1.0911236359717216, closeTo(thirdRootsOfZ[0].real, 1.0e-5));
    Assert.assertEquals(0.6299605249474365, closeTo(thirdRootsOfZ[0].imaginary, 1.0e-5));
    // test z_1
    Assert.assertEquals(-1.0911236359717216, closeTo(thirdRootsOfZ[1].real, 1.0e-5));
    Assert.assertEquals(0.6299605249474365, closeTo(thirdRootsOfZ[1].imaginary, 1.0e-5));
    // test z_2
    Assert.assertEquals(-2.3144374213981936E-16, closeTo(thirdRootsOfZ[2].real, 1.0e-5));
    Assert.assertEquals(-1.2599210498948732, closeTo(thirdRootsOfZ[2].imaginary, 1.0e-5));
  });

  /**
     * Test cornercases with NaN and Infinity.
     */

  test('NthRoot_cornercase_NAN_Inf', () {
    // NaN + finite -> NaN
    List<Complex> roots = oneNaN.nthRoot(3);
    Assert.assertEquals(1, roots.size());
    Assert.assertEquals(Complex.NaN, roots.get(0));

    roots = nanZero.nthRoot(3);
    Assert.assertEquals(1, roots.size());
    Assert.assertEquals(Complex.NaN, roots.get(0));

    // NaN + infinite -> NaN
    roots = nanInf.nthRoot(3);
    Assert.assertEquals(1, roots.size());
    Assert.assertEquals(Complex.NaN, roots.get(0));

    // finite + infinite -> Inf
    roots = oneInf.nthRoot(3);
    Assert.assertEquals(1, roots.size());
    Assert.assertEquals(Complex.INF, roots.get(0));

    // infinite + infinite -> Inf
    roots = negInfInf.nthRoot(3);
    Assert.assertEquals(1, roots.size());
    Assert.assertEquals(Complex.INF, roots.get(0));
  });

  /**
     * Test standard values
     */

  test('GetArgument', () {
    Complex z = new Complex(1, 0);
    Assert.assertEquals(0.0, closeTo(z.argument(), 1.0e-12));

    z = new Complex(1, 1);
    Assert.assertEquals(FastMath.PI / 4, closeTo(z.argument(), 1.0e-12));

    z = new Complex(0, 1);
    Assert.assertEquals(FastMath.PI / 2, closeTo(z.argument(), 1.0e-12));

    z = new Complex(-1, 1);
    Assert.assertEquals(3 * FastMath.PI / 4, closeTo(z.argument(), 1.0e-12));

    z = new Complex(-1, 0);
    Assert.assertEquals(FastMath.PI, closeTo(z.argument(), 1.0e-12));

    z = new Complex(-1, -1);
    Assert.assertEquals(-3 * FastMath.PI / 4, closeTo(z.argument(), 1.0e-12));

    z = new Complex(0, -1);
    Assert.assertEquals(-FastMath.PI / 2, closeTo(z.argument(), 1.0e-12));

    z = new Complex(1, -1);
    Assert.assertEquals(-FastMath.PI / 4, closeTo(z.argument(), 1.0e-12));

  });

  /**
     * Verify atan2-style handling of infinite parts
     */

  test('GetArgumentInf', () {
    Assert.assertEquals(FastMath.PI / 4, closeTo(infInf.argument(), 1.0e-12));
    Assert.assertEquals(FastMath.PI / 2, closeTo(oneInf.argument(), 1.0e-12));
    Assert.assertEquals(0.0, closeTo(infOne.argument(), 1.0e-12));
    Assert.assertEquals(FastMath.PI / 2, closeTo(zeroInf.argument(), 1.0e-12));
    Assert.assertEquals(0.0, closeTo(infZero.argument(), 1.0e-12));
    Assert.assertEquals(FastMath.PI, closeTo(negInfOne.argument(), 1.0e-12));
    Assert.assertEquals(-3.0 * FastMath.PI / 4, closeTo(negInfNegInf.argument(), 1.0e-12));
    Assert.assertEquals(-FastMath.PI / 2, closeTo(oneNegInf.argument(), 1.0e-12));
  });

  /**
     * Verify that either part NaN results in NaN
     */

  test('GetArgumentNaN', () {
    Assert.assertTrue(Double.isNaN(nanZero.argument));
    Assert.assertTrue(Double.isNaN(zeroNaN.argument));
    Assert.assertTrue(Double.isNaN(Complex.NaN.argument));
  });


  test('Serial', () {
    Complex z = new Complex(3.0, 4.0);
    Assert.assertEquals(z, TestUtils.serializeAndRecover(z));
    Complex ncmplx = TestUtils.serializeAndRecover(oneNaN) as Complex;
    Assert.assertEquals(nanZero, ncmplx);
    Assert.assertTrue(ncmplx.isNaN());
    Complex infcmplx = TestUtils.serializeAndRecover(infInf) as Complex;
    Assert.assertEquals(infInf, infcmplx);
    Assert.assertTrue(infcmplx.isInfinite());
    TestComplex tz = new TestComplex(3.0, 4.0);
    Assert.assertEquals(tz, TestUtils.serializeAndRecover(tz));
    TestComplex ntcmplx = TestUtils.serializeAndRecover(new TestComplex(oneNaN)) as TestComplex;
    Assert.assertEquals(nanZero, ntcmplx);
    Assert.assertTrue(ntcmplx.isNaN());
    TestComplex inftcmplx = TestUtils.serializeAndRecover(new TestComplex(infInf)) as TestComplex;
    Assert.assertEquals(infInf, inftcmplx);
    Assert.assertTrue(inftcmplx.isInfinite());
  });

}

/**
 * Class to test extending Complex
 */
class TestComplex extends Complex {

  TestComplex(double real, double imaginary) : super(real, imaginary);

  factory TestComplex.from(Complex other) {
    new TestComplex(other.real, other.imaginary);
  }

  @override
  TestComplex createComplex(double real, double imaginary) {
    return new TestComplex(real, imaginary);
  }

}
