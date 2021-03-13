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

import 'complex_matcher.dart';

void main() {
  const inf = double.infinity;
  const neginf = double.negativeInfinity;
  const nan = double.nan;
  const pi = math.pi;
  const oneInf = Complex(1, inf);
  const oneNegInf = Complex(1, neginf);
  const infOne = Complex(inf, 1);
  const infZero = Complex(inf, 0);
  const infNaN = Complex(inf, nan);
  const infNegInf = Complex(inf, neginf);
  const infInf = Complex(inf, inf);
  const negInfInf = Complex(neginf, inf);
  const negInfZero = Complex(neginf, 0);
  const negInfOne = Complex(neginf, 1);
  const negInfNaN = Complex(neginf, nan);
  const negInfNegInf = Complex(neginf, neginf);
  const oneNaN = Complex(1, nan);
  const zeroInf = Complex(0, inf);
  const zeroNaN = Complex(0, nan);
  const nanInf = Complex(nan, inf);
  const nanNegInf = Complex(nan, neginf);
  const nanZero = Complex(nan, 0);

  test('Constructor', () {
    const z = Complex(3.0, 4.0);
    expect(3.0, closeTo(z.real, 1.0e-5));
    expect(4.0, closeTo(z.imaginary, 1.0e-5));
  });

  test('ConstructorNaN', () {
    const z1 = Complex(3.0, double.nan);
    expect(z1.isNaN, isTrue);

    const z2 = Complex(double.nan, 4.0);
    expect(z2.isNaN, isTrue);

    const z3 = Complex(3.0, 4.0);
    expect(z3.isNaN, isFalse);
  });

  test('Abs', () {
    const z = Complex(3.0, 4.0);
    expect(5.0, closeTo(z.abs(), 1.0e-5));
  });

  test('AbsNaN', () {
    expect(Complex.nan.abs().isNaN, isTrue);
    const z = Complex(inf, nan);
    expect(z.abs().isNaN, isTrue);
  });

  test('AbsInfinite', () {
    const z1 = Complex(inf, 0);
    expect(inf, equals(z1.abs()));

    const z2 = Complex(0, neginf);
    expect(inf, equals(z2.abs()));

    const z3 = Complex(inf, neginf);
    expect(inf, equals(z3.abs()));
  });

  test('Add', () {
    const x = Complex(3.0, 4.0);
    const y = Complex(5.0, 6.0);
    final z = x + y;
    expect(8.0, closeTo(z.real, 1.0e-5));
    expect(10.0, closeTo(z.imaginary, 1.0e-5));
  });

  test('AddNaN', () {
    const x = Complex(3.0, 4.0);
    final z1 = x + Complex.nan;
    expect(Complex.nan, equals(z1));

    const z2 = Complex(1, nan);
    final w = x + z2;
    expect(Complex.nan, equals(w));
  });

  test('AddInf', () {
    const x = Complex(1, 1);
    const z = Complex(inf, 0);
    final w = x + z;
    expect(w.imaginary, equals(1));
    expect(inf, equals(w.real));

    const x1 = Complex(neginf, 0);
    expect((x1 + z).real.isNaN, isTrue);
  });

  test('ScalarAdd', () {
    const x = Complex(3.0, 4.0);
    const yDouble = 2.0;
    const yComplex = Complex(yDouble);
    expect(x + yComplex, equals(x + yDouble));
  });

  test('ScalarAddNaN', () {
    const x = Complex(3.0, 4.0);
    const yDouble = double.nan;
    const yComplex = Complex(yDouble);
    expect(x + yComplex, equals(x + yDouble));
  });

  test('ScalarAddInf', () {
    const x = Complex(1, 1);
    const yDouble = double.infinity;

    const yComplex = Complex(yDouble);
    expect(x + yComplex, x + yDouble);

    {
      const x = Complex(neginf, 0);
      expect(x + yComplex, x + yDouble);
    }
  });

  test('Conjugate', () {
    const x = Complex(3.0, 4.0);
    final z = x.conjugate();
    expect(3.0, closeTo(z.real, 1.0e-5));
    expect(-4.0, closeTo(z.imaginary, 1.0e-5));
  });

  test('ConjugateNaN', () {
    final z = Complex.nan.conjugate();
    expect(z.isNaN, isTrue);
  });

  test('ConjugateInfinite', () {
    {
      const z = Complex(0, inf);
      expect(neginf, equals(z.conjugate().imaginary));
    }
    {
      const z = Complex(0, neginf);
      expect(inf, equals(z.conjugate().imaginary));
    }
  });

  test('Divide', () {
    const x = Complex(3.0, 4.0);
    const y = Complex(5.0, 6.0);
    final z = x / y;
    expect(39.0 / 61.0, closeTo(z.real, 1.0e-5));
    expect(2.0 / 61.0, closeTo(z.imaginary, 1.0e-5));
  });

  test('DivideReal', () {
    const x = Complex(2, 3);
    const y = Complex(2, 0);
    expect(const Complex(1, 1.5), equals(x / y));
  });

  test('DivideImaginary', () {
    const x = Complex(2, 3);
    const y = Complex(0, 2);
    expect(const Complex(1.5, -1), equals(x / y));
  });

  test('DivideInf', () {
    const x = Complex(3, 4);

    {
      const w = Complex(neginf, inf);
      expect(x / w == Complex.zero, isTrue);

      final z = w / x;
      expect(z.real.isNaN, isTrue);
      expect(inf, equals(z.imaginary));
    }

    {
      const w = Complex(inf, inf);
      final z = w / x;
      expect(z.imaginary.isNaN, isTrue);
      expect(inf, equals(z.real));
    }

    {
      const w = Complex(1, inf);
      final z = w / w;
      expect(z.real.isNaN, isTrue);
      expect(z.imaginary.isNaN, isTrue);
    }
  });

  test('DivideZero', () {
    const x = Complex(3.0, 4.0);
    final z = x / Complex.zero;
    // expect(z, Complex.INF); // See MATH-657
    expect(z, equals(Complex.nan));
  });

  test('DivideZeroZero', () {
    const x = Complex(0.0, 0.0);
    final z = x / Complex.zero;
    expect(z, equals(Complex.nan));
  });

  test('DivideNaN', () {
    const x = Complex(3.0, 4.0);
    final z = x / Complex.nan;
    expect(z.isNaN, isTrue);
  });

  test('DivideNaNInf', () {
    {
      final z = oneInf / Complex.one;
      expect(z.real.isNaN, isTrue);
      expect(inf, equals(z.imaginary));
    }

    {
      final z = negInfNegInf / oneNaN;
      expect(z.real.isNaN, isTrue);
      expect(z.imaginary.isNaN, isTrue);
    }

    {
      final z = negInfInf / Complex.one;
      expect(z.real.isNaN, isTrue);
      expect(z.imaginary.isNaN, isTrue);
    }
  });

  test('ScalarDivide', () {
    const x = Complex(3.0, 4.0);
    const yDouble = 2.0;
    const yComplex = Complex(yDouble);
    expect(x / yComplex, equals(x / yDouble));
  });

  test('ScalarDivideNaN', () {
    const x = Complex(3.0, 4.0);
    const yDouble = double.nan;
    const yComplex = Complex(yDouble);
    expect(x / yComplex, equals(x / yDouble));
  });

  test('ScalarDivideInf', () {
    const x = Complex(1, 1);

    {
      const yDouble = double.infinity;
      const yComplex = Complex(yDouble);
      expect(x / yComplex, equals(x / yDouble));
    }

    {
      const yDouble = double.negativeInfinity;
      const yComplex = Complex(yDouble);
      expect(x / yComplex, equals(x / yDouble));
    }

    {
      const x = Complex(1, double.negativeInfinity);
      const yDouble = double.negativeInfinity;
      const yComplex = Complex(yDouble);
      expect(x / yComplex, equals(x / yDouble));
    }
  });

  test('ScalarDivideZero', () {
    const x = Complex(1, 1);
    expect(x / Complex.zero, equals(x / 0));
  });

  test('Reciprocal', () {
    const z = Complex(5.0, 6.0);
    final act = z.reciprocal();
    const expRe = 5.0 / 61.0;
    const expIm = -6.0 / 61.0;
    expect(expRe, closeTo(act.real, /*FastMath.ulp(expRe)*/ 1.0e-12));
    expect(expIm, closeTo(act.imaginary, /*FastMath.ulp(expIm)*/ 1.0e-12));
  });

  test('ReciprocalReal', () {
    const z = Complex(-2.0, 0.0);
    //expect(Complex.equals(Complex(-0.5, 0.0), z.reciprocal()), isTrue);
    expect(const Complex(-0.5, 0.0), equals(z.reciprocal()));
  });

  test('ReciprocalImaginary', () {
    const z = Complex(0.0, -2.0);
    expect(const Complex(0.0, 0.5), equals(z.reciprocal()));
  });

  test('ReciprocalInf', () {
    {
      const z = Complex(neginf, inf);
      expect(z.reciprocal() == Complex.zero, isTrue);
    }

    {
      final z = const Complex(1, inf).reciprocal();
      expect(z, equals(Complex.zero));
    }
  });

  test('ReciprocalZero', () {
    expect(Complex.zero.reciprocal(), equals(Complex.infinity));
  });

  test('ReciprocalNaN', () {
    expect(Complex.nan.reciprocal().isNaN, isTrue);
  });

  test('Multiply', () {
    const x = Complex(3.0, 4.0);
    const y = Complex(5.0, 6.0);
    final z = x * y;
    expect(-9.0, closeTo(z.real, 1.0e-5));
    expect(38.0, closeTo(z.imaginary, 1.0e-5));
  });

  test('MultiplyNaN', () {
    const x = Complex(3.0, 4.0);
    {
      final z = x * Complex.nan;
      expect(Complex.nan, equals(z));
    }
    {
      final z = Complex.nan * 5;
      expect(Complex.nan, equals(z));
    }
  });

  test('MultiplyInfInf', () {
    // expect(infInf.multiply(infInf).isNaN()); // MATH-620
    expect((infInf * infInf).isInfinite, isTrue);
  });

  test('MultiplyNaNInf', () {
    {
      const z = Complex(1, 1);
      final w = z * infOne;
      expect(w.real, equals(inf));
      expect(w.imaginary, equals(inf));
    }

    // [MATH-164]
    expect(const Complex(1, 0) * infInf == Complex.infinity, isTrue);
    expect(const Complex(-1, 0) * infInf == Complex.infinity, isTrue);
    expect(const Complex(1, 0) * negInfZero == Complex.infinity, isTrue);

    {
      final w = oneInf * oneNegInf;
      expect(w.real, equals(inf));
      expect(w.imaginary, equals(inf));
    }

    {
      final w = negInfNegInf * oneNaN;
      expect(w.real.isNaN, isTrue);
      expect(w.imaginary.isNaN, isTrue);
    }

    {
      const z = Complex(1, neginf);
      expect(Complex.infinity, equals(z * z));
    }
  });

  test('ScalarMultiply', () {
    const x = Complex(3.0, 4.0);
    const yDouble = 2.0;
    const yComplex = Complex(yDouble);
    expect(x * yComplex, equals(x * yDouble));
    const zInt = -5;
    final zComplex = Complex(zInt.toDouble());
    expect(x * zComplex, equals(x * zInt));
  });

  test('ScalarMultiplyNaN', () {
    const x = Complex(3.0, 4.0);
    const yDouble = double.nan;
    const yComplex = Complex(yDouble);
    expect(x * yComplex, equals(x * yDouble));
  });

  test('ScalarMultiplyInf', () {
    const x = Complex(1, 1);
    {
      const yDouble = double.infinity;
      const yComplex = Complex(yDouble);
      expect(x * yComplex, equals(x * yDouble));
    }

    {
      const yDouble = double.negativeInfinity;
      const yComplex = Complex(yDouble);
      expect(x * yComplex, equals(x * yDouble));
    }
  });

  test('Negate', () {
    const x = Complex(3.0, 4.0);
    final z = -x;
    expect(-3.0, closeTo(z.real, 1.0e-5));
    expect(-4.0, closeTo(z.imaginary, 1.0e-5));
  });

  test('NegateNaN', () {
    Complex z = -Complex.nan;
    expect(z.isNaN, isTrue);
  });

  test('Subtract', () {
    const x = Complex(3.0, 4.0);
    const y = Complex(5.0, 6.0);
    final z = x - y;
    expect(-2.0, closeTo(z.real, 1.0e-5));
    expect(-2.0, closeTo(z.imaginary, 1.0e-5));
  });

  test('SubtractNaN', () {
    const x = Complex(3.0, 4.0);
    {
      final z = x - Complex.nan;
      expect(Complex.nan, equals(z));
    }

    {
      const z = Complex(1, nan);
      final w = x - z;
      expect(Complex.nan, equals(w));
    }
  });

  test('SubtractInf', () {
    const z = Complex(neginf, 0);
    {
      const x = Complex(1, 1);

      final w = x - z;
      expect(w.imaginary, equals(1));
      expect(inf, equals(w.real));
    }

    {
      const x = Complex(neginf, 0);
      expect((x - z).real.isNaN, isTrue);
    }
  });

  test('ScalarSubtract', () {
    const x = Complex(3.0, 4.0);
    const yDouble = 2.0;
    const yComplex = Complex(yDouble);
    expect(x - yComplex, equals(x - yDouble));
  });

  test('ScalarSubtractNaN', () {
    const x = Complex(3.0, 4.0);
    const yDouble = double.nan;
    const yComplex = Complex(yDouble);
    expect(x - yComplex, equals(x - yDouble));
  });

  test('ScalarSubtractInf', () {
    const yDouble = double.infinity;
    const yComplex = Complex(yDouble);
    {
      const x = Complex(1, 1);

      expect(x - yComplex, equals(x - yDouble));
    }

    {
      const x = Complex(neginf, 0);
      expect(x - yComplex, equals(x - yDouble));
    }
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
    const x = Complex(3.0, 4.0);
    expect(x is TestComplex, isFalse);
  });

  test('EqualsSame', () {
    const x = Complex(3.0, 4.0);
    expect(x == x, isTrue);
  });

  test('FloatingPointEquals', () {
    const re = -3.21;
    const im = 456789e10;

    const x = Complex(re, im);
    const y = Complex(re, im);

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
    const x = Complex(3.0, 4.0);
    const y = Complex(3.0, 4.0);
    expect(x == y, isTrue);
  });

  test('EqualsRealDifference', () {
    const x = Complex(0.0, 0.0);
    const y = Complex(0.0 + double.minPositive, 0.0);
    expect(x == y, isFalse);
  });

  test('EqualsImaginaryDifference', () {
    const x = Complex(0.0, 0.0);
    const y = Complex(0.0, 0.0 + double.minPositive);
    expect(x == y, isFalse);
  });

  test('EqualsNaN', () {
    const realNaN = Complex(double.nan, 0.0);
    const imaginaryNaN = Complex(0.0, double.nan);
    Complex complexNaN = Complex.nan;
    expect(realNaN == imaginaryNaN, isTrue);
    expect(imaginaryNaN == complexNaN, isTrue);
    expect(realNaN == complexNaN, isTrue);
  });

  test('HashCode', () {
    const x = Complex(0.0, 0.0);
    {
      const y = Complex(0.0, 0.0 + double.minPositive);
      expect(x.hashCode == y.hashCode, isFalse);
    }
    {
      {
        const y = Complex(0.0 + double.minPositive, 0.0);
        expect(x.hashCode == y.hashCode, isFalse);
      }
      const realNaN = Complex(double.nan, 0.0);
      const imaginaryNaN = Complex(0.0, double.nan);
      expect(realNaN.hashCode, equals(imaginaryNaN.hashCode));
      expect(imaginaryNaN.hashCode, equals(Complex.nan.hashCode));
    }

    // TODO: review this section
    // MATH-1118
    // "equals" and "hashCode" must be compatible: if two objects have
    // different hash codes, "equals" must return false.
    //final msg = "'==' not compatible with 'hashCode'";

    //{
    //  const x = Complex(0.0, 0.0);
    //  const y = Complex(0.0, -0.0);
    //  expect(x.hashCode != y.hashCode, isTrue);
    //  expect(x == y, isFalse, reason: msg);
    //}
    //{
    //  const x = Complex(0.0, 0.0);
    //  const y = Complex(-0.0, 0.0);
    //  expect(x.hashCode != y.hashCode, isTrue);
    //  expect(x == y, isFalse, reason:msg);
    //}
  });

  test('Acos', () {
    const z = Complex(3, 4);
    const expected = Complex(0.936812, -2.30551);
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
    const z = Complex(3, 4);
    const expected = Complex(0.633984, 2.30551);
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
    const z = Complex(3, 4);
    const expected = Complex(1.44831, 0.158997);
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
    const z = Complex(3, 4);
    const expected = Complex(-27.03495, -3.851153);
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
    const z = Complex(3, 4);
    const expected = Complex(-6.58066, -7.58155);
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
    const z = Complex(3, 4);
    const expected = Complex(-13.12878, -15.20078);
    expect(expected, closeToZ(z.exp(), 1.0e-5));
    expect(Complex.one, closeToZ(Complex.zero.exp(), 10e-12));
    final iPi = Complex.i * const Complex(pi, 0);
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
    const z = Complex(3, 4);
    const expected = Complex(1.60944, 0.927295);
    expect(expected, closeToZ(z.log(), 1.0e-5));
  });

  test('LogNaN', () {
    expect(Complex.nan.log().isNaN, isTrue);
  });

  test('LogInf', () {
    expect(const Complex(inf, pi / 2), closeToZ(oneInf.log(), 10e-12));
    expect(const Complex(inf, -pi / 2), closeToZ(oneNegInf.log(), 10e-12));
    expect(infZero, closeToZ(infOne.log(), 10e-12));
    expect(const Complex(inf, pi), closeToZ(negInfOne.log(), 10e-12));
    expect(const Complex(inf, pi / 4), closeToZ(infInf.log(), 10e-12));
    expect(const Complex(inf, -pi / 4), closeToZ(infNegInf.log(), 10e-12));
    expect(const Complex(inf, 3 * pi / 4), closeToZ(negInfInf.log(), 10e-12));
    expect(
      const Complex(inf, -3 * pi / 4),
      closeToZ(negInfNegInf.log(), 10e-12),
    );
  });

  test('LogZero', () {
    expect(negInfZero, equals(Complex.zero.log()));
  });

  test('Pow', () {
    const x = Complex(3, 4);
    const y = Complex(5, 6);
    const expected = Complex(-1.860893, 11.83677);
    expect(expected, closeToZ(x.power(y), 1.0e-5));
  });

  test('PowNaNBase', () {
    const x = Complex(3, 4);
    expect(Complex.nan.power(x).isNaN, isTrue);
  });

  test('PowNaNExponent', () {
    const x = Complex(3, 4);
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
    expect(Complex.one, closeToZ((3.im - 1).power(Complex.zero), 10e-12));
  });

  test('ScalarPow', () {
    const x = Complex(3, 4);
    const yDouble = 5.0;
    const yComplex = Complex(yDouble);
    expect(x.power(yComplex), equals(x.pow(yDouble)));
  });

  test('ScalarPowNaNBase', () {
    const x = Complex.nan;
    const yDouble = 5.0;
    const yComplex = Complex(yDouble);
    expect(x.power(yComplex), equals(x.pow(yDouble)));
  });

  test('ScalarPowNaNExponent', () {
    const x = Complex(3, 4);
    const yDouble = double.nan;
    const yComplex = Complex(yDouble);
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
    expect(Complex.one, closeToZ(const Complex(-1, 3).pow(0.0), 10e-12));
  });

  test('Sin', () {
    const z = Complex(3, 4);
    const expected = Complex(3.853738, -27.01681);
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
    const z = Complex(3, 4);
    const expected = Complex(-6.54812, -7.61923);
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
    const z = Complex(3, 4);
    const expected = Complex(2, 1);
    expect(expected, closeToZ(z.sqrt(), 1.0e-5));
  });

  test('SqrtRealZero', () {
    const z = Complex(0.0, 4);
    const expected = Complex(1.41421, 1.41421);
    expect(expected, closeToZ(z.sqrt(), 1.0e-5));
  });

  test('SqrtRealNegative', () {
    const z = Complex(-3.0, 4);
    const expected = Complex(1, 2);
    expect(expected, closeToZ(z.sqrt(), 1.0e-5));
  });

  test('SqrtImaginaryZero', () {
    const z = Complex(-3.0, 0.0);
    final expected = Complex(0.0, math.sqrt(3));
    expect(z.sqrt(), expected);
  });

  test('SqrtImaginaryNegative', () {
    const z = Complex(-3.0, -4.0);
    const expected = Complex(1.0, -2.0);
    expect(expected, closeToZ(z.sqrt(), 1.0e-5));
  });

  test('SqrtPolar', () {
    var r = 1.0;
    for (var i = 0; i < 5; i++) {
      r += i;
      var theta = 0.0;
      for (var j = 0; j < 11; j++) {
        theta += pi / 12;
        final z = Complex.polar(r, theta);
        final sqrtz = Complex.polar(math.sqrt(r), theta / 2);
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
    const z = Complex(3, 4);
    const expected = Complex(4.08033, -2.94094);
    expect(expected, closeToZ(z.sqrt1z(), 1.0e-5));
  });

  test('Sqrt1zNaN', () {
    expect(Complex.nan.sqrt1z().isNaN, isTrue);
  });

  test('Tan', () {
    const z = Complex(3, 4);
    const expected = Complex(-0.000187346, 0.999356);
    expect(expected, closeToZ(z.tan(), 1.0e-5));
    /* Check that no overflow occurs (MATH-722) */

    {
      final actual = const Complex(3.0, 1E10).tan();
      const expected = Complex(0, 1);
      expect(expected, closeToZ(actual, 1.0e-5));
    }
    {
      final actual = const Complex(3.0, -1E10).tan();
      const expected = Complex(0, -1);
      expect(expected, closeToZ(actual, 1.0e-5));
    }
  });

  test('TanNaN', () {
    expect(Complex.nan.tan().isNaN, isTrue);
  });

  test('TanInf', () {
    expect(const Complex(0.0, 1.0), equals(oneInf.tan()));
    expect(const Complex(0.0, -1.0), equals(oneNegInf.tan()));
    expect(Complex.nan, equals(infOne.tan()));
    expect(Complex.nan, equals(negInfOne.tan()));
    expect(Complex.nan, equals(infInf.tan()));
    expect(Complex.nan, equals(infNegInf.tan()));
    expect(Complex.nan, equals(negInfInf.tan()));
    expect(Complex.nan, equals(negInfNegInf.tan()));
  });

  test('TanCritical', () {
    expect(infNaN, equals(const Complex(pi / 2, 0).tan()));
    expect(negInfNaN, equals(const Complex(-pi / 2, 0).tan()));
  });

  test('Tanh', () {
    const z = Complex(3, 4);
    const expected = Complex(1.00071, 0.00490826);
    expect(expected, closeToZ(z.tanh(), 1.0e-5));
    /* Check that no overflow occurs (MATH-722) */
    final actual = const Complex(1E10, 3.0).tanh();
    {
      const expected = Complex(1, 0);
      expect(expected, closeToZ(actual, 1.0e-5));
    }
    {
      final actual = const Complex(-1E10, 3.0).tanh();
      const expected = Complex(-1, 0);
      expect(expected, closeToZ(actual, 1.0e-5));
    }
  });

  test('TanhNaN', () {
    expect(Complex.nan.tanh().isNaN, isTrue);
  });

  test('TanhInf', () {
    expect(Complex.nan, equals(oneInf.tanh()));
    expect(Complex.nan, equals(oneNegInf.tanh()));
    expect(const Complex(1.0, 0.0), equals(infOne.tanh()));
    expect(const Complex(-1.0, 0.0), equals(negInfOne.tanh()));
    expect(Complex.nan, equals(infInf.tanh()));
    expect(Complex.nan, equals(infNegInf.tanh()));
    expect(Complex.nan, equals(negInfInf.tanh()));
    expect(Complex.nan, equals(negInfNegInf.tanh()));
  });

  test('TanhCritical', () {
    expect(nanInf, equals(const Complex(0, pi / 2).tanh()));
  });

  /// test issue MATH-221

  test('Math221', () {
    expect(const Complex(0, -1) == Complex.i * const Complex(-1, 0), isTrue);
  });

  /// Test: computing <b>third roots</b> of z.
  ///
  ///     z = -2 + 2 * i
  ///      => z_0 =  1      +          i
  ///      => z_1 = -1.3660 + 0.3660 * i
  ///      => z_2 =  0.3660 - 1.3660 * i

  test('NthRoot_normal_thirdRoot', () {
    // The complex number we want to compute all third-roots for.
    const z = Complex(-2, 2);
    // The List holding all third roots
    final thirdRootsOfZ = z.nthRoot(3); //.toArray(Complex[0]);
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
    const z = Complex(5, -2);
    // The List holding all fourth roots
    final fourthRootsOfZ = z.nthRoot(4); //.toArray(Complex[0]);
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
    // The number 8 has three third roots.
    // One we all already know is the number 2.
    // But there are two more complex roots.
    const z = Complex(8, 0);
    // The List holding all third roots
    final thirdRootsOfZ = z.nthRoot(3); //.toArray(Complex[0]);
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
    const z = Complex(0, 2);
    // The List holding all third roots
    final thirdRootsOfZ = z.nthRoot(3); //.toArray(Complex[0]);
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
    {
      final roots = oneNaN.nthRoot(3);
      expect(1, equals(roots.length));
      expect(Complex.nan, equals(roots[0]));
    }

    {
      final roots = nanZero.nthRoot(3);
      expect(1, equals(roots.length));
      expect(Complex.nan, equals(roots[0]));
    }

    // NaN + infinite -> NaN
    {
      final roots = nanInf.nthRoot(3);
      expect(1, equals(roots.length));
      expect(Complex.nan, equals(roots[0]));
    }

    // finite + infinite -> Inf
    {
      final roots = oneInf.nthRoot(3);
      expect(1, equals(roots.length));
      expect(Complex.infinity, equals(roots[0]));
    }

    // infinite + infinite -> Inf
    {
      final roots = negInfInf.nthRoot(3);
      expect(1, equals(roots.length));
      expect(Complex.infinity, equals(roots[0]));
    }
  });

  /// Test standard values

  test('argument', () {
    var z = const Complex(1, 0);
    expect(0.0, closeTo(z.argument(), 1.0e-12));

    z = const Complex(1, 1);
    expect(math.pi / 4, closeTo(z.argument(), 1.0e-12));

    z = const Complex(0, 1);
    expect(math.pi / 2, closeTo(z.argument(), 1.0e-12));

    z = const Complex(-1, 1);
    expect(3 * math.pi / 4, closeTo(z.argument(), 1.0e-12));

    z = const Complex(-1, 0);
    expect(math.pi, closeTo(z.argument(), 1.0e-12));

    z = const Complex(-1, -1);
    expect(-3 * math.pi / 4, closeTo(z.argument(), 1.0e-12));

    z = const Complex(0, -1);
    expect(-math.pi / 2, closeTo(z.argument(), 1.0e-12));

    z = const Complex(1, -1);
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
    final ntcmplx = TestUtils.serializeAndRecover(TestComplex.from(oneNaN)) as TestComplex;
    expect(nanZero, equals(ntcmplx));
    expect(ntcmplx.isNaN, isTrue);
    final inftcmplx = TestUtils.serializeAndRecover(TestComplex.from(infInf)) as TestComplex;
    expect(infInf, inftcmplx);
    expect(inftcmplx.isInfinite, isTrue);
  });*/
}
