// Copyright (c) 2005-2015 The Scalable Software Infrastructure Project.
// All rights reserved.
library complex.list;

import 'dart:typed_data' show Float64List;
import 'dart:collection' show ListBase;
import 'package:quiver/iterables.dart';

import 'complex.dart';

class ComplexList extends ListBase<Complex> {
  Float64List _value;

  ComplexList(int length) : _value = new Float64List(2 * length);

  factory ComplexList.from(Iterable<Complex> list) {
    var cl = new ComplexList(list.length);
    enumerate(list).forEach((iv) {
      cl[iv.index] = iv.value;
    });
    return cl;
  }

  factory ComplexList.polar(Iterable<double> r, Iterable<double> theta,
      [bool radians = true]) {
    return new ComplexList.from(zip([r, theta])
        .map((parts) => new Complex.polar(parts[0], parts[1], radians)));
  }

  void set length(int sz) {
    throw new UnsupportedError('fixed length');
  }

  int get length => _value.length ~/ 2;

  double re(int i) => _value[2 * i];

  double im(int i) => _value[2 * i + 1];

  setReal(int i, double v) => _value[2 * i] = v;
  setImag(int i, double v) => _value[2 * i + 1] = v;

  ComplexList operator *(ComplexList y) => copy()..pmul(y);

  ComplexList operator /(ComplexList y) => copy()..pdiv(y);

  ComplexList operator +(ComplexList y) {
    var z = new ComplexList(length);
    for (int i = 0; i < length; i++) {
      z.setReal(i, re(i) + y.re(i));
      z.setImag(i, im(i) + y.im(i));
    }
    return z;
  }

  ComplexList operator -(ComplexList y) {
    var z = new ComplexList(length);
    for (int i = 0; i < length; i++) {
      z.setReal(i, re(i) - y.re(i));
      z.setImag(i, im(i) - y.im(i));
    }
    return z;
  }

  ComplexList operator -() => copy()..scale(-Complex.ONE);

  Float64List real() {
    var l = new Float64List(length);
    for (int i = 0; i < length; i++) {
      l[i] = re(i);
    }
    return l;
  }

  Float64List imag() {
    var l = new Float64List(length);
    for (int i = 0; i < length; i++) {
      l[i] = im(i);
    }
    return l;
  }

  ComplexList conj() => new ComplexList.from(map((c) => c.conjugate()));

  Complex operator [](int i) {
//    print('warning: complex list array access');
    return new Complex(_value[2 * i], _value[2 * i + 1]);
  }

  void operator []=(int i, Complex value) {
    _value[2 * i] = value.real;
    _value[2 * i + 1] = value.imaginary;
  }

  ///     x <-> y
  static void swap(ComplexList x, ComplexList y) {
    for (var i = 0; i < x._value.length; i++) {
      var t = y._value[i];
      y._value[i] = x._value[i];
      x._value[i] = t;
    }
  }

  ///     y <- x
  ComplexList copy() {
    var y = new ComplexList(length);
    for (var i = 0; i < _value.length; i++) {
      y._value[i] = _value[i];
    }
    return y;
  }

  ///     y <- alpha * x + y
  static void axpy(Complex alpha, ComplexList x, ComplexList y) {
    for (var i = 0; i < x.length; i++) {
      var re = x.re(i) * alpha.real - x.im(i) * alpha.imaginary;
      var im = x.re(i) * alpha.imaginary + x.im(i) * alpha.real;

      y.setReal(i, re + y.re(i));
      y.setImag(i, im + y.im(i));
    }
  }

  ///     y <- x + alpha * y
  static void xpay(ComplexList x, Complex alpha, ComplexList y) {
    for (var i = 0; i < x.length; i++) {
      var re = y.re(i) * alpha.real - y.im(i) * alpha.imaginary;
      var im = y.re(i) * alpha.imaginary + y.im(i) * alpha.real;

      y.setReal(i, x.re(i) + re);
      y.setImag(i, x.im(i) + im);
    }
  }

  ///     z <- y + alpha * x
  static void axpyz(
      Complex alpha, ComplexList x, ComplexList y, ComplexList z) {
    for (var i = 0; i < x.length; i++) {
      z[i] = alpha * x[i] + y[i];
    }
  }

  ///     y <- alpha * x
  void scale(Complex alpha) {
    for (var i = 0; i < length; i++) {
      this[i] = alpha * this[i];
    }
  }

  ///     z_i <- x_i * y_i
  ComplexList pmul(ComplexList y) {
    ComplexList z = new ComplexList(length);
    for (var i = 0; i < length; i++) {
      z[i] = this[i] * y[i];
    }
    return z;
  }

  ///     z_i <- x_i / y_i
  ComplexList pdiv(ComplexList y) {
    ComplexList z = new ComplexList(length);
    for (var i = 0; i < length; i++) {
      z[i] = this[i] / y[i];
    }
    return z;
  }

  ///     x_i <- alpha
  void fill(Complex alpha) {
    for (var i = 0; i < length; i++) {
      this[i] = alpha;
    }
  }

  ///     x_i <- |x_i|
  ComplexList abs() {
    ComplexList y = new ComplexList(length);
    for (var i = 0; i < length; i++) {
      y[i] = new Complex(this[i].abs());
    }
    return y;
  }

  ComplexList argument() {
    ComplexList y = new ComplexList(length);
    for (var i = 0; i < length; i++) {
      y[i] = new Complex(this[i].argument());
    }
    return y;
  }

  ///     x_i <- 1 / x_i
  ComplexList reciprocal() {
    ComplexList y = new ComplexList(length);
    for (var i = 0; i < length; i++) {
      y[i] = this[i].reciprocal();
    }
    return y;
  }

  ///     x_i <- alpha + x_i
  void shift(Complex t) {
    for (var i = 0; i < length; i++) {
      this[i] = this[i] + t;
    }
  }

  ///     v <- x^T * y
  static Complex dot(ComplexList x, ComplexList y) {
    var value = 0.0;
    for (var i = 0; i < x.length; i++) {
      value = x[i] * y[i] + value;
    }
    return value;
  }

  ///     v <- ||x||_1
  double nrm1() {
    var t = 0.0;
    for (var i = 0; i < length; i++) {
      t += _sabs(this[i]);
    }
    return t;
  }

  static Complex _ssqrt(Complex s) {
    if (s is num) {
      return math.sqrt(s);
    } else if (s is Complex) {
      return s.sqrt();
    } else {
      return s.sqrt();
    }
  }

  ///     v <- ||x||_2
  double nrm2() {
    var t = 0.0;
    for (var i = 0; i < length; i++) {
      t += this[i] * this[i];
    }
    return _ssqrt(t);
  }

  ///     v <- ||x||_infinity
  double nrmi() {
    var t = 0.0;
    for (var i = 0; i < length; i++) {
      if (t < _sabs(this[i])) {
        t = _sabs(this[i]);
      }
    }
    return t;
  }

  ///     v <- sum x_i
  static Complex sum(ComplexList x) {
    var t = 0.0;
    for (var i = 0; i < x.length; i++) {
      t += x[i];
    }
    return t;
  }

  ///     c <- A * b
  static void matvec(ComplexList a, ComplexList x, ComplexList y, int op) {
    int n = x.length;

    /* y = A*x */

    if (op == INS_VALUE) {
      switch (n) {
        case 1:
          y[0] = a[0] * x[0];
          break;
        case 2:
          y[0] = a[0] * x[0] + a[2] * x[1];
          y[1] = a[1] * x[0] + a[3] * x[1];
          break;
        case 3:
          y[0] = a[0] * x[0] + a[3] * x[1] + a[6] * x[2];
          y[1] = a[1] * x[0] + a[4] * x[1] + a[7] * x[2];
          y[2] = a[2] * x[0] + a[5] * x[1] + a[8] * x[2];
          break;
        default:
          for (var i = 0; i < n; i++) {
            var t = 0.0;
            for (var j = 0; j < n; j++) {
              t += a[i + j * n] * x[j];
            }
            y[i] = t;
          }
          break;
      }
    } else if (op == SUB_VALUE) {
      switch (n) {
        case 1:
          y[0] -= a[0] * x[0];
          break;
        case 2:
          y[0] -= a[0] * x[0] + a[2] * x[1];
          y[1] -= a[1] * x[0] + a[3] * x[1];
          break;
        case 3:
          y[0] -= a[0] * x[0] + a[3] * x[1] + a[6] * x[2];
          y[1] -= a[1] * x[0] + a[4] * x[1] + a[7] * x[2];
          y[2] -= a[2] * x[0] + a[5] * x[1] + a[8] * x[2];
          break;
        default:
          for (var i = 0; i < n; i++) {
            var t = 0.0;
            for (var j = 0; j < n; j++) {
              t += a[i + j * n] * x[j];
            }
            y[i] -= t;
          }
          break;
      }
    } else {
      switch (n) {
        case 1:
          y[0] += a[0] * x[0];
          break;
        case 2:
          y[0] += a[0] * x[0] + a[2] * x[1];
          y[1] += a[1] * x[0] + a[3] * x[1];
          break;
        case 3:
          y[0] += a[0] * x[0] + a[3] * x[1] + a[6] * x[2];
          y[1] += a[1] * x[0] + a[4] * x[1] + a[7] * x[2];
          y[2] += a[2] * x[0] + a[5] * x[1] + a[8] * x[2];
          break;
        default:
          for (var i = 0; i < n; i++) {
            var t = 0.0;
            for (var j = 0; j < n; j++) {
              t += a[i + j * n] * x[j];
            }
            y[i] += t;
          }
          break;
      }
    }
  }

  ///     c <- A^T * b
  static void matvect(ComplexList a, ComplexList x, ComplexList y, int op) {
    int n = x.length;

    /* y = A*x */

    if (op == INS_VALUE) {
      switch (n) {
        case 1:
          y[0] = a[0] * x[0];
          break;
        case 2:
          y[0] = a[0] * x[0] + a[1] * x[1];
          y[1] = a[2] * x[0] + a[3] * x[1];
          break;
        case 3:
          y[0] = a[0] * x[0] + a[1] * x[1] + a[2] * x[2];
          y[1] = a[3] * x[0] + a[4] * x[1] + a[5] * x[2];
          y[2] = a[6] * x[0] + a[7] * x[1] + a[8] * x[2];
          break;
        default:
          for (var i = 0; i < n; i++) {
            t = 0.0;
            for (var j = 0; j < n; j++) {
              t += a[i * n + j] * x[j];
            }
            y[i] = t;
          }
          break;
      }
    } else if (op == SUB_VALUE) {
      switch (n) {
        case 1:
          y[0] -= a[0] * x[0];
          break;
        case 2:
          y[0] -= a[0] * x[0] + a[1] * x[1];
          y[1] -= a[2] * x[0] + a[3] * x[1];
          break;
        case 3:
          y[0] -= a[0] * x[0] + a[1] * x[1] + a[2] * x[2];
          y[1] -= a[3] * x[0] + a[4] * x[1] + a[5] * x[2];
          y[2] -= a[6] * x[0] + a[7] * x[1] + a[8] * x[2];
          break;
        default:
          for (var i = 0; i < n; i++) {
            var t = 0.0;
            for (var j = 0; j < n; j++) {
              t += a[i * n + j] * x[j];
            }
            y[i] -= t;
          }
          break;
      }
    } else {
      switch (n) {
        case 1:
          y[0] += a[0] * x[0];
          break;
        case 2:
          y[0] += a[0] * x[0] + a[1] * x[1];
          y[1] += a[2] * x[0] + a[3] * x[1];
          break;
        case 3:
          y[0] += a[0] * x[0] + a[1] * x[1] + a[2] * x[2];
          y[1] += a[3] * x[0] + a[4] * x[1] + a[5] * x[2];
          y[2] += a[6] * x[0] + a[7] * x[1] + a[8] * x[2];
          break;
        default:
          for (var i = 0; i < n; i++) {
            var t = 0.0;
            for (var j = 0; j < n; j++) {
              t += a[i * n + j] * x[j];
            }
            y[i] += t;
          }
          break;
      }
    }
  }

  ///     c <- A * b where A is not square
  static void matvec_ns(int m, int n, ComplexList a, int lda, ComplexList x,
      ComplexList y, int op) {
    int n = x.length;

    /* y = A*x */

    if (op == INS_VALUE) {
      for (var i = 0; i < m; i++) {
        var t = 0.0;
        for (var j = 0; j < n; j++) {
          t += a[i + j * lda] * x[j];
        }
        y[i] = t;
      }
    } else if (op == SUB_VALUE) {
      for (var i = 0; i < m; i++) {
        var t = 0.0;
        for (var j = 0; j < n; j++) {
          t += a[i + j * lda] * x[j];
        }
        y[i] -= t;
      }
    } else if (op == ADD_VALUE) {
      for (var i = 0; i < m; i++) {
        var t = 0.0;
        for (var j = 0; j < n; j++) {
          t += a[i + j * lda] * x[j];
        }
        y[i] += t;
      }
    } else {
      switch (n) {
        case 1:
          y[0] += a[0] * x[0];
          break;
        case 2:
          y[0] += a[0] * x[0] + a[2] * x[1];
          y[1] += a[1] * x[0] + a[3] * x[1];
          break;
        case 3:
          y[0] += a[0] * x[0] + a[3] * x[1] + a[6] * x[2];
          y[1] += a[1] * x[0] + a[4] * x[1] + a[7] * x[2];
          y[2] += a[2] * x[0] + a[5] * x[1] + a[8] * x[2];
          break;
        default:
          for (var i = 0; i < n; i++) {
            var t = 0.0;
            for (var j = 0; j < n; j++) {
              t += a[i + j * n] * x[j];
            }
            y[i] += t;
          }
          break;
      }
    }
  }

  ///     C <- A * B
  static void matmat(ComplexList a, ComplexList b, ComplexList c, int op) {
    int n = x.length;

    /* C = A*B */

    if (op == INS_VALUE) {
      switch (n) {
        case 1:
          c[0] = a[0] * b[0];
          break;
        case 2:
          c[0] = a[0] * b[0] + a[2] * b[1];
          c[1] = a[1] * b[0] + a[3] * b[1];
          c[2] = a[0] * b[2] + a[2] * b[3];
          c[3] = a[1] * b[2] + a[3] * b[3];
          break;
        case 3:
          c[0] = a[0] * b[0] + a[3] * b[1] + a[6] * b[2];
          c[1] = a[1] * b[0] + a[4] * b[1] + a[7] * b[2];
          c[2] = a[2] * b[0] + a[5] * b[1] + a[8] * b[2];
          c[3] = a[0] * b[3] + a[3] * b[4] + a[6] * b[5];
          c[4] = a[1] * b[3] + a[4] * b[4] + a[7] * b[5];
          c[5] = a[2] * b[3] + a[5] * b[4] + a[8] * b[5];
          c[6] = a[0] * b[6] + a[3] * b[7] + a[6] * b[8];
          c[7] = a[1] * b[6] + a[4] * b[7] + a[7] * b[8];
          c[8] = a[2] * b[6] + a[5] * b[7] + a[8] * b[8];
          break;
        default:
          for (var j = 0; j < n; j++) {
            for (var i = 0; i < n; i++) {
              c[i + j * n] = 0.0;
            }
            for (var l = 0; l < n; l++) {
              for (var i = 0; i < n; i++) {
                c[i + j * n] += a[i + l * n] * b[l + j * n];
              }
            }
          }
          break;
      }
    } else if (op == SUB_VALUE) {
      switch (n) {
        case 1:
          c[0] -= a[0] * b[0];
          break;
        case 2:
          c[0] -= a[0] * b[0] + a[2] * b[1];
          c[1] -= a[1] * b[0] + a[3] * b[1];
          c[2] -= a[0] * b[2] + a[2] * b[3];
          c[3] -= a[1] * b[2] + a[3] * b[3];
          break;
        case 3:
          c[0] -= a[0] * b[0] + a[3] * b[1] + a[6] * b[2];
          c[1] -= a[1] * b[0] + a[4] * b[1] + a[7] * b[2];
          c[2] -= a[2] * b[0] + a[5] * b[1] + a[8] * b[2];
          c[3] -= a[0] * b[3] + a[3] * b[4] + a[6] * b[5];
          c[4] -= a[1] * b[3] + a[4] * b[4] + a[7] * b[5];
          c[5] -= a[2] * b[3] + a[5] * b[4] + a[8] * b[5];
          c[6] -= a[0] * b[6] + a[3] * b[7] + a[6] * b[8];
          c[7] -= a[1] * b[6] + a[4] * b[7] + a[7] * b[8];
          c[8] -= a[2] * b[6] + a[5] * b[7] + a[8] * b[8];
          break;
        default:
          for (var j = 0; j < n; j++) {
            for (var l = 0; l < n; l++) {
              for (var i = 0; i < n; i++) {
                c[i + j * n] -= a[i + l * n] * b[l + j * n];
              }
            }
          }
          break;
      }
    } else {
      switch (n) {
        case 1:
          c[0] += a[0] * b[0];
          break;
        case 2:
          c[0] += a[0] * b[0] + a[2] * b[1];
          c[1] += a[1] * b[0] + a[3] * b[1];
          c[2] += a[0] * b[2] + a[2] * b[3];
          c[3] += a[1] * b[2] + a[3] * b[3];
          break;
        case 3:
          c[0] += a[0] * b[0] + a[3] * b[1] + a[6] * b[2];
          c[1] += a[1] * b[0] + a[4] * b[1] + a[7] * b[2];
          c[2] += a[2] * b[0] + a[5] * b[1] + a[8] * b[2];
          c[3] += a[0] * b[3] + a[3] * b[4] + a[6] * b[5];
          c[4] += a[1] * b[3] + a[4] * b[4] + a[7] * b[5];
          c[5] += a[2] * b[3] + a[5] * b[4] + a[8] * b[5];
          c[6] += a[0] * b[6] + a[3] * b[7] + a[6] * b[8];
          c[7] += a[1] * b[6] + a[4] * b[7] + a[7] * b[8];
          c[8] += a[2] * b[6] + a[5] * b[7] + a[8] * b[8];
          break;
        default:
          for (var j = 0; j < n; j++) {
            for (var l = 0; l < n; l++) {
              for (var i = 0; i < n; i++) {
                c[i + j * n] += a[i + l * n] * b[l + j * n];
              }
            }
          }
          break;
      }
    }
  }

  ///     C <- A * B where A and B are not square
  static void matmat_ns(int l, int m, int n, ComplexList a, int lda,
      ComplexList b, int ldb, ComplexList c, int ldc, int op) {
    int n = x.length;

    /* C = A*B */

    if (op == INS_VALUE) {
      for (var j = 0; j < m; j++) {
        for (var i = 0; i < l; i++) {
          c[i + j * ldc] = 0.0;
        }
        for (var k = 0; k < n; k++) {
          for (var i = 0; i < l; i++) {
            c[i + j * ldc] += a[i + k * lda] * b[k + j * ldb];
          }
        }
      }
    } else if (op == SUB_VALUE) {
      for (var j = 0; j < m; j++) {
        for (var k = 0; k < n; k++) {
          for (var i = 0; i < l; i++) {
            c[i + j * ldc] -= a[i + k * lda] * b[k + j * ldb];
          }
        }
      }
    } else {
      for (var j = 0; j < m; j++) {
        for (var k = 0; k < n; k++) {
          for (var i = 0; i < l; i++) {
            c[i + j * ldc] += a[i + k * lda] * b[k + j * ldb];
          }
        }
      }
    }
  }

  ///     A <- A^-1 with Gaussian elimination
  static void ge(ComplexList a) {
    int n = x.length;

    /* compute inverse matrix with Gaussian elimination */

    var lu = new ComplexList.from(a);
    for (var k = 0; k < n; k++) {
      lu[k + k * n] = 1.0 / lu[k + k * n];
      for (var i = k + 1; i < n; i++) {
        t = lu[i + k * n] * lu[k + k * n];
        for (var j = k + 1; j < n; j++) {
          lu[i + j * n] -= t * lu[k + j * n];
        }
        lu[i + k * n] = t;
      }
    }
    for (var k = 0; k < n; k++) {
      for (var i = 0; i < n; i++) {
        t = (i == k);
        for (var j = 0; j < i; j++) {
          t -= lu[i + j * n] * a[j + k * n];
        }
        a[i + k * n] = t;
      }
      for (var i = n - 1; i >= 0; i--) {
        t = a[i + k * n];
        for (var j = i + 1; j < n; j++) {
          t -= lu[i + j * n] * a[j + k * n];
        }
        a[k * n + i] = t * lu[i + i * n];
      }
    }
  }

  ///     x <- A^-1 b
  static void solve(
      ComplexList a, ComplexList b, ComplexList x, ComplexList w) {
    int n = x.length;

    for (i = 0; i < n * n; i++) w[i] = a[i];

    switch (n) {
      case 1:
        x[0] = b[0] / w[0];
        break;
      case 2:
        w[0] = 1.0 / w[0];
        w[1] *= w[0];
        w[3] -= w[1] * w[2];
        w[3] = 1.0 / w[3];
        /* forward sub */
        x[0] = b[0];
        x[1] = b[1] - w[1] * x[0];
        /* backward sub */
        x[1] *= w[3];
        x[0] -= w[2] * x[1];
        x[0] *= w[0];
        break;
      default:
        for (var k = 0; k < n; k++) {
          w[k + k * n] = 1.0 / w[k + k * n];
          for (var i = k + 1; i < n; i++) {
            t = w[i + k * n] * w[k + k * n];
            for (var j = k + 1; j < n; j++) {
              w[i + j * n] -= t * w[k + j * n];
            }
            w[i + k * n] = t;
          }
        }

        /* forward sub */
        for (var i = 0; i < n; i++) {
          x[i] = b[i];
          for (var j = 0; j < i; j++) {
            x[i] -= w[i + j * n] * x[j];
          }
        }
        /* backward sub */
        for (var i = n - 1; i >= 0; i--) {
          for (var j = i + 1; j < n; j++) {
            x[i] -= w[i + j * n] * x[j];
          }
          x[i] *= w[i + i * n];
        }
        break;
    }
  }

  ///     Q * R <- A with classical Gram-Schmidt
  static int cgs(ComplexList a, ComplexList q, ComplexList r) {
    int n = x.length;
    var tol = 1e-12;

    var a_k = new ComplexList(n);

    for (i = 0; i < n * n; i++) {
      q[i] = 0.0;
      r[i] = 0.0;
    }

    for (var k = 0; k < n; k++) {
      for (var i = 0; i < n; i++) {
        a_k[i] = a[i + k * n];
      }
      for (var j = 0; j < k; j++) {
        r[j + k * n] = 0;
        for (var i = 0; i < n; i++) {
          r[j + k * n] += q[i + j * n] * a[i + k * n];
        }
        for (var i = 0; i < n; i++) {
          a_k[i] -= r[j + k * n] * q[i + j * n];
        }
      }
      nrm2(n, /*&*/ a_k[0], /*&*/ nrm2);
      r[k + k * n] = nrm2;
      if (nrm2 < tol) break;
      for (var i = 0; i < n; i++) {
        q[i + k * n] = a_k[i] / nrm2;
      }
    }
  }

  ///     Q * R <- A with modified Gram-Schmidt
  static void mgs(ComplexList a, ComplexList q, ComplexList r) {
    int n = x.length;
    var tol = 1e-12;

    var a_j = new ComplexList(n);

    for (var i = 0; i < n * n; i++) {
      q[i] = 0.0;
      r[i] = 0.0;
    }

    for (var j = 0; j < n; j++) {
      for (var i = 0; i < n; i++) {
        a_j[i] = a[i + j * n];
      }
      nrm2(n, /*&*/ a_j[0], /*&*/ nrm2);
      r[j + j * n] = nrm2;
      for (var i = 0; i < n; i++) {
        if (nrm2 < tol) break;
        q[i + j * n] = a_j[i] / nrm2;
      }
      for (var k = j + 1; k < n; k++) {
        r[j + k * n] = 0;
        for (var i = 0; i < n; i++) {
          r[j + k * n] += q[i + j * n] * a[i + k * n];
        }
        for (var i = 0; i < n; i++) {
          a[i + k * n] -= r[j + k * n] * q[i + j * n];
        }
      }
    }
  }

  ///     QR algorithm
  static ComplexQR qr(ComplexList a, ComplexList q, ComplexList r) {
    int n = x.length;
    var maxiter = 100000;
    var tol = 1e-12;

    var a0 = new ComplexList(n * n);
    var iter = 0;
    while (iter < maxiter) {
      iter = iter + 1;
      cgs(n, a, q, r);
      for (var j = 0; j < n; j++) {
        for (var i = 0; i < n; i++) {
          a[i + j * n] = 0;
          for (var k = 0; k < n; k++) {
            a[i + j * n] += r[i + k * n] * q[k + j * n];
          }
        }
      }
      err = _ssqrt(a[1] * a[1]);
      if (err < tol) break;
    }

    return new ComplexQR(iter, err);
  }
}

class ComplexQR {
  final int qriter;
  final Complex qrerr;
  ComplexQR(this.qriter, this.qrerr);
}
