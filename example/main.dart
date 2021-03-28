import 'package:complex/complex.dart';

void main() {
  const z1 = Complex(1);
  const z2 = Complex(3, 4);

  print(z1.abs()); // 1.0
  print(z2.abs()); // 5.0
  print(z2.conjugate()); // (3.0, -4.0)
}
