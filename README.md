# Complex

A representation of a complex number, i.e. a number which has both a
real and an imaginary part.

Translated from
[The Apache Commons Mathematics Library](https://commons.apache.org/proper/commons-math/)
into [Dart](https://www.dartlang.org/).

Implementations of arithmetic operations handle `NaN` and
infinite values according to the rules for `double`, i.e.
`==` is an equivalence relation for all instances that have
a `NaN` in either real or imaginary part, e.g. the following are
considered equal:

- `1 + NaNi`
- `NaN + i`
- `NaN + NaNi`

Note that this is in contradiction with the IEEE-754 standard for floating
point numbers (according to which the test `x == x` must fail if
`x` is `NaN`).

## Example

```dart
const z1 = Complex(1);
const z2 = Complex(3, 4);
print(z1.abs()); // 1.0
print(z2.abs()); // 5.0
print(z2.conjugate()); // (3.0, -4.0)
```
