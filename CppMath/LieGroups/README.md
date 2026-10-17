# Frankel Lie groups — C++20 counterpart

Numerical and independently authored exact symbolic SO(3)/SU(2), alongside
[the idiomatic Rust implementation](../../Rust/lie-groups/README.md).

```cpp
#include "lie_groups.hpp"
using namespace mathphysics;
using SingleMatrix = Matrix<float, 3, 3>;
using DoubleMatrix = Matrix<double, 3, 3>;

auto m = SingleMatrix::identity();
auto r = So3Algebra<double>({0.0, 0.0, 1.5707963267948966}).exp();
auto y = r.apply({1.0, 0.0, 0.0}); // approximately (0,1,0)
auto u = Su2Algebra<double>({0.0, 0.0, 6.283185307179586}).exp();
// u is approximately -I, not +I; u.to_so3() is approximately I.
```

- General matrices use std::array, zero-based checked access and compile-time
  dimensions; no unchecked initializer-list copies.
- Group values **contain**, not inherit, matrices/quaternions. Private state
  prevents arbitrary matrix mutation or accidental determinant-minus-one values.
- Numerical float/double/long double are supported. The Rust side additionally
  supplies tested software binary16; native C++ half is not claimed here.
- Frankel E_i, hat/vee, cross, epsilon, both Pauli normalizations, exponentials,
  composition, inverses, checked imports and covering maps mirror the Rust API.
- Group/domain errors throw standard exceptions. Floating inverse/composition
  revalidate, rather than claiming exact closure at finite precision.
- include/symbolic.hpp adds arbitrary-precision cpp_rational coefficients,
  sparse real polynomials, a one-angle trigonometric quotient, and real
  stereographic SU(2)/SO(3) values. It uses installed Boost headers only.
- Rational matrix denominators are private and checked; move construction and
  move assignment leave source objects valid (zero/1 matrices or identity SU2),
  preventing a moved-from 0/0 false equality.

```bash
cmake -S . -B "$HOME/.cache/mathphysics-lie-groups/cpp-manual" -DCMAKE_BUILD_TYPE=Debug
cmake --build "$HOME/.cache/mathphysics-lie-groups/cpp-manual" -j2
ctest --test-dir "$HOME/.cache/mathphysics-lie-groups/cpp-manual" -V
```

The dependency-free test runner uses named test groups and enforced exceptions,
not NDEBUG-sensitive assert. It follows the existing per-component/convention
test style without requiring a GoogleTest download. CTest registers the main
contract suite and an epsilon TDD regression. The joint Rust validation script
also builds Release and ASan+UBSan configurations and compares both languages
against live SymPy/Sage.

See the Rust README/docs for derivations, sources, exact scope and evidence.
