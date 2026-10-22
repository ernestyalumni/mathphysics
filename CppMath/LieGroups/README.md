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

## Goldstein rotations (§§4.2–4.7)

`include/rotations.hpp` extends the same validated `So3<T>`/`Su2<T>` values,
matrix and algebra types (`T = float`, `double`, `long double`). Numerical
angles are radians; matrices act on column vectors with right-handed positive
angles. These are explicit products, rather than a claim that all engineering
software uses identical conventions:

| Constructor | Matrix | Chronological moving-axis sequence |
|---|---|---|
| `from_euler_zxz(phi,theta,psi)` | `Rz(phi) Rx(theta) Rz(psi)` | z, x', z'' |
| `from_euler_zyz(phi,theta,psi)` | `Rz(phi) Ry(theta) Rz(psi)` | z, y', z'' |
| `from_yaw_pitch_roll(psi,theta,phi)` | `Rz(psi) Ry(theta) Rx(phi)` | heading, pitch, roll |

The last product is a common aerospace/engineering ZYX convention; frame
handedness, NED/ENU axes and body-to-world versus world-to-body remain separate
choices. In the notes, `a_I=R a_B` is the passive coordinate conversion of a
fixed vector. The same numerical matrix gives an active rotation when both
input and output are expressed in one fixed frame. Goldstein's matrix is
`A_G=R.transpose()`, so its chronological passive Euler steps multiply in the
reverse displayed order, with opposite elementary angle signs.

```cpp
#include "rotations.hpp"
using namespace mathphysics;
auto r = from_euler_zxz(0.4, 0.7, -0.3);
auto e = to_euler_zxz(r); // first=phi, middle=theta, last=psi
auto q = quaternion_from_euler_zxz(0.4, 0.7, -0.3); // signed SU(2) lift
auto a = to_axis_angle(r); // principal angle [0,pi]
auto x = rotation_log(r);  // x.exp() reconstructs r, up to roundoff
auto ck = CayleyKlein<double>::from_quaternion(q);
auto a_goldstein = ck.goldstein_rotation(); // R(q)^T
```

`to_euler_zxz`, `to_euler_zyz` return the canonical middle angle in `[0,pi]`;
`to_yaw_pitch_roll` returns pitch in `[-pi/2,pi/2]`. Outer angles use `atan2`
in `[-pi,pi]`. At gimbal lock, `singular=true` and `last=0` selects one valid
representative; the original separate outer angles cannot be recovered.
The optional tolerance (default `32*epsilon(T)`) bounds the relevant sine or
cosine used to identify a singularity, so thresholding a near-singular value
introduces an error of that scale. Invalid/nonfinite tolerances are rejected.
Conversions are tested by reconstructing matrices, not comparing nonunique
angle triples. No globally continuous inverse parameterization exists.

`quaternion_from_rotation` uses the largest-component pivot and normalizes the
extracted quaternion to remove roundoff. It chooses `w>=0`; it does not preserve
the sign of a previously specified lift. Direct half-angle Euler constructors
and `CayleyKlein::from_quaternion` preserve that sign, and `q` and `-q` remain
distinct `Su2` values. `to_axis_angle` uses `atan2`, avoiding trace/`acos` loss
near identity and skew/`sin(phi)` division near a half-turn. At identity it
chooses axis `(1,0,0)` with `axis_line_unique=false`; otherwise the axis line is
unique. At `phi=pi`, either oriented axis sign gives the same rotation. The
principal logarithm is necessarily discontinuous there. A `So3` imported with
a loose tolerance can be normalized by extraction; no exact equality with an
arbitrary approximate imported matrix is claimed.

For scalar-first Hamilton `q=(w,x,y,z)`, the Goldstein Cayley–Klein dictionary is
`alpha=w+iz`, `beta=y+ix`, `gamma=-y+ix=-conj(beta)`,
`delta=w-iz=conj(alpha)`, and `|alpha|^2+|beta|^2=1`.
`K=U(q)^dagger` maps to `A_G=R(q)^T`.
`CayleyKlein::from_parameters` validates these constraints through the existing
checked SU(2) import; it rejects malformed conjugacy, non-unit and nonfinite
parameters. Here the dagger is the Hermitian matrix adjoint, distinct from
`adjoint_action(B,X)=hat(B*x)` for `X=hat(x)` in the Lie algebra.
`coordinate_change(B,A)` computes `B A B^T` for any finite 3x3 operator `A`
and validated rotation `B`, checking overflow/nonfinite results.

`cayley_hamilton_residual(M)` computes
`M^3-tr(M)M^2+c2(M)M-det(M)I` for a general 3x3 matrix over commuting scalars.
`c2` is the sum of principal 2x2 minors; no division or rotation assumption is
used. It works with the existing arbitrary-precision `Polynomial`, `Exact`
and one-angle `Trig` domains as well as floating point. The generic
`axis_rotation(axis,c,s)` and `cayley_klein_matrix(q)` also accept exact scalars;
these raw algebraic constructors do not assert unit-circle/unit-norm inputs.
For example, generic Euler expressions are simply products of three
`axis_rotation(axis,c,s)` matrices. The C++ symbolic layer still uses its
existing one-common-angle quotient; no independent multi-angle trigonometric
simplifier or general CAS is claimed.

The new CTest `goldstein_rotation_contracts` has ten named groups: all three
floating types, independently summed Taylor references and geometric Rodrigues,
Euler inverse branches/singularities, tiny/pi/large-angle cases, CK lift and
transpose dictionaries, conjugation/Ad, generic 9-symbol Cayley–Hamilton,
exact hat/CK/quaternion identities, and wrong-sign/order/domain controls.
`goldstein_rotation_interchange` emits 18 JSONL fixtures for independent Python
comparisons. Outputs and build trees belong outside the source checkout.

For ASan+UBSan under ptrace, use `ASAN_OPTIONS=detect_leaks=0 ctest ...`:
LeakSanitizer cannot operate in that environment. This disables leak checks,
not address/undefined-behavior instrumentation. It is an environment limitation,
not evidence that a failing sanitizer run passed.

## COM-centered 6DOF numerical dynamics

`include/rigid_body.hpp`, namespace `mathphysics::rigid_body`, adds a
double-precision `Model`, `State`, `Wrench`, `Derivative`, `StepOptions`, and
`StepResult`. The same model has a Rust implementation and independently
tested [Python dynamics](../../Python/rotations/DYNAMICS.md).

```cpp
#include "rigid_body.hpp"
using namespace mathphysics::rigid_body;
auto j = Mat::from_rows({{{2,.2,.1},{.2,3,.3},{.1,.3,4}}});
Model model(2.0, j);
State state; state.omega_body = {.3, -.4, .8};
auto derivative = model.rhs(state); // zero wrench by default
auto next = model.step(state, 1.0 / 240.0).state;
```

`position`, `velocity` are inertial COM coordinates; scalar-first Hamilton
`quaternion=(w,x,y,z)` maps body to inertial; `omega_body` and full symmetric
inertia use body axes at the COM. `Wrench` has `force_body`, `torque_body`,
and optional `force_inertial`, all zero by default. There is no implicit gravity
or damping. `kinetic_energy` is total kinetic energy; `momentum_inertial` is
spin angular momentum about the COM, not an orbital contribution.

The constructor checks positive finite mass and full symmetric SPD inertia
with physical principal-moment triangle inequalities (`1e-12` relative
roundoff). Accepted symmetry roundoff is averaged. Inertia solves use scaled
Cholesky factors; axes need not be principal. Quaternion inputs must satisfy
`abs(sum(q_i*q_i)-1)<=1e-10`; accepted roundoff is normalized, preserving sign.
Arbitrary nonunit quaternion data are rejected.

The RHS is `rdot=v`, `vdot=(R*force_body+force_inertial)/mass`,
`qdot=q*(0,omega_body)/2`, `J*omegadot=torque_body-omega_body cross J*omega_body`.
Steps solve implicit midpoint in `h=J*omega` with an analytic Jacobian and
pivoted Newton solve, then use the Cayley quaternion increment
`normalize((1,dt*omega_mid/2))`. Body force uses the quaternion midpoint;
position uses midpoint velocity. Wrench inputs remain constant during a step.
Defaults: 20 Newton iterations, residual infinity norm bounded by
`1e-12*(1+|h0|inf+dt*|torque_body|inf)`. `StepResult` reports iterations and
residual. Finite `dt>0` has no arbitrary cap; invalid input, singular Jacobian,
nonconvergence, or nonfinite arithmetic throw. Large converged steps may still
be inaccurate. Output quaternion normalization only corrects roundoff.
Newton options require tolerance in `(0,1e-4]` and iterations in `[1,100]`.

In exact arithmetic free steps preserve rotational energy, `|h|`, and `R*h`.
Forced steps satisfy the discrete body-torque work identity; no exact inertial
midpoint torque quadrature is claimed. This is a second-order integrator.
`rigid_body_contracts` has eight groups covering exact translation, sphere and
principal rotations, full-inertia conservation, intermediate-axis flips,
exact symmetric-top attitude and convergence, forcing/RHS signs, opposite
lifts/discrete work, and invalid inputs/solver failure.
`rigid_body_interchange` emits eight JSONL fixtures for independent Python
DOP853 and Rust comparisons. Build artifacts stay outside the checkout.
