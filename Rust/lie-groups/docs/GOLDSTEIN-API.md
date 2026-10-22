# Goldstein rotation APIs

`goldstein` extends the existing checked `So3<T>` and `Su2<T>` types. All matrices
act on column vectors. `R` is active right-handed, equivalently the body-to-inertial
component map; Goldstein's passive inertial-to-body matrix is `A_G = R^T`.

All Euler arrays use **product order**:

| Convention | Input/output angles | Matrix |
|---|---|---|
| `Zxz` | `[phi, theta, psi]` | `Rz(phi) Rx(theta) Rz(psi)` |
| `Zyz` | `[phi, theta, psi]` | `Rz(phi) Ry(theta) Rz(psi)` |
| `Zyx` | `[yaw psi, pitch theta, roll phi]` | `Rz(psi) Ry(theta) Rx(phi)` |

```rust
use lie_groups::goldstein::{from_euler, extract_euler, axis_angle, EulerConvention};
let r = from_euler(EulerConvention::Zyx, [0.4_f64, 0.2, -0.3]).unwrap();
let angles = extract_euler(&r, EulerConvention::Zyx);
assert!(!angles.singular);
let finite_rotation = axis_angle(&r).unwrap();
```

`from_euler`, `euler_lift` and `CayleyKlein<T>` support the existing `Real`
backends (`f64`, `f32`, software binary16 `Half`). Numerical membership is checked;
operations may return errors on excessive drift, and do not silently normalize.
Inverse Euler/axis-angle/lift extraction is specifically `f64`.

`extract_euler` gives proper-Euler middle angle `[0,pi]` or ZYX pitch
`[-pi/2,pi/2]`. At poles it fixes the last product angle to zero and marks
`singular`; pole classification uses 64 machine epsilons. Near-pole angle values
are ill-conditioned even when reconstruction of the matrix remains accurate.
`axis_angle` uses a stable quaternion pivot and atan2, including near zero/pi;
its principal angle is `[0,pi]`, with arbitrary axis `[1,0,0]` at the identity.
`lift_rotation` revalidates at default tolerance without normalization; loosely
validated imported matrices can therefore produce an error. It chooses a lift
with nonnegative scalar part; this choice cannot
recover an input lift's winding/sign. `euler_lift` instead composes the
half-angle lifts, and `CayleyKlein` retains `q` versus `-q` as distinct values.

For scalar-first Hamilton `q=[w,x,y,z]`, Goldstein's lift is
`K=U(q)^dagger=[[w+iz,y+ix],[-y+ix,w-iz]]`; its SO(3) action is `R(q)^T`.
`CayleyKlein::try_from_parameters` checks all four complex parameters.
The generic `cayley_klein_parameters` and `cayley_klein_matrix` functions instead
return exact polynomial expressions with unit-norm/real-symbol assumptions
explicitly left to the caller.

`axis_rotation`, `euler_matrix`, `rodrigues`, `orthogonal_conjugate` and
`cayley_hamilton_residual` are scalar-generic raw matrix functions. Rodrigues
requires a unit axis and a unit-circle cosine/sine pair. Orthogonal conjugation
uses `B A B^T`; the vector dictionary `Ad_B(hat(w))=hat(Bw)` requires **proper**
`B`. Cayley–Hamilton is the characteristic-polynomial identity for any 3x3
matrix, distinct from the Cayley–Klein representation.

`multitrig::MultiTrig<N>` implements the exact coefficient ring
`Q[s_i,c_i]/(s_i^2+c_i^2-1)` with at most eight independent real angles. It permits
three-angle Euler expressions and half-to-full quaternion checks using explicit
double-angle substitutions. It has no general assumptions engine, inverse
trigonometry, addition-identity search, parser, or CAS integration. Normal forms
have cosine exponents 0/1, total sine degree at most 64 and at most 4096 terms;
resource exhaustion and inherited checked-i128 overflow panic. Rational
evaluation rejects assignments off the unit circles. Display output uses fixed
`s0,c0,...` symbols for trusted-producer external CAS interchange.

Run from this crate, with build artifacts outside the repository:

```sh
CARGO_TARGET_DIR=/tmp/mathphysics-rust-goldstein cargo test --offline
CARGO_TARGET_DIR=/tmp/mathphysics-rust-goldstein cargo run --quiet --offline --example goldstein_fixtures
CARGO_TARGET_DIR=/tmp/mathphysics-rust-goldstein cargo run --quiet --offline --example goldstein_fixtures -- --symbolic
```

The fixture example produces NDJSON: 24 numerical edge/regular records, or six
exact full/half-angle records under `--symbolic`. These are mathematical
verification fixtures; names and field values do not certify arbitrary inputs.
