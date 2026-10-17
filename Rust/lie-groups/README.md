# SO(3), so(3), SU(2), su(2) — numerical and symbolic Rust

An independently authored, scalar-generic toolkit following Frankel §§19.1a–b.
The original C++/Rust experiments, book corpus and study notes remain unchanged.
Matching C++20 implementation: [../../CppMath/LieGroups](../../CppMath/LieGroups).

## Compile-time scalar and dimension selection

```rust
use lie_groups::{Matrix, So3, So3Algebra, Half};

type SingleMatrix = Matrix<f32, 3, 3>;
type DoubleMatrix = Matrix<f64, 3, 3>;
type HalfMatrix = Matrix<Half, 3, 3>;

let m = SingleMatrix::identity();
let r: So3<f64> = So3Algebra::new([0.0, 0.0, std::f64::consts::FRAC_PI_2])
    .exp().unwrap();
let y = r.apply(&[1.0, 0.0, 0.0]).unwrap(); // approximately [0,1,0]
```

Rust generics and const generics provide the C++ template mechanism:
Matrix<T,R,C> is monomorphized, and incompatible multiplication dimensions
do not compile. Nested arrays avoid unstable generic-const R*C expressions.
Real/RationalScalar/Scalar traits express required operations; there is no
dynamic dispatch here. See the compiled [example](examples/frankel.rs).

Composition is preferable to matrix inheritance **in both languages**:
a rotation is not closed under matrix addition, arbitrary mutation or scaling.
So3 contains a private matrix; Su2 contains a private unit quaternion.
Read-only accessors expose representations without breaking their invariants.
General matrices support addition, scaling, products and indexing normally.

### Half precision: what is actually supported

Stable rustc 1.97.1 rejects native f16 (E0658; reproduced locally).
Half is our two-byte IEEE binary16 storage type, not an alias for f32.
It provides direct f64 conversion with nearest/ties-even rounding and
f32-promoted software arithmetic/transcendentals, rounded back to binary16.

All 65,536 encodings are tested: non-NaN encodings round-trip, NaNs remain NaNs;
ties, subnormals, signed zero and overflow have explicit controls.
A 128-fixture rotation comparison to f64 at the same represented inputs
measured maximum component error about 0.00142 (test bound 0.01).
This is not native f16 acceleration or a performance benchmark.
Add an in-crate scalar implementation or a local newtype for an external half
type; Rust orphan rules prevent arbitrary downstream foreign-trait/foreign-type
implementations.

## Numerical APIs

- Matrix<T,R,C>: row-major, zero-based, stack arrays; transpose, trace,
  dimension-checked product, vector application, 2x2/3x3 determinant.
- Complex<T>: explicit real and imaginary components, arithmetic and adjoint.
- So3Algebra<T>: coordinates, Frankel basis, hat/vee, bracket.
- So3<T:Real>: checked matrix construction, axis-angle, exponential,
  composition, checked inverse, vector action.
- Su2Algebra<T>: coordinates in T_i=-i sigma_i/2, basis, checked matrix
  conversion, cross-product bracket and differential map to so(3).
- Su2<T:Real>: checked unit quaternion/matrix construction, exponential,
  full-Pauli-over-i exponential, composition, inverse, opposite lift, SO(3) cover.
- levi_civita<const N>(indices): const-evaluable alternating symbol with
  explicit out-of-range errors. The 3D symbol is not Lorentz/spinor epsilon.

Numerical values are **approximations**, not exact group-membership proofs.
Construction checks finite entries, orthogonality/unitarity and determinant+1.
Default entrywise tolerance is 32*scalar epsilon; supplied tolerances must be
finite, nonnegative and<1/8. Composition and SO(3) inversion revalidate at the
default tolerance and can fail on drift. No silent normalization or projection
onto the group occurs. SU(2) inversion/antipodes preserve its quaternion norm
by sign changes. Small-angle series avoid cancellation; large-angle error is
still limited by scalar precision and library trigonometric argument reduction.

## Exact symbolic APIs: three different encodings

1. **Cadabra-like algebra coefficients:** So3Algebra<Poly> and
   Su2Algebra<Poly> hold generic linear combinations. Bilinearity, epsilon,
   brackets, cross products and Jacobi are actual polynomial calculations.
   The compact real form requires real coordinates. Generic complex coefficients
   instead complexify the algebra; the type name alone does not impose a real
   form. Polynomial adjoints conjugate coefficients while fixing formal variables
   (interpreted as real here), not arbitrary complex-valued assignments.
2. **SymPy-like explicit symbolic matrices:** Trig is the exact quotient
   Q[s,c]/(s²+c²-1), with normal form a(s)+c*b(s). It implements a **single
   common formal angle**, not a general trigonometric simplifier.
   so3_exponential accepts an exact rational unit axis and represents
   exp(theta*hat(n)); su2_exponential retains argument_scale=1 for
   exp(theta*(-i n.sigma)), versus 1/2 for exp(theta*(-i n.sigma/2)).
   **Keep the argument metadata:** the formal matrices alone look identical.
3. **Sage-like constrained symbolic group values:** SymbolicSu2 uses
   q=±(1-|t|²,2t)/(1+|t|²) over real symbols, with sealed construction and
   exact rational-function matrices. Two charts cover SU(2). SymbolicSo3
   projects a lift and identifies opposite lifts; one chart covers SO(3).
   These are stereographic **quaternion** parameters (t=n*tan(theta/4)),
   not the usual matrix Cayley chart that misses pi rotations.

Q/C/Poly come from the local qft-symbolics crate: checked i128 coefficients,
degree/term budgets, no dependency downloads. Overflow/resource exhaustion
can panic; CLI boundaries never promote a panic to PASS. RationalMatrix
denominators are private/nonzero polynomials; evaluation poles are errors,
and equality is on the common nonzero domain. Chart parameters and symbols
are interpreted as real; complex-coefficient inputs/evaluations are rejected.
Exact evaluation currently accepts rational values, not arbitrary radicals;
formal real variables are more general than those finite witnesses.

This is worth implementing as a focused physics tool. It is **not** a Rust
translation of all SymPy/Sage/Cadabra: no general parser, Jordan exponential,
Gröbner engine, assumptions engine, logarithm branch selection or all Lie types.

## Conventions that the tests enforce

Column vectors, active right-handed rotations; code epsilon(0,1,2)=+1 corresponds
to book epsilon123=+1. Hat(v)w=v cross w, so positiveZ quarter-turn sendsX toY.

    [E_i,E_j] = epsilon_ijk E_k
    [sigma_i,sigma_j] = 2 i epsilon_ijk sigma_k
    T_i = -i sigma_i/2: [T_i,T_j] = epsilon_ijk T_k
    F_i = -i sigma_i:   [F_i,F_j] = 2 epsilon_ijk F_k

The Pauli matrices are Hermitian and determinant-1: they are neither su(2)
generators as a real anti-Hermitian algebra nor SU(2) elements themselves.
The SU(2) half-angle exponential maps to the SO(3) full-angle exponential.
A 2pi spinor turn gives-I, a 4pi turn gives+I, and q/-q project to the same
rotation without being the same SU(2) value.

See [derivation source](docs/frankel-lie-groups.tex),
[source/CAS audit](docs/SOURCE-AUDIT.md), and [verification](docs/VERIFICATION.md).

## Run

```bash
export CARGO_TARGET_DIR="$HOME/.cache/mathphysics-lie-groups/target"
cargo test --offline
cargo run --offline --example frankel
cargo run --offline -- verify
cargo run --offline -- verify --wrong-sign # FAIL, exit 1
cargo run --offline -- symbolic
bash validation/run.sh
```

Full validation uses the **existing** SageMath/Cadabra images and locally cloned
SymPy, with Docker --network none / --pull=never. It includes Rust debug/release,
C++ debug/release/ASan+UBSan, exact SymPy/Sage checks of both implementations,
Cadabra indexed epsilon reductions, independent numerical Taylor references,
LaTeX compilation and the old QFT kernel's complete regression suite.
Artifacts stay under ~/.cache/mathphysics-lie-groups/.

CLI exit 0 means PASS or COMPUTED (not interchangeable), 1 means FAIL, 2 means
INCONCLUSIVE. The fixed-symbol interchange files are trusted producer outputs;
CAS expression evaluators are **not safe arbitrary-model-input parsers**.

## Licensing

New code is covered by the repository MIT license. All algorithms are independently
derived from mathematics, not translations of Sage/Cadabra GPL source.
SymPy(BSD), Sage(GPL), Cadabra(GPL) run as external verification tools.
C++ exact coefficients use installed Boost.Multiprecision headers (Boost license);
no third-party source is vendored. Public Sage source was inspected by fetching
only four selected files plus its license to external cache, not cloning Sage.
