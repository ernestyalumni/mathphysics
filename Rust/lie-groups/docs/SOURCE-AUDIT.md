# Source audit — 2026-10-17

No originals, CAS source, Propulsion code or study notes were modified.
New implementations are independent mathematical code, not source translations.

## Frankel and current notes

Original: user-supplied *The Geometry of Physics*,3rd edition PDF.
Canonical existing corpus: Frankel-GeometryOfPhysics-3e/book.md and pages/.
Read chapter 19 pages 0555–0563; directly rendered and visually checked original
PDF 555–556 (printed 492–493) for Eqs. 19.1–19.6. E_i signs, cross-product action,
Pauli entries and both generator normalizations agree with the source image.

The parsed 19.3 has E_i omega^j; the original image has E_j omega^j.
We follow the original, without silently treating OCR as textbook authority.
No fresh whole-book OCR or book-wide source certification was performed.
Rendered review images/hashes remain in the external source-audit cache.

Existing mathphysics notes inspected:
- documents/mono-notes/topics/quaternions-spin.tex:166–270,390–477,636–723.
- topics/su2-so3-covering.tex:44–167,477–576,644–722.
- topics/rigid-body-foundations.tex:21–43.

Those notes sometimes call sigma/(2i) **F**. The new API avoids that collision
by naming the normalization HalfPauliOverI or PauliOverI; this project's
explanations call the corresponding generators T and F respectively.

## Ernest's existing implementations

Propulsion checkout inspected at be47b68cd9919d31a3b684ab9e0dee94044501db:
- Cosmos/Source/Algebra/Modules/Vectors/Mat3x3.h:24–55:
  template<Field>, nine row-major entries, but a virtual constructor and an
  undefined return type T prevent this header compiling as written. Its
  unbounded initializer-list copy is unsafe. We did not change it.
- Stunticons/wildrider/src/algebra/rings/rotations.rs:4–157:
  generic row-major arithmetic, but Rotation3 accepts arbitrary matrices and
  arbitrary mutation. The new general matrix test retains its 9-entry fixture.
- Stunticons/wildrider/src/algebra/modules/vectors/quaternions.rs:39–134:
  right-handed cross product, but one-based vector indexing differs from the
  neighboring zero-based matrix.
- Cosmos/QuaternionConventionLab/{cpp/quaternion.hpp,rust/src/lib.rs}:
  scalar-first Hamilton quaternions and active body/local-to-world convention.
- CUDACFD/module.cpp:107–158: float skew generators and weighted hat matrix.
- T1000/Source/OrbitalMechanics/StateVector.hpp:61–75: right-handed cross product.
- Cosmos/Source/UnitTests/Algebra/Modules/Vectors/Vector3_tests.cpp:
  named GoogleTest groups, per-component checks and type-trait assertions.

mathphysics/CppMath has algebra/ring/GCD experiments, not an established general
dense matrix framework. The new CppMath/LieGroups module is separate.

## Cadabra: indexed formal tensors, not inferred group values

Checkout67cdebce798dbb2e230ae3246221ca3302a03129:
- tests/field_theory.cdb:151–162: three-dimensional EpsilonTensor(delta=delta),
  epsilon_to_delta and expand_delta.
- tests/components.cdb:270–274,841–845: epsilon components.
- tests/substitute.cdb:168,182: SelfNonCommuting indexed objects.
- core/modules/Lie.cc: representation/singlet calculations through LiE, not a
  continuously parameterized group-element constructor.

The new live Cadabra script declares noncommuting E_a and explicit epsilon
structure constants; it actually reduces epsilon contractions and Jacobi,
with a wrong-sign control left nonzero. Cadabra does not deduce a Lie algebra
from the name E or automatically know that an arbitrary matrix belongs to SO3.

## SymPy: explicit symbolic matrices and evaluated representations

The user's checkout is nested: repos/sympy/sympy, Git
50e68debd256f0ff52b1856bba54868e5499448e, version 1.15.0.dev.
Our oracle imports this checkout read-only in the existing Sage container
(to reuse its installed mpmath; no package was installed).

Files relative to that checkout:
- sympy/functions/special/tensor_functions.py:19–59: LeviCivita evaluates
  alternating symbols; its numeric integer cases are not a metric volume form.
- sympy/physics/matrices.py:9–43: actual Pauli matrices, msigma(1..3).
- sympy/physics/quantum/spin.py:434–743: formal Rotation with evaluated Wigner
  representations, passive z-y-z convention. j=1 is a spherical-basis
  representation, not already a Cartesian real 3x3 matrix.
- sympy/algebras/quaternion.py:658–667,1194–1201,1283–1303:
  quaternion normalization and actual rotation-matrix conversion.
- sympy/matrices/matrixbase.py:4936–4969: matrix exponentials via Jordan forms,
  potentially NotImplemented for difficult symbolic inputs.
- sympy/liealgebras/: root/Cartan/Weyl classification is not automatically a
  compact-real SU2/SO3 element API.

The oracle actually executes LeviCivita, msigma, Matrix.exp and Quaternion
conversion, comparing both Rust and C++ matrices exactly.

## Sage: inspected selected files, and executed installed Sage 10.8

Only four source files plus COPYING.txt were fetched, not a huge clone.
Public source pin0835358cba62dc19fccdff54f32e52fded9f22b0:
https://github.com/sagemath/sage/tree/0835358cba62dc19fccdff54f32e52fded9f22b0

- src/sage/groups/matrix_gps/orthogonal.py:147–222,338–373,472–545:
  ring/form-selected orthogonal groups; default formI; generic check X F X^T=F.
- src/sage/groups/matrix_gps/unitary.py:91–146,267–301,386–437:
  specified Hermitian form, coefficient conjugation, unitary/determinant gates.
  A finite-field SU is not the physical continuous SU(2).
- src/sage/algebras/lie_algebras/classical_lie_algebra.py:557–705,1082–1139:
  special-linear/split orthogonal matrix models versus separate compact forms.
  Cartan type A1 does not by itself select compact su2 over the reals.
- src/sage/manifolds/differentiable/metric.py:1681–1858:
  metric volume form, with ordered component sqrt(abs(det(g))) in a positive
  frame. The metric and orientation matter; a bare alternating symbol is not
  the whole tensor in curvilinear coordinates.

The live Sage oracle constructs an alternating form and Euclidean metric volume
form; an explicit compact structure-constant Lie algebra with negative-definite
Killing matrix; exact matrices and their conjugation covering; native group
positive/negative witnesses.

**Concrete validator trap:** the inspected orthogonal factory passes a false
special flag in its generic path. In the installed Sage 10.8 image,
SO(3,QQ)(diag(-1,1,1)) actually succeeds. This was reproduced, not inferred only
from source. Our oracle independently requires determinant+1 and records that
diagnostic. We neither fixed Sage nor treated its class name as ground truth.
This is a scoped observation of those artifacts, not a claim about every Sage
backend/version. Native SU2 correctly rejects the Pauli negative witness.

The downloaded source pin and installed image are separate artifacts; build
equivalence was not established. Sage image:
sha256:aeef59a8c17212357aeb83bd7f1cfac43bb1c789ae3fa2d887227c763bdfcac3.

## What was worth implementing

Keep generic dense matrices, constrained groups, algebra coordinates and symbolic
expressions as separate layers. Reuse simple finite mathematical encodings
rather than translate an entire CAS or confuse abstract classification with a
chosen real form. All mathematical results here are known, not novel physics.

Licenses: new code MIT; SymPy BSD-3-Clause; Sage/Cadabra GPL; Boost.Multiprecision
uses the Boost Software License. No GPL implementation/test text is copied.
