# Verification contract and receipt

Run bash validation/run.sh from this crate. The latest full-run location is
recorded below; source hashes bind the actual tested tree.

## Goldstein extension — 2026-10-21

The [Goldstein verification receipt](GOLDSTEIN-VERIFICATION.md) records the
new Euler/axis-angle/Cayley–Klein/characteristic-polynomial coverage and
completed independent Python/SymPy/Cadabra and full native/QFT regressions.
The original October 17 receipt below is historical; new counts and domains
are recorded in the extension receipt and [API guide](GOLDSTEIN-API.md).

## Acceptance gates

- Rust 27 integration tests and 4 doctests (three compile-fail, one positive),
  debug and release; formatting; Clippy all targets with warnings denied.
- Half checks all 65,536 encodings, exact conversion boundaries/ties, and 128
  rotation witnesses against f64 at identical represented inputs.
- C++ 13 named runtime test groups plus epsilon regression, in Debug, Release
  and ASan+UBSan core-test configurations. Assertions are not disabled by NDEBUG.
- All 9 basis commutators for SO3, Pauli/T/F normalization distinctions, symbolic
  cross/hat/bracket/Jacobi, exact traces, matrix/algebra import rejection.
- Symbolic closed exponentials obey orthogonality/unitarity and determinant 1
  in the one-angle quotient. Half-angle metadata and covering are tested.
- Arbitrary-real-parameter rational charts have exact unit norm, correct
  determinants, inverse/covering identities; q/-q are distinct in SU2 but not SO3.
- Exact exported Rust/C++ matrices compared with actual SymPy and actual Sage;
  Cadabra reduces indexed epsilon contractions/Jacobi, rejects a wrong sign.
- 25 f64 cases, no skipped inputs: Rust/C++ parity and independent 80-term
  matrix-exponential references, tolerance 1e-12.
- CLI verified/wrong/unsupported verdicts and exit 0/1/2; examples compile/run.
- New derivation LaTeX compiles externally; old qft-symbolics full regression
  suite runs, including its SMGA/Python/Cadabra checks.

## TDD and review evidence

Initial epsilon implementation returned 0. Rust's 3 initial tests genuinely
failed, and the initial C++ epsilon regression failed. Red receipts are in
~/.cache/mathphysics-lie-groups/probes/{rust-red,cpp-red}.log.
Later tests were written alongside each bounded implementation stage.

Read-only source/math scouts inspected existing implementations and CAS code.
Two broad reviews exceeded the child turn limit (not successful reviews).
Two narrowed reviews completed and identified:
1. **Unchecked floating inverse:** bounding M^T M-I entrywise does not bound
   M M^T-I identically. A supplied counterexample passed construction at 0.12
   but its transpose failed even the maximum allowed tolerance. Rust inverse
   now returns Result and revalidates; C++ revalidates/throws. Regression added.
2. **C++ moved-from 0/0:** private constructors do not protect denominators when
   default moves empty the underlying polynomial maps. Explicit moves now
   leave zero/1 rational matrices and identity SU2 values behind. Move
   construction/assignment regression tests cover both.

The first review counterexample was reproduced as a failing C++ test before
the fix (probes/review-red.log). Symbolic zero denominators are also rejected
directly. These corrections matter more than a large raw test count.

## Interpretation

Exact identities mean zero residuals in the stated real-symbol/polynomial or
single-angle quotient domains. Numerical agreement is not an all-input proof.
Generic floating values are approximate group representatives at finite tolerance;
operations may reject accumulated drift. Work budgets are not OS memory bounds.

Rust coefficients inherit checked i128 overflow/resource panics; CLI catches
these as INCONCLUSIVE. C++ coefficients are arbitrary-precision rationals but
polynomial degree/term bounds still throw. No general CAS or parser is claimed.
CAS evaluators consume only fixed-symbol trusted producer outputs, not arbitrary
model-proposed strings.

A native Sage SO class accepts a reflection in the observed generic QQ path;
our oracle explicitly checks determinant independently. See SOURCE-AUDIT.md.

## Completed run — 2026-10-17

Evidence: ~/.cache/mathphysics-lie-groups/runs/20261017T163207Z/

- All gates above passed. SymPy 1.15.0.dev: 80 checks; Sage 10.8: 86 checks,
  each across both Rust and C++ producer records. Cadabra: three exact zero
  reductions and the wrong-sign residual 4*delta left nonzero.
- All 25 f64 cases passed, no exclusions; largest component/reference/parity
  residual 5.551115123125783e-16 against the stated 1e-12 bound.
- Rust binary16 128-fixture maximum matrix error versus same-input f64:
  0.0014193515119980127, below the 0.01 regression bound.
- Original QFT suite: 55 integration tests plus two compile-fail doctests;
  13 selected textbook entries; all 200 SMGA fixtures and live Cadabra checks.
- New derivation compiled with two LaTeX passes: no final undefined references,
  overfull boxes or LaTeX warnings.
- All 51 source/test/oracle/LaTeX/build-input hash entries matched after the
  run. Base HEAD was 14f43b1 plus the then-uncommitted feature tree; hashing,
  rather than base HEAD alone, identifies the code that was tested.
- rustc/cargo 1.97.1; GCC 11.4.0; CMake 3.31.2. Image IDs and source checkout
  commit are recorded in provenance.txt. Images ran without network.
- Only public Sage source files were fetched during the read-only audit;
  no new packages, cloned Sage repository, dependency downloads or background
  service were introduced.
- Existing unrelated mono-notes PDF/build-report modifications remain unstaged.
