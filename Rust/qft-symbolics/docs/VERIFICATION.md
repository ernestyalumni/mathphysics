# Verification receipt — 2026-10-16

This is a curated summary, not a replacement for the runnable checks.

## Reproduce

From `Rust/qft-symbolics/`:

```bash
bash validation/run.sh
```

First-slice full run after review corrections:
`~/.cache/mathphysics-rust/runs/20261016T151728Z/`.
Latest two-slice full run:
`~/.cache/mathphysics-rust/runs/20261017T032014Z/` (UTC timestamp;
local session date October 16). Tested base HEAD: 39d96d3, with this second
slice still uncommitted. Source/test/oracle/LaTeX hashes bind the dirty source.
All gates passed. Generated LaTeX PDF/auxiliary files are in that run's
external latex/ directory, not committed.
Each run stores commands' outputs, git SHA/branch/dirty status, tool versions,
source hashes and Cadabra image ID **outside the repository**.

At that first-slice execution, base HEAD was
`f8c3cf2c8ff358c52922ad8645c5f337aaafb40b`,
branch `feat/rust-qft-verifier`, with the new source uncommitted.
The per-source hashes in `provenance.txt` bind that tested dirty source;
a later commit identifier alone must not be mistaken for the runtime source ID.

Environment: rustc/cargo 1.97.1; Python 3.13.3; NumPy 2.5.1.
Existing Cadabra image ID:
`sha256: 48e6f90b4716e1e107deebe52da8b83e572b6ef6df36a1eabcf8fc72ba78ef50`.
The inspected Cadabra **source checkout** and this prebuilt image are separate
artifacts; their build equivalence was not established.

## Outcomes

- Formatting: PASS.
- Clippy, all targets, warnings denied: PASS.
- Debug: **54 integration tests +2 compile-fail doctests**, all PASS.
- Release: same tests, all PASS; overflow checks enabled.
- Corpus crosswalk: 13 textbook entries point to actual listed tests.
- CLI: exact zero-residual PASS; deliberately wrong sign gives FAIL/exit 1;
  invalid domain/over-budget gives INCONCLUSIVE/exit 2.
- Exact SMGA recurrence/formula: **all 200 required fixtures** PASS,
  n=3..10, seeds 0..24, **133 zero and 67 nonzero**, no exclusions.
- Python parity: both old NumPy paths match those 200 Rust cases,
  absolute tolerance 1e-10.
- Live Cadabra: Schouten residual 0, wrong-sign residual nonzero,
  joined/canonicalised Clifford anticommutator residual 0.
  Real transformations were performed, not merely expression parsing.
- Typed tensor layer: 20 tests, including12 computed shared residuals. Live
  Cadabra agrees on all 12: nine zero, three known-nonzero controls.
- Formal-d47.2: two exact symbolic normal-word identities and independent
  Cadabra joins/contractions; wrong-sign control nonzero. Exhaustive explicit
  matrix checks cover498 assignments of lengths 0..4 in d=1..4 subspaces.
- Fully evaluated 50.4: exact polynomial identities for 16 independent vector
  components and for massless bracket products; epsilon sign, proper boost
  and parity controls. See TEXTBOOK-EXTENDED.md.
- Angular50.2a: general written factorization derivation and 48 exact rational
  executable witnesses, direct spinor-phase checks and invalid-input controls.

## Review corrections before the final run

A read-only reviewer identified two acceptance risks:

1. Public Ratio fields allowed bypassing the constructor with a zero
   denominator, potentially producing a zero cross-multiplication residual
   for undefined 0/0. Fields are now private; read-only accessors preserve
   the nonzero-polynomial invariant. Two compile-fail doctests prevent
   direct construction or denominator mutation.
2. Demo/test logic originally counted any SMGA error as an excluded wall.
   It now requires **every one of the 200 specified fixtures** to succeed;
   any domain error is a failure/inconclusive result, not a skipped point.
   Intentional singular/wall cases remain separate negative tests.

These fixes do not retroactively strengthen old scripts' verification.

Second-slice read-only reviews found no concrete algebraic false-PASS defect.
Corrections from those reviews:
- Charge symmetry variants before eagerly allocating/cloning a large orbit;
  rank-eight symmetric/exhausted-budget regression added.
- Replace a repetitive component sampling pattern with all 498 assignments.
- Reject out-of-range angular components before squaring (including a huge
  invalid-angle regression), document checked-i128 arithmetic limitations,
  and compare phi/star directly so a phase mistake cannot hide in momentum.

## Scope boundaries

- Known results only; no novelty claim or held-out agent benchmark.
- Polynomial zero residuals are exact over their explicit domain.
- Finite component identities are restricted to the chosen representation.
- Formal-d47.2 now has symbolic and Cadabra verification, not a general
  dimensional-regularization or d-dimensional gamma5 scheme.
- 50.4 includes the fully evaluated invariant/epsilon expression.50.2a has
  an explicit derivation and exact witnesses; 50.2b remains absent.
- Tensor algebra is commuting, mono-term and bounded; no multi-term Bianchi
  reduction, Grassmann/gamma factor ordering or spinor epsilon lowering.
- Work budgets are not OS memory limits. Legacy coefficient/Clifford overflow
  still panics, caught by the CLI rather than uniformly structured in the API.
- 60.2a has a five-leg symbolic expansion plus an explicit all-n linearity
  derivation, not an arbitrary constrained-ideal simplifier.
- Parke–Taylor symmetry/weight tests do not prove a Yang–Mills amplitude.
- The SMGA recurrence is a special half-collinear (2, 2) construction, **not**
  ordinary Lorentzian gluon recursion. 200 cases are not an all-n proof.
- No full Python/Cadabra parity; no complete textbook chapter is marked done.
- Existing old ch60 script was separately observed to fail at its reference
  spinor Dirac equation; see MIGRATION.md.

## Inputs left untouched

The user's original textbook LaTeX, all old Python scripts, local Cadabra
checkout and mono-notes PDFs/build report were not modified by this work.
No external publishing, push, model API call, dependency download or system
installation was performed.
