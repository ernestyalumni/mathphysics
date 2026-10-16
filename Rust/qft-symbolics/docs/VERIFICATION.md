# Verification receipt — 2026-10-16

This is a curated summary, not a replacement for the runnable checks.

## Reproduce

From `Rust/qft-symbolics/`:

```bash
bash validation/run.sh
```

Latest full run after review corrections:
`~/.cache/mathphysics-rust/runs/20261016T151728Z/`.
Earlier 15:12 and 15:15 runs are superseded by this one.
Each run stores commands' outputs, git SHA/branch/dirty status, tool versions,
source hashes and Cadabra image ID **outside the repository**.

At execution, base HEAD was
`f8c3cf2c8ff358c52922ad8645c5f337aaafb40b`,
branch `feat/rust-qft-verifier`, with the new source uncommitted.
The per-source hashes in `provenance.txt` bind that tested dirty source;
a later commit identifier alone must not be mistaken for the runtime source ID.

Environment: rustc/cargo 1.97.1; Python 3.13.3; NumPy 2.5.1.
Existing Cadabra image ID:
`sha256:48e6f90b4716e1e107deebe52da8b83e572b6ef6df36a1eabcf8fc72ba78ef50`.
The inspected Cadabra **source checkout** and this prebuilt image are separate
artifacts; their build equivalence was not established.

## Outcomes

- Formatting: PASS.
- Clippy, all targets, warnings denied: PASS.
- Debug: **25 integration tests + 2 compile-fail doctests**, all PASS.
- Release: same tests, all PASS; overflow checks enabled.
- Corpus crosswalk: 12 textbook entries point to actual listed tests.
- CLI: exact zero-residual PASS; deliberately wrong sign gives FAIL/exit 1;
  invalid domain/over-budget gives INCONCLUSIVE/exit 2.
- Exact SMGA recurrence/formula: **all 200 required fixtures** PASS,
  n=3..10, seeds 0..24, **133 zero and 67 nonzero**, no exclusions.
- Python parity: both old NumPy paths match those 200 Rust cases,
  absolute tolerance 1e-10.
- Live Cadabra: Schouten residual 0, wrong-sign residual nonzero,
  joined/canonicalised Clifford anticommutator residual 0.
  Real transformations were performed, not merely expression parsing.

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

## Scope boundaries

- Known results only; no novelty claim or held-out agent benchmark.
- Polynomial zero residuals are exact over their explicit domain.
- Finite component identities are restricted to the chosen representation.
- 47.2 is only tested at d=4, not dimensionally continued.
- 50.4 verifies the chiral trace identity but does not yet provide its fully
  evaluated invariant/epsilon expression.
- 60.2a has a five-leg symbolic expansion plus an explicit all-n linearity
  derivation, not an arbitrary constrained-ideal simplifier.
- Parke–Taylor symmetry/weight tests do not prove a Yang–Mills amplitude.
- The SMGA recurrence is a special half-collinear (2,2) construction, **not**
  ordinary Lorentzian gluon recursion. 200 cases are not an all-n proof.
- No full Python/Cadabra parity; no complete textbook chapter is marked done.
- Existing old ch60 script was separately observed to fail at its reference
  spinor Dirac equation; see MIGRATION.md.

## Inputs left untouched

The user's original textbook LaTeX, all old Python scripts, local Cadabra
checkout and mono-notes PDFs/build report were not modified by this work.
No external publishing, push, model API call, dependency download or system
installation was performed.
