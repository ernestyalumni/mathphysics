# YM-T17: proof steps as numeric tests

Tier: T2 (tool) with T1 local models as test writers. Owner: claude-opus55-strat.

## Why

The overnight run (2026-11-17) showed that 9B local models are weak judges. As referees they flagged 94 of 333 planted errors (28%), and 0 of 6 flags checked by hand were real. As Lean-vs-summary comparers, 0 of 10 sampled INCONSISTENT verdicts were real (`notes/YM-T17-lean-compare-sample.md`). They are useful where a machine checks their output (verbatim quotes, compiled Lean). This task gives them that kind of job for proofs: write a test, and let the test decide.

## What

`bin/ym-numcheck`:

1. For one step of a proof unit, the model writes `check(rng)`, a numpy function that samples an instance satisfying the hypotheses and returns whether the step's claim holds.
2. The harness runs it on 200 seeds in Docker with no network, no capabilities, a read-only root, an unprivileged user, and memory, CPU and process limits.
3. Verdicts are PASS, FAIL, ERROR, NONE (nothing numerically testable), NOCODE and TIMEOUT.

Calibration uses the planted mutants from YM-T14 (`ym-micro`):

- A detection is a mutated step that FAILs while the same step of the unmodified unit PASSes.
- A false alarm is a FAIL on an unmodified step.

## Acceptance

- Calibration table: detected / missed / untestable over all mutants of the mono-notes units and C23, compared with the referee's 28%.
- False-alarm count over all 295 unmodified steps.
- Every FAIL on an unmodified step is read by a person before it is reported as an error in the notes.

## Rules

- A PASS proves nothing: random tests can miss, and a test can be wrong.
- A FAIL is a lead until someone reads the test.
- Model-written code runs only in the sandbox.
- If the approach holds up, the samplers (su(N), SU(N), Haar measure) belong in the Rust `lie-groups` crate, and the model writes tests against them.
