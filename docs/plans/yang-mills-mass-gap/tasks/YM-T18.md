# YM-T18: verified samplers for numeric tests, and agreement across samples

Tier: T2 (tool). Owner: claude-opus55-strat. Branch: `ym/YM-T18-claude-opus55-samplers`. Follows YM-T17.

## Why

In YM-T17 most wrong tests went wrong in building the instance, not in stating the claim: a wrong su(N) construction, a spectral measure that breaks the hypotheses, a lattice built wrongly. And a single test per step gave about one false alarm per detection. Two remedies:

1. Give the test writer instances from code that is itself tested, so it only writes the claim.
2. Call a step FAIL only when several independently written tests agree.

## What

- `bin/ymlib.py`: samplers for su(N), u(N), Hermitian matrices, Haar U(N) and SU(N), the generators T_a, a matrix exponential, probability vectors, a finite Hamiltonian with vacuum and gap together with its spectral measure, and positive semidefinite matrices.
- `bin/ymlib_selftest.py`, run with `ym-numcheck selftest` inside the sandbox:
  - exact identities, including Gell-Mann orthonormality, completeness and the Casimir;
  - Haar moments with known values, among them E (tr U)^N = 1 on SU(N) against 0 on U(N), which tells the two measures apart.
- `ym-numcheck run --ymlib` mounts the library and lists it in the prompt; `--steps FILE` tests an explicit list of steps; `score --vote K` calls a step FAIL only when at least K seeds' tests fail.

## Acceptance

- `ym-numcheck selftest` passes in both sandbox images (numpy 2.4.1; numpy 1.26.4 with scipy 1.11.4).
- A calibration with `--ymlib`, compared with the same model and seeds without it (`notes/YM-T17-day-runs.md`).

## Rules

As YM-T17: a PASS proves nothing; a FAIL is a lead until a person reads the test; model code runs only in the sandbox. A sampler that the self-test does not cover is not added to the prompt.
