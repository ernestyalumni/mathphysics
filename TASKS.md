# Current Execution Queue — TMP / MHV Discovery Engine

**Last updated:** 2026-07-25
**Deadline driving this queue:** YC application deadline 2026-07-27, 8:00 PM PT
**Detailed task packets:** [`tasks/README.md`](tasks/README.md)

## Do now, in order

1. ~~**TASK-000 — baseline reproducibility**~~ DONE 2026-07-25. `uv run pytest` now works
   from a fresh checkout with no ad hoc flags: focused suite 43 passed, full suite 64 passed
   (the previously-noted WIP failures are gone — committed in `619b505`).

2. ~~**TASK-009 — two-minute verification demo**~~ DONE 2026-07-25. See
   [`amplitudes/DEMO.md`](amplitudes/DEMO.md) and
   `amplitudes/scripts/09_verification_demo.py`.

3. **TASK-002 / TASK-003 — durable result reports**
   Convert the already-working checks into human- and machine-readable artifacts with seeds, versions, and explicit scope.

4. **Application review**
   `yc/` contains public-safe product language and video rehearsals. Authenticated form state and private founder/company facts stay outside this public research repo.

## Parallel only if it cannot delay the critical path

- **TASK-005 — verifier wrapper:** useful if it can be demonstrated quickly with a known-true and known-false claim.
- **TASK-008 — public research/product README:** polish after the demo behavior is frozen.

## After submission

- **TASK-004:** amplitude consistency checks.
- **TASK-006:** adjacent-region exploration.
- **TASK-010:** spinor-helicity eval suite v0.
- **TASK-007:** NMHV pattern search.

## Current evidence

```text
Focused amplitude suite: 43 passed, no warnings
Full suite:              64 passed, 0 failures
Script 08:               OVERALL: PASS
Demo (script 09):        correct candidate PASS n=3..10; broken candidate correctly FAILs
```

Reproduction command, verified 2026-07-25 (no ad hoc `--with`/`--no-project` flags needed
anymore — see `tasks/TASK-000-baseline-reproducibility.md`):

```bash
cd amplitudes/MHVamplitudes
uv run pytest -q tests/unit_tests/amplitudes tests/integration_tests/test_smga_reproduction.py
```

Two-minute demo command, from repo root:

```bash
uv run --project amplitudes/MHVamplitudes python amplitudes/scripts/09_verification_demo.py
```

## Agent rules

- Claim one task at a time using the protocol in `tasks/README.md`.
- Never commit on `master`.
- Do not stage the pre-existing edits in `spinors/invariants.py` or `spinors/massless_spinors.py` unless Ernest explicitly assigns that work.
- A generated report must include command, git SHA, branch, seed, backend, and pass/fail criteria.
- Publishing, arXiv submission, application submission, and external messages require Ernest's approval.
