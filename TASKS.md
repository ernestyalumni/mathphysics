# Current Execution Queue — TMP / MHV Discovery Engine

**Last updated:** 2026-07-19
**Deadline driving this queue:** application-ready demo by 2026-07-25
**Detailed task packets:** [`tasks/README.md`](tasks/README.md)

## Do now, in order

1. **TASK-000 — baseline reproducibility**
   Make the focused MHV test command work from a fresh checkout and document the two unrelated WIP failures separately.

2. **TASK-009 — two-minute verification demo**
   Produce one command and one static report that show closed form vs independent recurrence through `n = 10`.

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
Focused amplitude suite: 43 passed, 1 marker warning
Full suite:              62 passed, 2 failures in uncommitted massless-spinor WIP
Script 08:               OVERALL: PASS
```

Reproduction command used on 2026-07-19:

```bash
uv run --no-project --with numpy env \
  PYTHONPATH=amplitudes/MHVamplitudes/src \
  python amplitudes/scripts/08_smga_stripped_amplitudes.py
```

Focused test command used on 2026-07-19:

```bash
cd amplitudes/MHVamplitudes
uv run --no-project --with numpy --with pytest env PYTHONPATH=src \
  pytest -q tests/unit_tests/amplitudes \
  tests/integration_tests/test_smga_reproduction.py
```

## Agent rules

- Claim one task at a time using the protocol in `tasks/README.md`.
- Never commit on `master`.
- Do not stage the pre-existing edits in `spinors/invariants.py` or `spinors/massless_spinors.py` unless Ernest explicitly assigns that work.
- A generated report must include command, git SHA, branch, seed, backend, and pass/fail criteria.
- Publishing, arXiv submission, application submission, and external messages require Ernest's approval.
