# TASK-000 — Reproducible package and test entry point

Status: DONE (2026-07-25)
Priority: P0
Depends on: —
Estimated size: 1–2 hours
Branch suggestion: `feat/task-000-reproducible-baseline`

## Goal

Make one short, documented command run the focused MHV verifier from a fresh checkout. The current physics path works, but package metadata and dependency declarations make the default developer experience brittle.

## Baseline to preserve

On 2026-07-19:

- the focused amplitude suite passed 43 tests;
- `amplitudes/scripts/08_smga_stripped_amplitudes.py` reported `OVERALL: PASS`;
- the full suite had 62 passes and 2 failures from pre-existing uncommitted `massless_spinors.py` work;
- `uv run pytest` initially failed because `pyproject.toml` points to a missing package `README.md`.

Do not “fix” or stage the uncommitted spinor files as part of this task.

## Deliverables

1. `amplitudes/MHVamplitudes/README.md` explaining package scope, truth labels, and quick start.
2. Correct runtime/test dependency declarations so the quick-start command does not require ad hoc `--with` flags.
3. A registered `slow` pytest marker, or removal of the warning by an equally explicit mechanism.
4. Updated `amplitudes/MHVamplitudes/TESTING.md` with:
   - focused MHV command;
   - full-suite command;
   - the known unrelated WIP failure note;
   - Python version used.

## Definition of done

From `amplitudes/MHVamplitudes/` in a fresh environment:

```bash
uv run pytest -q tests/unit_tests/amplitudes \
  tests/integration_tests/test_smga_reproduction.py
```

passes with no unknown-marker warning. Record the exact output and date in `TESTING.md`.

## Out of scope

- Repairing `massless_spinors.py`.
- Editing the physics recurrence or closed-form implementation.
- Deleting the untracked `Cadabra2/` tree or build artifacts.

## Resolution (2026-07-25)

- Root cause: `pyproject.toml` declared `readme = "README.md"` but no package-level
  `README.md` existed, which made `uv`'s editable-install build fail outright (not just
  require extra flags). Added `amplitudes/MHVamplitudes/README.md`.
- Added `numpy` to `[project] dependencies` and `pytest` to a `dev` `[dependency-groups]`
  entry (uv includes the `dev` group by default), so `uv run pytest ...` needs no `--with` /
  `--no-project` flags.
- Registered the `slow` marker via `[tool.pytest.ini_options]` in `pyproject.toml`.
- Rewrote `TESTING.md` with the focused command, full-suite command, and the current (clean)
  state.
- Verified from a fresh checkout (`.venv`, `uv.lock`, `.pytest_cache` deleted first):
  focused suite `43 passed`, full suite `64 passed` — the previously-noted 2 WIP failures in
  `massless_spinors.py` no longer exist; that work was committed in `619b505` (2026-07-20).
