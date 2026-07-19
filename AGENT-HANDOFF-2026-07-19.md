# Agent Handoff — Verified MHV Demo Sprint

**Date:** 2026-07-19
**Branch:** `feat/yc-mhv-execution-2026-07`
**Start here:** `AGENTS.md` → `RESEARCH-PLAN.md` → `TASKS.md` → `tasks/README.md`

## Mission for this sprint

Turn the existing MHV/SMGA reproduction into a credible, repeatable demonstration of the core product loop:

```text
candidate expression → independent computation → verifier verdict → persisted evidence
```

This sprint packages known results. It does not claim original physics.

## What already works

- Closed form: `amplitudes/MHVamplitudes/src/mhvamplitudes/amplitudes/smga.py`
- Specialized recurrence: `amplitudes/MHVamplitudes/src/mhvamplitudes/amplitudes/berends_giele.py`
- Region sampler: `amplitudes/MHVamplitudes/src/mhvamplitudes/kinematics/phase_space.py`
- Integration comparison through `n = 10`: `amplitudes/MHVamplitudes/tests/integration_tests/test_smga_reproduction.py`
- Standalone verification narrative: `amplitudes/scripts/08_smga_stripped_amplitudes.py`

On 2026-07-19, the focused amplitude suite passed 43 tests and Script 08 reported `OVERALL: PASS`.

## Known worktree state—preserve it

These changes existed before this handoff and are not part of the demo sprint:

- modified `amplitudes/MHVamplitudes/src/mhvamplitudes/spinors/invariants.py`
- modified `amplitudes/MHVamplitudes/src/mhvamplitudes/spinors/massless_spinors.py`
- untracked `.DS_Store` files
- untracked `Cadabra2/` tree containing build/dependency artifacts

Do not delete, rewrite, stage, or “clean up” those paths. The two full-suite failures currently come from the uncommitted massless-spinor implementation calling `GammaMatrices()` without its required argument.

## First task for a new agent

Claim `tasks/TASK-000-baseline-reproducibility.md`. If it is already claimed, claim `TASK-009` or a dependency-free documentation task. Do not start `TASK-006` or `TASK-007` during the application critical path.

## Acceptance standard

A demo is ready only when a reviewer can:

1. clone or open the repo;
2. run one documented command;
3. see the two computation paths and their agreement;
4. inspect a persisted report containing scope and run metadata;
5. understand that this is a verified reproduction, not a novelty claim.

## Coordination

- Change task `Status:` before starting and when finishing.
- Keep each commit limited to one task packet.
- Put public research/demo artifacts here. The `yc/` directory is a sanitized companion only.
- Keep authenticated form links, personal priorities, equity/employment details, users/revenue facts, and other founder-sensitive material out of this repo.
