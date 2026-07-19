# tasks/ — Agent-Sized Task Board for the TMP Discovery Engine

**Mission:** AI-driven discovery engine for theoretical & mathematical physics. First target
area: MHV / single-minus gluon scattering amplitudes, benchmarked against arXiv:2602.12176v2
(SMGA). See `../TMP-MISSION.md` and `../amplitudes/DISCOVERY-PLAN.md` for the big picture.

This directory breaks the mission into tasks small enough for **one agent session each**
(2–6 focused hours). Each `TASK-NNN-*.md` file is self-contained: an agent should be able to
read only that file (plus the files it points to) and produce the deliverable.

## Protocol for agents

1. **Branch:** never work on `master`. Use `feat/task-NNN-<slug>` or continue an existing
   feature branch. Verify with `git branch --show-current` before every commit.
2. **Claim:** set `Status:` in the task file to `IN PROGRESS (<agent>, <date>)` in your first
   commit; set to `DONE (<date>)` with a pointer to the deliverable when finished.
3. **Verify before claiming done:** every task has a "Definition of done" section with a
   runnable check. If the check doesn't pass, the task is not done — report what failed.
4. **Environment:** Python via `uv venv --python 3.13` (never system Python 3.9). The
   numeric package lives at `amplitudes/MHVamplitudes/` (`pyproject.toml`, `src/`, `tests/`).
5. **Conventions:** Srednicki mostly-plus metric. Cadabra2-first for symbolic claims where
   feasible; SymPy is acceptable and is the default for SMGA numerics (decided 2026-04-23).
   LaTeX is the canonical writeup artifact; PDFs are derived.
6. **No overclaiming:** results are "verified" only if a script in the repo reproduces them.
   "Conjectured" and "verified" must never be conflated in any writeup.

## Task index

| Task | Title | Depends on | Priority |
|------|-------|-----------|----------|
| [TASK-001](TASK-001-berends-giele-numeric.md) | Berends-Giele numeric recursion, hardened + tested | — | P0 |
| [TASK-002](TASK-002-smga-reproduction.md) | Reproduce SMGA stripped amplitudes A_123 … A_123456 | 001 | P0 |
| [TASK-003](TASK-003-eq16-verification.md) | Verify SMGA Eq.(16) for n = 3…10 | 001 | P0 |
| [TASK-004](TASK-004-consistency-checks.md) | Cyclicity, KK relations, soft & collinear limits | 002 | P1 |
| [TASK-005](TASK-005-cadabra2-verifier.md) | Cadabra2 verifier wrapper + known-answer tests | — | P1 |
| [TASK-006](TASK-006-beyond-r1-exploration.md) | Explore kinematic regions beyond R₁ (discovery) | 002, 003 | P1 |
| [TASK-007](TASK-007-nmhv-patterns.md) | NMHV / multi-minus pattern search (discovery) | 002 | P2 |
| [TASK-008](TASK-008-paper-skeleton.md) | arXiv paper skeleton wired to computed results | 002, 003 | P1 |
| [TASK-009](TASK-009-demo-dashboard.md) | Live demo artifact: verification dashboard | 003 | P0 (YC) |
| [TASK-010](TASK-010-eval-suite.md) | Physics-hallucination eval suite v0 | 005 | P1 |

P0 tasks are on the critical path for the **YC application demo (deadline 2026-07-27)** —
see `../yc/README.md`.

## Relationship to older boards

`../TASKS.md`, `../RESEARCH-PLAN.md`, and `../amplitudes/TASKS.md` (talk prep, April 2026)
predate this board. Where they conflict, **this board wins**; the old files remain as context.
