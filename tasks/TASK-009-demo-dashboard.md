# TASK-009 — Two-minute verification demo

Status: DONE (2026-07-25)
Priority: P0
Depends on: TASK-000; consumes reports from TASK-002/TASK-003 when available
Estimated size: 3–5 hours
Branch suggestion: `feat/task-009-verification-demo`

## Goal

Create a demo that makes the central product behavior obvious in under two minutes:

```text
known candidate formula
        ↓
independent specialized recurrence
        ↓
comparison across multiplicities and seeded points
        ↓
PASS / FAIL / INCONCLUSIVE ledger with residuals
```

The demo must remain useful if no web UI is built. Prefer a crisp CLI and static HTML/Markdown report over a fragile dashboard.

## Required demo path

1. One command from the repo root.
2. Print the candidate formula and provenance.
3. Run or load a fresh comparison through `n = 10`.
4. Show at least one deliberately wrong candidate returning `FAIL`; a verifier that only displays green checks is not persuasive.
5. Save a report containing timestamp, git SHA, branch, backend, seeds, number of points, and residual summary.
6. End with an explicit label: “verified reproduction of known result; no novelty claim.”

## Deliverables

- `DEMO.md` — presenter and reviewer instructions.
- A thin CLI or script that orchestrates existing code without duplicating physics logic.
- `amplitudes/results/demo-latest.md` or static HTML generated from structured results.
- A 90-second screen-recording shot list; application-sensitive narration belongs in the private application repo.

## Definition of done

A new reviewer can run the command and see both a passing known expression and a failing perturbed expression. Total interactive time is under two minutes on the development Mac.

## Anti-overbuild rule

Do not build authentication, accounts, databases, deployment, or a general dashboard during this task. A local demo is enough for the current milestone.

## Resolution (2026-07-25)

- `amplitudes/scripts/09_verification_demo.py` — thin CLI orchestrating the existing
  `mhvamplitudes.amplitudes.smga` / `berends_giele` / `kinematics.phase_space` modules; no
  physics logic duplicated. Command: `uv run --project amplitudes/MHVamplitudes python
  amplitudes/scripts/09_verification_demo.py` from repo root. Runtime ~0.13s.
- Correct-candidate path: closed form vs. recursion agree exactly (`0.00e+00` residual) for
  `n = 3..10`, seeded and reproducible.
- Deliberate-failure path: a "broken" candidate (wrong normalization exponent, `2**(n-1)`
  instead of `2**(n-2)`, built from the real `sg_ij`/`sg_i_set` primitives) is run for
  `n = 4, 5, 6` and correctly FAILs every case — restricted to small `n` because the true
  amplitude becomes increasingly likely to be exactly 0 at larger `n` (more sign factors in
  the product), which would spuriously "pass" a broken candidate too; this is a demo-design
  choice, not a change to any physics code.
- Report: `amplitudes/results/demo-latest.md`, regenerated every run, with timestamp, git
  SHA, branch, backend versions, seeds, trial counts, and a residual table.
- `amplitudes/DEMO.md` — presenter/reviewer instructions.
- `yc/demo-video-script.md` updated to match the real command and output (previously blocked
  on this task).
