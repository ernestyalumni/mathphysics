# TASK-009 — Two-minute verification demo

Status: TODO
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
