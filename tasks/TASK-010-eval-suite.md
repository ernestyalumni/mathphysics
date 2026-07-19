# TASK-010 — Spinor-helicity eval suite v0

Status: TODO AFTER DEMO
Priority: P1
Depends on: TASK-005
Estimated size: 4–6 hours
Branch suggestion: `feat/task-010-spinor-eval-v0`

## Goal

Build the smallest benchmark that measures whether an AI model can produce machine-verifiable spinor-helicity reasoning rather than plausible text.

## Scope

Exactly 10 public problems:

- 2 antisymmetry/sign-convention checks;
- 2 Schouten-identity manipulations;
- 2 little-group weight checks;
- 2 momentum-conservation identities;
- 1 soft-limit problem;
- 1 deliberately underspecified problem where admitting uncertainty is correct.

Keep any future held-out set outside public prompts; do not pretend these first 10 are a statistically meaningful leaderboard.

## Deliverables

- `evals/spinor-helicity-v0/README.md`
- `evals/spinor-helicity-v0/problems.jsonl`
- `evals/spinor-helicity-v0/verifier.py`
- `evals/spinor-helicity-v0/harness.py`
- verifier unit tests containing known-pass, known-fail, and inconclusive cases

## Scoring

Report:

- exact verifier passes;
- confident wrong answers;
- admitted unknowns;
- inconclusive verifier outcomes;
- prompt/model/verifier versions.

## Definition of done

The verifier tests pass locally, a mock model adapter can exercise the harness end to end, and the README states that no model baseline has been measured unless one actually ran.
