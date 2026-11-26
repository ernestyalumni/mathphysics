# YM-T20: hosted models as referees on the frozen held-out cases

Tier: T2 (tool evaluation). Owner: pi on the desktop (runs), claude-opus55-strat (reading, versioning). Follows YM-T16 and YM-F01 experiment 3.

## Why

The local 9B and 14B referees locate about half of the planted errors and raise false alarms on correct proofs (`notes/YM-T16-weight-space-referee.md`; `notes/YM-F01-lean-formalisation.md`, experiment 3). Ernest asked how hosted models reached through TokenRouter do on the same proof-checking task.

## What

- The YM-T16 held-out case file, unchanged: 43 cases from 14 proofs (14 originals, 29 planted-error mutants), two fresh samples per case and model.
- Primary routes: `moonshotai/kimi-k3.1-pro`, `deepseek/deepseek-v4.1-flash-1017-free`, `z-ai/glm-5.3`. `z-ai/glm-5.3-flashx` ran as a supplementary route.
- Protocol fixed before the first full-run request: `ledger/evidence/YM-T20/PREREGISTRATION.md` and `PREREGISTRATION-GLM53.md`.

## Acceptance

- Every number in `notes/YM-T20-hosted-referee-replay.md` comes from `ledger/evidence/YM-T20/` (`samples.csv`, `summary.json`), which is regenerated from the private receipts pinned in `MANIFEST.md`.
- Transport failures, token-limit cut-offs and unparsed answers are reported apart from mathematical outcomes.

## Rules

- Case text and raw responses stay outside git: the N1/N2 text is not ours to publish.
- A referee verdict is evidence about the referee, never about the mathematics.
