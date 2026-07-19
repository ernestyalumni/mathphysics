# TASK-010 — Physics-hallucination eval suite v0

Status: TODO
Priority: P1 (this is the most commercially legible artifact — see `../yc/README.md`)
Depends on: TASK-005 (verifier wrapper); TASK-004 helps
Estimated size: 1–2 sessions
Branch suggestion: `feat/task-010-eval-suite`

## Goal

v0 of an eval suite that measures whether an LLM produces *silently wrong* physics —
wrong signs, index placement, tensor structure, dropped terms, invalid limits — with every
item machine-gradable by the symbolic/numeric verifiers in this repo. This is the
"frontier labs pay for hard verified evals" wedge in the TMP mission.

## Design

- Item format (one YAML/JSON per item in `evals/items/`):
  `{id, domain, prompt, ground_truth_expr, verifier: {backend: sympy|cadabra2|numeric, script},
  distractor_notes, difficulty, source: path-or-citation}`.
- Grading: model answer → parsed expression → verifier computes (answer − truth) and
  simplifies / evaluates numerically at random points → pass/fail. Parsing failures are
  scored separately from wrong-physics failures.
- **Every ground truth must itself be verified** by a repo script before an item ships.

## v0 scope: 25 items

- 8 × spinor-helicity identities and ⟨⟩[] manipulations (source: Srednicki 48–50 +
  existing `amplitudes/scripts/03_spinor_product_identities.py`).
- 6 × amplitude values / limits (Parke-Taylor cases, soft/collinear leading behavior).
- 6 × sign/convention traps (mostly-plus vs mostly-minus translations, conjugation,
  ordering reversals) — the classic silent-failure zone.
- 5 × "is this derivation valid?" items with a planted flawed step; ground truth is the
  location of the flaw, gradable because the flawed identity fails the verifier.

## Deliverables

- `evals/README.md` (format, grading protocol, how to add items, honesty rules —
  no items whose ground truth wasn't machine-verified).
- `evals/items/` with 25 verified items.
- `evals/run_eval.py --model <name>` — runs via the Claude API (default
  `claude-sonnet-5`; get current model ids from the `claude-api` skill when implementing),
  emits per-domain pass rates to `evals/results/<model>-<date>.json`.
- A baseline run on at least one model, committed.

## Definition of done

All 25 ground truths pass their own verifiers; one baseline model run committed; README
documents the scoring so a stranger could re-run it.
