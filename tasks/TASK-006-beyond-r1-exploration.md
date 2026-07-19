# TASK-006 — Explore kinematic regions beyond R₁ (discovery)

Status: TODO
Priority: P1 (first genuine-discovery task)
Depends on: TASK-002, TASK-003
Estimated size: 1–2 sessions; open-ended — timebox each session
Branch suggestion: `feat/task-006-beyond-r1`

## Goal

SMGA established the single-minus closed form in one region, R₁. The obvious open move —
flagged in `../amplitudes/DISCOVERY-PLAN.md` Phase 3 — is: what happens in the other
half-collinear regions R_k, and can a closed form be found there?

## Method (conjecture → verify loop)

1. Enumerate the analogous regions (permute which subsets go collinear; follow the paper's
   region-labeling). Implement samplers per region with asserted defining inequalities.
2. Evaluate single-minus amplitudes numerically (TASK-001 engine) across each region for
   n = 4…8; tabulate values in the paper's stripped normalization.
3. Pattern-hunt: does an Eq.(16)-like product of sign functions fit? Fit candidate closed
   forms (products over sign-function combinations; small rational prefactors), then test
   each candidate at ≥ 100 fresh points per (region, n).
4. Every surviving conjecture gets: exact statement in LaTeX, verification script, and a
   `conjectured — numerically verified at N points, not proven` label. No stronger claims.

## Deliverables

- `amplitudes/scripts/10_beyond_r1.py` (sampler + evaluator + fitter + verifier).
- `amplitudes/results/beyond_r1.md` — per-region findings, including negative results
  ("no product formula of form X fits region R₂" is a real result; record it).
- If anything survives: a new section drafted in `amplitudes/07-open-directions.tex`.

## Definition of done

The script runs end-to-end and the results file exists with at least regions adjacent to R₁
covered for n = 4…6. Discovery is not required for "done" — an honest map of what was tried
and what fit/failed is the deliverable.
