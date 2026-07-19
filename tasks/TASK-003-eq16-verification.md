# TASK-003 — Verify SMGA Eq.(16) for n = 3…10

Status: TODO
Priority: P0 (feeds the YC demo, TASK-009)
Depends on: TASK-001
Estimated size: one session (2–4 h)
Branch suggestion: `feat/task-003-eq16`

## Goal

Independently verify the SMGA closed-form conjecture

    A_{1...n}|_{R₁} = (1/2^{n-2}) ∏_{m=2}^{n-1} (sg_{m,m+1} + sg_{1,2...m})

against direct Berends-Giele evaluation for n = 3…10, at many random points inside R₁.
This is exactly the "AI conjectured, symbolic/numeric engine verified" loop the discovery
engine is about — the paper's formula was conjectured by an LLM; we check it harder and
further (they may not have pushed to n = 10; going further is a talking point).

## Inputs

- SMGA paper source (see TASK-002 for path), especially the definitions of the sign
  functions sg_{i,j} and sg_{1,2...m} and the boundaries of R₁.
- Berends-Giele engine from TASK-001.

## Steps

1. Implement the sign functions exactly per the paper, with docstring citations.
2. Sampler for random kinematic points strictly inside R₁ (rejection-sample; assert the
   defining inequalities hold).
3. For n = 3…10, ≥ 100 points each: compare Eq.(16) to direct evaluation; record max
   relative error per n.
4. Also probe the boundary of R₁ (points approaching the defining surfaces) and record
   behavior — any structured deviation is a potential discovery lead; log it, don't hide it.
5. Emit machine-readable results: `amplitudes/results/eq16_verification.json`
   (per-n: points, max_rel_err, boundary notes) + human summary
   `amplitudes/results/eq16_verification.md`.

## Definition of done

`uv run python amplitudes/scripts/09_eq16_verification.py` (new script) exits 0 and writes
both result files; max relative error < 1e-8 for all n in the interior, or the deviation is
documented with reproduction seeds.

## Discovery hooks (report, don't chase yet)

- Does the formula degrade or stay exact as n grows?
- Anything special at n where new sign-function combinations first appear?
