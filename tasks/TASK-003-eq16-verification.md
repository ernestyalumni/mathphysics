# TASK-003 — Verify SMGA Eq.(16) for n = 3…10

Status: PARTIAL — n=3..10 integration tests exist; durable report and broader sampling are missing
Priority: P0 (feeds the YC demo, TASK-009)
Depends on: TASK-000
Estimated size: one session (2–4 h)
Branch suggestion: `feat/task-003-eq16`

## Goal

Compare the SMGA closed-form expression

    A_{1...n}|_{R₁} = (1/2^{n-2}) ∏_{m=2}^{n-1} (sg_{m,m+1} + sg_{1,2...m})

against a separately implemented version of the paper's Berends–Giele recurrence for
n = 3…10 at many random points inside R₁.
This is a useful engineering version of the "AI conjectured, verifier checked" loop. The
paper gives and proves an all-n formula, so finite tests through n = 10 do **not** extend
the physics result. They validate our implementation and create a reproducible regression
artifact for the product demo.

## Inputs

- SMGA paper source (see TASK-002 for path), especially the definitions of the sign
  functions sg_{i,j} and sg_{1,2...m} and the boundaries of R₁.
- Existing specialized recurrence and `tests/integration_tests/test_smga_reproduction.py`.

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

## Diagnostic hooks (report, don't market as discoveries)

- Does runtime or numerical behavior degrade with n?
- Do chamber-boundary samples expose implementation ambiguity or invalid-input handling?
- Does a deliberately perturbed expression reliably produce a nonzero residual?
