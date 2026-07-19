# TASK-006 — Map half-collinear chambers beyond the R₁ decay channel

Status: TODO AFTER SUBMISSION
Priority: P1 (first bounded discovery experiment)
Depends on: TASK-002, TASK-003, TASK-004
Estimated size: 1–2 time-boxed sessions
Branch suggestion: `experiment/task-006-half-collinear-chambers`

## Goal

The paper's `R_k` regions obtained by choosing a different single ingoing particle are related to `R₁` by cyclic relabeling; treating those as a new result would be wrong. This task identifies allowed half-collinear chambers that are not merely cyclic images of the single-decay channel, then maps the existing recurrence there.

The paper also says a longer general formula outside `R₁` will appear elsewhere, so checking the literature and subsequent author work is part of the task.

## Method

1. Re-read the paper's chamber and region definitions. Record which cases are cyclic images of `R₁` and exclude them from novelty claims.
2. Search subsequent primary literature and author releases for the promised general formula.
3. Define a bounded family of other allowed sign/chamber configurations and implement samplers with asserted inequalities.
4. Evaluate single-minus amplitudes for `n = 4..8` in the paper's stripped normalization.
5. Predeclare a small candidate-expression family, fit on training seeds, and test survivors on at least 100 fresh points per chamber and multiplicity.
6. Label every survivor `conjectured — numerically verified at N held-out points, not proven`.

## Deliverables

- `amplitudes/experiments/half_collinear_chambers/README.md` with scope, literature status, candidate family, and stop condition.
- One sampler/evaluator/verifier script.
- `amplitudes/results/beyond_r1.md` with:
  - cases excluded as cyclic relabelings;
  - cases already covered by later literature;
  - positive and negative fit results;
  - seeds, tolerances, and run metadata.
- A discovery log only if a real proposer/critic/verifier loop ran.

## Definition of done

The script runs end to end and the report covers at least two nontrivial chamber families for `n = 4..6`. Discovery is not required: an honest map of what was known, equivalent, tried, fitted, and rejected is the deliverable.

## Stop condition

Stop after two sessions or 20 candidate families, whichever comes first. Re-scope before continuing.
