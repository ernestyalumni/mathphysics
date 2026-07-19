# TASK-007 — NMHV / multi-minus pattern search

Status: DEFERRED UNTIL AFTER 2026-07-27
Priority: P2
Depends on: TASK-002, TASK-004
Estimated size: multiple time-boxed sessions
Branch suggestion: `experiment/task-007-nmhv-patterns`

## Goal

Determine whether the verifier-and-pattern-search workflow transfers from the single-minus half-collinear problem to the smallest nontrivial multi-minus configurations. A useful negative map is acceptable; novelty is not required.

## First bounded experiment

1. Choose exactly one helicity family and multiplicities `n = 5, 6`.
2. Write down conventions and the known reference expression before generating candidates.
3. Generate reproducible kinematic samples with held-out seeds.
4. Fit only a predeclared, small candidate-expression class.
5. Test survivors on held-out points and against soft/collinear limits.
6. Run a literature search before using the word “new.”

## Deliverables

- `amplitudes/experiments/nmhv/README.md` with hypothesis class and stop condition.
- One runnable experiment script.
- `amplitudes/results/nmhv-first-map.md`, including failed candidate families.
- A discovery log only if an actual proposer/critic/verifier loop ran.

## Definition of done

The experiment can be reproduced from a fixed command, all seeds and tolerances are logged, and the report clearly separates known reference results, fitted conjectures, held-out verification, and unresolved cases.

## Stop condition

Stop after two sessions or 20 candidate families, whichever comes first. Re-scope before spending more time.
