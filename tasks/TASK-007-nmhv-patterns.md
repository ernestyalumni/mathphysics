# TASK-007 — NMHV / multi-minus pattern search (discovery)

Status: TODO
Priority: P2
Depends on: TASK-002 (engine + regime machinery)
Estimated size: 1–2 sessions; timebox
Branch suggestion: `feat/task-007-nmhv`

## Goal

Look for closed-form or structured behavior of NMHV (three-minus) and multi-minus
amplitudes in half-collinear regimes — the natural "one step harder" ladder after the
single-minus results. Also the first place graviton generalization could enter (the SMGA
paper mentions it directly).

## Method

1. In the R₁ regime machinery from TASK-002/003, evaluate NMHV configurations for
   n = 5…8 (all inequivalent minus-position placements, using cyclic/reflection symmetry
   to reduce the set).
2. Normalize the same way SMGA strips their amplitudes; tabulate.
3. Pattern-hunt with the same candidate family as TASK-006 (sign-function products,
   ratios of sign-function products, low-degree rational functions of the regime
   parameters). Fit on one point set, verify on a disjoint set.
4. Optional stretch: repeat one small case (n = 5) for gravitons via KLT/double-copy from
   the gluon engine, if the gluon side is solid.

## Deliverables

- `amplitudes/scripts/11_nmhv_patterns.py`
- `amplitudes/results/nmhv_patterns.md` — tables + candidate fits + verification stats,
  negative results included.

## Definition of done

n = 5, 6 NMHV tables exist with at least one candidate-form fit attempted and its
verification outcome recorded. Anything promising gets escalated into its own task file
(create `TASK-0XX-*.md` and link it from `tasks/README.md`).
