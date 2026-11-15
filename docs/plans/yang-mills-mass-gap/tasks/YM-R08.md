# YM-R08 — The infrared end of the Lean: units IV.1, V.1, V.3

Tier T2–T3 · Workstream A/C · Depends on: none · Blocks: M2, M4, YM-G05 · Review: T3 · Budget: two sessions

## Inputs
`ARM/outline/units/IV-1.md`, `V-1.md`, `V-3.md`; `ARM/outline/work-items.md`; the Lean files they link, especially `LEAN/YM/EveryPhase/Decay.lean` and `LEAN/YM/Phase/Mixture.lean`.

## Do
1. For each unit: the list of hypotheses it consumes from upstream (as named Lean structures or statements) and the list of conclusions it exports.
2. Quote the exact formal statement of the polynomial low-energy bound (`Decay.lean:173–192`) and of vacuum uniqueness in every phase.
3. Identify every place where the argument would **change** if an exponential clustering bound were available as an extra hypothesis: which lemmas become unnecessary, which conclusions strengthen (X a point? full-sequence convergence?).
4. Write the hypothetical interface: a Lean-style structure `InfraredInput` (statement only, not compiled) whose fields are the extra bounds, such that item 3's strengthened conclusions would follow.

## Deliverable
`PLAN/notes/YM-R08-infrared-interface.md`.

## Acceptance
Item 1 is a two-column table per unit with declaration names. Item 4 type-checks informally: every identifier it mentions exists in the pinned source (give path:line) or is declared new.
