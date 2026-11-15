# YM-R05 — The unit-scale system: define the class 𝒞(β̄) from the Lean

Tier T3 · Workstream C · Depends on: none (YM-O01 helps but is not required; read source text) · Blocks: YM-G01, G0 · Review: T3, plus Lean-literate T2 · Budget: two sessions

## Goal
Write down exactly what the ultraviolet construction leaves at level 0: the density, the norms and bounds its pieces satisfy, its symmetries, and whether it is reflection positive. This is the object every infrared conjecture is about.

## Inputs
`ARM/outline/units/III-9.md`, `III-10.md`, `III-12.md`, `III-14.md`, `II-2.md`; `ARM/glossary.md`; `ARM/data/key-statements.tsv`; then the Lean files those outlines link (`LEAN/YM/RG/Activities.lean`, `PathActions.lean`, `History.lean`, `AnalyticBridge/VERIFIED_RESET_STEP/Evaluation.lean`, and what they import).

## Do
1. Give the level-m density ρ_m with every term defined: 𝒦, ℒ, ℰ, 𝒢, the vacuum constant, the small-field indicators.
2. List the bounds available at level 0, each with its Lean declaration name, path and line: regular activity norm, forced activity norm, coupling band, locality/decay of activities, analyticity.
3. Symmetries of ρ_0: gauge, translations, hyperoctahedral, reflections.
4. Reflection positivity: the outlines say it is proved for the **original** Wilson measure, not for an effective action. Confirm from the Lean. Say what, if anything, is known about ρ_0 in this respect.
5. Exact preservation of the partition function: find the declaration that states ∫ρ_m = Z_Wilson and quote it.
6. Propose a definition of 𝒞(β̄): the weakest set of properties from items 1–5 that the construction provides. Mark which properties are kernel-checked and which you inferred from prose.

## Deliverable
`PLAN/notes/YM-R05-unit-scale-class.md`.

## Acceptance
- Every bound in item 2 has `LEAN@3456e5e:path:line`.
- Item 4 ends in one of: "reflection positive (declaration …)", "not claimed", "unknown".
- The proposed definition distinguishes kernel-checked from inferred properties line by line.

## Traps
Names like `Treasury` or `Jewel` mean nothing. Use `ARM/glossary.md`, then unfold. The outline for III.13 describes an *alternative* family not used by the main theorem; do not mix its properties in.
