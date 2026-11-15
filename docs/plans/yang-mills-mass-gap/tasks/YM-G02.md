# YM-G02 — Defect ⇒ gap lemma for reflection-positive lattice gauge theory

Tier T3 · Workstream C · Depends on: YM-R02 · Blocks: milestone G2 · Review: T3 + red team · Budget: two sessions

## Goal
A written proof, and a Lean-style statement, of: for the Wilson lattice gauge theory on a torus (and for which members of 𝒞(β̄)?), a bound on a partition-function defect on one box implies a lower bound on the transfer-matrix gap uniform in larger spatial boxes.

## Do
1. Check the elementary inequality of `PLAN/03-MASS-GAP-STRATEGY.md` §4 line by line. It claims: for positive trace-class T with eigenvalues λ_0 ≥ λ_1 ≥ … ≥ 0, D(n) := (Tr Tⁿ)²/Tr T²ⁿ − 1 ≥ Σ_{i≥1}(λ_i/λ_0)ⁿ. Confirm or correct. Test numerically on random positive matrices (Python is fine; keep the script as evidence).
2. Establish what is needed to apply it: existence, positivity and self-adjointness of the transfer operator for the Wilson action on a periodic lattice. Cite a source you opened.
3. Port N2 §2.2: the step that makes the estimate uniform as the spatial volume grows. Write the gauge-theory version. Mark every step that does not transfer with `GAP:`.
4. State the result in Lean syntax against Mathlib (not compiled is acceptable; say so).

## Deliverable
`PLAN/notes/YM-G02-defect-gap-lemma.md`; evidence script under `PLAN/ledger/evidence/YM-G02/`.

## Acceptance
Item 1 has a verdict and a script whose output is quoted. Item 3 either is a complete argument or lists its gaps explicitly. The label on the final lemma is `conjectured` unless a human or Lean has checked it.
