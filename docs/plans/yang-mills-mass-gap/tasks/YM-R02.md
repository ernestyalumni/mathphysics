# YM-R02 — Digest N2 §2: from a partition-function defect to a full transfer gap

Tier T3 · Workstream C · Depends on: none · Blocks: YM-G02, M1 · Review: T3 · Budget: one session

## Goal
Extract the finite-volume criterion that turns "a small partition-function defect on one box" into "a spectral estimate uniform in larger circumferences".

## Inputs
`N2/build/main.tex`, `N2/build/sections/finite-volume.tex` (check the file list for the right section file), `N2/paper.pdf` pp. 9–12; §3.12 for the Hilbert-space construction it refers to.

## Do
1. Digest template for §2.1 (transfer eigenvalues and a measure of excited-state weight) and §2.2 (increasing the spatial volume).
2. State the defect quantity exactly as defined, and the theorem(s) with all constants.
3. List every property of the model the argument uses: reflection positivity in which directions, symmetry between time and space periods, positivity of the transfer operator, compactness, anything else.
4. Compare with the elementary inequality in `PLAN/03-MASS-GAP-STRATEGY.md` §4 (D(n) ≥ Σ (λ_i/λ_0)^n). Say whether N2 uses it, a stronger form, or something different, and how §2.2 removes the volume dependence.
5. Section "Port to lattice gauge theory": for each property in item 3, say whether 4D Wilson lattice gauge theory on a torus has it, with a reference you actually opened (Seiler 1982 or Osterwalder–Seiler if available locally; otherwise mark `(from memory, unverified)`).

## Deliverable
`PLAN/notes/YM-R02-n2-defect-to-gap.md`.

## Acceptance
- The criterion is stated so that someone could implement a numerical test of it (YM-N03 will).
- Item 5 has a yes / no / unknown for every property.
