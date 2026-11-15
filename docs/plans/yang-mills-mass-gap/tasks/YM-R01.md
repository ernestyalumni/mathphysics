# YM-R01 — Digest N1: exponential decay in 2D O(n) models

Tier T3 · Workstream C · Depends on: none · Blocks: YM-G03, YM-G04, G0 · Review: T3, different family · Budget: one session

## Goal
A complete working digest of the 20-page paper that carries the infrared idea we want to port to gauge theory.

## Inputs
`N1/build/main.tex` and `N1/build/sections/{rotation,annulus,clusters}.tex`; `N1/paper.pdf` for page numbers.

## Do
1. Follow the digest template in `tasks/README.md` for §2 (fourth-order rotation estimate), §3 (annulus), §4 (sign clusters), and the reduction from n to 3 components.
2. Reproduce in full: the definitions of H(A), W(A) (eq. 2.1), the second-variation identity, eqs. (2.8), (2.10), (2.11), (2.12), and the statement of Proposition 2.2.
3. For Proposition 2.2 state exactly: which terms form the nonnegative square, which are the "controlled local terms", and what bounds them.
4. For §3 state exactly: the geometry of the annulus, what is pinned, the test-function family, where log N and (log N)² come from, and the final inequality for Z(π)/Z(0).
5. For §4 state exactly: the identity between the crossing probability and 1 − Z(π)/Z(0), which correlation inequalities are used (Ginibre, FKG, Campbell–Chayes) and what property of the model each needs.
6. Last section, "What is two-dimensional, what is spin-specific": list every step that uses planarity, every step that uses that the target is a sphere, every step that uses ferromagnetism.

## Deliverable
`PLAN/notes/YM-R01-n1-digest.md`.

## Acceptance
- Every numbered equation used is reproduced and matches the LaTeX source character for character in content.
- A reader could re-derive Prop. 2.2 for O(3) from the digest alone.
- The last section has at least one entry under each of its three headings, each tied to an equation or lemma number.

## Do not
Judge whether the paper is correct. Speculate about gauge theory beyond the last section's list.
