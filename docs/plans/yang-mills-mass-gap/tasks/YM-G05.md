# YM-G05 — Warm-up lemma: clustering of the mixture collapses the phases

Tier T3 (proof) + T2 (Lean statement) · Workstream C · Depends on: YM-R08 · Review: T3 · Budget: one session

Assessment §6.2 says: if S^mix clusters, "a short calculation shows that ∫|S^ξ(A) − S^mix(A)|² dν = 0. Continuity and full support would then identify every component."

## Do
1. Write that calculation in full, with the exact form of clustering it needs and the exact properties of the phase decomposition it uses (cite the Lean declarations for S^mix = ∫S^ξ dν, continuity in ξ, full support of ν).
2. State it as a Lean theorem over the existing hierarchy structures. If YM-O01 is done, try to prove it in Lean in a scratch project that imports the pinned build. Report honestly how far it got.

## Why
It is the smallest true statement that goes beyond the submission, and it exercises the whole pipeline: read, state, prove, formalize, review.

## Deliverable
`PLAN/notes/YM-G05-clustering-collapses-phases.md`; any Lean under `REPOS/yang-mills-work/ymext/` (outside the pinned repos).
