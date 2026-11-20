# YM-F01 — Lean formalisation of the conditional-gap chain

Tier: T2. Claimed by Claude Opus 5.5 on `ym/YM-F01-claude-opus55-lean`. Request: `request-2026-11-20-lean-formalisation.txt`.

## Goal
1. Choose the steps of `documents/mono-notes/topics/yang-mills-lattice.tex` that the conditional lattice gap depends on.
2. Formalise them in Lean 4 with Mathlib, run in Docker.
3. Test whether a local Lean prover on the RTX 3060 can prove them, with Lean as the only judge.

## Done when
- `formal/check.sh --strict` passes: no `sorry`, standard axioms only.
- The note `notes/YM-F01-lean-formalisation.md` lists every difference between the Lean statements and the source.
- The prover run, pre-registered in that note, is reported.

## Out of scope
- The trace formula and the transfer operator, which need measure theory on the configuration space.
- The defect recurrence, Conjecture `ym:conj:recurrence`, which is open.
- The complexification step of `casimir-sun`.
