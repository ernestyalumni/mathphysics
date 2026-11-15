# YM-G03 — The gauge fourth variation: is there a positive square?

Tier T3 · Workstream C · Depends on: YM-R01 · Blocks: milestone G3, route R1 · Review: T3 + red team. **Run twice independently** by two model families.

## Goal
Decide, by exact computation, whether the structure of N1 Proposition 2.2 survives when spins on sites become SU(2) elements on links and edge twists become plaquette twists.

## Setup
Finite 4D (start with 2D and 3D) periodic or open lattice, links U_ℓ ∈ SU(2), action S = β Σ_p (1 − ½ Re tr U_p), optional plaquette strengths b_p ∈ [0, β]. For a field B = (B_p) of su(2) elements define H(B) = Σ_p b_p ½ Re tr( U_p (exp(B_p) − 1) ) and W(B) = ⟨exp H(B)⟩, the ratio of twisted to untwisted partition functions. (Fix the exact form after YM-R01; the point is that it is the direct analogue of N1 eq. (2.1).) Twists that are lattice coboundaries of link fields, B = dA with the appropriate parallel transport, are changes of variables and give W = 1: this is the source of Ward identities.

## Do
1. **Phase A (symbolic, Python/SymPy oracle allowed).** Derive DᵏH for k = 1..4 at B = 0. Derive the second-variation Ward identity. Then the fourth variation along two colour directions X_a, X_b in the pattern (a, b, a, b), as in N1 eqs. (2.10)–(2.12). Identify the commutator field and whether the mixed fourth-order coefficient cancels as in (2.12).
2. Write the summed fourth derivative over a family of test fields as (candidate square) + (remainder). State exactly what is and is not sign-definite.
3. **Phase B (exact numerics).** On the smallest lattices where the quantities are nontrivial, evaluate both sides by exact or high-precision integration over SU(2) (character expansion, or quadrature with error control). Tabulate. A single clean counterexample to sign-definiteness is a valid, valuable result.
4. Power counting: for test families supported at dyadic frequencies between 1 and N inside a container at scale N, estimate the size of the commutator response and of the remainder. Does 4D give (log N)² against log N, as 2D does? Show the computation; do not assert.

## Deliverable
`PLAN/notes/YM-G03-gauge-fourth-variation-<session>.md`; scripts under `PLAN/ledger/evidence/YM-G03/`. When unit T11 exists, the Rust version replaces the oracle for item 1.

## Acceptance
- The O(3) case is reproduced first with the same code and matches N1 (2.10)–(2.12). No gauge result is accepted before this passes.
- Item 2 ends with one of: "square structure present (formula)", "absent (counterexample)", "undecided (what blocks)".

## Traps
Gauge transformations make many variations trivially zero; quotient them out before looking for structure. Ordering of non-commuting factors around a plaquette matters; state your orientation and base-point convention once in `PLAN/notes/conventions.md`.
