# YM-T02 — Compact simple Lie algebras: structure constants, Killing form, Casimir constant (Rust)

Tier T2 · Workstream B · Depends on: PDD intent approved · Crate: `Rust/lie-groups` · Review: T2

## Behaviour wanted
For su(N), so(N), sp(N) in their defining matrix representations with exact rational (Gaussian-rational) arithmetic: an explicit basis, structure constants, the Killing form, an invariant form ⟨X,Y⟩ = −(1/m) Re tr(XY), and the constant C with tr(ad X ad Y) = −C ⟨X,Y⟩. A check that the algebra is simple is not required; a check that the Killing form is a multiple of the invariant form is.

## Fixtures
- su(2): dim 3, C = 8. su(3): dim 8, C = 18. General su(N): C = 2N² (test N = 2..6).
- Oracle: `PLAN/ledger/evidence/check_oneloop.py` second half; an independent hand value for so(3) ≅ su(2) under its own normalization (derive and state it in the test comment).
- Negative control: a non-simple algebra (u(2) or so(4)) must fail the "Killing ∝ invariant form" check.

## Process and acceptance
As YM-T01. Reuse existing matrix and exact-scalar types in `lie-groups`; do not add a second matrix type.
