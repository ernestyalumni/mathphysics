# YM-T11 — Calculus on words of group elements (Rust)

Tier T3 · Workstream B · Depends on: YM-R01; PDD intent approved · Crate `Rust/qft-symbolics` (extends the ordered-algebra module) · Review: T3

## Behaviour wanted
Expressions that are traces of ordered products of group-valued symbols and exponentials of Lie-algebra-valued symbols, e.g. tr(U₁ exp(B₁) U₂ exp(B₂) …). Operations: nth directional derivative at B = 0 along given algebra directions (exact, via the ordered exponential series); substitution of su(2) relations (X_a X_b = −δ_ab/4·1 + ½ ε_abc X_c in a stated normalization) and so(3) relations; collection of terms by commutator structure.

## Fixtures
- O(3): reproduce N1 eqs. (2.10), (2.11), (2.12), including the cancellation of the mixed fourth-order coefficient for the pattern (a, b, a, b). Take the equations from `N1/build/sections/rotation.tex`.
- d/dt tr(U exp(tX)) at 0 = tr(UX); second derivative = tr(UX²); for su(2) with |X| = 1 in the stated normalization this is −¼ tr U.
- Negative control: swapping the order of two non-commuting factors must change the result.

This is the unit that lets YM-G03 be done in our own tools instead of by hand.
