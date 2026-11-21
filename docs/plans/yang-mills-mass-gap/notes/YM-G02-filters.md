# The conditional-gap chain against the two filters

Label: unreviewed (claude-opus55-strat, 2026-11-21). Objects: `ym:prop:purity`, `ym:prop:defect`, `ym:prop:conditional-gap`, `ym:conj:recurrence` and `ym:conj:defect` in `documents/mono-notes/topics/yang-mills-lattice.tex`.

Two filters apply to any proposed gap argument:
- it must fail for 4D U(1), which has a massless Coulomb phase;
- it must use some input that distinguishes θ = 0 from θ = π. At θ = π, SU(2) Yang–Mills has a mixed anomaly between time reversal and the 1-form symmetry, so it cannot be trivially gapped (Gaiotto–Kapustin–Komargodski–Seiberg 2017).

## θ = π: passed, through positivity

`ym:prop:purity` assumes a self-adjoint trace-class T with eigenvalues λ₁ ≥ λ₂ ≥ ⋯ ≥ 0, and `ym:prop:defect` assumes T ≥ 0 on every spatial torus. Both come from reflection positivity of the θ = 0 Wilson action.

A θ term makes the Euclidean weight e^{iθQ} non-positive; at θ = π it is ±1. T is then no longer nonnegative, so r_j ∈ [0, 1], P ≤ 1 and the telescope fail. The chain thus rests on an input that fails at θ = π, which is what the filter requires. Any lemma that replaces positivity must keep this property.

## U(1): the non-abelian input must sit in the conjectures

The proved part of the chain (purity, telescope, conditional gap) holds verbatim for compact U(1), whose Wilson transfer operator is nonnegative too. Since 4D U(1) has no gap at weak coupling, at least one of the two conjectures must fail for U(1):
- `ym:conj:defect`: Δ_β(n(β)) ≤ η* at a bounded physical scale;
- `ym:conj:recurrence`: Δ(n) ≤ 1/(2c₁) ⟹ Δ(2n) ≤ c₁Δ(n)² with c₁ depending only on G.

Measured so far (`notes/YM-N03-defect-measurement.md`):
- SU(2), β = 1.8: Δ(2) ≈ 26, Δ(4) ≈ 0. The defect collapses once the torus exceeds the deconfinement length.
- U(1), n = 2: in the Coulomb phase the finite-size integrand stays large (17.8 at β = 1.6) and Δ(2) ≈ 80–88.

The expected outcome is that U(1)'s Δ(n) never becomes small at weak coupling. Then the recurrence's premise is never met and `ym:conj:defect` is the conjecture that fails. A proof would have to use non-abelian structure exactly there, in why small defect is reached at a bounded physical scale. U(1) n = 4 (8⁴ scans) is running.

## Testing the recurrence itself

A direct test needs Δ(8) near the n = 4 threshold, where Δ(4) ≈ 0.3 at β ≈ 2.0. The conjecture then predicts Δ(8) ≤ 0.09 c₁. β-integration would need about 67 GPU-hours.

The block-gluing estimator (`Rust/lattice-gauge/src/glued.rs`) interpolates between the (2n)⁴ torus and 2⁴ independent n⁴ tori over the boundary plaquettes only, so the integrand is measured in one ensemble per λ. It reproduces the exact 2D square defect, and in 4D at n = 2 it agrees with β-integration (0.193(81) against 0.187(23) at β = 0.6; 1.41(12) against 1.456(31) at β = 1.0). n = 4 at β = 1.8 and 2.0 is running as a precision test before n = 8.
