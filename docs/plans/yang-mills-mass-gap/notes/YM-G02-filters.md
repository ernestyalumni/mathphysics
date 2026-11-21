# The conditional-gap chain against the two filters

Label: unreviewed (claude-opus55-strat, 2026-11-21). Objects: `ym:prop:purity`, `ym:prop:defect`, `ym:prop:conditional-gap`, `ym:conj:recurrence` and `ym:conj:defect` in `documents/mono-notes/topics/yang-mills-lattice.tex`.

Two filters apply to any proposed gap argument:
- it must fail for 4D U(1), which has a massless Coulomb phase;
- it must use some input that distinguishes θ = 0 from θ = π. At θ = π, SU(2) Yang–Mills has a mixed anomaly between time reversal and the 1-form symmetry, so it cannot be trivially gapped (Gaiotto–Kapustin–Komargodski–Seiberg 2017).

## θ = π: what the filter does and does not require

`ym:prop:purity` assumes a self-adjoint trace-class T with eigenvalues λ₁ ≥ λ₂ ≥ ⋯ ≥ 0, and `ym:prop:defect` assumes T ≥ 0 on every rectangular torus. Both come from reflection positivity of the θ = 0 Wilson action. A θ term makes the Euclidean weight e^{iθQ} non-positive (±1 at θ = π), so the chain does not apply at θ ≠ 0.

The anomaly at θ = π excludes only a *unique, symmetry-preserving, trivially gapped* vacuum. Two time-reversed gapped vacua are allowed. So the filter constrains arguments that also fix the symmetry realization and rule out degeneracy and topological order (YM-O03 §5). An earlier version of this note stated it more broadly.

## U(1): the non-abelian input must sit in the seed

The proved part of the chain holds verbatim for compact U(1), whose Wilson transfer operators are nonnegative too. Since YM-O03 §1, so does the recurrence: Δ(2n) ≤ 3740 Δ(n)² whenever Δ(n) ≤ 1/7480, from positivity alone. Checked by claude-opus55-strat and claude-opus55-q7nd; the β = 0 rank-one case is in F01's RankOne.lean.

Since 4D U(1) has no gap at weak coupling, it follows that **no seed Δ_β(n) ≤ 1/7480 can exist at any n for U(1) in its Coulomb phase**. Otherwise the theorem would give a gap. So `ym:conj:defect`, the seed at a bounded physical size, is exactly where any proof must use non-abelian structure.

The measured U(1) defect at n = 2 grows into the hundreds across β ≈ 0.9 and stays at about 80–88 in the Coulomb phase. U(1) n = 4 is running. These comparisons locate where seeds fail; they do not test the (now proved) quadratic law.

## What remains numerically

- The seed threshold 1/7480 is far below the value 1/4 used earlier as a cut.
- From the strong-coupling formula, the seed holds at n = 2 for β ≲ 0.10 and at n = 4 for β ≲ 1.36. That is the strong-coupling regime, where a gap is already known.
- Whether a seed exists at a bounded physical size as β → ∞ is the open problem, and present lattice sizes cannot settle it.
- The block-gluing estimator (`Rust/lattice-gauge/src/glued.rs`) remains a validated tool for defect differences, but it is no longer needed for the recurrence.
