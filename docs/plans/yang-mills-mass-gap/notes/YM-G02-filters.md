# The conditional-gap chain against the two filters

Label: unreviewed (claude-opus55-strat, 2026-11-21). Objects: `ym:prop:purity`, `ym:prop:defect`, `ym:prop:conditional-gap`, `ym:conj:recurrence` and `ym:conj:defect` in `documents/mono-notes/topics/yang-mills-lattice.tex`.

Two filters apply to any proposed gap argument:
- it must fail for 4D U(1), which has a massless Coulomb phase;
- it must use some input that distinguishes θ = 0 from θ = π. At θ = π, SU(2) Yang–Mills has a mixed anomaly between time reversal and the 1-form symmetry, so it cannot be trivially gapped (Gaiotto–Kapustin–Komargodski–Seiberg 2017).

## θ = π: what the filter does and does not require

`ym:prop:purity` assumes a self-adjoint trace-class T with eigenvalues λ₁ ≥ λ₂ ≥ ⋯ ≥ 0, and `ym:prop:defect` assumes T ≥ 0 on every rectangular torus. Both come from reflection positivity of the θ = 0 Wilson action. A θ term makes the Euclidean weight e^{iθQ} non-positive (±1 at θ = π), so the chain does not apply at θ ≠ 0.

The anomaly at θ = π excludes only a *unique, symmetry-preserving, trivially gapped* vacuum. Two time-reversed gapped vacua are allowed. So the filter constrains arguments that also fix the symmetry realization and rule out degeneracy and topological order (YM-O03 §5). An earlier version of this note stated it more broadly.

## U(1): the non-abelian input must sit in the seed

The proved part of the chain holds verbatim for compact U(1), whose Wilson transfer operators are nonnegative too. Since YM-O03 §1, so does the recurrence: Δ(2n) ≤ 3740 Δ(n)² whenever Δ(n) ≤ 1/7480, from positivity alone. Checked by claude-opus55-strat and claude-opus55-q7nd, and proved in Lean (F01 `Rectangle.lean`, 2ed351b, standard axioms only: `rectangle_recurrence`, `box_recurrence`, `conditional_gap_of_boxes`). The β = 0 rank-one case is in `RankOne.lean`. The only inputs outside Lean are positivity of the four-direction rectangular transfer operators (the C23 identification, which gives nonnegative doubling costs and the purity estimate) and the seed.

A seed would give λ₂/λ₁ ≤ 2^{−1/n} on the dyadic boxes at that β. To conclude that **no seed can exist at any n for U(1) in its Coulomb phase**, one needs the bridge to the infinite-volume state: vacuum local states along these boxes, a local limit, and a tested observable with nonzero gapless photon spectral weight, in the same boundary conditions, sector and action as a masslessness theorem (for example Fröhlich–Spencer). With that bridge the contrapositive holds; it is not a consequence of the recurrence alone or of two data points (YM-O04 §2). So `ym:conj:defect` is where a proof must use non-abelian structure.

Measured: the U(1) defect at n = 2 rises through β ≈ 0.9 (where its integrand peaks at about 249) and stays at about 80–88 in the Coulomb phase. At n = 4 the increments from β = 0.6 are 70.8 at β = 1.2 and 80.3 at β = 1.6. These comparisons locate where seeds fail; they do not test the (now proved) quadratic law.

## What remains numerically

- The seed threshold 1/7480 is far below the value 1/4 used earlier as a cut.
- The truncated surface formula crosses the threshold at β ≈ 0.097 (n = 2) and 1.355 (n = 4); these are not certified. A certified seed comes from Δ_β(n) ≤ 192 β n⁴: at n = 2, any β ≤ 1/(7480·3072).
- Whether a seed exists at a bounded physical size as β → ∞ is the open problem, and present lattice sizes cannot settle it.
- The block-gluing estimator (`Rust/lattice-gauge/src/glued.rs`) remains a validated tool for defect differences, but it is no longer needed for the recurrence.
