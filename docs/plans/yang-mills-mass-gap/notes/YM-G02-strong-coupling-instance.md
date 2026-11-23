# YM-G02: a first end-to-end instance of the conditional-gap chain, at strong coupling

Label: unreviewed (claude-opus55-strat, 2026-11-21). It follows YM-O04 §6.1. It proves nothing about the continuum, and its strong-coupling conclusion is weaker than what Osterwalder–Seiler (1978) already gives by cluster expansion. Its purpose is to close every interface of the chain at one coupling, so that the chain can be machine-checked end to end. Lean pointers are to F01 at 2ed351b (`docs/plans/yang-mills-mass-gap/formal/YMFormal/`).

## The chain and its inputs

`conditional_gap_of_boxes` (`Rectangle.lean`) takes the function F(x) = log Z(x) on boxes x = (x₀, x₁, x₂, x₃) with positive sides, and four hypotheses:
- `nonneg`: every doubling cost 2F(x) − F(x with side i doubled) is ≥ 0;
- `estimate`: a doubling cost p ≤ 1/4 at x bounds the cost p′ at the doubled box by 4p²;
- `from_defect`: at side N = 2ᵏn, r_k^N ≤ exp(Δ(N)) − 1, where r_k = λ₂/λ₁ for the transfer operator of the spatial torus of side N;
- `seed`: Δ(n) = 16F(n,n,n,n) − F(2n,2n,2n,2n) ≤ 1/7480.

Its conclusion: r_k^n ≤ 1/2 for every k. That is, λ₂/λ₁ ≤ 2^{−1/n}, a gap of at least (log 2)/n in lattice units, on every spatial torus of side 2ᵏn, at the given β.

The first three hypotheses follow from one structural fact (H) about the Wilson weight. The seed is a separate estimate. The sections below give each input.

## 1. The seed: an elementary bound

Let Z_β(x) = ∫ ∏_ℓ dU_ℓ exp(β Σ_p a₀(U_p)), with normalized Haar measure and a₀ = ½ tr. The sum runs over the P(x) = 6 x₀x₁x₂x₃ plaquettes. Since |a₀| ≤ 1, each integrand lies between e^{−βP} and e^{βP}, and so do the integrals:

  −βP(x) ≤ log Z_β(x) ≤ βP(x).

With P(2n,…) = 16 P(n,…) = 96 n⁴:

  Δ_β(n) ≤ 16 βP(n) + βP(2n) = 32 β P(n) = 192 β n⁴.

**Corollary (seed).** If 0 ≤ β ≤ 1/(7480 · 192 · n⁴), then Δ_β(n) ≤ 1/7480. For n = 2 this is β ≤ 1/22 978 560 ≈ 4.35 × 10⁻⁸. For n = 1 it is β ≤ 1/1 436 160 ≈ 6.96 × 10⁻⁷, if the conventions for a side-1 torus are accepted: there U_μ(x + μ̂) = U_μ(x), and the plaquettes in a plane containing μ are commutators.

This is the same for any compact gauge group with a₀ real and |a₀| ≤ 1 (U(1), Z₂, SU(N) with a₀ = Re tr/N). It needs no cluster expansion. Lean needs only monotonicity of the integral and of the logarithm.

## 2. The structural fact (H) and the three spectral hypotheses

**(H)** For every box x with positive sides and every direction i, Z(x) = Tr T^{x_i}. Here T = T_{x̂ᵢ} is an operator built from the other three sides only, independent of x_i, nonzero, self-adjoint, positive and trace class.

Given (H), the three hypotheses are the content of `Purity.lean` and of YM-O03 §1, YM-O04 §2:
- Write λ₁ ≥ λ₂ ≥ ⋯ ≥ 0 for the eigenvalues with multiplicity, and s_t = Σ_{j≥2} (λ_j/λ₁)^t for the excited tail.
- The doubling cost in direction i at time extent t is p = 2 log Tr T^t − log Tr T^{2t} = 2 log(1 + s_t) − log(1 + s_{2t}) ≥ 0 (`nonneg`).
- The tail satisfies s_t ≤ e^p − 1, and the doubled cost satisfies p′ ≤ 2 e^{2p} p² ≤ 4p² for p ≤ 1/4 (`estimate`).
- At side N, r^N ≤ s_N ≤ exp(cost) − 1. The cost of doubling the time extent N → 2N on the N⁴ box is the first term of the defect telescope, so it is ≤ Δ(N), because all terms are ≥ 0 (`cost_le_defect_binary`; `from_defect`).

### Proof of (H) for the Wilson action (standard; Lüscher, Commun. Math. Phys. 54 (1977) 283; Osterwalder–Seiler, Ann. Phys. 110 (1978) 440)

By the hypercubic symmetry of the Wilson action, it suffices to take i = 0 as the time direction.
- **Spatial configurations.** Let 𝒰 = G^{E_s} be the spatial links of one time slice, with product Haar measure, and let ℋ = L²(𝒰).
- **Temporal gauge.** Fix the temporal links to 1 on every time slice except the last; this is a change of variables. The remaining temporal links at the last slice act on the spatial links of slice 0 as a gauge transformation g. Integrating over g gives the projection P = ∫ dg R(g) onto gauge-invariant functions, where R(g) is the unitary action of gauge transformations on ℋ.
- **The trace formula.** Write V(U) = exp(β Σ_{spatial p} a₀(U_p)) for the spatial plaquettes of one slice, and K(U, U′) = ∏_{ℓ ∈ E_s} k(U_ℓ U′_ℓ^{−1}) with k(g) = exp(β a₀(g)) for the temporal plaquettes. Then Z(x) = Tr(P T̃^{x₀}) with T̃ = V^{1/2} K V^{1/2}. Since P commutes with V and K, Z(x) = Tr T^{x₀} with T = P T̃ P on Pℋ.
- **Positivity, the elementary way.** For SU(2), identify U with a unit vector u ∈ S³ ⊂ ℝ⁴ (quaternions). Then a₀(U U′^{−1}) = ½ tr(U U′^†) = u · u′, the Euclidean dot product. The kernel (u, u′) ↦ u · u′ is a Gram kernel, hence positive semidefinite. By Schur's product theorem, each power (u · u′)^m is positive semidefinite. Therefore k(U U′^{−1}) = exp(β u · u′) = Σ_m β^m (u · u′)^m / m! is positive semidefinite for β ≥ 0.
  - The same holds for U(1) (cos(θ − θ′) is the dot product of points on S¹) and for Z₂ (σσ′).
  - Tensor products of positive semidefinite kernels are positive semidefinite, so K ≥ 0. Hence T̃ = V^{1/2} K V^{1/2} ≥ 0, with V bounded and positive. This avoids Peter–Weyl.
- **Trace class.** k is a continuous positive semidefinite kernel on a compact space with a finite measure. By Mercer's theorem, the integral operator is trace class with trace ∫ k(u, u) du = e^β per link. A finite tensor product of trace-class operators is trace class, and multiplication by the bounded V^{1/2} preserves this. T ≠ 0 because Tr T = Z > 0.

So (H) holds for the Wilson action of SU(2), U(1) and Z₂ at every β ≥ 0, every positive side, and every direction. For general SU(N), a₀(UU′^{−1}) = Re tr(U U′^†)/N is again a real inner product (Frobenius) on the unitary matrices, and the same argument works.

### Formalization note

For Z₂ every space is finite-dimensional: ℋ = ℝ^{2^{|E_s|}}, T is a positive semidefinite matrix, and (H) is a finite-sum identity. **Z₂ lattice gauge theory is therefore a target for a fully formal end-to-end theorem in Lean with no outside hypotheses.**
- Statement: for 0 ≤ β ≤ 1/(7480·192·n⁴) and every k, the Z₂ transfer matrix on the spatial torus of side 2ᵏn has λ₂/λ₁ ≤ 2^{−1/n}.
- What it needs: matrices, Schur's product theorem (in mathlib for `Matrix.PosSemidef` Hadamard products, to be checked), the trace formula as a finite sum, and the existing F01 chain.
- What it isn't: Z₂ is not Yang–Mills. But every interface is the same, and the SU(2) version differs only in replacing matrices by trace-class integral operators (Mercer).

## 3. From the finite-box gap to an infinite-volume statement (sketch)

YM-O04 §2 lists this step as a separate hypothesis. Here is a proposed proof, to be reviewed.

- **Cylinder states.** Fix the spatial side N = 2ᵏn and let the time extent go to infinity. Since λ₁ > λ₂, the torus expectations converge to the vacuum state ω_N of T_N. Let θ be the reflection in the time-zero plane of sites. Reflection positivity gives, for every bounded F supported at times in [0, s], an Osterwalder–Schrader vector ψ_F with ⟨ψ_F, ψ_G⟩ = ω_N(θF̄ · G) and ⟨ψ_F, (T_N/λ₁)^t ψ_G⟩ = ω_N(θF̄ · τ_t G), where τ_t is translation by t in time. For F, G supported in [0, s] and t ≥ 0:

  |ω_N(θF̄ · τ_t G) − ω_N(θF̄) ω_N(G)| ≤ (λ₂/λ₁)^t ‖F‖_∞ ‖G‖_∞ ≤ 2^{−t/n} ‖F‖_∞ ‖G‖_∞.

  Proof: the vacuum term of the spectral decomposition of T_N/λ₁ is the product of expectations. Apply Cauchy–Schwarz on the orthogonal complement of the vacuum, with ‖ψ_F‖² = ω_N(θF̄ · F) ≤ ‖F‖²_∞.
  - The constant 2^{−t/n} is uniform in k: this is where the chain enters.
- **Local limit.** The configuration space G^{E(ℤ⁴)} is compact, so the states ω_{2ᵏn}, viewed as measures on a neighbourhood of the origin that eventually fits inside each torus, have weakly convergent subsequences (Prokhorov, or Tychonoff and Riesz). For fixed local F, G and t, the inequality passes to every limit point ω.
- **Reconstruction.** ω is translation invariant and reflection positive (both are closed conditions), so the Osterwalder–Seiler construction gives a Hilbert space ℋ_ω, a vacuum Ω, and a positive contraction T_ω for one time step. By construction, the vectors ψ_F for local F are dense.
- **Gap.** For local F orthogonal to Ω, the spectral measure μ_F of T_ω satisfies ∫ λ^t dμ_F = ⟨ψ_F, T_ω^t ψ_F⟩ ≤ 2^{−t/n} ‖F‖²_∞ for all t. So μ_F is supported in [0, 2^{−1/n}]. By density, the spectrum of T_ω on Ω^⊥ lies in [0, 2^{−1/n}]. In particular the vacuum is unique, and the mass gap is m ≥ (log 2)/n in lattice units.
  - The vacuum subtraction is exactly the first term above, so there is no separate cyclicity issue for the vectors used.
  - What remains to be checked is that the vectors ψ_F with F supported in single slabs exhaust the reconstruction. In the Osterwalder–Seiler construction they do.
- **Which state.** Every limit point has the gap. At such small β the infinite-volume Gibbs state is unique (Dobrushin's criterion), but that is not needed for the gap statement.

**Points for review.**
- (a) Convergence of the cylinder limit, using only λ₁ > λ₂ and trace class.
- (b) Passing the inequality to limit points for observables whose support needs the time slab inside the torus.
- (c) The density claim in the reconstruction.

## 4. What the end-to-end instance proves and what it does not

- **Proves (once §3 is reviewed):** for SU(2), U(1) and Z₂ Wilson lattice gauge theory at 0 ≤ β ≤ 1/(7480·192·n⁴), the infinite-volume lattice theory along dyadic tori has a mass gap m ≥ (log 2)/n in lattice units.
  - For U(1) this is the confined strong-coupling phase. It does not contradict the Coulomb phase at large β.
  - The Z₂ finite-box part can be fully formal.
- **Does not prove:** anything at weak coupling, anything uniform as a → 0, or anything in physical units. Osterwalder–Seiler's cluster expansion already gives a much larger gap at small β (to leading order m ≈ −4 log u, with u = I₂(β)/I₁(β) ≈ β/4).

## 5. The uniform seed and the OpenAI construction: three audit questions

The precise target (YM-N03, "The precise seed target"): there is L < ∞ such that, along the scaling trajectory and sufficiently far toward the continuum, some n has a(β) n ≤ L and Δ_β(n) ≤ 1/7480, with the same bare parameters on both boxes.

`01-WHAT-OPENAI-PROVED.md` §3 says the construction's renormalization group stops at unit scale. The bare coupling is tuned so that the level-0 coupling is one fixed large number β̄, for every N and every volume index p. Three questions decide whether that construction can supply or test a seed:
1. **Is the bare tuning independent of the volume index p?** No. The handbook (`yang-mills-lattice.tex`, "What a proof of the gap requires", item 1) records that the bare coupling β_{N,p} depends on p through the condition β_{N,p,0} = β̄, controlled only to o(log(1/a_N)). The defect compares two volumes at the same bare parameters, so any seed taken from that construction must be stable under the difference β_{N,p} − β_{N,p′}.
2. **Is the level-0 effective theory reflection positive in all four directions, so that (H) holds for it?** The chain needs (H) for whatever theory the seed is stated in. The renormalized actions are not Wilson actions.
3. **What is the error on log Z at unit scale?** Volume terms cancel in Δ, but an error of order the volume on each box is useless (YM-O04 §6.2). The seed needs a remainder below 1/7480 in absolute terms.

Even with positive answers, the seed for the level-0 theory is a statement at the fixed **weak** effective coupling β̄. The elementary bound of §1 does not apply there. That is the open infrared problem. The finite-size criteria surveyed in YM-G06 (Dobrushin–Shlosman, Martinelli–Olivieri) have the same shape: a finite-box condition implying decay. For spin systems they need boxes large compared with the correlation length, which is the same obstruction.

## Next

1. Review of §3 (Codex), then a decision on the Z₂ fully formal instance (F01 owner).
2. The three audit questions in §5, against the OpenAI construction's Lean sources (YM-A tasks).
