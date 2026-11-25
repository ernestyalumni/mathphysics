# YM-G02: a first end-to-end instance of the conditional-gap chain, at strong coupling

Label: reviewed by Codex (YM-O05, `notes/YM-O05-codex-review-2026-11-22.md` §1–3); corrections applied 2026-11-23. It follows YM-O04 §6.1. It proves nothing about the continuum, and its strong-coupling conclusion is weaker than what Osterwalder–Seiler (1978) already gives by cluster expansion. Its purpose is to close every interface of the chain at one coupling, so that the chain can be machine-checked end to end. Lean pointers are to F01 at da4d739 (`docs/plans/yang-mills-mass-gap/formal/YMFormal/`).

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

**Corollary (seed).** If 0 ≤ β ≤ 1/(7480 · 192 · n⁴), then Δ_β(n) ≤ 1/7480. For n = 2 this is β ≤ 1/22 978 560 ≈ 4.35 × 10⁻⁸. For n = 1 it is β ≤ 1/1 436 160 ≈ 6.96 × 10⁻⁷. Sides of length 1 are allowed, with links kept as indexed objects (self-links and repeated incidences are not merged into a simple graph).
- On the all-one torus a plaquette is a commutator.
- On a rectangle where only side μ has length 1, the plaquette in the plane (μ, ν) is U_μ(x) U_ν(x) U_μ(x + ν̂)^{−1} U_ν(x)^{−1}, and its two μ-links need not coincide.

The bound holds in all these cases, including the comparison of side 1 with side 2. Boundedness gives this seed, but not transfer positivity for arbitrary bounded actions.

This is the same for any compact gauge group with a₀ real and |a₀| ≤ 1 (U(1), Z₂, SU(N) with a₀ = Re tr/N). It needs no cluster expansion. Lean needs only monotonicity of the integral and of the logarithm.

## 2. The structural fact (H) and the three spectral hypotheses

**(H)** For every box x with positive sides and every direction i, Z(x) = Tr T^{x_i}. Here T = T_{x̂ᵢ} is an operator built from the other three sides only, independent of x_i, nonzero, self-adjoint, positive and trace class.

Given (H), the three hypotheses are the content of `Purity.lean` and of YM-O03 §1, YM-O04 §2:
- Write λ₁ ≥ λ₂ ≥ ⋯ ≥ 0 for the eigenvalues with multiplicity, and s_t = Σ_{j≥2} (λ_j/λ₁)^t for the excited tail.
- The doubling cost in direction i at time extent t is p = 2 log Tr T^t − log Tr T^{2t} = 2 log(1 + s_t) − log(1 + s_{2t}) ≥ 0 (`nonneg`).
- The tail satisfies s_t ≤ e^p − 1, and the doubled cost satisfies p′ ≤ 2 e^{2p} p² ≤ 4p² for p ≤ 1/4 (`estimate`).
- At side N = 2ᵏn, take time as the **first** doubled direction of the telescope, so its spatial torus has side N. In Δ(N) = 8p₁ + 4p₂ + 2p₃ + p₄ all costs are ≥ 0, so 8p₁ ≤ Δ(N), and r^N ≤ s_N ≤ e^{p₁} − 1 ≤ e^{Δ(N)} − 1 (`from_defect`, including k = 0 on spatial side n).
  - Eigenvalues are ordered with multiplicity and normalized by the positive top eigenvalue; a missing second eigenvalue is 0.
  - From r_k^n ≤ 1/2, take logarithms only when r_k > 0. In the rank-one case the ratio is 0 and the gap is +∞.

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
- What it needs: matrices, the trace formula as a finite sum, and the existing F01 chain.
- Positivity is simpler than Schur here. With normalized counting measure, the one-link matrix ½[[e^β, e^{−β}], [e^{−β}, e^β]] has eigenvalues cosh β (constant vector) and sinh β (sign vector), both ≥ 0 for β ≥ 0.
- The steps: tensor over links; conjugate by √V; average over site signs; prove that the periodic finite sum equals the trace on the invariant subspace; apply the spectral inequalities.
- Include the normalized counting factors, side-1 self-links, and β = 0.
- The theorem quantifies over slice sizes, so it is not a single finite check.
- What it isn't: Z₂ is not Yang–Mills. For SU(2) the finite box additionally needs Haar changes of variables, measurable insertions, spectral summability and trace identities (Mercer), and the infinite-volume step of §3 is not finite-dimensional.

## 3. From the finite-box gap to an infinite-volume statement

This is the argument as specified by YM-O05 §2. Fix the dyadic sequence N = 2ᵏn and the given β. Write q = 2^{−1/n}.

In Lean (F02, `LocalLimit`, workspace-2c), steps 1–2 supply the field `finite_bound` (|corr k t| ≤ normSq · r_k^t), and step 3 supplies `tendsto`. `LocalLimit.decay` is the resulting bound on the limit correlation. Steps 4–5, the reconstruction and its gap, are outside that structure.

**Observable sector.** Use gauge-invariant, bounded, continuous local functions, supported in time slabs of any finite width at nonnegative times. Use the site reflection θ in the time-zero plane, with reflected links correctly oriented and inverted, and with no link of a positive-time observable crossing to negative times.
- Alternatively, work with the full link algebra and the operator A = P T̃ P on the **full** slice space (extra zero spectrum), sewing without a projection at time 0.
- Projecting the shared time-0 boundary twice is wrong: a gauge-variant function of the time-0 links can have nonzero OS norm but zero gauge average.

1. **Cylinder limit at fixed N.** Let S = T_N/λ₁, let h be the normalized top vector, and let Q = |h⟩⟨h|. Then ‖S^j − Q‖₁ = Σ_{l>1} (λ_l/λ₁)^j ≤ r^{j−1} Σ_{l>1} λ_l/λ₁ → 0, and Tr S^j → 1.
   - An observable F in an s-step slab enters the cyclic integral through an inserted slab operator B_F. Its kernel is bounded in absolute value by ‖F‖_∞ times the positive s-step kernel, so B_F is bounded.
   - The long complementary power converges in trace norm, so normalized expectations converge as the time extent grows.
   - This uses the local Wilson kernel, not only partition-function traces.
2. **OS vectors at fixed N.** Integrate the positive half-cylinder conditioned on the time-0 spatial links. The future boundary vector of F is the normalized inserted slab applied to h (for F = 1 it is h). Sewing it to the reflected past gives the OS inner product, and t empty layers in between give S^t. By Cauchy–Schwarz on h^⊥:

   |ω_N(θF̄ · τ_t G) − ω_N(θF̄) ω_N(G)| ≤ r^t ‖F‖_∞ ‖G‖_∞ ≤ q^t ‖F‖_∞ ‖G‖_∞,

   including overlapping boundaries at t = 0. Gauge-invariant insertions keep these vectors in Pℋ.
3. **Spatial limit.** Take the cylinder limit first. Then lift the spatial configurations periodically to the infinite lattice.
   - G is compact metrizable and there are countably many links, so the states have weakly convergent subsequences.
   - For fixed supports and fixed t, the translated and reflected products fit inside every sufficiently large N.
   - Weak convergence passes expectations of **continuous** local functions and of their products, not of arbitrary bounded measurable events. Translation invariance, site reflection positivity and the bound with q^t survive in every limit point ω.
   - These are limit points of these cylinders. They are not a classification of boundary conditions, and not a continuum state.
4. **Reconstruction and one-step positivity.** The OS classes of the local functions, over all slab widths, span a dense subspace of ℋ_ω by construction. Centring F ↦ F − ω(F) gives a dense subset of Ω^⊥.
   - Positivity of one time step does not follow from translation invariance and site reflection positivity alone. Pass the finite-cylinder inequality ω_N(θF̄ · τ₁F) ≥ 0, together with the contraction bounds, to the limit on the same dense algebra. Then 0 ≤ S_ω ≤ 1.
5. **Gap.** For a centred F and any b > q: b^t μ_F([b, 1]) ≤ ⟨ψ_F, S_ω^t ψ_F⟩ ≤ C_F q^t for all t. So the spectral projection of S_ω on [b, 1] annihilates a dense subset of Ω^⊥, and by continuity all of it.
   - The fixed vectors of S_ω are the vacuum line: the vacuum is unique **in this reconstruction**. That is distinct from uniqueness of the infinite-volume Gibbs state, which is not used.
   - The transfer gap is at least (log 2)/n in lattice units.

## 4. What the end-to-end instance proves and what it does not

- **Proves, with §3 as specified:** for SU(2), U(1) and Z₂ Wilson lattice gauge theory at 0 ≤ β ≤ 1/(7480·192·n⁴), every limit point of the vacuum states along the dyadic tori reconstructs to a lattice theory whose transfer operator has a unique vacuum and a gap of at least (log 2)/n in lattice units.
  - For U(1) this is the confined strong-coupling phase. It does not contradict the Coulomb phase at large β.
  - The Z₂ finite-box part can be fully formal.
- **Does not prove:** anything at weak coupling, anything uniform as a → 0, or anything in physical units. Osterwalder–Seiler's cluster expansion already gives a much larger gap at small β (to leading order m ≈ −4 log u, with u = I₂(β)/I₁(β) ≈ β/4).

## 5. The uniform seed and the OpenAI construction: three audit questions

The precise target (YM-N03, "The precise seed target"): there is L < ∞ such that, along the scaling trajectory and sufficiently far toward the continuum, some n has a(β) n ≤ L and Δ_β(n) ≤ 1/7480, with the same bare parameters on both boxes.

`01-WHAT-OPENAI-PROVED.md` §3 says the construction's renormalization group stops at unit scale. The bare coupling is tuned so that the level-0 coupling is one fixed large number β̄, for every N and every volume index p. Three questions decide whether that construction can supply or test a seed:
1. **Is the bare tuning independent of the volume index p?** The construction does not establish it.
   - The public object is two-indexed: `PublicWilson.lean` has `beta : ℕ → ℕ → ℝ`, and `ChallengeSpec.lean` asks for a uniform slope, not equality in p.
   - The informalization fixes the coarse endpoint independently of p and shoots uniquely for each (N, p). No theorem identifies two bare tunings, though none shows that particular values differ.
   - A seed taken from the construction needs matched boxes and bare parameters, or control of retuning. Since |∂_β log Z| ≤ P for the shifted weight, changing β on the doubled box alone costs up to P(2n)|δβ|; leading-log agreement is not enough. (YM-O05 §3; the handbook records the same dependence.)
2. **Is the level-0 effective theory reflection positive in every direction, so that (H) holds for it?**
   - This is required only to apply the recurrence to the effective actions themselves. That would also need common rectangular time-independent trace kernels, not merely reflection positivity on periodic cubes.
   - It is not required if the renormalization group is used only to bound the **fine** Wilson seed, for which (H) holds by §2.
   - No proof of coarse (H) is in the cited result.
3. **What is the error on log Z at unit scale?**
   - The normalized block step preserves the total unnormalized mass exactly on its own trajectory, but the normalized Wilson probability alone loses the scalar.
   - A seed needs the scalar normalizations kept, β matched, and 16 log Z_n − log Z_{2n} bounded, including non-extensive terms, with margin below 1/7480.
   - Two extensive error bars, or small polymer activities without a finite-volume remainder, do not do this.

Even with positive answers, the seed for the level-0 theory is a statement at the fixed **weak** effective coupling β̄. The elementary bound of §1 does not apply there. That is the open infrared problem. The finite-size criteria surveyed in YM-G06 (Dobrushin–Shlosman, Martinelli–Olivieri) have the same shape: a finite-box condition implying decay. For spin systems they need boxes large compared with the correlation length, which is the same obstruction.

## Next

1. The Z₂ finite-box instance is under way in F02 (workspace-2c), approved by Ernest.
2. A necessary condition for any seed (YM-O05 §5): at side n, n · m_box ≥ −log(e^{1/59840} − 1) = 10.9994, where m_box is the true first excitation on the spatial torus of side n.
