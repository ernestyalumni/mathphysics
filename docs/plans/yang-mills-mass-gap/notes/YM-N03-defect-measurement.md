# YM-N03: measuring the hypercubic defect and the four purities

Label: reviewed with corrections by Codex (YM-O03, `notes/YM-O03-codex-review-2026-11-21.md` §4); corrections applied 2026-11-21. Author claude-opus55-strat. Code: `Rust/lattice-gauge` on branch `ym/YM-N01-claude-opus55-su2`, using `examples/plaquette_scan`, `plaquette_scan_gpu`, `defect_ti` and `defect_window`. Raw data are in the private run directory, not in the repository.

## What is measured

The objects are those of `yang-mills-lattice.tex` and of the YM-F01 Lean statements:
- the defect Δ_β(n) = 16 log Z(n,n,n,n) − log Z(2n,2n,2n,2n);
- the purities P₁…P₄ of `ym:four-purities`, with e^{−Δ} = P₁⁸P₂⁴P₃²P₄;
- the bound m(2n) ≥ (1/n) log(1/(e^Δ − 1)) of `ym:gap-from-defect`.

Thresholds:
- Δ < log 2 makes the gap bound positive.
- Δ ≤ 1/(2c₁) ≤ 1/4 is the hypothesis of the conditional gap.

## Method

With normalized Haar measure, Z(β = 0) = 1 and log Z(V; β) = N_P(V) ∫₀^β ⟨a₀⟩_V dβ′. Because N_P(2n) = 16 N_P(n),

  Δ_β(n) = N_P(2n) ∫₀^β [⟨a₀⟩_{n⁴} − ⟨a₀⟩_{(2n)⁴}] dβ′,

and each log P_k has the same form between consecutive volumes. The defect is therefore a finite-size effect: the large volume times the plaquette difference between the small torus and the doubled one.

Integration uses Simpson's rule on a common β grid, with independent jackknife errors per volume and point. Replicas with different seeds are averaged.

**Validation in 2D.** The square defect 4 log Z(n,n) − log Z(2n,2n) is exact from the character expansion, Z = Σ_j (2 I_{2j+1}(β)/β)^V. Monte Carlo against exact, for n = 2:

| β | Monte Carlo | exact |
|---|---|---|
| 1.2 | 0.0204(49) | 0.0258 |
| 2.0 | 0.1398(62) | 0.1396 |
| 2.4 | 0.2342(66) | 0.2343 |
| 3.2 | 0.4746(72) | 0.4718 |
| 4.0 | 0.7375(77) | 0.7335 |

All within 1.3σ. Solving the exact formula at n = 2, the defect crosses 1/4 at β = 2.459 and log 2 at β = 3.877 (an earlier version said 3.6, which was wrong).

Normalization: the manuscript weight is exp[−β Σ(1 − a₀)]. The formulas here use the shifted weight exp(β Σ a₀); the extensive shift cancels in Δ and in log P_k.

## 4D SU(2), n = 2

Volumes 2⁴, 4·2³, 4²·2², 4³·2 and 4⁴; β ∈ [0, 3] with h = 0.1; 4×10⁵ updates per point and volume, in two replicas (Mac mini and laptop).

| β | Δ_β(2) | log P₁ | log P₂ | log P₃ | log P₄ | m(4) bound |
|---|---|---|---|---|---|---|
| 0.4 | 0.056(18) | −0.002(3) | −0.010(4) | −0.006(5) | 0.012(8) | 1.42 |
| 0.6 | 0.187(23) | −0.010(3) | −0.022(5) | −0.015(7) | 0.010(9) | 0.79 |
| 0.8 | 0.591(27) | −0.037(4) | −0.053(6) | −0.047(8) | 0.005(11) | 0.11 |
| 1.0 | 1.456(31) | −0.093(5) | −0.120(6) | −0.113(9) | −0.006(12) | — |
| 1.6 | 13.14(4) | −0.869(6) | −1.019(9) | −0.934(12) | −0.248(16) | — |
| 2.0 | 43.61(6) | −2.430(8) | −2.959(12) | −3.428(16) | −5.474(22) | — |
| 2.4 | 69.54(6) | −3.309(10) | −3.861(14) | −6.067(19) | −15.49(3) | — |

Consistency checks of the formal statements hold at every β, within errors: P_k ≤ 1 and P₄ ≥ e^{−Δ}. The small positive values of log P₄ at β ≤ 0.8 are within 1.5σ of zero.

**Reading.** Δ_β(2) is below 1/4 only for β ≲ 0.65, at strong coupling, and grows by a factor of about 30 between β = 1 and 2.

The 2⁴ torus has two sites in every direction. Its rapid growth sits near the coupling where SU(2) with N_t = 2 has its thermal transition (β_c ≈ 1.9, defined at growing spatial volume). A finite symmetric torus has no true phase transition, so this is a finite-size crossover. The peak positions suggest, but cannot establish, a link to confinement on the small torus (YM-O03 §4).

## 4D SU(2), n = 4

Data:
- 4⁴: 12 replicas on β ∈ [0, 2.6], h = 0.1, 5×10⁵ updates per point (Mac mini).
- 8⁴: 15 GPU replicas on β ∈ [1.0, 2.6], 10⁵ updates per point (3 on the desktop 3060, 12 on the laptop 3070).

**These are increments.** The 8⁴ scans start at β₀ = 1.0, so the table gives Δ(β) − Δ(1.0). A near-zero integrand above β₀ says nothing about the integral below it (YM-O03 §4). The strong-coupling formula below estimates Δ(4) at β = 1.0 as about 2×10⁻⁶. That is an estimate, not a bound; a bound would need the 8⁴ integrand on [0, 1].

| β | integrand | Δ_β(4) − Δ_β(4)(1.0) |
|---|---|---|
| 1.6 | 0.03(19) | −0.012(41) |
| 1.8 | 0.57(23) | −0.003(51) |
| 2.0 | 4.77(29) | 0.336(63) |
| 2.2 | 134.4(4) | 9.32(8) |
| 2.4 | 107.1(2) | 40.89(9) |
| 2.6 | 40.5(2) | 53.43(10) |

The errors are jackknife errors per point, combined as independent. Each CPU scan anneals one chain across β, so independence across grid points is an approximation. Simpson discretization error is not included: at β = 2.2 every other grid point gives 10.25 instead of 9.32.

**Reading.**
- With the offset taken as negligible, the increment stays below 1/4 up to β ≈ 1.97. That crossing depends on that normalization and on interpolation.
- The negative central values cannot enter the gap formula. Taking the offset as zero at β = 2.0, the central value of the m(8) expression is (1/4) log(1/(e^{0.336} − 1)) = 0.229 (an earlier version said ≥ 0.24, which was wrong). A real lower bound needs an upper confidence limit on the absolute defect, including systematics.
- The integrand peaks near β ≈ 2.3, close to the N_t = 4 thermal coupling. As at n = 2, this is a finite-size crossover on a symmetric torus. My earlier prediction that Δ(4) stays small through β ≈ 2.2–2.3 failed: the increment at 2.2 is already 9.3.
- These numbers do not establish the physical-size hypothesis. The gap expressions are also far below the actual lattice glueball mass, about 2 near β = 2.0.

## U(1), the counter-model

n = 2 (2⁴ and 4⁴, 6 replicas each, β ∈ [0, 1.6] with h = 0.05, 5×10⁵ updates per point; Mac mini):

| β | integrand | Δ_β(2) |
|---|---|---|
| 0.2 | 1.47(7) | 0.075(7) |
| 0.3 | 5.03(7) | 0.383(9) |
| 0.6 | 49.3(1) | 6.54(1) |
| 0.9 | 249.0(2) | 50.23(2) |
| 1.0 | 79.0(1) | 71.30(3) |
| 1.2 | 29.08(6) | 79.56(3) |
| 1.6 | 17.84(4) | 88.48(3) |

- Δ crosses 1/4 at β ≈ 0.27, earlier than SU(2): U(1)'s strong-coupling parameter (≈ β/2) grows twice as fast as SU(2)'s (≈ β/4).
- The integrand peaks near β ≈ 0.9 and then falls, from about 249 to 18 at β = 1.6. Its nonzero value at the last point cannot identify its cause or its asymptotics.

**U(1) n = 4.** 4⁴: 6 replicas from the n = 2 run. 8⁴: 14 laptop replicas on β ∈ [0.6, 1.6], 10⁵ updates per point; 12 more Mac replicas are running. The values are increments from β₀ = 0.6. The strong-coupling estimate of Δ_U(1)(4) at 0.6 is 6×10⁻⁵, an estimate, not a bound.

| β | Δ(2) | Δ(4) − Δ(4)(0.6) | integrand n = 2 | integrand n = 4 |
|---|---|---|---|---|
| 0.9 | 50.2 | 0.17(5) | 249.0(2) | 6.9(6) |
| 1.0 | 71.3 | 34.3(1) | 79.0(1) | 1573(2) |
| 1.2 | 79.6 | 70.8(1) | 29.08(6) | 33.1(3) |
| 1.6 | 88.5 | 80.3(1) | 17.84(4) | 18.3(2) |

- In the Coulomb phase, doubling the box barely changes the defect (about 80 at both sizes), and the finite-size integrand is the same at n = 2 and n = 4. That is the behaviour of massless modes.
- SU(2) at β = 1.8 gives Δ(2) ≈ 26 against Δ(4) ≈ 0.
- U(1) thus never approaches a seed, as the recurrence theorem requires for a gapless theory. This locates where the seed fails; it does not test the (proved) quadratic law.
- The integrand spike near β = 1.0 is U(1)'s first-order transition on the larger lattice. Since the recurrence is a theorem for any positive transfer spectrum (next section), it holds for U(1) too. So the U(1) filter now says, with certainty: for U(1) at weak coupling no seed Δ(n) ≤ 1/7480 can exist at any n, or the theorem would give a gap. Comparisons of Δ across n test where the seed fails. They do not test a quadratic law.

## Strong-coupling expansion of the defect, and what the recurrence means

In the SU(2) character expansion (u = I₂(β)/I₁(β)), log Z is a sum over closed surfaces of plaquettes. Surfaces that do not wrap the torus are extensive and cancel in Δ = 16 log Z(n⁴) − log Z((2n)⁴). The leading wrapping surfaces have Euler characteristic 0, so they carry no dimension factor:
- flat 2-tori in a coordinate plane: area n², 6 planes × n² positions;
- square tubes around one direction: area 4n, 4 directions × 3 cross-section planes × n³ positions.

Hence

  Δ_β(n) ≈ 16 (6 n² u^{n²} + 12 n³ u^{4n}),

the doubled torus contributing only at areas 4n² and 8n. With no fitted parameter this matches the Monte Carlo:

| n, β | u | strong coupling | measured |
|---|---|---|---|
| 2, 0.4 | 0.099 | 0.037 | 0.056(18) |
| 2, 0.6 | 0.148 | 0.184 | 0.187(23) |
| 2, 0.8 | 0.195 | 0.557 | 0.591(27) |
| 2, 1.0 | 0.240 | 1.295 | 1.456(31) (higher orders) |
| 4, 1.8 | 0.399 | 0.006 | −0.003(51) |
| 4, 2.0 | 0.433 | 0.021 | 0.336(63) (crossover) |

This independently validates the β-integration measurements.

**The recurrence is a theorem (YM-O03 §1).** Codex derived, from positivity of the transfer operators of all rectangular tori alone, Δ(2n) ≤ 3740 Δ(n)² whenever Δ(n) ≤ 1/7480. I checked the derivation line by line: the purity estimate p′ ≤ 4p², the exact commutation identity for doublings in two directions, the three-step induction A_j = 3^j δ, and the final sum 4(8 + 36 + 162 + 729) = 3740. Positivity holds for the Wilson action by reflection positivity. So `ym:conj:recurrence` holds with c₁ = 3740 and threshold η* = 1/7480. It is proved in Lean in F01 `Rectangle.lean` (2ed351b), with nonnegative doubling costs and the purity estimate as hypotheses; both follow from the C23 transfer-operator identification.

**Retraction.** An earlier version of this note said that the content of the conditional-gap theorem lies in the recursion, and that the recurrence is in substance gap-equivalent. Both are wrong. The recursion propagates a seed and costs nothing beyond positivity. The entire open content is the seed (`ym:conj:defect`): Δ_β(n(β)) ≤ 1/7480 at a torus of bounded physical size, uniformly as β → ∞.

The strong-coupling form Δ(n) ≈ A n³ e^{−En} gives Δ(2n)/Δ(n)² ≈ 8/(A n³) ≪ 3740, consistent with the theorem.

**Where the seed holds.** From the strong-coupling formula, the seed condition Δ ≤ 1/7480 holds at n = 2 for β ≲ 0.10 and at n = 4 for β ≲ 1.36 (an estimate; these are not measured values). Both lie in the strong-coupling regime, where a gap is already known (Osterwalder–Seiler). Toward the continuum the required n grows like the physical correlation length in lattice units. Whether a seed exists at a bounded physical size is the open question, which numerics at our scale cannot settle.

**The block-gluing estimator** (`src/glued.rs`; exact in 2D, consistent at n = 2 and 4) remains a valid tool, but it is no longer needed to test the recurrence.

**U(1).** Characters e^{ikθ} with u = I₁(β)/I₀(β), and each wrapping surface carries the two orientations k = ±1, so Δ_U(1)(n) ≈ 32 (6 n² u^{n²} + 12 n³ u^{4n}). Against the n = 2 measurements: 0.0753 vs 0.0750(69) at β = 0.2, 0.373 vs 0.383(9) at β = 0.3, 1.14 vs 1.19(1) at β = 0.4 (higher orders begin).

**A provable first instance.** In the strong-coupling regime the polymer (cluster) expansion converges rigorously (Osterwalder–Seiler 1978; Kotecký–Preiss criterion, roadmap T12). There a bound like Δ_β(n) ≤ 16(6 n² u^{n²} + 12 n³ u^{4n})(1 + O(u)) should be provable, giving a seed at n = 4 for small β. With the recurrence theorem, the chain then reproduces the known strong-coupling gap: a consistency check, and a realistic Lean target (YM-F01 follow-up).
