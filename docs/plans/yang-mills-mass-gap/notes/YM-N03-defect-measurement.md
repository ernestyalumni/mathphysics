# YM-N03: measuring the hypercubic defect and the four purities

Label: unreviewed (claude-opus55-strat, 2026-11-21). Code: `Rust/lattice-gauge` on branch `ym/YM-N01-claude-opus55-su2`, using `examples/plaquette_scan`, `plaquette_scan_gpu`, `defect_ti` and `defect_window`. Raw data are in the private run directory, not in the repository.

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

All within 1.3σ. In 2D at n = 2 the defect crosses 1/4 near β ≈ 2.45 and log 2 near β ≈ 3.6.

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

**Reading.** Δ_β(2) is below 1/4 only for β ≲ 0.65, at strong coupling. It then grows by two orders of magnitude between β ≈ 1 and 2.
- A 2⁴ torus has two sites in every direction, so it is a lattice at temperature 1/(2a) in each direction.
- SU(2) deconfines at N_t = 2 near β ≈ 1.9, and the small torus is then in the wrong phase.
- The defect measures that transition, not the bulk mass gap.

**Consequence for the conditional gap.** Δ_β(n) can be small only once the n⁴ torus is in the confined phase, roughly n > 1/(a T_c). With T_c/√σ ≈ 0.69 and Teper's a√σ (hep-th/9812187, Table 10), that means n ≳ 4 at β = 2.3 and n ≳ 8 at β = 2.5: a fixed physical size that grows like 1/a toward the continuum.

The hypothesis "Δ small at some scale" is therefore about as strong as finite-size confinement at that scale, and the content of the conditional-gap theorem lies in the recursion from n to 2n (`ym:conj:recurrence` in the YM-F01 note). The test of this reading is n = 4: the prediction is that Δ_β(4) stays small up to β ≈ 2.2–2.3, near the N_t = 4 deconfinement coupling β_c ≈ 2.30.

## 4D SU(2), n = 4

Data:
- 4⁴: 12 replicas on β ∈ [0, 2.6], h = 0.1, 5×10⁵ updates per point (Mac mini).
- 8⁴: 15 GPU replicas on β ∈ [1.0, 2.6], 10⁵ updates per point (3 on the desktop 3060, 12 on the laptop 3070).

The integral starts at β₀ = 1.0. The integrand N_P(8⁴)(⟨a₀⟩₄ − ⟨a₀⟩₈) is consistent with zero for β = 1.0–1.7 (errors 0.15–0.20), which justifies Δ(β₀) ≈ 0.

| β | integrand | Δ_β(4) | m(8) bound |
|---|---|---|---|
| 1.6 | 0.03(19) | −0.012(41) | — |
| 1.8 | 0.57(23) | −0.003(51) | positive (Δ ≈ 0) |
| 2.0 | 4.77(29) | 0.336(63) | ≥ 0.24 |
| 2.2 | 134.4(4) | 9.32(8) | — |
| 2.4 | 107.1(2) | 40.89(9) | — |
| 2.6 | 40.5(2) | 53.43(10) | — |

**Reading.**
- Δ_β(4) < 1/4 for β ≲ 1.97, against β ≲ 0.65 at n = 2. The threshold moves to weaker coupling as n grows.
- The integrand peaks at β ≈ 2.3, the deconfinement coupling of an N_t = 4 lattice: the defect jumps where the small torus deconfines.
- The thresholds lie below those deconfinement couplings (β_c ≈ 1.88 for N_t = 2, 2.30 for N_t = 4). A confined small torus is necessary but not sufficient.
- In physical units (a√σ ≈ 0.69 near β = 2.0, √σ ≈ 440 MeV), the hypothesis holds at n = 4 for a torus about 1.2 fm across.
- At fixed β the defect collapses with n: at β = 1.8, Δ(2) ≈ 26 but Δ(4) ≈ 0.
- The gap bounds are far below the actual lattice glueball mass, about 2 near β = 2.0.
- Whether β*(n) keeps pace with the continuum limit, n·a(β*(n)) bounded, needs n = 8 (8⁴ against 16⁴). That is beyond the present statistics: the small volume dominates the error, which grows like N_P(2n)/√N_P(n).

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
- The integrand peaks near β ≈ 0.9, the deconfinement of the tiny torus.
- In the Coulomb phase the integrand stays large (17.8 at β = 1.6) instead of decaying: the massless photon keeps the finite-size effects.

The decisive comparison is the n-dependence at fixed β. In SU(2) at β = 1.8 the defect collapses from Δ(2) ≈ 26 to Δ(4) ≈ 0. In U(1)'s Coulomb phase a massless theory should give Δ(4) ≳ Δ(2) (zero-mode terms ~ log n in log Z). If instead it collapsed, the criterion would not be detecting a gap. U(1) n = 4 (8⁴ scans on the laptop) is running. n = 2 on β ∈ [0, 1.6], and an 8⁴ scan for n = 4; U(1) has a massless Coulomb phase above β ≈ 1.0. Expectation: in the Coulomb phase Δ does not become small as n grows. If it did, the criterion would not be detecting a gap.

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

**The recurrence.** For large n the tubes dominate, so Δ(n) ≈ A n³ e^{−En} with A = 192 and E = −4 log u. Then Δ(2n)/Δ(n)² = 8/(A n³), small, and `ym:conj:recurrence` holds with room to spare. The same holds for any Δ(n) of the form A n³ e^{−En} with E > 0, whatever E is.

The recurrence is therefore the statement that finite-volume corrections decay exponentially in the torus size, with constants uniform as β grows: in substance, a gap at the doubling scale. In a massless theory (U(1)'s Coulomb phase) finite-size effects decay only as powers, so the ratio grows with n and the recurrence fails.

The conditional-gap theorem thus reformulates the problem as a renormalization-group-type inequality rather than reducing it. A proof of the recurrence uniform in β is of the same nature and difficulty as the gap.

**Numerical outlook.** A direct test at n = 8 (8⁴ against 16⁴) is out of reach with our hardware:
- β-integration would need about 67 GPU-hours;
- the block-gluing estimator (`src/glued.rs`; exact in 2D, consistent at n = 2 and 4: −1.4(1.0) at β 1.8, 2.0(1.2) at β 2.0) is at most 2–3× cheaper;
- multihit and a GPU path would still leave roughly 60–250 GPU-hours.
