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

## 4D SU(2), n = 4 (preliminary, 2026-11-21 07:00Z)

Data so far:
- 4⁴: 12 Mac replicas plus the laptop replicas, β ∈ [0, 2.6], h = 0.1, 5×10⁵ updates per point.
- 8⁴: GPU replicas (3 on the desktop 3060, 2–3 on the laptop 3070) on β ∈ [1.0, 2.6], 10⁵ updates per point.

The integral starts at β₀ = 1.0. The integrand N_P(8⁴)(⟨a₀⟩₄ − ⟨a₀⟩₈) is consistent with zero for β = 1.0–1.7 (|value| < 1σ, errors about 0.17–0.24), which justifies Δ(β₀) ≈ 0.

| β | integrand | Δ_β(4) − Δ(1.0) | m(8) bound |
|---|---|---|---|
| 1.6 | −0.13(22) | 0.013(47) | — |
| 1.8 | 0.87(31) | 0.023(59) | ≥ 0.94 |
| 2.0 | 4.56(40) | 0.30(8) | ≥ 0.26 |
| 2.1 | 32.6(4) | (Δ explodes) | — |

**Reading.**
- At n = 4 the defect stays below 1/4 up to β ≈ 1.95, against β ≈ 0.65 at n = 2. The threshold moves to weaker coupling as n grows, the predicted trend.
- The thresholds lie below the deconfinement couplings of the small tori (β_c ≈ 1.88 for N_t = 2, ≈ 2.30 for N_t = 4). A confined small torus is necessary but not sufficient: the defect is extensive and picks up finite-size effects earlier.
- In physical units, using a√σ ≈ 0.69 near β = 2.0 and √σ ≈ 440 MeV, the hypothesis holds at n = 4 for a torus about 1.2 fm across.
- The resulting bounds (m(8) ≥ 0.94 at β = 1.8) are far below the actual lattice glueball mass, about 2 near β = 2.0.
- Whether the threshold β*(n) keeps pace with the continuum limit, n·a(β*(n)) bounded as n grows, needs n = 8 (8⁴ against 16⁴). That is beyond the present statistics: the error is dominated by the small volume and grows like N_P(2n)/√N_P(n).

## U(1), the counter-model

Pending. n = 2 on β ∈ [0, 1.6], and an 8⁴ scan for n = 4; U(1) has a massless Coulomb phase above β ≈ 1.0. Expectation: in the Coulomb phase Δ does not become small as n grows. If it did, the criterion would not be detecting a gap.
