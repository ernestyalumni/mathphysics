import YMFormal.Rectangle

/-!
# The strong-coupling seed and the infinite-volume hypothesis (F02, option C)

Ernest, 2026-11-22, approved option C: "Seed bound plus the explicit infinite-volume hypothesis
only. A small step, and it adds the hypothesis Codex asked to have stated."

**Seed (Codex, review YM-O04).**
- With a probability measure (normalized Haar) and plaquette terms `|a_p| ≤ 1`, a box with `P`
  plaquettes has `e^{-βP} ≤ Z ≤ e^{βP}`, so `|log Z| ≤ βP`.
- A 4-dimensional periodic box with sides `x` has `P(x) = 6 ∏ x_k` plaquettes, so `P(2n) = 16 P(n)`
  and `Δ(n) = 16 log Z(n) - log Z(2n) ≤ 32 β P(n) = 192 β n⁴`.
- Shifting the action by a constant per plaquette, as in `β Σ_p (1 - a_p)`, changes `log Z` by a
  multiple of `P` and leaves `Δ` unchanged, because `P(2n) = 16 P(n)`.

**Infinite volume (Codex, review YM-O04 §2).** The finite-volume chain bounds `λ₂/λ₁` on the dyadic
boxes. Passing to an infinite-volume state is a separate hypothesis. `LocalLimit` states it
explicitly. It is assumed, not proved.

The fields:
- `corr k t`: the vacuum time correlation of a bounded, centred local observable `A` on the `k`-th
  box;
- `finite_bound`: `|corr k t| ≤ ‖A‖² (λ₂/λ₁)ᵗ`, from the spectral theorem on each box;
- `limit`, `tendsto`: the infinite-volume correlation is the pointwise limit at each fixed `t`.

**Not covered.**
- The gap of the whole reconstructed sector, which needs cyclic or dense states and control of the
  vacuum subtraction.
- The identification of the limiting state with that of a particular theorem.
-/

namespace YMFormal

open Real MeasureTheory Filter Topology

/-- `|log Z| ≤ β P` for `Z = ∫ exp (β S) dμ` with `μ` a probability measure and `|S| ≤ P`. -/
theorem log_partition_bound {Ω : Type*} [MeasurableSpace Ω] (μ : Measure Ω) [IsProbabilityMeasure μ]
    {S : Ω → ℝ} (hS : Measurable S) {P β : ℝ} (hβ : 0 ≤ β) (hbound : ∀ ω, |S ω| ≤ P) :
    |log (∫ ω, exp (β * S ω) ∂μ)| ≤ β * P := by
  have hmeas : AEStronglyMeasurable (fun ω => exp (β * S ω)) μ :=
    (measurable_exp.comp (hS.const_mul β)).aestronglyMeasurable
  have hpt : ∀ ω, exp (-(β * P)) ≤ exp (β * S ω) ∧ exp (β * S ω) ≤ exp (β * P) := by
    intro ω
    have h := abs_le.mp (hbound ω)
    constructor
    · exact exp_le_exp.mpr (by nlinarith [h.1])
    · exact exp_le_exp.mpr (by nlinarith [h.2])
  have hint : Integrable (fun ω => exp (β * S ω)) μ :=
    (integrable_const (exp (β * P))).mono' hmeas
      (ae_of_all _ fun ω => by rw [Real.norm_eq_abs, abs_of_pos (exp_pos _)]; exact (hpt ω).2)
  have hlo : exp (-(β * P)) ≤ ∫ ω, exp (β * S ω) ∂μ := by
    have := integral_mono (integrable_const (exp (-(β * P)))) hint (fun ω => (hpt ω).1)
    simpa using this
  have hhi : ∫ ω, exp (β * S ω) ∂μ ≤ exp (β * P) := by
    have := integral_mono hint (integrable_const (exp (β * P))) (fun ω => (hpt ω).2)
    simpa using this
  have hZ : 0 < ∫ ω, exp (β * S ω) ∂μ := lt_of_lt_of_le (exp_pos _) hlo
  rw [abs_le]
  constructor
  · have := log_le_log (exp_pos _) hlo
    rw [log_exp] at this
    linarith
  · have := log_le_log hZ hhi
    rwa [log_exp] at this

/-- The Wilson form: a sum of `P` plaquette terms, each of absolute value at most 1. -/
theorem log_partition_bound_plaquettes {Ω ι : Type*} [MeasurableSpace Ω] (μ : Measure Ω)
    [IsProbabilityMeasure μ] (plaquettes : Finset ι) {a : ι → Ω → ℝ} (ha : ∀ p, Measurable (a p))
    (hbound : ∀ p ω, |a p ω| ≤ 1) {β : ℝ} (hβ : 0 ≤ β) :
    |log (∫ ω, exp (β * ∑ p ∈ plaquettes, a p ω) ∂μ)| ≤ β * plaquettes.card := by
  refine log_partition_bound μ (Finset.measurable_sum _ fun p _ => ha p) hβ (fun ω => ?_)
  calc |∑ p ∈ plaquettes, a p ω| ≤ ∑ p ∈ plaquettes, |a p ω| := Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ _p ∈ plaquettes, (1 : ℝ) := Finset.sum_le_sum fun p _ => hbound p ω
    _ = plaquettes.card := by simp

/-- **The seed bound.** If `|log Z(x)| ≤ β · 6 ∏ x_k` on boxes with positive sides, then
`Δ(n) ≤ 192 β n⁴`. -/
theorem seed_of_bound (F : (Fin 4 → ℕ) → ℝ) {β : ℝ} {n : ℕ} (hn : 1 ≤ n)
    (hF : ∀ x : Fin 4 → ℕ, (∀ k, 0 < x k) → |F x| ≤ β * (6 * ∏ k, (x k : ℝ))) :
    boxDefect F n ≤ 192 * β * n ^ 4 := by
  have h1 := hF (fun _ => n) (fun _ => by omega)
  have h2 := hF (fun _ => 2 * n) (fun _ => by omega)
  simp only [Finset.prod_const, Finset.card_univ, Fintype.card_fin, Nat.cast_mul, Nat.cast_ofNat]
    at h1 h2
  unfold boxDefect
  have a1 := (abs_le.mp h1).1
  have a2 := (abs_le.mp h2).1
  have b1 := (abs_le.mp h1).2
  nlinarith

/-- For a probability measure and a measurable weight `lo ≤ w ≤ hi` with `0 < lo`, the integral `Z`
satisfies `log lo ≤ log Z ≤ log hi`. -/
theorem log_integral_bounds {Ω : Type*} [MeasurableSpace Ω] (μ : Measure Ω) [IsProbabilityMeasure μ]
    {w : Ω → ℝ} (hw : Measurable w) {lo hi : ℝ} (hlo : 0 < lo) (hb : ∀ ω, lo ≤ w ω ∧ w ω ≤ hi) :
    log lo ≤ log (∫ ω, w ω ∂μ) ∧ log (∫ ω, w ω ∂μ) ≤ log hi := by
  have hint : Integrable w μ :=
    (integrable_const hi).mono' hw.aestronglyMeasurable (ae_of_all _ fun ω => by
      rw [Real.norm_eq_abs, abs_of_pos (hlo.trans_le (hb ω).1)]; exact (hb ω).2)
  have h1 : lo ≤ ∫ ω, w ω ∂μ := by
    simpa using integral_mono (integrable_const lo) hint (fun ω => (hb ω).1)
  have h2 : ∫ ω, w ω ∂μ ≤ hi := by
    simpa using integral_mono hint (integrable_const hi) (fun ω => (hb ω).2)
  exact ⟨log_le_log hlo h1, log_le_log (hlo.trans_le h1) h2⟩

/-- **The handbook convention.** The handbook's weight is `e^{-β Σ_p (1 - a_p)}`, with
`a_p = Re tr U_p / m ∈ [-1, 1]`. Then `-2βP ≤ log Z ≤ 0`. -/
theorem log_partition_bound_wilson {Ω ι : Type*} [MeasurableSpace Ω] (μ : Measure Ω)
    [IsProbabilityMeasure μ] (plaquettes : Finset ι) {a : ι → Ω → ℝ} (ha : ∀ p, Measurable (a p))
    (hbound : ∀ p ω, |a p ω| ≤ 1) {β : ℝ} (hβ : 0 ≤ β) :
    -(2 * β * plaquettes.card) ≤ log (∫ ω, exp (-(β * ∑ p ∈ plaquettes, (1 - a p ω))) ∂μ) ∧
      log (∫ ω, exp (-(β * ∑ p ∈ plaquettes, (1 - a p ω))) ∂μ) ≤ 0 := by
  have hmeas : Measurable fun ω => exp (-(β * ∑ p ∈ plaquettes, (1 - a p ω))) :=
    measurable_exp.comp ((Finset.measurable_sum _ fun p _ => measurable_const.sub (ha p)).const_mul β).neg
  have hsum : ∀ ω, 0 ≤ ∑ p ∈ plaquettes, (1 - a p ω) ∧
      ∑ p ∈ plaquettes, (1 - a p ω) ≤ 2 * plaquettes.card := by
    intro ω
    constructor
    · exact Finset.sum_nonneg fun p _ => by linarith [(abs_le.mp (hbound p ω)).2]
    · calc ∑ p ∈ plaquettes, (1 - a p ω) ≤ ∑ _p ∈ plaquettes, (2 : ℝ) :=
            Finset.sum_le_sum fun p _ => by linarith [(abs_le.mp (hbound p ω)).1]
        _ = 2 * plaquettes.card := by simp [mul_comm]
  have h := log_integral_bounds μ hmeas (exp_pos (-(2 * β * plaquettes.card))) (fun ω => by
    constructor
    · exact exp_le_exp.mpr (by nlinarith [(hsum ω).2])
    · rw [exp_le_one_iff]; nlinarith [(hsum ω).1])
  rw [log_exp, log_one] at h
  exact h

/-- **The seed in the handbook convention.** If `-2β · 6∏x_k ≤ log Z(x) ≤ 0` on boxes with positive
sides, then `Δ(n) ≤ 0 + 2β P(2n) = 192 β n⁴`. -/
theorem seed_of_wilson_bound (F : (Fin 4 → ℕ) → ℝ) {β : ℝ} {n : ℕ} (hn : 1 ≤ n)
    (hF : ∀ x : Fin 4 → ℕ, (∀ k, 0 < x k) →
      -(2 * β * (6 * ∏ k, (x k : ℝ))) ≤ F x ∧ F x ≤ 0) :
    boxDefect F n ≤ 192 * β * n ^ 4 := by
  have h1 := hF (fun _ => n) (fun _ => by omega)
  have h2 := hF (fun _ => 2 * n) (fun _ => by omega)
  simp only [Finset.prod_const, Finset.card_univ, Fintype.card_fin, Nat.cast_mul, Nat.cast_ofNat]
    at h1 h2
  unfold boxDefect
  nlinarith [h1.2, h2.1]

/-- At strong coupling the seed holds: `Δ(n) ≤ 192 β n⁴` and `β ≤ 1/(7480 · 192 · n⁴)` give
`Δ(n) ≤ 1/7480`. -/
theorem seed_strong_coupling (F : (Fin 4 → ℕ) → ℝ) {β : ℝ} {n : ℕ} (hn : 1 ≤ n)
    (hβ : β ≤ 1 / (7480 * 192 * (n : ℝ) ^ 4)) (hb : boxDefect F n ≤ 192 * β * n ^ 4) :
    boxDefect F n ≤ 1 / 7480 := by
  have hn4 : (0 : ℝ) < (n : ℝ) ^ 4 := by
    have : (0 : ℝ) < n := by exact_mod_cast hn
    positivity
  have : 192 * β * (n : ℝ) ^ 4 ≤ 1 / 7480 := by
    have h := mul_le_mul_of_nonneg_right hβ (le_of_lt (by positivity : (0 : ℝ) < 192 * (n : ℝ) ^ 4))
    have e : 1 / (7480 * 192 * (n : ℝ) ^ 4) * (192 * (n : ℝ) ^ 4) = 1 / 7480 := by
      field_simp
    nlinarith
  linarith

/-- **The chain at strong coupling, in the handbook convention.** For
`0 ≤ β ≤ 1/(7480 · 192 · n⁴)`, the inputs are:
- the plaquette bound `-2β · 6∏x_k ≤ log Z(x) ≤ 0` (`log_partition_bound_wilson`);
- the spectral hypotheses `nonneg`, `estimate` and the ratio bound at each scale.

The conclusion: `(λ₂/λ₁)ⁿ ≤ 1/2` on every dyadic box `2ᵏ n`. -/
theorem conditional_gap_strong_coupling (F : (Fin 4 → ℕ) → ℝ) {β : ℝ} {n : ℕ} (hn : 1 ≤ n)
    (hβ : β ≤ 1 / (7480 * 192 * (n : ℝ) ^ 4))
    (hF : ∀ x : Fin 4 → ℕ, (∀ k, 0 < x k) →
      -(2 * β * (6 * ∏ k, (x k : ℝ))) ≤ F x ∧ F x ≤ 0)
    {r : ℕ → ℝ} (hr : ∀ k, 0 ≤ r k)
    (nonneg : ∀ i x, (∀ k, 0 < x k) → 0 ≤ sideCost F i x)
    (estimate : ∀ i x, (∀ k, 0 < x k) →
      sideCost F i x ≤ 1 / 4 → sideCost F i (sideDouble i x) ≤ 4 * sideCost F i x ^ 2)
    (from_defect : ∀ k, r k ^ (2 ^ k * n) ≤ exp (boxDefect F (2 ^ k * n)) - 1) :
    ∀ k, r k ^ n ≤ 1 / 2 :=
  conditional_gap_of_boxes F hn hr nonneg estimate
    (seed_strong_coupling F hn hβ (seed_of_wilson_bound F hn hF)) from_defect

/-- **The infinite-volume hypothesis of review YM-O04 §2, stated explicitly.** It is assumed, not
proved. For a bounded, centred local observable `A` with `‖A‖² = normSq`, `corr k t` is its vacuum
time correlation on the `k`-th dyadic box. The spectral theorem on that box gives `finite_bound`, with
`r k = λ₂/λ₁` there. The infinite-volume correlation `limit t` is the pointwise limit along the boxes
at each fixed time. -/
structure LocalLimit where
  normSq : ℝ
  normSq_nonneg : 0 ≤ normSq
  r : ℕ → ℝ
  r_nonneg : ∀ k, 0 ≤ r k
  corr : ℕ → ℕ → ℝ
  finite_bound : ∀ k t, |corr k t| ≤ normSq * r k ^ t
  limit : ℕ → ℝ
  tendsto : ∀ t, Tendsto (fun k => corr k t) atTop (𝓝 (limit t))

/-- From the finite-volume gap `(λ₂/λ₁)ⁿ ≤ 1/2` on every box, the infinite-volume correlation decays
at rate `log 2 / n`: `|C(t)| ≤ ‖A‖² 2^{-t/n}`. -/
theorem LocalLimit.decay (L : LocalLimit) {n : ℕ} (hn : 1 ≤ n) (gap : ∀ k, L.r k ^ n ≤ 1 / 2) (t : ℕ) :
    |L.limit t| ≤ L.normSq * (1 / 2 : ℝ) ^ ((t : ℝ) / n) := by
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  -- each box: `r ≤ (1/2)^(1/n)`, hence `rᵗ ≤ (1/2)^(t/n)`
  have hr : ∀ k, L.r k ^ t ≤ (1 / 2 : ℝ) ^ ((t : ℝ) / n) := by
    intro k
    have h0 := L.r_nonneg k
    have hroot : L.r k ≤ (1 / 2 : ℝ) ^ ((1 : ℝ) / n) := by
      have := Real.rpow_le_rpow (pow_nonneg h0 n) (gap k) (by positivity : (0 : ℝ) ≤ 1 / n)
      rwa [← Real.rpow_natCast, ← Real.rpow_mul h0, mul_one_div_cancel hn'.ne', Real.rpow_one] at this
    calc L.r k ^ t ≤ ((1 / 2 : ℝ) ^ ((1 : ℝ) / n)) ^ t := pow_le_pow_left₀ h0 hroot t
      _ = (1 / 2 : ℝ) ^ ((t : ℝ) / n) := by
        rw [← Real.rpow_natCast, ← Real.rpow_mul (by norm_num)]
        ring_nf
  have hbound : ∀ k, |L.corr k t| ≤ L.normSq * (1 / 2 : ℝ) ^ ((t : ℝ) / n) := fun k =>
    (L.finite_bound k t).trans (mul_le_mul_of_nonneg_left (hr k) L.normSq_nonneg)
  exact le_of_tendsto' ((L.tendsto t).abs) hbound

end YMFormal
