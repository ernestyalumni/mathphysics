import YMFormal.Defect
import YMFormal.ConditionalGap

/-!
# The gap chain in ratio form, including rank one

Codex's review YM-O03 (2026-11-21, §1, "Endpoint correction") points out that at `β = 0` the
transfer operator has rank one: `λ₂ = 0`, `𝖯 = 1` and `Δ = 0`, so the logarithmic statements
`ym:gap-from-purity`, `ym:gap-from-defect` and `ym:prop:conditional-gap` divide by zero or take the
log of zero. The Lean versions avoid this by hypothesis (`0 < lam 1`, `0 < Δ`, `0 < δ k`), which
excludes rank one.

Here the same chain is stated for the ratio `r = λ₂/λ₁ ∈ [0, 1]` with no positivity hypotheses.
Rank one is the case `r = 0`, an infinite gap. The logarithmic forms are the corollaries for
`r > 0`, since `𝗆 = -log r`.
-/

namespace YMFormal

open Real

/-- `ym:gap-from-purity` in ratio form: `(λ₂/λ₁)ⁿ ≤ 𝖯⁻¹ - 1`, with no assumption `λ₂ > 0`. -/
theorem ratio_pow_le_impurity {lam : ℕ → ℝ} (h : Spectrum lam) {n : ℕ} (hn : 1 ≤ n) :
    (lam 1 / lam 0) ^ n ≤ 1 / purity lam n - 1 := by
  have hunit := ratio_mem_unit h n
  calc (lam 1 / lam 0) ^ n = ratio lam n 0 := rfl
    _ ≤ excited lam n := (summable_ratio h hn).le_tsum 0 (fun i _ => (hunit i).1)
    _ ≤ 1 / purity lam n - 1 := (purity_bounds h hn).2.2.1

/-- `ym:gap-from-defect` in ratio form: if `e^{-Δ} ≤ 𝖯(n)` then `(λ₂/λ₁)ⁿ ≤ e^Δ - 1`.
No assumption `Δ > 0` or `λ₂ > 0`; `Δ = 0` forces `λ₂ = 0`. -/
theorem ratio_pow_le_defect {lam : ℕ → ℝ} (h : Spectrum lam) {n : ℕ} (hn : 1 ≤ n) {Δ : ℝ}
    (hP : exp (-Δ) ≤ purity lam n) : (lam 1 / lam 0) ^ n ≤ exp Δ - 1 := by
  have hPpos : 0 < purity lam n := lt_of_lt_of_le (exp_pos _) hP
  have hinv : 1 / purity lam n ≤ exp Δ := by
    rw [div_le_iff₀ hPpos]
    calc (1 : ℝ) = exp Δ * exp (-Δ) := by rw [← exp_add]; simp
      _ ≤ exp Δ * purity lam n := mul_le_mul_of_nonneg_left hP (exp_pos _).le
  linarith [ratio_pow_le_impurity h hn]

/-- `ym:prop:conditional-gap` in ratio form. `r k` stands for `λ₂/λ₁` at the `k`-th doubling.
Defects may vanish (`0 ≤ δ k`); the conclusion `(r k)ⁿ ≤ 1/2` is `𝗆 ≥ log 2 / n` when `r k > 0`
and an infinite gap when `r k = 0`. It needs no `n ≥ 1` either. -/
theorem conditional_gap_ratio {c₁ : ℝ} (hc : 2 ≤ c₁) {n : ℕ} {δ r : ℕ → ℝ}
    (hδ : ∀ k, 0 ≤ δ k) (hr : ∀ k, 0 ≤ r k)
    (recurrence : ∀ k, δ k ≤ 1 / (2 * c₁) → δ (k + 1) ≤ c₁ * δ k ^ 2)
    (start : δ 0 ≤ 1 / (2 * c₁))
    (from_defect : ∀ k, r k ^ (2 ^ k * n) ≤ exp (δ k) - 1) :
    ∀ k, r k ^ n ≤ 1 / 2 := by
  intro k
  have hc₀ : 0 < c₁ := by linarith
  have hdouble := defect_double_exponential hc hδ recurrence start k
  have hhalf : (0 : ℝ) ≤ (1 / 2) ^ 2 ^ k := by positivity
  have hle1 : (1 / 2 : ℝ) ^ 2 ^ k ≤ 1 := pow_le_one₀ (by norm_num) (by norm_num)
  -- δ k ≤ (1/2)^(2^k) / c₁ ≤ (1/2)^(2^k) / 2 ≤ 1
  have hδk : δ k ≤ (1 / 2) ^ 2 ^ k / c₁ := by rw [le_div_iff₀ hc₀]; linarith
  have hδk' : δ k ≤ (1 / 2) ^ 2 ^ k / 2 :=
    hδk.trans (div_le_div_of_nonneg_left hhalf (by norm_num) hc)
  have hexp : exp (δ k) - 1 ≤ (1 / 2) ^ 2 ^ k := by
    have := exp_sub_one_le_two_mul (hδ k) (by linarith)
    linarith
  have hpow : (r k ^ n) ^ 2 ^ k ≤ (1 / 2) ^ 2 ^ k := by
    rw [← pow_mul, mul_comm]
    exact (from_defect k).trans hexp
  exact (pow_le_pow_iff_left₀ (pow_nonneg (hr k) n) (by norm_num) (pow_ne_zero _ two_ne_zero)).mp hpow

/-- The rank-one endpoint (`β = 0`): a zero defect forces `λ₂ = 0`, an infinite gap. -/
theorem rank_one_of_zero_defect {lam : ℕ → ℝ} (h : Spectrum lam) {n : ℕ} (hn : 1 ≤ n)
    (hP : exp (-0) ≤ purity lam n) : lam 1 = 0 := by
  have hle := ratio_pow_le_defect h hn hP
  simp only [exp_zero, sub_self] at hle
  have hr : 0 ≤ lam 1 / lam 0 := div_nonneg (h.nonneg 1) h.top_pos.le
  have hz : (lam 1 / lam 0) ^ n = 0 := le_antisymm hle (pow_nonneg hr n)
  have := pow_eq_zero_iff (n := n) (by omega) |>.mp hz
  rcases div_eq_zero_iff.mp this with h1 | h0
  · exact h1
  · exact absurd h0 h.top_pos.ne'

end YMFormal
