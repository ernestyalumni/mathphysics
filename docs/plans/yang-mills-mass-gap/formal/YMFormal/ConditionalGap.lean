import Mathlib

/-!
# Conditional lattice gap (`ym:prop:conditional-gap`)

Source: `documents/mono-notes/topics/yang-mills-lattice.tex`, Conjecture `ym:conj:recurrence`
and Proposition `ym:prop:conditional-gap`.

`δ k` is the defect `Δ_β(2ᵏ n)` and `gap k` is the transfer gap `𝗆(2^{k+1} n)`. Two inputs are
hypotheses:
- `recurrence` is Conjecture `ym:conj:recurrence`, which is open.
- `from_defect` is Eq. `ym:gap-from-defect` at side `2ᵏ n`, formalised as `gap_from_defect`
  in `Defect.lean`.

The defects are assumed positive. A zero defect means a purity equal to 1, so the
second eigenvalue is 0 and the gap is infinite.
-/

open Real

namespace YMFormal

/-- `e^x - 1 ≤ 2x` for `0 ≤ x ≤ 1`, the convexity step of the source's proof. -/
theorem exp_sub_one_le_two_mul {x : ℝ} (h₀ : 0 ≤ x) (h₁ : x ≤ 1) : exp x - 1 ≤ 2 * x := by
  have h := Real.exp_bound (x := x) (by rw [abs_of_nonneg h₀]; exact h₁) (n := 1) one_pos
  norm_num [abs_of_nonneg h₀] at h
  linarith [le_abs_self (exp x - 1)]

/-- The induction of the source's proof: `c₁ δ_k ≤ 2^{-2^k}`. -/
theorem defect_double_exponential {c₁ : ℝ} (hc : 2 ≤ c₁) {δ : ℕ → ℝ}
    (hδ : ∀ k, 0 ≤ δ k)
    (recurrence : ∀ k, δ k ≤ 1 / (2 * c₁) → δ (k + 1) ≤ c₁ * δ k ^ 2)
    (start : δ 0 ≤ 1 / (2 * c₁)) :
    ∀ k, c₁ * δ k ≤ (1 / 2) ^ (2 ^ k) := by
  intro k
  have hc₀ : 0 < c₁ := by linarith
  induction k with
  | zero =>
    have h := mul_le_mul_of_nonneg_left start hc₀.le
    rw [show c₁ * (1 / (2 * c₁)) = 1 / 2 by field_simp] at h
    simpa using h
  | succ k ih =>
    have hhalf : (1 / 2 : ℝ) ^ (2 ^ k) ≤ 1 / 2 := by
      calc (1 / 2 : ℝ) ^ (2 ^ k) ≤ (1 / 2) ^ 1 :=
            pow_le_pow_of_le_one (by norm_num) (by norm_num) Nat.one_le_two_pow
        _ = 1 / 2 := by norm_num
    have hk : δ k ≤ 1 / (2 * c₁) := by
      rw [le_div_iff₀ (by linarith)]
      nlinarith
    calc c₁ * δ (k + 1) ≤ c₁ * (c₁ * δ k ^ 2) := mul_le_mul_of_nonneg_left (recurrence k hk) hc₀.le
      _ = (c₁ * δ k) ^ 2 := by ring
      _ ≤ ((1 / 2) ^ (2 ^ k)) ^ 2 := pow_le_pow_left₀ (mul_nonneg hc₀.le (hδ k)) ih 2
      _ = (1 / 2) ^ (2 ^ (k + 1)) := by rw [← pow_mul, pow_succ]

/-- Proposition `ym:prop:conditional-gap`, Eq. `ym:conditional-gap`. -/
theorem conditional_gap {c₁ : ℝ} (hc : 2 ≤ c₁) {n : ℕ} (hn : 1 ≤ n) {δ gap : ℕ → ℝ}
    (hδ : ∀ k, 0 < δ k)
    (recurrence : ∀ k, δ k ≤ 1 / (2 * c₁) → δ (k + 1) ≤ c₁ * δ k ^ 2)
    (start : δ 0 ≤ 1 / (2 * c₁))
    (from_defect : ∀ k, (1 / (2 ^ k * n)) * log (1 / (exp (δ k) - 1)) ≤ gap k) :
    ∀ k, log 2 / n ≤ gap k := by
  intro k
  have hc₀ : 0 < c₁ := by linarith
  have hd := defect_double_exponential hc (fun k => (hδ k).le) recurrence start k
  have hpos := hδ k
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  have hle : δ k ≤ 1 := by
    have : (1 / 2 : ℝ) ^ (2 ^ k) ≤ 1 := pow_le_one₀ (by norm_num) (by norm_num)
    nlinarith
  have hexp := exp_sub_one_le_two_mul hpos.le hle
  have hexp_pos : 0 < exp (δ k) - 1 := by linarith [Real.add_one_lt_exp hpos.ne']
  have h₁ : log (1 / (2 * δ k)) ≤ log (1 / (exp (δ k) - 1)) :=
    log_le_log (by positivity) (one_div_le_one_div_of_le hexp_pos hexp)
  have h₂ : c₁ * 2 ^ (2 ^ k) / 2 ≤ 1 / (2 * δ k) := by
    have hprod : c₁ * δ k * 2 ^ (2 ^ k) ≤ 1 := by
      calc c₁ * δ k * 2 ^ (2 ^ k) ≤ (1 / 2) ^ (2 ^ k) * 2 ^ (2 ^ k) :=
            mul_le_mul_of_nonneg_right hd (by positivity)
        _ = 1 := by rw [← mul_pow]; norm_num
    rw [div_le_div_iff₀ (by norm_num) (by positivity)]
    nlinarith
  have h₃ : log (c₁ * 2 ^ (2 ^ k) / 2) = log (c₁ / 2) + 2 ^ k * log 2 := by
    rw [show c₁ * 2 ^ (2 ^ k) / 2 = (c₁ / 2) * 2 ^ (2 ^ k) by ring,
      log_mul (by positivity) (by positivity), log_pow]
    push_cast
    ring
  have h₄ : 0 ≤ log (c₁ / 2) := log_nonneg (by linarith)
  have h₅ := log_le_log (by positivity) h₂
  have hL : 2 ^ k * log 2 ≤ log (1 / (exp (δ k) - 1)) := by linarith
  calc log 2 / n = 1 / (2 ^ k * n) * (2 ^ k * log 2) := by field_simp
    _ ≤ 1 / (2 ^ k * n) * log (1 / (exp (δ k) - 1)) :=
        mul_le_mul_of_nonneg_left hL (by positivity)
    _ ≤ gap k := from_defect k

end YMFormal
