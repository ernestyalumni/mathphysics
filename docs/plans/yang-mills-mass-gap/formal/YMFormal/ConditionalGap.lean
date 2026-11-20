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
  sorry

/-- The induction of the source's proof: `c₁ δ_k ≤ 2^{-2^k}`. -/
theorem defect_double_exponential {c₁ : ℝ} (hc : 2 ≤ c₁) {δ : ℕ → ℝ}
    (hδ : ∀ k, 0 ≤ δ k)
    (recurrence : ∀ k, δ k ≤ 1 / (2 * c₁) → δ (k + 1) ≤ c₁ * δ k ^ 2)
    (start : δ 0 ≤ 1 / (2 * c₁)) :
    ∀ k, c₁ * δ k ≤ (1 / 2) ^ (2 ^ k) := by
  sorry

/-- Proposition `ym:prop:conditional-gap`, Eq. `ym:conditional-gap`. -/
theorem conditional_gap {c₁ : ℝ} (hc : 2 ≤ c₁) {n : ℕ} (hn : 1 ≤ n) {δ gap : ℕ → ℝ}
    (hδ : ∀ k, 0 < δ k)
    (recurrence : ∀ k, δ k ≤ 1 / (2 * c₁) → δ (k + 1) ≤ c₁ * δ k ^ 2)
    (start : δ 0 ≤ 1 / (2 * c₁))
    (from_defect : ∀ k, (1 / (2 ^ k * n)) * log (1 / (exp (δ k) - 1)) ≤ gap k) :
    ∀ k, log 2 / n ≤ gap k := by
  sorry

end YMFormal
