import Mathlib

/-!
# Purity of the trace (`ym:prop:purity`)

Source: `documents/mono-notes/topics/yang-mills-lattice.tex`, Proposition `ym:prop:purity`
and Eq. `ym:gap-from-purity`.

The transfer operator is represented only by its eigenvalues: `lam 0 ≥ lam 1 ≥ ⋯ ≥ 0`,
`lam 0 > 0`, summable (trace class). Indices are shifted by one from the source:
`lam 0` is `λ₁` and `lam 1` is `λ₂`.
That these are the eigenvalues of the transfer operator `𝖳`, and that
`Z(n) = Tr 𝖳ⁿ = ∑ λⱼⁿ`, is not formalised here.
-/

open Real

namespace YMFormal

/-- `r_{j+2}` of the source: `(λ_{j+2}/λ₁)ⁿ`, for the excited states only. -/
noncomputable def ratio (lam : ℕ → ℝ) (n j : ℕ) : ℝ := (lam (j + 1) / lam 0) ^ n

/-- `s_n = ∑_{j ≥ 2} r_j`. -/
noncomputable def excited (lam : ℕ → ℝ) (n : ℕ) : ℝ := ∑' j, ratio lam n j

/-- `q_n = ∑_{j ≥ 2} r_j²`. -/
noncomputable def excitedSquares (lam : ℕ → ℝ) (n : ℕ) : ℝ := ∑' j, ratio lam n j ^ 2

/-- `Z(n) = ∑_j λ_jⁿ`. -/
noncomputable def partition (lam : ℕ → ℝ) (n : ℕ) : ℝ := ∑' j, lam j ^ n

/-- `𝖯(n) = Z(2n)/Z(n)²`. -/
noncomputable def purity (lam : ℕ → ℝ) (n : ℕ) : ℝ := partition lam (2 * n) / partition lam n ^ 2

/-- The hypotheses of `ym:prop:purity` on the spectrum. -/
structure Spectrum (lam : ℕ → ℝ) : Prop where
  top_pos : 0 < lam 0
  nonneg : ∀ j, 0 ≤ lam j
  le_top : ∀ j, lam j ≤ lam 0
  summable : Summable lam

/-- Algebraic core of the purity bounds: only `0 ≤ q ≤ s` is used. -/
theorem purity_bounds_of_sums {s q : ℝ} (hq : 0 ≤ q) (hqs : q ≤ s) :
    (1 + q) / (1 + s) ^ 2 ≤ 1 ∧
    s ≤ 1 / ((1 + q) / (1 + s) ^ 2) - 1 ∧
    1 / ((1 + q) / (1 + s) ^ 2) - 1 ≤ 2 * s + s ^ 2 := by
  sorry

/-- Each excited ratio lies in `[0, 1]`. -/
theorem ratio_mem_unit {lam : ℕ → ℝ} (h : Spectrum lam) (n j : ℕ) :
    0 ≤ ratio lam n j ∧ ratio lam n j ≤ 1 := by
  sorry

/-- The excited ratios are summable for `n ≥ 1` (trace class). -/
theorem summable_ratio {lam : ℕ → ℝ} (h : Spectrum lam) {n : ℕ} (hn : 1 ≤ n) :
    Summable (ratio lam n) := by
  sorry

/-- `Z(n) = λ₁ⁿ (1 + s_n)`. -/
theorem partition_eq {lam : ℕ → ℝ} (h : Spectrum lam) {n : ℕ} (hn : 1 ≤ n) :
    partition lam n = lam 0 ^ n * (1 + excited lam n) := by
  sorry

/-- `s_{2n} = q_n`. -/
theorem excited_double (lam : ℕ → ℝ) (n : ℕ) : excited lam (2 * n) = excitedSquares lam n := by
  sorry

/-- `ym:prop:purity`, Eq. `ym:purity-bounds`. -/
theorem purity_bounds {lam : ℕ → ℝ} (h : Spectrum lam) {n : ℕ} (hn : 1 ≤ n) :
    purity lam n = (1 + excitedSquares lam n) / (1 + excited lam n) ^ 2 ∧
    purity lam n ≤ 1 ∧
    excited lam n ≤ 1 / purity lam n - 1 ∧
    1 / purity lam n - 1 ≤ 2 * excited lam n + excited lam n ^ 2 ∧
    excited lam (2 * n) ≤ excited lam n ^ 2 := by
  sorry

/-- Eq. `ym:gap-from-purity`, first part: `e^{-n𝗆} = r₂ ≤ s_n ≤ 𝖯(n)⁻¹ - 1`,
with `𝗆 = log(λ₁/λ₂)`; it needs `λ₂ > 0` (otherwise the gap is infinite). -/
theorem exp_neg_gap_le {lam : ℕ → ℝ} (h : Spectrum lam) (h₂ : 0 < lam 1) {n : ℕ} (hn : 1 ≤ n) :
    exp (-(n * log (lam 0 / lam 1))) ≤ 1 / purity lam n - 1 := by
  sorry

/-- Eq. `ym:gap-from-purity`, second part. The source states it when `𝖯(n) > 1/2`, which
makes the bound positive; the inequality itself holds whenever `λ₂ > 0`. -/
theorem gap_from_purity {lam : ℕ → ℝ} (h : Spectrum lam) (h₂ : 0 < lam 1) {n : ℕ} (hn : 1 ≤ n) :
    (1 / n) * log (1 / (1 / purity lam n - 1)) ≤ log (lam 0 / lam 1) := by
  sorry

end YMFormal
