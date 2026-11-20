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
  have hs : 0 ≤ s := hq.trans hqs
  have h1q : 0 < 1 + q := by linarith
  have h1s : 0 < 1 + s := by linarith
  have hpos : 0 < (1 + s) ^ 2 := by positivity
  refine ⟨?_, ?_, ?_⟩
  · rw [div_le_one hpos]; nlinarith
  · rw [one_div_div, le_sub_iff_add_le, le_div_iff₀ h1q]
    nlinarith [mul_le_mul_of_nonneg_left hqs h1s.le]
  · rw [one_div_div, sub_le_iff_le_add, div_le_iff₀ h1q]
    nlinarith [mul_nonneg hpos.le hq]

/-- Each excited ratio lies in `[0, 1]`. -/
theorem ratio_mem_unit {lam : ℕ → ℝ} (h : Spectrum lam) (n j : ℕ) :
    0 ≤ ratio lam n j ∧ ratio lam n j ≤ 1 := by
  have hb₀ : 0 ≤ lam (j + 1) / lam 0 := div_nonneg (h.nonneg _) h.top_pos.le
  have hb₁ : lam (j + 1) / lam 0 ≤ 1 := (div_le_one h.top_pos).mpr (h.le_top _)
  exact ⟨pow_nonneg hb₀ n, pow_le_one₀ hb₀ hb₁⟩

/-- The excited ratios are summable for `n ≥ 1` (trace class). -/
theorem summable_ratio {lam : ℕ → ℝ} (h : Spectrum lam) {n : ℕ} (hn : 1 ≤ n) :
    Summable (ratio lam n) := by
  have hshift : Summable (fun j => lam (j + 1) / lam 0) :=
    ((summable_nat_add_iff 1).mpr h.summable).div_const _
  refine Summable.of_nonneg_of_le (fun j => (ratio_mem_unit h n j).1) (fun j => ?_) hshift
  exact pow_le_of_le_one (div_nonneg (h.nonneg _) h.top_pos.le)
    ((div_le_one h.top_pos).mpr (h.le_top _)) (by omega)

/-- `Z(n) = λ₁ⁿ (1 + s_n)`. -/
theorem partition_eq {lam : ℕ → ℝ} (h : Spectrum lam) {n : ℕ} (hn : 1 ≤ n) :
    partition lam n = lam 0 ^ n * (1 + excited lam n) := by
  have htop : lam 0 ^ n ≠ 0 := pow_ne_zero _ h.top_pos.ne'
  have hterm : ∀ j, lam j ^ n = lam 0 ^ n * (lam j / lam 0) ^ n := fun j => by
    rw [div_pow, mul_div_cancel₀ _ htop]
  have hsum : Summable (fun j => lam j ^ n) := by
    have hratio : Summable (fun j => (lam j / lam 0) ^ n) := by
      have hlam : Summable (fun j => lam j / lam 0) := h.summable.div_const _
      refine Summable.of_nonneg_of_le (fun j => pow_nonneg (div_nonneg (h.nonneg _) h.top_pos.le) n)
        (fun j => ?_) hlam
      exact pow_le_of_le_one (div_nonneg (h.nonneg _) h.top_pos.le)
        ((div_le_one h.top_pos).mpr (h.le_top _)) (by omega)
    simpa [← hterm] using hratio.mul_left (lam 0 ^ n)
  unfold partition excited ratio
  rw [hsum.tsum_eq_zero_add]
  simp_rw [fun j => hterm (j + 1), tsum_mul_left]
  ring

/-- `s_{2n} = q_n`. -/
theorem excited_double (lam : ℕ → ℝ) (n : ℕ) : excited lam (2 * n) = excitedSquares lam n := by
  unfold excited excitedSquares ratio
  congr 1
  funext j
  rw [mul_comm, pow_mul]

/-- `ym:prop:purity`, Eq. `ym:purity-bounds`. -/
theorem purity_bounds {lam : ℕ → ℝ} (h : Spectrum lam) {n : ℕ} (hn : 1 ≤ n) :
    purity lam n = (1 + excitedSquares lam n) / (1 + excited lam n) ^ 2 ∧
    purity lam n ≤ 1 ∧
    excited lam n ≤ 1 / purity lam n - 1 ∧
    1 / purity lam n - 1 ≤ 2 * excited lam n + excited lam n ^ 2 ∧
    excited lam (2 * n) ≤ excited lam n ^ 2 := by
  have h2n : 1 ≤ 2 * n := by omega
  have hsum := summable_ratio h hn
  have hunit := ratio_mem_unit h n
  have hsq : Summable (fun j => ratio lam n j ^ 2) :=
    Summable.of_nonneg_of_le (fun j => sq_nonneg _)
      (fun j => by nlinarith [hunit j]) hsum
  have hq : 0 ≤ excitedSquares lam n := tsum_nonneg (fun j => sq_nonneg _)
  have hs : 0 ≤ excited lam n := tsum_nonneg (fun j => (hunit j).1)
  have hqs : excitedSquares lam n ≤ excited lam n :=
    Summable.tsum_le_tsum (fun j => by nlinarith [hunit j]) hsq hsum
  have hformula : purity lam n = (1 + excitedSquares lam n) / (1 + excited lam n) ^ 2 := by
    unfold purity
    rw [partition_eq h h2n, partition_eq h hn, excited_double]
    have htop : lam 0 ^ n ≠ 0 := pow_ne_zero _ h.top_pos.ne'
    have h1s : (1 + excited lam n) ≠ 0 := by linarith
    rw [pow_mul', mul_pow]
    field_simp
  obtain ⟨b₁, b₂, b₃⟩ := purity_bounds_of_sums hq hqs
  refine ⟨hformula, hformula ▸ b₁, hformula ▸ b₂, hformula ▸ b₃, ?_⟩
  rw [excited_double]
  have hle : ∀ j, ratio lam n j ≤ excited lam n := fun j =>
    hsum.le_tsum j (fun i _ => (hunit i).1)
  calc excitedSquares lam n = ∑' j, ratio lam n j ^ 2 := rfl
    _ ≤ ∑' j, ratio lam n j * excited lam n :=
        Summable.tsum_le_tsum (fun j => by rw [sq]; exact mul_le_mul_of_nonneg_left (hle j) (hunit j).1)
          hsq (hsum.mul_right _)
    _ = excited lam n ^ 2 := by rw [tsum_mul_right, sq]; rfl

/-- Eq. `ym:gap-from-purity`, first part: `e^{-n𝗆} = r₂ ≤ s_n ≤ 𝖯(n)⁻¹ - 1`,
with `𝗆 = log(λ₁/λ₂)`; it needs `λ₂ > 0` (otherwise the gap is infinite). -/
theorem exp_neg_gap_le {lam : ℕ → ℝ} (h : Spectrum lam) (h₂ : 0 < lam 1) {n : ℕ} (hn : 1 ≤ n) :
    exp (-(n * log (lam 0 / lam 1))) ≤ 1 / purity lam n - 1 := by
  have hpos : 0 < lam 1 / lam 0 := div_pos h₂ h.top_pos
  have hexp : exp (-(n * log (lam 0 / lam 1))) = ratio lam n 0 := by
    unfold ratio
    rw [← Real.rpow_natCast, Real.rpow_def_of_pos hpos, Real.log_div h₂.ne' h.top_pos.ne',
      Real.log_div h.top_pos.ne' h₂.ne']
    congr 1
    ring
  rw [hexp]
  have hunit := ratio_mem_unit h n
  calc ratio lam n 0 ≤ excited lam n := (summable_ratio h hn).le_tsum 0 (fun i _ => (hunit i).1)
    _ ≤ 1 / purity lam n - 1 := (purity_bounds h hn).2.2.1

/-- Eq. `ym:gap-from-purity`, second part. The source states it when `𝖯(n) > 1/2`, which
makes the bound positive; the inequality itself holds whenever `λ₂ > 0`. -/
theorem gap_from_purity {lam : ℕ → ℝ} (h : Spectrum lam) (h₂ : 0 < lam 1) {n : ℕ} (hn : 1 ≤ n) :
    (1 / n) * log (1 / (1 / purity lam n - 1)) ≤ log (lam 0 / lam 1) := by
  have hle := exp_neg_gap_le h h₂ hn
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  have hx : 0 < 1 / purity lam n - 1 := lt_of_lt_of_le (exp_pos _) hle
  have hlog : -(n * log (lam 0 / lam 1)) ≤ log (1 / purity lam n - 1) :=
    (Real.le_log_iff_exp_le hx).mpr hle
  rw [one_div (1 / purity lam n - 1), Real.log_inv, one_div_mul_eq_div, div_le_iff₀ hn']
  linarith

end YMFormal
