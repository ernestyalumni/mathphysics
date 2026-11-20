import YMFormal.Purity

/-!
# Defect and gap on the doubled torus (`ym:prop:defect`)

Source: `documents/mono-notes/topics/yang-mills-lattice.tex`, Definition `ym:def:defect`,
Proposition `ym:prop:defect`.

`z₁, …, z₅` stand for `Z(n,n,n,n)`, `Z(2n,n,n,n)`, `Z(2n,2n,n,n)`, `Z(2n,2n,2n,n)` and
`Z(2n,2n,2n,2n)`. That each ratio `𝖯_k` is the purity of a transfer operator is the
first sentence of the source's proof. It rests on the symmetry of `Z` and the trace formula
`ym:trace-formula`, and enters here as hypotheses: `0 < 𝖯_k ≤ 1`, and `𝖯₄ = purity lam n`
for the spectrum `lam` on the spatial torus `(2n,2n,2n)`.
-/

open Real

namespace YMFormal

/-- `Δ_β(n) = 16 log Z(n,n,n,n) - log Z(2n,2n,2n,2n)`, Eq. `ym:defect`. -/
noncomputable def defect (z₁ z₅ : ℝ) : ℝ := 16 * log z₁ - log z₅

/-- Eq. `ym:defect-telescope`, the product of the four purities. -/
theorem defect_telescope {z₁ z₂ z₃ z₄ z₅ : ℝ} (h₁ : 0 < z₁) (h₂ : 0 < z₂) (h₃ : 0 < z₃)
    (h₄ : 0 < z₄) (h₅ : 0 < z₅) :
    exp (-defect z₁ z₅) =
      (z₂ / z₁ ^ 2) ^ 8 * (z₃ / z₂ ^ 2) ^ 4 * (z₄ / z₃ ^ 2) ^ 2 * (z₅ / z₄ ^ 2) := by
  have hexp : exp (-defect z₁ z₅) = z₅ / z₁ ^ 16 := by
    unfold defect
    rw [neg_sub, exp_sub, exp_log h₅, show (16 : ℝ) * log z₁ = log (z₁ ^ 16) by
      rw [Real.log_pow]; norm_num, exp_log (by positivity)]
  rw [hexp]
  field_simp

/-- Eq. `ym:defect-telescope`: `Δ_β(n) ≥ 0` and `𝖯₄ ≥ e^{-Δ_β(n)}`. -/
theorem defect_nonneg_and_last {z₁ z₂ z₃ z₄ z₅ : ℝ} (h₁ : 0 < z₁) (h₂ : 0 < z₂) (h₃ : 0 < z₃)
    (h₄ : 0 < z₄) (h₅ : 0 < z₅) (hP₁ : z₂ / z₁ ^ 2 ≤ 1) (hP₂ : z₃ / z₂ ^ 2 ≤ 1)
    (hP₃ : z₄ / z₃ ^ 2 ≤ 1) (hP₄ : z₅ / z₄ ^ 2 ≤ 1) :
    0 ≤ defect z₁ z₅ ∧ exp (-defect z₁ z₅) ≤ z₅ / z₄ ^ 2 := by
  have htel := defect_telescope h₁ h₂ h₃ h₄ h₅
  have p₁ : 0 < z₂ / z₁ ^ 2 := by positivity
  have p₂ : 0 < z₃ / z₂ ^ 2 := by positivity
  have p₃ : 0 < z₄ / z₃ ^ 2 := by positivity
  have p₄ : 0 < z₅ / z₄ ^ 2 := by positivity
  have hfront : (z₂ / z₁ ^ 2) ^ 8 * (z₃ / z₂ ^ 2) ^ 4 * (z₄ / z₃ ^ 2) ^ 2 ≤ 1 := by
    have a := pow_le_one₀ p₁.le hP₁ (n := 8)
    have b := pow_le_one₀ p₂.le hP₂ (n := 4)
    have c := pow_le_one₀ p₃.le hP₃ (n := 2)
    have ab : (z₂ / z₁ ^ 2) ^ 8 * (z₃ / z₂ ^ 2) ^ 4 ≤ 1 := by
      simpa using mul_le_mul a b (by positivity) zero_le_one
    simpa using mul_le_mul ab c (by positivity) zero_le_one
  have hlast : exp (-defect z₁ z₅) ≤ z₅ / z₄ ^ 2 := by
    rw [htel]
    calc (z₂ / z₁ ^ 2) ^ 8 * (z₃ / z₂ ^ 2) ^ 4 * (z₄ / z₃ ^ 2) ^ 2 * (z₅ / z₄ ^ 2)
        ≤ 1 * (z₅ / z₄ ^ 2) := mul_le_mul_of_nonneg_right hfront p₄.le
      _ = z₅ / z₄ ^ 2 := one_mul _
  refine ⟨?_, hlast⟩
  have : exp (-defect z₁ z₅) ≤ 1 := hlast.trans hP₄
  have := Real.exp_le_one_iff.mp this
  linarith

/-- Eq. `ym:gap-from-defect`. `lam` is the spectrum of the transfer operator on the spatial
torus `(2n,2n,2n)` and `purity lam n = 𝖯₄ ≥ e^{-Δ}`. The source assumes `Δ < log 2`, which makes
the bound positive; the inequality itself needs only `Δ > 0` and `λ₂ > 0`. -/
theorem gap_from_defect {lam : ℕ → ℝ} (h : Spectrum lam) (h₂ : 0 < lam 1) {n : ℕ} (hn : 1 ≤ n)
    {Δ : ℝ} (hΔ : 0 < Δ) (hP : exp (-Δ) ≤ purity lam n) :
    (1 / n) * log (1 / (exp Δ - 1)) ≤ log (lam 0 / lam 1) := by
  have hle := exp_neg_gap_le h h₂ hn
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  have hPpos : 0 < purity lam n := lt_of_lt_of_le (exp_pos _) hP
  have hinv : 1 / purity lam n ≤ exp Δ := by
    rw [div_le_iff₀ hPpos]
    calc (1 : ℝ) = exp Δ * exp (-Δ) := by rw [← exp_add]; simp
      _ ≤ exp Δ * purity lam n := mul_le_mul_of_nonneg_left hP (exp_pos _).le
  have hx : 0 < exp Δ - 1 := by linarith [Real.add_one_lt_exp hΔ.ne']
  have hle' : exp (-(n * log (lam 0 / lam 1))) ≤ exp Δ - 1 := by linarith
  have hlog : -(n * log (lam 0 / lam 1)) ≤ log (exp Δ - 1) := (Real.le_log_iff_exp_le hx).mpr hle'
  rw [one_div (exp Δ - 1), Real.log_inv, one_div_mul_eq_div, div_le_iff₀ hn']
  linarith

end YMFormal
