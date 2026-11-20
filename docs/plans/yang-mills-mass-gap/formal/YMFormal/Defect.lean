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
  sorry

/-- Eq. `ym:defect-telescope`: `Δ_β(n) ≥ 0` and `𝖯₄ ≥ e^{-Δ_β(n)}`. -/
theorem defect_nonneg_and_last {z₁ z₂ z₃ z₄ z₅ : ℝ} (h₁ : 0 < z₁) (h₂ : 0 < z₂) (h₃ : 0 < z₃)
    (h₄ : 0 < z₄) (h₅ : 0 < z₅) (hP₁ : z₂ / z₁ ^ 2 ≤ 1) (hP₂ : z₃ / z₂ ^ 2 ≤ 1)
    (hP₃ : z₄ / z₃ ^ 2 ≤ 1) (hP₄ : z₅ / z₄ ^ 2 ≤ 1) :
    0 ≤ defect z₁ z₅ ∧ exp (-defect z₁ z₅) ≤ z₅ / z₄ ^ 2 := by
  sorry

/-- Eq. `ym:gap-from-defect`. `lam` is the spectrum of the transfer operator on the spatial
torus `(2n,2n,2n)` and `purity lam n = 𝖯₄ ≥ e^{-Δ}`. The source assumes `Δ < log 2`, which makes
the bound positive; the inequality itself needs only `Δ > 0` and `λ₂ > 0`. -/
theorem gap_from_defect {lam : ℕ → ℝ} (h : Spectrum lam) (h₂ : 0 < lam 1) {n : ℕ} (hn : 1 ≤ n)
    {Δ : ℝ} (hΔ : 0 < Δ) (hP : exp (-Δ) ≤ purity lam n) :
    (1 / n) * log (1 / (exp Δ - 1)) ≤ log (lam 0 / lam 1) := by
  sorry

end YMFormal
