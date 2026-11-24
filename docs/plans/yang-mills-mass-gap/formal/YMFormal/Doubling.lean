import YMFormal.Rectangle

/-!
# Doubling costs from a transfer spectrum (handbook: `ym:cost-nonneg`, Lemma `ym:lem:one-doubling`)

The rectangle recurrence takes two hypotheses on the doubling costs: they are nonnegative
(`nonneg`) and they satisfy Lemma `ym:lem:one-doubling` (`estimate`). Here both are proved for
`Z(t) = ∑_j λ_jᵗ` with `Spectrum lam`, so that once a box's `Z` is the trace of a power of a
nonnegative transfer operator they are no longer assumptions.

With `p(t) = 2 log Z(t) - log Z(2t) = -log 𝖯(t)`:
- `doubling_cost_nonneg`: `0 ≤ p(t)`;
- `doubling_cost_estimate`: `p(t) ≤ 1/4 → p(2t) ≤ 4 p(t)²`.
-/

namespace YMFormal

open Real

/-- `Z(t) > 0` for `t ≥ 1`. -/
theorem partition_pos {lam : ℕ → ℝ} (h : Spectrum lam) {t : ℕ} (ht : 1 ≤ t) :
    0 < partition lam t := by
  rw [partition_eq h ht]
  have hs : 0 ≤ excited lam t :=
    tsum_nonneg fun j => (ratio_mem_unit h t j).1
  have := pow_pos h.top_pos t
  positivity

theorem purity_pos {lam : ℕ → ℝ} (h : Spectrum lam) {t : ℕ} (ht : 1 ≤ t) : 0 < purity lam t := by
  unfold purity
  have h1 := partition_pos h ht
  have h2 := partition_pos h (by omega : 1 ≤ 2 * t)
  positivity

/-- The doubling cost `p(t) = 2 log Z(t) - log Z(2t)`. -/
noncomputable def doublingCost (lam : ℕ → ℝ) (t : ℕ) : ℝ :=
  2 * log (partition lam t) - log (partition lam (2 * t))

theorem doublingCost_eq {lam : ℕ → ℝ} (h : Spectrum lam) {t : ℕ} (ht : 1 ≤ t) :
    doublingCost lam t = -log (purity lam t) := by
  unfold doublingCost purity
  have h1 := partition_pos h ht
  have h2 := partition_pos h (by omega : 1 ≤ 2 * t)
  rw [log_div h2.ne' (by positivity), log_pow]
  push_cast
  ring

/-- `ym:cost-nonneg`. -/
theorem doubling_cost_nonneg {lam : ℕ → ℝ} (h : Spectrum lam) {t : ℕ} (ht : 1 ≤ t) :
    0 ≤ doublingCost lam t := by
  rw [doublingCost_eq h ht]
  have hP := (purity_bounds h ht).2.1
  have := log_nonpos (purity_pos h ht).le hP
  linarith

/-- `exp (1/2) < 2`. -/
theorem exp_half_lt_two : exp (1 / 2) < 2 := by
  have h := exp_one_lt_d9
  have hsq : exp (1 / 2) ^ 2 = exp 1 := by rw [← exp_nat_mul]; norm_num
  nlinarith [exp_pos (1 / 2)]

/-- Lemma `ym:lem:one-doubling`: `p(t) ≤ 1/4 → p(2t) ≤ 4 p(t)²`. -/
theorem doubling_cost_estimate {lam : ℕ → ℝ} (h : Spectrum lam) {t : ℕ} (ht : 1 ≤ t)
    (hp : doublingCost lam t ≤ 1 / 4) :
    doublingCost lam (2 * t) ≤ 4 * doublingCost lam t ^ 2 := by
  have ht2 : 1 ≤ 2 * t := by omega
  set p := doublingCost lam t with hpdef
  have hp0 : 0 ≤ p := doubling_cost_nonneg h ht
  have hb := purity_bounds h ht
  have hb2 := purity_bounds h ht2
  set s := excited lam t
  set s2 := excited lam (2 * t)
  have hs0 : 0 ≤ s := tsum_nonneg fun j => (ratio_mem_unit h t j).1
  have hs20 : 0 ≤ s2 := tsum_nonneg fun j => (ratio_mem_unit h (2 * t) j).1
  -- `s ≤ 1/𝖯 - 1 = e^p - 1 ≤ p e^p`
  have hPp : 1 / purity lam t = exp p := by
    rw [hpdef, doublingCost_eq h ht, exp_neg, exp_log (purity_pos h ht), one_div]
  have hs : s ≤ p * exp p := by
    have h1 : s ≤ exp p - 1 := by rw [← hPp]; exact hb.2.2.1
    have h2 : exp p - 1 ≤ p * exp p := by
      have := add_one_le_exp (-p)
      have he := exp_pos p
      have : exp (-p) * exp p = 1 := by rw [← exp_add]; simp
      nlinarith
    linarith
  -- `s₂ₜ ≤ s²`
  have hs2 : s2 ≤ s ^ 2 := hb.2.2.2.2
  -- `p(2t) = log (1/𝖯(2t)) ≤ log ((1 + s₂ₜ)²) ≤ 2 s₂ₜ`
  have hq : doublingCost lam (2 * t) ≤ 2 * s2 := by
    rw [doublingCost_eq h ht2]
    have hP2 := purity_pos h ht2
    have hinv : 1 / purity lam (2 * t) ≤ (1 + s2) ^ 2 := by have := hb2.2.2.2.1; nlinarith
    have hlog : log (1 / purity lam (2 * t)) ≤ log ((1 + s2) ^ 2) :=
      log_le_log (by positivity) hinv
    rw [one_div, log_inv, log_pow] at hlog
    have hl1 : log (1 + s2) ≤ s2 := by
      have := add_one_le_exp s2
      have := log_le_log (by positivity) this
      rwa [log_exp, add_comm] at this
    push_cast at hlog
    linarith
  -- `e^{2p} ≤ e^{1/2} < 2`
  have he : exp p ^ 2 ≤ 2 := by
    have : exp p ^ 2 = exp (2 * p) := by rw [← exp_nat_mul]; norm_num
    rw [this]
    exact (exp_le_exp.mpr (by linarith)).trans exp_half_lt_two.le
  have hsq : s ^ 2 ≤ (p * exp p) ^ 2 := pow_le_pow_left₀ hs0 hs 2
  calc doublingCost lam (2 * t) ≤ 2 * s2 := hq
    _ ≤ 2 * s ^ 2 := by linarith
    _ ≤ 2 * (p * exp p) ^ 2 := by linarith
    _ = 2 * p ^ 2 * exp p ^ 2 := by ring
    _ ≤ 2 * p ^ 2 * 2 := by have := sq_nonneg p; nlinarith
    _ = 4 * p ^ 2 := by ring

end YMFormal
