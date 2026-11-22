import YMFormal.RankOne

/-!
# The rectangle recurrence (Codex, review YM-O03 §1)

Codex's derivation, 2026-11-21: for the four-direction lattice partition function,
`Δ(2n) ≤ 3740 Δ(n)²` whenever `Δ(n) ≤ 1/7480`. This is the form of the hypothesis `recurrence` of
`conditional_gap` with `c₁ = 3740`, so it would discharge Conjecture `ym:conj:recurrence`.

**Setting.** A box has four side lengths. Index the boxes `n · 2^e` by exponent vectors
`e : Fin 4 → ℕ` and write `G e` for `log Z` of that box. Doubling side `i` adds one to `e i`.
- `cost G i e = 2 G(e) - G(Dᵢ e)` is the log impurity of doubling side `i` (`pᵢ` in the review);
- `cubeDefect G c = 16 G(c,c,c,c) - G(c+1,c+1,c+1,c+1)`; at `c = 0` this is `Δ(n)`, at `c = 1` it is `Δ(2n)`.

**Hypotheses, not proved here.** They are where the physics enters:
- `nonneg`: every cost is nonnegative (each ratio is a purity, at most 1);
- `estimate`: Estimate 1 of the review, `pᵢ ≤ 1/4 → pᵢ(Dᵢ e) ≤ 4 pᵢ²`.

Both follow from `purity_bounds` once `Z` of every box is the trace of a power of a positive
transfer operator in each of the four directions (reflection positivity, C23). That identification
is not formalised.

**Proved here.** The rest of the derivation: the commutation identity, the path sums, the induction
on the number of long sides, the constant 3740, and the chain into `conditional_gap_ratio`.
-/

namespace YMFormal

open Real Finset

/-- Exponent vectors: `e` stands for the box with sides `n · 2^(e k)`. -/
abbrev Exps := Fin 4 → ℕ

/-- Doubling side `i`. -/
def dbl (i : Fin 4) (e : Exps) : Exps := e + Pi.single i 1

/-- The cost of doubling side `i`: `pᵢ(e) = 2 G(e) - G(Dᵢ e)`. -/
def cost (G : Exps → ℝ) (i : Fin 4) (e : Exps) : ℝ := 2 * G e - G (dbl i e)

/-- `16 G(c,c,c,c) - G(c+1,…)`: `Δ(n)` at `c = 0`, `Δ(2n)` at `c = 1`. -/
def cubeDefect (G : Exps → ℝ) (c : ℕ) : ℝ := 16 * G (fun _ => c) - G (fun _ => c + 1)

theorem dbl_comm (i j : Fin 4) (e : Exps) : dbl i (dbl j e) = dbl j (dbl i e) := by
  unfold dbl; abel

/-- Estimate 2 of the review: doublings in two directions commute. -/
theorem cost_swap (G : Exps → ℝ) (i j : Fin 4) (e : Exps) :
    cost G i (dbl j e) = 2 * cost G i e - 2 * cost G j e + cost G j (dbl i e) := by
  unfold cost; rw [dbl_comm i j e]; ring

theorem dbl_apply (i : Fin 4) (e : Exps) (k : Fin 4) : dbl i e k = e k + if k = i then 1 else 0 := by
  simp [dbl, Pi.single_apply]

theorem sum_dbl (i : Fin 4) (e : Exps) : ∑ k, dbl i e k = (∑ k, e k) + 1 := by
  simp [dbl, Finset.sum_add_distrib]

/-! ### The base case: boxes with sides in `{n, 2n}` -/

/-- A potential along doublings: `φ(e) = 16 G(e) / 2^(Σ e)`. -/
noncomputable def pot (G : Exps → ℝ) (e : Exps) : ℝ := 16 * G e / 2 ^ (∑ k, e k)

theorem pot_drop (G : Exps → ℝ) (i : Fin 4) (e : Exps) :
    pot G e - pot G (dbl i e) = 8 / 2 ^ (∑ k, e k) * cost G i e := by
  unfold pot cost
  rw [sum_dbl, pow_succ]
  field_simp
  ring

def Binary (e : Exps) : Prop := ∀ k, e k ≤ 1

theorem binary_sum_le {e : Exps} (h : Binary e) : ∑ k, e k ≤ 4 := by
  calc ∑ k, e k ≤ ∑ _k : Fin 4, 1 := Finset.sum_le_sum (fun k _ => h k)
    _ = 4 := by simp

/-- Undoing a doubling: if `e k = 1` then `e = Dₖ (e with side k set to 0)`. -/
theorem dbl_update_zero {e : Exps} {k : Fin 4} (hk : e k = 1) :
    dbl k (Function.update e k 0) = e := by
  funext x
  rw [dbl_apply]
  by_cases hx : x = k
  · subst hx; simp [hk]
  · simp [hx]

theorem pot_down_from_zero (G : Exps → ℝ) (nonneg : ∀ i e, 0 ≤ cost G i e) :
    ∀ m : ℕ, ∀ e : Exps, Binary e → ∑ k, e k = m → pot G e ≤ pot G 0 := by
  intro m
  induction m with
  | zero =>
    intro e _ hs
    have : e = 0 := by
      funext k
      exact (Finset.sum_eq_zero_iff.mp hs) k (Finset.mem_univ k)
    rw [this]
  | succ m ih =>
    intro e he hs
    obtain ⟨k, -, hk⟩ : ∃ k ∈ Finset.univ, e k ≠ 0 := by
      by_contra hcon
      push_neg at hcon
      rw [Finset.sum_eq_zero hcon] at hs
      omega
    have hk1 : e k = 1 := by have := he k; omega
    set e' := Function.update e k 0 with he'
    have hde : dbl k e' = e := dbl_update_zero hk1
    have hb' : Binary e' := by
      intro x; by_cases hx : x = k
      · subst hx; simp [e']
      · simp [e', hx, he x]
    have hs' : ∑ x, e' x = m := by
      have := sum_dbl k e'
      rw [hde, hs] at this
      omega
    have hdrop := pot_drop G k e'
    rw [hde] at hdrop
    have hpos : 0 ≤ 8 / 2 ^ (∑ x, e' x) * cost G k e' := by
      have := nonneg k e'; positivity
    have := ih e' hb' hs'
    linarith

theorem pot_up_to_one (G : Exps → ℝ) (nonneg : ∀ i e, 0 ≤ cost G i e) :
    ∀ m : ℕ, ∀ e : Exps, Binary e → (∑ k, e k) + m = 4 → pot G (fun _ => 1) ≤ pot G e := by
  intro m
  induction m with
  | zero =>
    intro e he hs
    have : e = fun _ => 1 := by
      funext k
      by_contra hne
      have hk0 : e k = 0 := by have := he k; omega
      have hlt : ∑ x, e x < ∑ _x : Fin 4, 1 := by
        apply Finset.sum_lt_sum (fun x _ => he x)
        exact ⟨k, Finset.mem_univ k, by omega⟩
      simp at hlt
      omega
    rw [this]
  | succ m ih =>
    intro e he hs
    obtain ⟨k, hk⟩ : ∃ k, e k = 0 := by
      by_contra hcon
      push_neg at hcon
      have : ∀ x, e x = 1 := fun x => by have := he x; have := hcon x; omega
      have : ∑ x, e x = 4 := by simp [this]
      omega
    have hb' : Binary (dbl k e) := by
      intro x; rw [dbl_apply]; by_cases hx : x = k
      · subst hx; simp [hk]
      · simp [hx, he x]
    have hs' : (∑ x, dbl k e x) + m = 4 := by rw [sum_dbl]; omega
    have hdrop := pot_drop G k e
    have hpos : 0 ≤ 8 / 2 ^ (∑ x, e x) * cost G k e := by
      have := nonneg k e; positivity
    have := ih (dbl k e) hb' hs'
    linarith

theorem pot_zero (G : Exps → ℝ) : pot G 0 = 16 * G 0 := by simp [pot]

theorem pot_one (G : Exps → ℝ) : pot G (fun _ => 1) = G (fun _ => 1) := by
  simp [pot]; norm_num

theorem cubeDefect_zero_eq (G : Exps → ℝ) :
    cubeDefect G 0 = pot G 0 - pot G (fun _ => 1) := by
  rw [pot_zero, pot_one]; rfl

theorem cubeDefect_nonneg (G : Exps → ℝ) (nonneg : ∀ i e, 0 ≤ cost G i e) :
    0 ≤ cubeDefect G 0 := by
  rw [cubeDefect_zero_eq]
  have := pot_up_to_one G nonneg 4 0 (fun _ => by simp) (by simp)
  linarith

/-- Every doubling inside the first cube costs at most `Δ(n)`. -/
theorem cost_le_defect_binary (G : Exps → ℝ) (nonneg : ∀ i e, 0 ≤ cost G i e)
    {i : Fin 4} {e : Exps} (he : Binary e) (hi : e i = 0) : cost G i e ≤ cubeDefect G 0 := by
  rw [cubeDefect_zero_eq]
  have hb' : Binary (dbl i e) := by
    intro x; rw [dbl_apply]; by_cases hx : x = i
    · subst hx; simp [hi]
    · simp [hx, he x]
  have hs : ∑ k, e k ≤ 3 := by have := binary_sum_le hb'; rw [sum_dbl] at this; omega
  have h1 := pot_down_from_zero G nonneg _ e he rfl
  have h2 := pot_up_to_one G nonneg (4 - ∑ k, dbl i e k) (dbl i e) hb'
    (by have := binary_sum_le hb'; omega)
  have hdrop := pot_drop G i e
  have hw : 1 ≤ 8 / (2 : ℝ) ^ (∑ k, e k) := by
    rw [le_div_iff₀ (by positivity)]
    have : (2 : ℝ) ^ (∑ k, e k) ≤ 2 ^ 3 := pow_le_pow_right₀ (by norm_num) hs
    norm_num at this ⊢; linarith
  have hc := nonneg i e
  nlinarith

/-! ### The induction on the number of sides of length `4n` -/

/-- The number of sides of length `4n` (exponent 2). -/
def twos (e : Exps) : ℕ := (Finset.univ.filter (fun k => e k = 2)).card

def Ternary (e : Exps) : Prop := ∀ k, e k ≤ 2

theorem twos_update_one {e : Exps} {h : Fin 4} (hh : e h = 2) :
    twos (Function.update e h 1) + 1 = twos e := by
  unfold twos
  have : Finset.univ.filter (fun k => Function.update e h 1 k = 2) =
      (Finset.univ.filter (fun k => e k = 2)).erase h := by
    ext k
    by_cases hk : k = h
    · subst hk; simp
    · simp [hk]
  rw [this, Finset.card_erase_of_mem (by simp [hh])]
  have : 0 < (Finset.univ.filter (fun k => e k = 2)).card :=
    Finset.card_pos.mpr ⟨h, by simp [hh]⟩
  omega

/-- Codex's claim: `pᵢ(x) ≤ 3^j Δ(n)` when side `i` has length `n`, the others have lengths
`n`, `2n` or `4n`, and `j` of them have length `4n` (`j ≤ 3`, `9 Δ(n) ≤ 1/4`). -/
theorem cost_le_pow_three (G : Exps → ℝ) (nonneg : ∀ i e, 0 ≤ cost G i e)
    (estimate : ∀ i e, cost G i e ≤ 1 / 4 → cost G i (dbl i e) ≤ 4 * cost G i e ^ 2)
    (small : 9 * cubeDefect G 0 ≤ 1 / 4) :
    ∀ j ≤ 3, ∀ i e, e i = 0 → Ternary e → twos e = j → cost G i e ≤ 3 ^ j * cubeDefect G 0 := by
  have hδ := cubeDefect_nonneg G nonneg
  intro j
  induction j with
  | zero =>
    intro _ i e hi he ht
    have hb : Binary e := by
      intro k
      have hk2 : e k ≠ 2 := by
        intro h2
        have : k ∈ Finset.univ.filter (fun k => e k = 2) := by simp [h2]
        unfold twos at ht
        rw [Finset.card_eq_zero] at ht
        simp [ht] at this
      have := he k; omega
    simpa using cost_le_defect_binary G nonneg hb hi
  | succ j ih =>
    intro hj i e hi he ht
    obtain ⟨h, hh⟩ : ∃ h, e h = 2 := by
      have : 0 < twos e := by omega
      obtain ⟨h, hmem⟩ := Finset.card_pos.mp this
      exact ⟨h, by simpa using hmem⟩
    have hhi : h ≠ i := by rintro rfl; omega
    -- `e = D_h y` with `y h = 1`
    set y := Function.update e h 1 with hy
    have hye : dbl h y = e := by
      funext x; rw [dbl_apply]; by_cases hx : x = h
      · subst hx; simp [y, hh]
      · simp [y, hx]
    have hyi : y i = 0 := by simp [y, Ne.symm hhi, hi]
    have hyt : Ternary y := by
      intro x; by_cases hx : x = h
      · subst hx; simp [y]
      · simp [y, hx, he x]
    have hytw : twos y = j := by have := twos_update_one hh; rw [← hy] at this; omega
    -- `D_i y = D_h z` with `z h = 0`
    set z := Function.update (dbl i y) h 0 with hz
    have hdy : dbl i y h = 1 := by rw [dbl_apply]; simp [y, hhi]
    have hze : dbl h z = dbl i y := dbl_update_zero hdy
    have hzh : z h = 0 := by simp [z]
    have hzt : Ternary z := by
      intro x; by_cases hx : x = h
      · subst hx; simp [z]
      · simp only [z, Function.update_of_ne hx, dbl_apply]
        by_cases hxi : x = i
        · subst hxi; simp [hyi]
        · simp [hxi, hyt x]
    have hztw : twos z = j := by
      rw [← hytw]
      unfold twos
      congr 1
      ext x
      simp only [Finset.mem_filter, Finset.mem_univ, true_and]
      by_cases hx : x = h
      · subst hx; simp [z, y]
      · simp only [z, y, Function.update_of_ne hx, dbl_apply]
        by_cases hxi : x = i
        · subst hxi; simp [hi]
        · simp [hxi]
    have hA : 3 ^ j * cubeDefect G 0 ≤ 1 / 4 := by
      have : (3 : ℝ) ^ j ≤ 9 := by
        calc (3 : ℝ) ^ j ≤ 3 ^ 2 := pow_le_pow_right₀ (by norm_num) (by omega)
          _ = 9 := by norm_num
      nlinarith
    have h1 := ih (by omega) i y hyi hyt hytw
    have h2 := ih (by omega) h z hzh hzt hztw
    have hz0 := nonneg h z
    have h3 : cost G h (dbl i y) ≤ 4 * cost G h z ^ 2 := by
      rw [← hze]; exact estimate h z (h2.trans hA)
    have hsq : cost G h z ^ 2 ≤ (3 ^ j * cubeDefect G 0) ^ 2 := pow_le_pow_left₀ hz0 h2 2
    have hy0 := nonneg h y
    have hswap := cost_swap G i h y
    rw [hye] at hswap
    have hAnn : 0 ≤ 3 ^ j * cubeDefect G 0 := by positivity
    rw [pow_succ]
    nlinarith

/-! ### The recurrence -/

/-- **Rectangle recurrence.** `Δ(2n) ≤ 3740 Δ(n)²` whenever `Δ(n) ≤ 1/7480`. -/
theorem rectangle_recurrence (G : Exps → ℝ) (nonneg : ∀ i e, 0 ≤ cost G i e)
    (estimate : ∀ i e, cost G i e ≤ 1 / 4 → cost G i (dbl i e) ≤ 4 * cost G i e ^ 2)
    (small : cubeDefect G 0 ≤ 1 / 7480) :
    cubeDefect G 1 ≤ 3740 * cubeDefect G 0 ^ 2 := by
  have hδ := cubeDefect_nonneg G nonneg
  have h9 : 9 * cubeDefect G 0 ≤ 1 / 4 := by linarith
  have key := cost_le_pow_three G nonneg estimate h9
  set δ := cubeDefect G 0
  -- the path (1,1,1,1) → (2,1,1,1) → (2,2,1,1) → (2,2,2,1) → (2,2,2,2)
  let v0 : Exps := fun _ => 1
  let w0 : Exps := ![0, 1, 1, 1]
  let w1 : Exps := ![2, 0, 1, 1]
  let w2 : Exps := ![2, 2, 0, 1]
  let w3 : Exps := ![2, 2, 2, 0]
  have e0 : dbl 0 w0 = v0 := by decide
  have e1 : dbl 1 w1 = dbl 0 v0 := by decide
  have e2 : dbl 2 w2 = dbl 1 (dbl 0 v0) := by decide
  have e3 : dbl 3 w3 = dbl 2 (dbl 1 (dbl 0 v0)) := by decide
  have e4 : dbl 3 (dbl 2 (dbl 1 (dbl 0 v0))) = fun _ => 1 + 1 := by decide
  have hpath : cubeDefect G 1 = 8 * cost G 0 v0 + 4 * cost G 1 (dbl 0 v0) +
      2 * cost G 2 (dbl 1 (dbl 0 v0)) + cost G 3 (dbl 2 (dbl 1 (dbl 0 v0))) := by
    unfold cubeDefect cost
    rw [e4]
    ring
  -- each step is the doubling of a side of length 2n, bounded by Estimate 1
  have step : ∀ (k : Fin 4) (w : Exps) (j : ℕ), j ≤ 3 → w k = 0 → Ternary w → twos w = j →
      cost G k (dbl k w) ≤ 4 * (3 ^ j * δ) ^ 2 := by
    intro k w j hj hk hw ht
    have hc := key j hj k w hk hw ht
    have hA : 3 ^ j * δ ≤ 1 / 4 := by
      have : (3 : ℝ) ^ j ≤ 27 := by
        calc (3 : ℝ) ^ j ≤ 3 ^ 3 := pow_le_pow_right₀ (by norm_num) hj
          _ = 27 := by norm_num
      nlinarith
    have hsq := pow_le_pow_left₀ (nonneg k w) hc 2
    have := estimate k w (hc.trans hA)
    nlinarith
  have s0 := step 0 w0 0 (by norm_num) (by decide) (by unfold Ternary; decide) (by decide)
  have s1 := step 1 w1 1 (by norm_num) (by decide) (by unfold Ternary; decide) (by decide)
  have s2 := step 2 w2 2 (by norm_num) (by decide) (by unfold Ternary; decide) (by decide)
  have s3 := step 3 w3 3 (by norm_num) (by decide) (by unfold Ternary; decide) (by decide)
  rw [e0] at s0; rw [e1] at s1; rw [e2] at s2; rw [e3] at s3
  rw [hpath]
  nlinarith

/-! ### In terms of side lengths, and into the conditional gap -/

/-- Doubling side `i` of a box with side lengths `x`. -/
def sideDouble (i : Fin 4) (x : Fin 4 → ℕ) : Fin 4 → ℕ := Function.update x i (2 * x i)

/-- `pᵢ(x) = 2 log Z(x) - log Z(Dᵢ x)` for `F = log Z`. -/
def sideCost (F : (Fin 4 → ℕ) → ℝ) (i : Fin 4) (x : Fin 4 → ℕ) : ℝ :=
  2 * F x - F (sideDouble i x)

/-- `Δ(n) = 16 log Z(n,n,n,n) - log Z(2n,2n,2n,2n)`. -/
def boxDefect (F : (Fin 4 → ℕ) → ℝ) (n : ℕ) : ℝ := 16 * F (fun _ => n) - F (fun _ => 2 * n)

/-- The exponent picture of `F` at base side `n`. -/
def atScale (F : (Fin 4 → ℕ) → ℝ) (n : ℕ) : Exps → ℝ := fun e => F (fun k => n * 2 ^ e k)

theorem atScale_dbl (n : ℕ) (i : Fin 4) (e : Exps) :
    (fun k => n * 2 ^ (dbl i e k)) = sideDouble i (fun k => n * 2 ^ e k) := by
  funext k
  rw [dbl_apply, sideDouble]
  by_cases hk : k = i
  · subst hk; simp [pow_succ]; ring
  · simp [hk]

theorem cost_atScale (F : (Fin 4 → ℕ) → ℝ) (n : ℕ) (i : Fin 4) (e : Exps) :
    cost (atScale F n) i e = sideCost F i (fun k => n * 2 ^ e k) := by
  unfold cost sideCost atScale; rw [atScale_dbl]

theorem cubeDefect_atScale (F : (Fin 4 → ℕ) → ℝ) (n : ℕ) (c : ℕ) :
    cubeDefect (atScale F n) c = boxDefect F (n * 2 ^ c) := by
  unfold cubeDefect boxDefect atScale
  congr 2
  funext k; rw [pow_succ]; ring

/-- The rectangle recurrence for side lengths. The hypotheses are required only for boxes with
positive sides. -/
theorem box_recurrence (F : (Fin 4 → ℕ) → ℝ) {n : ℕ} (hn : 1 ≤ n)
    (nonneg : ∀ i x, (∀ k, 0 < x k) → 0 ≤ sideCost F i x)
    (estimate : ∀ i x, (∀ k, 0 < x k) →
      sideCost F i x ≤ 1 / 4 → sideCost F i (sideDouble i x) ≤ 4 * sideCost F i x ^ 2)
    (small : boxDefect F n ≤ 1 / 7480) :
    boxDefect F (2 * n) ≤ 3740 * boxDefect F n ^ 2 := by
  have hpos : ∀ e : Exps, ∀ k, 0 < n * 2 ^ e k := fun e k => by positivity
  have h := rectangle_recurrence (atScale F n)
    (fun i e => by rw [cost_atScale]; exact nonneg i _ (hpos e))
    (fun i e he => by
      simp only [cost_atScale] at he ⊢
      rw [atScale_dbl]
      exact estimate i _ (hpos e) he)
    (by rw [cubeDefect_atScale]; simpa using small)
  rw [cubeDefect_atScale, cubeDefect_atScale] at h
  simpa [mul_comm] using h

theorem boxDefect_nonneg (F : (Fin 4 → ℕ) → ℝ) {n : ℕ} (hn : 1 ≤ n)
    (nonneg : ∀ i x, (∀ k, 0 < x k) → 0 ≤ sideCost F i x) : 0 ≤ boxDefect F n := by
  have hpos : ∀ e : Exps, ∀ k, 0 < n * 2 ^ e k := fun e k => by positivity
  have := cubeDefect_nonneg (atScale F n) (fun i e => by rw [cost_atScale]; exact nonneg i _ (hpos e))
  rw [cubeDefect_atScale] at this
  simpa using this

/-- **The conditional gap with the recurrence discharged.** With the rectangle recurrence the only
remaining inputs are the spectral hypotheses (`nonneg`, `estimate`, `from_defect`) and the seed
`Δ(n) ≤ 1/7480`. `r k` is `λ₂/λ₁` at side `2ᵏ n`; `(r k)ⁿ ≤ 1/2` is `𝗆 ≥ log 2 / n`. -/
theorem conditional_gap_of_boxes (F : (Fin 4 → ℕ) → ℝ) {n : ℕ} (hn : 1 ≤ n) {r : ℕ → ℝ}
    (hr : ∀ k, 0 ≤ r k)
    (nonneg : ∀ i x, (∀ k, 0 < x k) → 0 ≤ sideCost F i x)
    (estimate : ∀ i x, (∀ k, 0 < x k) →
      sideCost F i x ≤ 1 / 4 → sideCost F i (sideDouble i x) ≤ 4 * sideCost F i x ^ 2)
    (seed : boxDefect F n ≤ 1 / 7480)
    (from_defect : ∀ k, r k ^ (2 ^ k * n) ≤ exp (boxDefect F (2 ^ k * n)) - 1) :
    ∀ k, r k ^ n ≤ 1 / 2 := by
  have hk1 : ∀ k, 1 ≤ 2 ^ k * n := fun k => by
    have : 1 ≤ 2 ^ k := Nat.one_le_two_pow
    nlinarith
  refine conditional_gap_ratio (c₁ := 3740) (by norm_num) (δ := fun k => boxDefect F (2 ^ k * n))
    (fun k => boxDefect_nonneg F (hk1 k) nonneg) hr ?_
    (by show boxDefect F (2 ^ 0 * n) ≤ 1 / (2 * 3740); rw [pow_zero, one_mul]; linarith) from_defect
  intro k hk
  have := box_recurrence F (hk1 k) nonneg estimate (by norm_num at hk ⊢; linarith)
  simpa [pow_succ, mul_assoc, mul_comm, mul_left_comm] using this

end YMFormal
