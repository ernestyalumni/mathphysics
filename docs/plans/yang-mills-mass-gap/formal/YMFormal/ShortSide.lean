import YMFormal.Z2

/-!
# The short-side recurrence (YM-F03)

The rectangle recurrence (`Rectangle.lean`) is driven by the whole defect `Δ(n)`. Here the same
commutation identity is run inside the class of boxes with exactly one side of length `n` and all
other sides `2n` or `4n`. That class is closed under the induction, so the seed is the single
doubling cost

  `κᵢ(n) = 2 log Z(n; 2n, 2n, 2n) - log Z(2n; 2n, 2n, 2n)`   (short side in direction `i`),

which is at most `Δ(n)`. With the sharper one-doubling estimate `p(2t) ≤ (29/50) p(t)²` for
`p(t) ≤ 1/9`, this gives `κ(2n) ≤ 40 κ(n)²` for `κ(n) ≤ 1/80`.
-/

namespace YMFormal

open Real Finset

/-! ### The sharper one-doubling estimate -/

/-- `s ≤ (y/2)(1 + y²)` with `y = 1/𝖯 - 1 < 1`: from `1/𝖯 = (1+s)²/(1+q)` and `q ≤ s²`. -/
theorem excited_le_half {s q y : ℝ} (hs : 0 ≤ s) (hq : 0 ≤ q) (hqs : q ≤ s ^ 2)
    (hy : y = (1 + s) ^ 2 / (1 + q) - 1) (hsy : s ≤ y) :
    s ≤ y / 2 * (1 + y ^ 2) := by
  have h1 : (1 + s) ^ 2 / (1 + s ^ 2) ≤ (1 + s) ^ 2 / (1 + q) :=
    div_le_div_of_nonneg_left (by positivity) (by positivity) (by linarith)
  have h2 : 2 * s / (1 + s ^ 2) ≤ y := by
    have : (1 + s) ^ 2 / (1 + s ^ 2) = 1 + 2 * s / (1 + s ^ 2) := by field_simp; ring
    rw [hy]; linarith
  rw [div_le_iff₀ (by positivity)] at h2
  have hy0 : 0 ≤ y := hs.trans hsy
  have h3 : s ^ 2 ≤ y ^ 2 := pow_le_pow_left₀ hs hsy 2
  nlinarith

/-- `exp p - 1 ≤ (10577/10000) p` for `0 ≤ p ≤ 1/9`. -/
theorem exp_sub_one_le_ninth {p : ℝ} (h0 : 0 ≤ p) (h1 : p ≤ 1 / 9) :
    exp p - 1 ≤ 10578 / 10000 * p := by
  have h := Real.exp_bound' h0 (by linarith : p ≤ 1) (n := 4) (by norm_num)
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial] at h
  norm_num at h
  have p2 : p ^ 2 ≤ 1 / 9 * p := by nlinarith
  have p3 : p ^ 3 ≤ 1 / 81 * p := by nlinarith
  have p4 : p ^ 4 ≤ 1 / 729 * p := by nlinarith [mul_nonneg h0 (sq_nonneg p)]
  nlinarith

/-- `p(t) ≤ 1/9 → p(2t) ≤ (29/50) p(t)²`. -/
theorem doubling_cost_sharp {lam : ℕ → ℝ} (h : Spectrum lam) {t : ℕ} (ht : 1 ≤ t)
    (hp : doublingCost lam t ≤ 1 / 9) :
    doublingCost lam (2 * t) ≤ 29 / 50 * doublingCost lam t ^ 2 := by
  have ht2 : 1 ≤ 2 * t := by omega
  set p := doublingCost lam t with hpdef
  have hp0 : 0 ≤ p := doubling_cost_nonneg h ht
  have hb := purity_bounds h ht
  have hb2 := purity_bounds h ht2
  set s := excited lam t
  set s2 := excited lam (2 * t)
  have hs0 : 0 ≤ s := tsum_nonneg fun j => (ratio_mem_unit h t j).1
  have hs20 : 0 ≤ s2 := tsum_nonneg fun j => (ratio_mem_unit h (2 * t) j).1
  have hPp : 1 / purity lam t = exp p := by
    rw [hpdef, doublingCost_eq h ht, exp_neg, exp_log (purity_pos h ht), one_div]
  set y := exp p - 1 with hydef
  have hsy : s ≤ y := by rw [hydef, ← hPp]; exact hb.2.2.1
  have hsq : s2 ≤ s ^ 2 := hb.2.2.2.2
  have hq : excitedSquares lam t = s2 := (excited_double lam t).symm
  have hyform : y = (1 + s) ^ 2 / (1 + s2) - 1 := by
    rw [hydef, ← hPp, hb.1, hq, one_div_div]
  have hhalf : s ≤ y / 2 * (1 + y ^ 2) := excited_le_half hs0 hs20 hsq hyform hsy
  have hy0 : 0 ≤ y := hs0.trans hsy
  have hyp : y ≤ 10578 / 10000 * p := exp_sub_one_le_ninth hp0 hp
  have hy1 : y ≤ 1176 / 10000 := by linarith
  have hy2 : y ^ 2 ≤ 1384 / 100000 := by nlinarith
  have hs' : s ≤ 5363 / 10000 * p := by
    have h1 : y / 2 * (1 + y ^ 2) ≤ y / 2 * (101384 / 100000) :=
      mul_le_mul_of_nonneg_left (by linarith) (by linarith)
    linarith
  have hcost : doublingCost lam (2 * t) ≤ 2 * s2 := by
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
  have hss : s ^ 2 ≤ (5363 / 10000 * p) ^ 2 := pow_le_pow_left₀ hs0 hs' 2
  nlinarith [sq_nonneg p]

/-! ### The induction inside the one-short-side class -/

/-- Side `i` short (exponent 0), every other side with exponent 1 or 2. -/
def OneShort (i : Fin 4) (e : Exps) : Prop := e i = 0 ∧ ∀ k, k ≠ i → e k = 1 ∨ e k = 2

/-- The box with one short side `i` and all others doubled. -/
def shortE (i : Fin 4) : Exps := fun k => if k = i then 0 else 1

theorem oneShort_twos_zero {i : Fin 4} {e : Exps} (he : OneShort i e) (ht : twos e = 0) :
    e = shortE i := by
  funext k
  by_cases hk : k = i
  · subst hk; simp [shortE, he.1]
  · simp only [shortE, hk, if_false]
    rcases he.2 k hk with h | h
    · exact h
    · exfalso
      have : k ∈ Finset.univ.filter (fun k => e k = 2) := by simp [h]
      unfold twos at ht
      rw [Finset.card_eq_zero] at ht
      simp [ht] at this

/-- One step of the induction: the bound `A` on all one-short-side boxes with `j` sides `4n` gives
`2A + (29/50) A²` on those with `j + 1`. -/
theorem oneShort_step (G : Exps → ℝ) (nonneg : ∀ i e, 0 ≤ cost G i e)
    (estimate : ∀ i e, cost G i e ≤ 1 / 9 → cost G i (dbl i e) ≤ 29 / 50 * cost G i e ^ 2)
    {A : ℝ} (hA : A ≤ 1 / 9) {j : ℕ}
    (ih : ∀ i e, OneShort i e → twos e = j → cost G i e ≤ A) :
    ∀ i e, OneShort i e → twos e = j + 1 → cost G i e ≤ 2 * A + 29 / 50 * A ^ 2 := by
  intro i e he ht
  obtain ⟨h, hh⟩ : ∃ h, e h = 2 := by
    have : 0 < twos e := by omega
    obtain ⟨h, hmem⟩ := Finset.card_pos.mp this
    exact ⟨h, by simpa using hmem⟩
  have hhi : h ≠ i := by rintro rfl; have := he.1; omega
  set y := Function.update e h 1 with hy
  have hye : dbl h y = e := by
    funext x; rw [dbl_apply]; by_cases hx : x = h
    · subst hx; simp [y, hh]
    · simp [y, hx]
  have hyi : y i = 0 := by simp [y, Ne.symm hhi, he.1]
  have hyS : OneShort i y := by
    refine ⟨hyi, fun k hk => ?_⟩
    by_cases hx : k = h
    · subst hx; simp [y]
    · simpa [y, hx] using he.2 k hk
  have hytw : twos y = j := by have := twos_update_one hh; rw [← hy] at this; omega
  set z := Function.update (dbl i y) h 0 with hz
  have hdy : dbl i y h = 1 := by rw [dbl_apply]; simp [y, hhi]
  have hze : dbl h z = dbl i y := dbl_update_zero hdy
  have hzS : OneShort h z := by
    refine ⟨by simp [z], fun k hk => ?_⟩
    simp only [z, Function.update_of_ne hk, dbl_apply]
    by_cases hxi : k = i
    · subst hxi; simp [hyi]
    · simpa [hxi] using hyS.2 k hxi
  have hztw : twos z = j := by
    rw [← hytw]
    unfold twos
    congr 1
    ext x
    simp only [Finset.mem_filter, Finset.mem_univ, true_and]
    by_cases hx : x = h
    · subst hx; simp [z, y]
    · simp only [z, Function.update_of_ne hx, dbl_apply]
      by_cases hxi : x = i
      · subst hxi; simp [hyi]
      · simp [hxi]
  have h1 := ih i y hyS hytw
  have h2 := ih h z hzS hztw
  have hz0 := nonneg h z
  have h3 : cost G h (dbl i y) ≤ 29 / 50 * cost G h z ^ 2 := by
    rw [← hze]; exact estimate h z (h2.trans hA)
  have hsq : cost G h z ^ 2 ≤ A ^ 2 := pow_le_pow_left₀ hz0 h2 2
  have hy0 := nonneg h y
  have hswap := cost_swap G i h y
  rw [hye] at hswap
  linarith

/-- **The short-side recurrence.** If every one-short-side doubling cost is at most `κ ≤ 1/80`, then
the same costs one scale up are at most `40 κ²`. -/
theorem short_side_recurrence (G : Exps → ℝ) (nonneg : ∀ i e, 0 ≤ cost G i e)
    (estimate : ∀ i e, cost G i e ≤ 1 / 9 → cost G i (dbl i e) ≤ 29 / 50 * cost G i e ^ 2)
    {κ : ℝ} (small : κ ≤ 1 / 80) (seed : ∀ i, cost G i (shortE i) ≤ κ) (i : Fin 4) :
    cost G i (shortE i + 1) ≤ 40 * κ ^ 2 := by
  have hκ : 0 ≤ κ := (nonneg 0 _).trans (seed 0)
  have b0 : ∀ i e, OneShort i e → twos e = 0 → cost G i e ≤ κ := by
    intro i e he ht; rw [oneShort_twos_zero he ht]; exact seed i
  have b1 : ∀ i e, OneShort i e → twos e = 1 → cost G i e ≤ 201 / 100 * κ := fun i e he ht =>
    (oneShort_step G nonneg estimate (by linarith) b0 i e he ht).trans (by nlinarith)
  have b2 : ∀ i e, OneShort i e → twos e = 2 → cost G i e ≤ 405 / 100 * κ := fun i e he ht =>
    (oneShort_step G nonneg estimate (by linarith) b1 i e he ht).trans (by nlinarith)
  have b3 : ∀ i e, OneShort i e → twos e = 3 → cost G i e ≤ 822 / 100 * κ := fun i e he ht =>
    (oneShort_step G nonneg estimate (by linarith) b2 i e he ht).trans (by nlinarith)
  set w : Exps := fun k => if k = i then 0 else 2 with hw
  have hwS : OneShort i w := ⟨by simp [w], fun k hk => Or.inr (by simp [w, hk])⟩
  have hwt : twos w = 3 := by
    have hall : ∀ i : Fin 4, (Finset.univ.filter
        (fun k : Fin 4 => (if k = i then 0 else 2 : ℕ) = 2)).card = 3 := by decide
    exact hall i
  have hc := b3 i w hwS hwt
  have hdw : dbl i w = shortE i + 1 := by
    funext k; rw [dbl_apply]; by_cases hk : k = i
    · subst hk; simp [w, shortE]
    · simp [w, shortE, hk]
  have he := estimate i w (hc.trans (by linarith))
  rw [hdw] at he
  have hsq := pow_le_pow_left₀ (nonneg i w) hc 2
  nlinarith [sq_nonneg κ]

/-! ### In terms of side lengths, and into the conditional gap -/

/-- The box with side `n` in direction `i` and sides `2n` in the others. -/
def shortBox (n : ℕ) (i : Fin 4) : Fin 4 → ℕ := fun k => if k = i then n else 2 * n

/-- `κᵢ(n) = 2 log Z(n; 2n, 2n, 2n) - log Z(2n, 2n, 2n, 2n)`, short side in direction `i`. -/
def shortCost (F : (Fin 4 → ℕ) → ℝ) (n : ℕ) (i : Fin 4) : ℝ := sideCost F i (shortBox n i)

theorem shortCost_eq (F : (Fin 4 → ℕ) → ℝ) (n : ℕ) (i : Fin 4) :
    shortCost F n i = 2 * F (shortBox n i) - F (fun _ => 2 * n) := by
  unfold shortCost sideCost
  congr 2
  funext k
  by_cases hk : k = i
  · subst hk; simp [sideDouble, shortBox]
  · simp [sideDouble, shortBox, hk]

theorem atScale_shortE (F : (Fin 4 → ℕ) → ℝ) (n : ℕ) (i : Fin 4) :
    cost (atScale F n) i (shortE i) = shortCost F n i := by
  rw [cost_atScale, shortCost]
  congr 1
  funext k
  by_cases hk : k = i
  · subst hk; simp [shortE, shortBox]
  · simp [shortE, shortBox, hk, mul_comm]

theorem atScale_shortE_succ (F : (Fin 4 → ℕ) → ℝ) (n : ℕ) (i : Fin 4) :
    cost (atScale F n) i (shortE i + 1) = shortCost F (2 * n) i := by
  rw [cost_atScale, shortCost]
  congr 1
  funext k
  by_cases hk : k = i
  · subst hk; simp [shortE, shortBox, mul_comm]
  · simp [shortE, shortBox, hk]; ring

theorem box_short_recurrence (F : (Fin 4 → ℕ) → ℝ) {n : ℕ} (hn : 1 ≤ n)
    (nonneg : ∀ i x, (∀ k, 0 < x k) → 0 ≤ sideCost F i x)
    (estimate : ∀ i x, (∀ k, 0 < x k) →
      sideCost F i x ≤ 1 / 9 → sideCost F i (sideDouble i x) ≤ 29 / 50 * sideCost F i x ^ 2)
    {κ : ℝ} (small : κ ≤ 1 / 80) (seed : ∀ i, shortCost F n i ≤ κ) (i : Fin 4) :
    shortCost F (2 * n) i ≤ 40 * κ ^ 2 := by
  have hpos : ∀ e : Exps, ∀ k, 0 < n * 2 ^ e k := fun e k => by positivity
  have h := short_side_recurrence (atScale F n)
    (fun i e => by rw [cost_atScale]; exact nonneg i _ (hpos e))
    (fun i e he => by
      simp only [cost_atScale] at he ⊢
      rw [atScale_dbl]
      exact estimate i _ (hpos e) he)
    small (fun i => by rw [atScale_shortE]; exact seed i) i
  rwa [atScale_shortE_succ] at h

/-- The short-side cost is one term of the hypercubic defect. -/
theorem shortCost_le_boxDefect (F : (Fin 4 → ℕ) → ℝ) {n : ℕ} (hn : 1 ≤ n)
    (nonneg : ∀ i x, (∀ k, 0 < x k) → 0 ≤ sideCost F i x) (i : Fin 4) :
    shortCost F n i ≤ boxDefect F n := by
  have hpos : ∀ e : Exps, ∀ k, 0 < n * 2 ^ e k := fun e k => by positivity
  have h := cost_le_defect_binary (atScale F n)
    (fun j e => by rw [cost_atScale]; exact nonneg j _ (hpos e))
    (i := i) (e := shortE i) (fun k => by unfold shortE; split_ifs <;> omega) (by simp [shortE])
  rw [atScale_shortE, cubeDefect_atScale] at h
  simpa using h

/-- **The conditional gap from one purity.** Inputs: the spectral hypotheses in every direction and the
seed `κᵢ(n) ≤ 1/80` for each `i`. `r k` is `λ₂/λ₁` on a transverse cube of side `2ᵏ · 2n`. -/
theorem conditional_gap_short (F : (Fin 4 → ℕ) → ℝ) {n : ℕ} (hn : 1 ≤ n) {r : ℕ → ℝ} (i : Fin 4)
    (hr : ∀ k, 0 ≤ r k)
    (nonneg : ∀ i x, (∀ k, 0 < x k) → 0 ≤ sideCost F i x)
    (estimate : ∀ i x, (∀ k, 0 < x k) →
      sideCost F i x ≤ 1 / 9 → sideCost F i (sideDouble i x) ≤ 29 / 50 * sideCost F i x ^ 2)
    (seed : ∀ i, shortCost F n i ≤ 1 / 80)
    (from_cost : ∀ k, r k ^ (2 ^ k * n) ≤ exp (shortCost F (2 ^ k * n) i) - 1) :
    ∀ k, r k ^ n ≤ 1 / 2 := by
  have hk1 : ∀ k, 1 ≤ 2 ^ k * n := fun k => by
    have : 1 ≤ 2 ^ k := Nat.one_le_two_pow
    nlinarith
  let δ : ℕ → ℝ := fun k => Finset.univ.sup' Finset.univ_nonempty (fun i => shortCost F (2 ^ k * n) i)
  have hδle : ∀ k i, shortCost F (2 ^ k * n) i ≤ δ k := fun k i =>
    Finset.le_sup' (fun i => shortCost F (2 ^ k * n) i) (Finset.mem_univ i)
  have hbox : ∀ k i, 0 ≤ shortCost F (2 ^ k * n) i := fun k i =>
    nonneg i _ (fun l => by unfold shortBox; have := hk1 k; split_ifs <;> omega)
  refine conditional_gap_ratio (c₁ := 40) (by norm_num) (δ := δ)
    (fun k => (hbox k 0).trans (hδle k 0)) hr ?_ ?_
    (fun k => (from_cost k).trans (by have := exp_le_exp.mpr (hδle k i); linarith))
  · intro k hk
    refine Finset.sup'_le _ _ fun i _ => ?_
    have := box_short_recurrence F (hk1 k) nonneg estimate (κ := δ k) (by norm_num at hk ⊢; linarith)
      (hδle k) i
    simpa [pow_succ, mul_assoc, mul_comm, mul_left_comm] using this
  · refine Finset.sup'_le _ _ fun i _ => ?_
    have := seed i
    simpa using (by norm_num; linarith : shortCost F n i ≤ 1 / (2 * 40))

/-! ### Z₂ lattice gauge theory -/

theorem z2_estimate_sharp {β : ℝ} (hβ : 0 ≤ β) (i : Fin 4) (x : Fin 4 → ℕ) (hx : ∀ k, 0 < x k)
    (h : sideCost (z2F β) i x ≤ 1 / 9) :
    sideCost (z2F β) i (sideDouble i x) ≤ 29 / 50 * sideCost (z2F β) i x ^ 2 := by
  have hd : sideDouble i x i = 2 * x i := by simp [sideDouble]
  have hx' : 0 < sideDouble i x i := by rw [hd]; have := hx i; omega
  rw [sideCost_eq i hβ _ hx', transverse_sideDouble, hd]
  rw [sideCost_eq i hβ x (hx i)] at h ⊢
  exact doubling_cost_sharp (z2Spec_spectrum i hβ _) (hx i) h

/-- `(λ₂/λ₁)ⁿ ≤ e^{κᵢ(n)} - 1` on the transverse cube of side `2n`. -/
theorem z2_from_shortCost {β : ℝ} (hβ : 0 ≤ β) (i : Fin 4) {n : ℕ} (hn : 1 ≤ n) :
    z2Ratio hβ i (2 * n) ^ n ≤ exp (shortCost (z2F β) n i) - 1 := by
  have ht : transverse i (shortBox n i) = fun _ => 2 * n := by
    funext j; simp [transverse, shortBox, j.2]
  have hside : shortBox n i i = n := by simp [shortBox]
  have hS := z2Spec_spectrum i hβ (fun _ => 2 * n)
  have hcost : shortCost (z2F β) n i = -log (purity (z2Spec i hβ (fun _ => 2 * n)) n) := by
    unfold shortCost
    rw [sideCost_eq i hβ _ (by rw [hside]; omega), ht, hside]
    exact doublingCost_eq hS hn
  refine ratio_pow_le_defect hS hn ?_
  rw [← exp_log (purity_pos hS hn)]
  exact exp_le_exp.mpr (by linarith)

/-- **Z₂, from the short-side seed alone.** If `κᵢ(n) ≤ 1/80` in every direction then
`(λ₂/λ₁)ⁿ ≤ 1/2` on every transverse cube of side `2ᵏ · 2n`. -/
theorem z2_gap_of_short_seed {β : ℝ} (hβ : 0 ≤ β) {n : ℕ} (hn : 1 ≤ n)
    (seed : ∀ i, shortCost (z2F β) n i ≤ 1 / 80) (i : Fin 4) (k : ℕ) :
    z2Ratio hβ i (2 * (2 ^ k * n)) ^ n ≤ 1 / 2 :=
  conditional_gap_short (z2F β) hn (r := fun k => z2Ratio hβ i (2 * (2 ^ k * n))) i
    (fun _ => z2Ratio_nonneg hβ i _) (z2_nonneg hβ) (z2_estimate_sharp hβ) seed
    (fun k => z2_from_shortCost hβ i (Nat.one_le_iff_ne_zero.mpr (by positivity))) k

/-- The elementary seed with the new threshold: `0 ≤ β ≤ 1/(80 · 192 · n⁴)`. -/
theorem z2_gap_short {β : ℝ} (hβ : 0 ≤ β) {n : ℕ} (hn : 1 ≤ n)
    (hβn : β ≤ 1 / (80 * 192 * (n : ℝ) ^ 4)) (i : Fin 4) (k : ℕ) :
    z2Ratio hβ i (2 * (2 ^ k * n)) ^ n ≤ 1 / 2 := by
  refine z2_gap_of_short_seed hβ hn (fun i => ?_) i k
  have h1 := shortCost_le_boxDefect (z2F β) hn (z2_nonneg hβ) i
  have h2 := seed_of_wilson_bound (z2F β) hn (fun x _ => z2F_bounds hβ x)
  have hn4 : (0 : ℝ) < (n : ℝ) ^ 4 := by
    have : (0 : ℝ) < n := by exact_mod_cast hn
    positivity
  rw [le_div_iff₀ (by positivity)] at hβn
  linarith

end YMFormal
