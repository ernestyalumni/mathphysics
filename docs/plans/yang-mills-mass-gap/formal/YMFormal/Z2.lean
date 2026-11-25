import YMFormal.Layer
import YMFormal.Seed

/-!
# Z₂ lattice gauge theory on a periodic four-dimensional box

A box has sides `x : Fin 4 → ℕ`. Sites are `s : (k : Fin 4) → Fin (x k)`; a link is a site and a
direction, so that with a side of length 1 distinct links stay distinct. The energy is
`∑_s ∑_{μ<ν} (1 - U_p)` and `Z` is normalized by the number of configurations (the handbook's
convention with normalized Haar measure). `Z` is defined from the lattice alone.

`z2Z_eq_partition`: for every direction `i`, `Z(x) = ∑_j λ_jᵗ` with `t = x i`, where the `λ_j` are the
eigenvalues of a nonnegative transfer matrix that depends only on the other three sides.

`z2_gap` is the end-to-end theorem. It has no hypothesis beyond the range of `β`.
-/

namespace YMFormal

open Real Finset

/-! ### The lattice -/

abbrev Site (x : Fin 4 → ℕ) := (k : Fin 4) → Fin (x k)

/-- One step in direction `μ`, periodic. -/
def step {x : Fin 4 → ℕ} (μ : Fin 4) (s : Site x) : Site x :=
  Function.update s μ (finRotate (x μ) (s μ))

abbrev Config (x : Fin 4 → ℕ) := Site x → Fin 4 → ℤˣ

/-- The plaquette at `s` in the plane `(μ, ν)`. -/
def plaq {x : Fin 4 → ℕ} (U : Config x) (s : Site x) (μ ν : Fin 4) : ℤˣ :=
  U s μ * U (step μ s) ν * U (step ν s) μ * U s ν

/-- The Wilson energy. -/
def energy {x : Fin 4 → ℕ} (U : Config x) : ℝ :=
  ∑ s, ∑ μ, ∑ ν, if μ < ν then 1 - sgn (plaq U s μ ν) else 0

/-- The partition function. -/
noncomputable def z2Z (β : ℝ) (x : Fin 4 → ℕ) : ℝ :=
  (Fintype.card (Config x) : ℝ)⁻¹ * ∑ U : Config x, exp (-(β * energy U))

noncomputable def z2F (β : ℝ) (x : Fin 4 → ℕ) : ℝ := log (z2Z β x)

/-! ### Elementary bounds -/

theorem card_site (x : Fin 4 → ℕ) : (Fintype.card (Site x) : ℝ) = ∏ k, (x k : ℝ) := by
  simp [Site]

theorem energy_nonneg {x : Fin 4 → ℕ} (U : Config x) : 0 ≤ energy U := by
  refine Finset.sum_nonneg fun s _ => Finset.sum_nonneg fun μ _ => Finset.sum_nonneg fun ν _ => ?_
  split_ifs
  · rcases sgn_cases (plaq U s μ ν) with h | h <;> rw [h] <;> norm_num
  · exact le_refl 0

theorem energy_le {x : Fin 4 → ℕ} (U : Config x) : energy U ≤ 2 * (6 * ∏ k, (x k : ℝ)) := by
  have h : ∀ s : Site x, (∑ μ, ∑ ν, if μ < ν then 1 - sgn (plaq U s μ ν) else 0) ≤ 12 := by
    intro s
    have hb : ∀ μ ν : Fin 4, (if μ < ν then 1 - sgn (plaq U s μ ν) else 0) ≤
        (if μ < ν then (2 : ℝ) else 0) := by
      intro μ ν
      split_ifs
      · rcases sgn_cases (plaq U s μ ν) with h | h <;> rw [h] <;> norm_num
      · exact le_refl 0
    calc _ ≤ ∑ μ : Fin 4, ∑ ν : Fin 4, (if μ < ν then (2 : ℝ) else 0) :=
          Finset.sum_le_sum fun μ _ => Finset.sum_le_sum fun ν _ => hb μ ν
      _ = 12 := by simp [Fin.sum_univ_four]; norm_num
  calc energy U ≤ ∑ _s : Site x, (12 : ℝ) := Finset.sum_le_sum fun s _ => h s
    _ = 2 * (6 * ∏ k, (x k : ℝ)) := by
      rw [Finset.sum_const, Finset.card_univ, nsmul_eq_mul, card_site]; ring

theorem z2Z_bounds {β : ℝ} (hβ : 0 ≤ β) (x : Fin 4 → ℕ) :
    exp (-(2 * β * (6 * ∏ k, (x k : ℝ)))) ≤ z2Z β x ∧ z2Z β x ≤ 1 := by
  have hc : (0 : ℝ) < Fintype.card (Config x) := by exact_mod_cast Fintype.card_pos
  unfold z2Z
  constructor
  · rw [inv_mul_eq_div, le_div_iff₀ hc]
    calc _ = ∑ _U : Config x, exp (-(2 * β * (6 * ∏ k, (x k : ℝ)))) := by
          rw [Finset.sum_const, Finset.card_univ, nsmul_eq_mul]; ring
      _ ≤ _ := Finset.sum_le_sum fun U _ => exp_le_exp.mpr (by nlinarith [energy_le U])
  · rw [inv_mul_eq_div, div_le_one hc]
    calc _ ≤ ∑ _U : Config x, (1 : ℝ) := Finset.sum_le_sum fun U _ => by
          rw [exp_le_one_iff]; nlinarith [energy_nonneg U]
      _ = _ := by simp

theorem z2F_bounds {β : ℝ} (hβ : 0 ≤ β) (x : Fin 4 → ℕ) :
    -(2 * β * (6 * ∏ k, (x k : ℝ))) ≤ z2F β x ∧ z2F β x ≤ 0 := by
  obtain ⟨h1, h2⟩ := z2Z_bounds hβ x
  have hpos : 0 < z2Z β x := lt_of_lt_of_le (exp_pos _) h1
  refine ⟨?_, log_nonpos hpos.le h2⟩
  have := log_le_log (exp_pos _) h1
  rwa [log_exp] at this

/-! ### Slicing in direction `i` -/

section Slice
variable (i : Fin 4)

/-- The transverse directions. -/
abbrev Dir := {j : Fin 4 // j ≠ i}

/-- The sites of a slice with transverse sides `y`. -/
abbrev SliceSite (y : Dir i → ℕ) := (j : Dir i) → Fin (y j)

variable {i}

def sliceStep {y : Dir i → ℕ} (j : Dir i) (a : SliceSite i y) : SliceSite i y :=
  Function.update a j (finRotate (y j) (a j))

/-- The spatial energy of a slice. -/
def sliceE {y : Dir i → ℕ} (σ : SliceSite i y → Dir i → ℤˣ) : ℝ :=
  ∑ a, ∑ j : Dir i, ∑ k : Dir i, if j.1 < k.1 then
    1 - sgn (σ a j * σ (sliceStep j a) k * σ (sliceStep k a) j * σ a k) else 0

variable (i)

/-- The transverse sides of a box. -/
def transverse (x : Fin 4 → ℕ) : Dir i → ℕ := fun j => x j.1

/-- A site is a time and a slice site. -/
def siteEquiv (x : Fin 4 → ℕ) : Site x ≃ Fin (x i) × SliceSite i (transverse i x) :=
  Equiv.piSplitAt i (fun k => Fin (x k))

theorem siteEquiv_fst (x : Fin 4 → ℕ) (s : Site x) : (siteEquiv i x s).1 = s i := rfl
theorem siteEquiv_snd (x : Fin 4 → ℕ) (s : Site x) (j : Dir i) : (siteEquiv i x s).2 j = s j.1 := rfl

theorem siteEquiv_step_time (x : Fin 4 → ℕ) (s : Site x) :
    siteEquiv i x (step i s) = (finRotate (x i) (siteEquiv i x s).1, (siteEquiv i x s).2) := by
  refine Prod.ext ?_ ?_
  · simp [siteEquiv_fst, step]
  · funext j
    rw [siteEquiv_snd]
    simp only [siteEquiv_snd, step]
    rw [Function.update_of_ne j.2]

theorem siteEquiv_step_space (x : Fin 4 → ℕ) (s : Site x) (j : Dir i) :
    siteEquiv i x (step j.1 s) = ((siteEquiv i x s).1, sliceStep j (siteEquiv i x s).2) := by
  refine Prod.ext ?_ ?_
  · rw [siteEquiv_fst]
    simp only [siteEquiv_fst, step]
    rw [Function.update_of_ne (Ne.symm j.2)]
  · funext k
    rw [siteEquiv_snd]
    by_cases hk : k = j
    · subst hk
      exact (Function.update_self k.1 _ s).trans
        (Function.update_self k (finRotate (transverse i x k) ((siteEquiv i x s).2 k)) (siteEquiv i x s).2).symm
    · have : k.1 ≠ j.1 := fun h => hk (Subtype.ext h)
      simp only [step, sliceStep]
      rw [Function.update_of_ne this, Function.update_of_ne hk, siteEquiv_snd]

/-- A configuration is the spatial and the temporal links of every slice. -/
def configEquiv (x : Fin 4 → ℕ) : Config x ≃
    (Fin (x i) → SliceSite i (transverse i x) → Dir i → ℤˣ) × (Fin (x i) → SliceSite i (transverse i x) → ℤˣ) where
  toFun U := (fun t a j => U ((siteEquiv i x).symm (t, a)) j.1, fun t a => U ((siteEquiv i x).symm (t, a)) i)
  invFun p := fun s μ => if h : μ = i then p.2 (siteEquiv i x s).1 (siteEquiv i x s).2
    else p.1 (siteEquiv i x s).1 (siteEquiv i x s).2 ⟨μ, h⟩
  left_inv U := by
    funext s μ
    by_cases h : μ = i
    · subst h; simp
    · simp [h]
  right_inv p := by
    refine Prod.ext ?_ ?_
    · funext t a j
      simp [j.2]
    · funext t a
      simp

/-- Splitting the six planes into those inside a slice and those containing `i`. -/
theorem planes_split (f : Fin 4 → Fin 4 → ℝ) (hf : ∀ μ ν, f μ ν = f ν μ) :
    (∑ μ, ∑ ν, if μ < ν then f μ ν else 0) =
      (∑ j : Dir i, ∑ k : Dir i, if j.1 < k.1 then f j.1 k.1 else 0) + ∑ j : Dir i, f j.1 i := by
  have sub : ∀ g : Fin 4 → ℝ, ∑ j : Dir i, g j.1 = ∑ j, if j ≠ i then g j else 0 := by
    intro g
    rw [← Finset.sum_filter]
    exact (Finset.sum_subtype _ (by simp) g).symm
  have h1 : (∑ j : Dir i, ∑ k : Dir i, if j.1 < k.1 then f j.1 k.1 else 0) =
      ∑ j, if j ≠ i then (∑ k, if k ≠ i then (if j < k then f j k else 0) else 0) else 0 := by
    rw [← sub (fun j => ∑ k, if k ≠ i then (if j < k then f j k else 0) else 0)]
    exact Finset.sum_congr rfl fun j _ => sub (fun k => if j.1 < k then f j.1 k else 0)
  rw [h1, sub (fun j => f j i)]
  have h10 := hf 1 0; have h20 := hf 2 0; have h30 := hf 3 0
  have h21 := hf 2 1; have h31 := hf 3 1; have h32 := hf 3 2
  fin_cases i <;> simp [Fin.sum_univ_four] <;> linarith

theorem plaq_symm {x : Fin 4 → ℕ} (U : Config x) (s : Site x) (μ ν : Fin 4) :
    plaq U s μ ν = plaq U s ν μ := by
  unfold plaq
  ac_rfl

theorem energy_slice (x : Fin 4 → ℕ) (U : Config x) :
    energy U = ∑ t, sliceE ((configEquiv i x U).1 t) +
      ∑ t, tempE (fun j a => sliceStep j a) ((configEquiv i x U).1 t) ((configEquiv i x U).2 t)
        ((configEquiv i x U).1 (finRotate (x i) t)) := by
  have htime : ∀ t a, step i ((siteEquiv i x).symm (t, a)) =
      (siteEquiv i x).symm (finRotate (x i) t, a) := by
    intro t a
    rw [Equiv.eq_symm_apply, siteEquiv_step_time]; simp
  have hspace : ∀ t a (j : Dir i), step j.1 ((siteEquiv i x).symm (t, a)) =
      (siteEquiv i x).symm (t, sliceStep j a) := by
    intro t a j
    rw [Equiv.eq_symm_apply, siteEquiv_step_space]; simp
  unfold energy
  rw [← (siteEquiv i x).symm.sum_comp, Fintype.sum_prod_type, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun t _ => ?_
  unfold sliceE tempE
  rw [← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun a _ => ?_
  rw [planes_split i (fun μ ν => 1 - sgn (plaq U ((siteEquiv i x).symm (t, a)) μ ν))
    (fun μ ν => by rw [plaq_symm])]
  congr 1
  · refine Finset.sum_congr rfl fun j _ => Finset.sum_congr rfl fun k _ => ?_
    simp only [plaq, hspace, configEquiv, Equiv.coe_fn_mk]
  · refine Finset.sum_congr rfl fun j _ => ?_
    simp only [plaq, hspace, htime, configEquiv, Equiv.coe_fn_mk]

instance instA (t : ℕ) (y : Dir i → ℕ) : Fintype (Fin t → SliceSite i y → Dir i → ℤˣ) := inferInstance
instance instB (t : ℕ) (y : Dir i → ℕ) : Fintype (Fin t → SliceSite i y → ℤˣ) := inferInstance

/-- The eigenvalues of the transfer matrix in direction `i`, for transverse sides `y`. -/
noncomputable def z2Spec {β : ℝ} (hβ : 0 ≤ β) (y : Dir i → ℕ) : ℕ → ℝ :=
  layerSpec (fun j a => sliceStep (i := i) (y := y) j a) sliceE hβ

theorem z2Spec_spectrum {β : ℝ} (hβ : 0 ≤ β) (y : Dir i → ℕ) : Spectrum (z2Spec i hβ y) :=
  layerSpec_spectrum _ _ hβ

/-- **Trace formula in direction `i`.** -/
theorem z2Z_eq_partition {β : ℝ} (hβ : 0 ≤ β) (x : Fin 4 → ℕ) (hx : 0 < x i) :
    z2Z β x = partition (z2Spec i hβ (transverse i x)) (x i) := by
  unfold z2Spec
  rw [← layerZ_eq_partition _ _ hβ hx, z2Z, layerZ]
  congr 1
  · rw [Fintype.card_congr (configEquiv i x), Fintype.card_prod]
    congr 1
    unfold layerCard
    rw [mul_pow, Fintype.card_fun (α := Fin (x i)), Fintype.card_fun (α := Fin (x i)), Fintype.card_fin]
    push_cast
    rfl
  · rw [← (configEquiv i x).symm.sum_comp, Fintype.sum_prod_type]
    refine Finset.sum_congr rfl fun σ _ => Finset.sum_congr rfl fun τ _ => ?_
    rw [energy_slice i x, Equiv.apply_symm_apply]

theorem transverse_sideDouble (x : Fin 4 → ℕ) : transverse i (sideDouble i x) = transverse i x := by
  funext j
  simp [transverse, sideDouble, Function.update_of_ne j.2]

theorem sideCost_eq {β : ℝ} (hβ : 0 ≤ β) (x : Fin 4 → ℕ) (hx : 0 < x i) :
    sideCost (z2F β) i x = doublingCost (z2Spec i hβ (transverse i x)) (x i) := by
  have hd : sideDouble i x i = 2 * x i := by simp [sideDouble]
  unfold sideCost doublingCost z2F
  rw [z2Z_eq_partition i hβ x hx, z2Z_eq_partition i hβ (sideDouble i x) (by rw [hd]; omega),
    transverse_sideDouble, hd]

end Slice

/-! ### The hypotheses of the chain, discharged -/

theorem z2_nonneg {β : ℝ} (hβ : 0 ≤ β) (i : Fin 4) (x : Fin 4 → ℕ) (hx : ∀ k, 0 < x k) :
    0 ≤ sideCost (z2F β) i x := by
  rw [sideCost_eq i hβ x (hx i)]
  exact doubling_cost_nonneg (z2Spec_spectrum i hβ _) (hx i)

theorem z2_estimate {β : ℝ} (hβ : 0 ≤ β) (i : Fin 4) (x : Fin 4 → ℕ) (hx : ∀ k, 0 < x k)
    (h : sideCost (z2F β) i x ≤ 1 / 4) :
    sideCost (z2F β) i (sideDouble i x) ≤ 4 * sideCost (z2F β) i x ^ 2 := by
  have hd : sideDouble i x i = 2 * x i := by simp [sideDouble]
  have hx' : 0 < sideDouble i x i := by rw [hd]; have := hx i; omega
  rw [sideCost_eq i hβ _ hx', transverse_sideDouble, hd]
  rw [sideCost_eq i hβ x (hx i)] at h ⊢
  exact doubling_cost_estimate (z2Spec_spectrum i hβ _) (hx i) h

/-- `λ₂/λ₁` for the transfer matrix in direction `i` on the transverse cube of side `N`. -/
noncomputable def z2Ratio {β : ℝ} (hβ : 0 ≤ β) (i : Fin 4) (N : ℕ) : ℝ :=
  z2Spec i hβ (fun _ => N) 1 / z2Spec i hβ (fun _ => N) 0

theorem z2Ratio_nonneg {β : ℝ} (hβ : 0 ≤ β) (i : Fin 4) (N : ℕ) : 0 ≤ z2Ratio hβ i N :=
  div_nonneg ((z2Spec_spectrum i hβ _).nonneg 1) (z2Spec_spectrum i hβ _).top_pos.le

/-- `(λ₂/λ₁)ᴺ ≤ e^{Δ(N)} - 1` on the transverse cube of side `N`. -/
theorem z2_from_defect {β : ℝ} (hβ : 0 ≤ β) (i : Fin 4) {N : ℕ} (hN : 1 ≤ N) :
    z2Ratio hβ i N ^ N ≤ exp (boxDefect (z2F β) N) - 1 := by
  have hS := z2Spec_spectrum i hβ (fun _ => N)
  have hpos : ∀ e : Exps, ∀ k, 0 < N * 2 ^ e k := fun e k => by positivity
  have hle : sideCost (z2F β) i (fun _ => N) ≤ boxDefect (z2F β) N := by
    have h := cost_le_defect_binary (atScale (z2F β) N)
      (fun j e => by rw [cost_atScale]; exact z2_nonneg hβ j _ (hpos e))
      (i := i) (e := fun _ => 0) (fun _ => by simp) rfl
    rw [cost_atScale, cubeDefect_atScale] at h
    simpa using h
  have hcost : sideCost (z2F β) i (fun _ => N) = -log (purity (z2Spec i hβ (fun _ => N)) N) := by
    rw [sideCost_eq i hβ (fun _ => N) hN]
    exact doublingCost_eq hS hN
  refine ratio_pow_le_defect hS hN ?_
  rw [← exp_log (purity_pos hS hN)]
  exact exp_le_exp.mpr (by linarith)

/-- **Z₂ lattice gauge theory, end to end.** For `0 ≤ β ≤ 1/(7480 · 192 · n⁴)`, in every direction
`i` and for every `k`, the transfer matrix on the transverse cube of side `2ᵏ n` satisfies
`(λ₂/λ₁)ⁿ ≤ 1/2`. -/
theorem z2_gap {β : ℝ} (hβ : 0 ≤ β) {n : ℕ} (hn : 1 ≤ n)
    (hβn : β ≤ 1 / (7480 * 192 * (n : ℝ) ^ 4)) (i : Fin 4) (k : ℕ) :
    z2Ratio hβ i (2 ^ k * n) ^ n ≤ 1 / 2 :=
  conditional_gap_strong_coupling (z2F β) hn hβn (fun x _ => z2F_bounds hβ x)
    (r := fun k => z2Ratio hβ i (2 ^ k * n)) (fun _ => z2Ratio_nonneg hβ i _)
    (z2_nonneg hβ) (z2_estimate hβ)
    (fun k => z2_from_defect hβ i (Nat.one_le_iff_ne_zero.mpr (by positivity))) k

/-- The gap form: if `λ₂ > 0` then `-log (λ₂/λ₁) ≥ log 2 / n`. If `λ₂ = 0` the gap is infinite. -/
theorem z2_gap_log {β : ℝ} (hβ : 0 ≤ β) {n : ℕ} (hn : 1 ≤ n)
    (hβn : β ≤ 1 / (7480 * 192 * (n : ℝ) ^ 4)) (i : Fin 4) (k : ℕ)
    (hr : 0 < z2Ratio hβ i (2 ^ k * n)) : log 2 / n ≤ -log (z2Ratio hβ i (2 ^ k * n)) := by
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  have h := log_le_log (pow_pos hr n) (z2_gap hβ hn hβn i k)
  rw [log_pow, one_div, log_inv] at h
  rw [div_le_iff₀ hn']
  linarith

/-- Sanity check on the definition: `Z = 1` on the box with all sides 1, because every plaquette is
trivial there. -/
theorem z2Z_one (β : ℝ) : z2Z β (fun _ => 1) = 1 := by
  have hE : ∀ U : Config (fun _ => 1), energy U = 0 := by
    intro U
    refine Finset.sum_eq_zero fun s _ => Finset.sum_eq_zero fun μ _ => Finset.sum_eq_zero fun ν _ => ?_
    have hstep : ∀ ρ, step ρ s = s := fun ρ => Subsingleton.elim _ _
    have hp : plaq U s μ ν = 1 := by
      unfold plaq
      rw [hstep, hstep]
      have : U s μ * U s ν * U s μ * U s ν = (U s μ * U s μ) * (U s ν * U s ν) := by ac_rfl
      rw [this, Int.units_mul_self, Int.units_mul_self, one_mul]
    simp [hp]
  have hc : (Fintype.card (Config (fun _ => 1)) : ℝ) ≠ 0 := by
    exact_mod_cast Fintype.card_pos.ne'
  unfold z2Z
  rw [Finset.sum_congr rfl fun U _ => by rw [hE U, mul_zero, neg_zero, exp_zero]]
  rw [Finset.sum_const, Finset.card_univ, nsmul_eq_mul, mul_one, inv_mul_cancel₀ hc]

end YMFormal
