import YMFormal.MatrixSpectrum

/-!
# A layered `±1` system and its transfer matrix

`Λ` is the set of sites of one time slice, `D` the transverse directions, `shift j` the step in
direction `j`. A configuration has spatial links `σ t : Λ → D → ℤˣ` and temporal links
`τ t : Λ → ℤˣ` for `t : Fin T`, with time periodic through `finRotate`. The spatial energy of a slice
is an arbitrary function `Es`. Proved here:
- `layerZ_eq_trace`: the partition function is `Tr Wᵀ` for the transfer matrix `W`;
- `transfer_posSemidef`: `W` is positive semidefinite for `β ≥ 0`, from the expansion
  `e^{βuv} = cosh β + u v sinh β` for `u, v = ±1`;
- `transfer_trace_pos`.
-/

namespace YMFormal

open Real Finset Matrix

/-- A sign as a real number. -/
def sgn (u : ℤˣ) : ℝ := ((u : ℤ) : ℝ)

@[simp] theorem sgn_one : sgn 1 = 1 := by simp [sgn]
@[simp] theorem sgn_mul (u v : ℤˣ) : sgn (u * v) = sgn u * sgn v := by simp [sgn]
theorem sgn_cases (u : ℤˣ) : sgn u = 1 ∨ sgn u = -1 := by
  rcases Int.units_eq_one_or u with h | h <;> simp [sgn, h]
theorem sgn_mul_self (u : ℤˣ) : sgn u * sgn u = 1 := by
  rcases sgn_cases u with h | h <;> simp [h]
theorem sgn_prod {ι : Type*} (s : Finset ι) (f : ι → ℤˣ) : sgn (∏ i ∈ s, f i) = ∏ i ∈ s, sgn (f i) := by
  classical
  induction s using Finset.induction_on with
  | empty => simp
  | insert a s ha ih => rw [Finset.prod_insert ha, Finset.prod_insert ha, sgn_mul, ih]

/-! ### Closed walks and traces of powers -/

section Walk
variable {S : Type*} [Fintype S] [DecidableEq S]

theorem walk_sum (W : Matrix S S ℝ) (n : ℕ) : ∀ g : S → S → ℝ,
    ∑ σ : Fin (n + 1) → S, (∏ j : Fin n, W (σ j.castSucc) (σ j.succ)) * g (σ (Fin.last n)) (σ 0) =
      ∑ a, ∑ b, (W ^ n) a b * g b a := by
  induction n with
  | zero =>
    intro g
    rw [← (Equiv.funUnique (Fin 1) S).symm.sum_comp]
    simp [Matrix.one_apply]
  | succ n ih =>
    intro g
    rw [← (Fin.consEquiv (fun _ : Fin (n + 2) => S)).sum_comp, Fintype.sum_prod_type]
    refine Finset.sum_congr rfl fun s _ => ?_
    have h := ih (fun b a => W s a * g b s)
    have e : ∀ σ : Fin (n + 1) → S,
        (∏ j : Fin (n + 1), W ((Fin.consEquiv (fun _ : Fin (n + 2) => S)) (s, σ) j.castSucc)
          ((Fin.consEquiv (fun _ : Fin (n + 2) => S)) (s, σ) j.succ)) *
          g ((Fin.consEquiv (fun _ : Fin (n + 2) => S)) (s, σ) (Fin.last (n + 1)))
            ((Fin.consEquiv (fun _ : Fin (n + 2) => S)) (s, σ) 0) =
        (∏ j : Fin n, W (σ j.castSucc) (σ j.succ)) * (W s (σ 0) * g (σ (Fin.last n)) s) := by
      intro σ
      rw [Fin.prod_univ_succ]
      simp only [Fin.consEquiv_apply, Fin.castSucc_zero, Fin.cons_zero, Fin.cons_succ,
        Fin.castSucc_succ, ← Fin.succ_last]
      ring
    rw [Finset.sum_congr rfl fun σ _ => e σ, h]
    simp only [pow_succ' W n, Matrix.mul_apply, Finset.sum_mul]
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun b _ => Finset.sum_congr rfl fun a _ => ?_
    ring

/-- A sum over closed walks is the trace of a power. -/
theorem closed_walk_sum (W : Matrix S S ℝ) {t : ℕ} (ht : 1 ≤ t) :
    ∑ σ : Fin t → S, ∏ j, W (σ j) (σ (finRotate t j)) = (W ^ t).trace := by
  obtain ⟨n, rfl⟩ : ∃ n, t = n + 1 := ⟨t - 1, by omega⟩
  have h := walk_sum W n (fun b a => W b a)
  have e : ∀ σ : Fin (n + 1) → S, ∏ j, W (σ j) (σ (finRotate (n + 1) j)) =
      (∏ j : Fin n, W (σ j.castSucc) (σ j.succ)) * W (σ (Fin.last n)) (σ 0) := by
    intro σ
    rw [Fin.prod_univ_castSucc, finRotate_last]
    congr 1
    refine Finset.prod_congr rfl fun j _ => ?_
    congr 2
    ext
    rw [coe_finRotate_of_ne_last (by simpa using (Fin.castSucc_lt_last j).ne)]
    simp
  rw [Finset.sum_congr rfl fun σ _ => e σ, h, pow_succ]
  simp [Matrix.trace, Matrix.mul_apply]

end Walk

/-! ### The layered system -/

section Layer
variable {Λ D : Type*} [Fintype Λ] [DecidableEq Λ] [Fintype D] [DecidableEq D]

/-- The temporal energy between two slices. -/
def tempE (shift : D → Λ → Λ) (σ : Λ → D → ℤˣ) (τ : Λ → ℤˣ) (σ' : Λ → D → ℤˣ) : ℝ :=
  ∑ a, ∑ j, (1 - sgn (σ a j * τ (shift j a) * σ' a j * τ a))

/-- Number of configurations of one layer, spatial and temporal. -/
def layerCard (Λ D : Type*) [Fintype Λ] [DecidableEq Λ] [Fintype D] [DecidableEq D] : ℝ :=
  Fintype.card (Λ → D → ℤˣ) * Fintype.card (Λ → ℤˣ)

theorem layerCard_pos : 0 < layerCard Λ D := by
  unfold layerCard
  have h1 : 0 < Fintype.card (Λ → D → ℤˣ) := Fintype.card_pos
  have h2 : 0 < Fintype.card (Λ → ℤˣ) := Fintype.card_pos
  positivity

/-- The transfer matrix. -/
noncomputable def transfer (shift : D → Λ → Λ) (Es : (Λ → D → ℤˣ) → ℝ) (β : ℝ) :
    Matrix (Λ → D → ℤˣ) (Λ → D → ℤˣ) ℝ := fun σ σ' =>
  (layerCard Λ D)⁻¹ * (exp (-(β * Es σ) / 2) * exp (-(β * Es σ') / 2)) *
    ∑ τ : Λ → ℤˣ, exp (-(β * tempE shift σ τ σ'))

/-- The partition function of the layered system, with normalized counting measure. -/
noncomputable def layerZ (shift : D → Λ → Λ) (Es : (Λ → D → ℤˣ) → ℝ) (β : ℝ) (t : ℕ) : ℝ :=
  ((layerCard Λ D) ^ t)⁻¹ * ∑ σ : Fin t → Λ → D → ℤˣ, ∑ τ : Fin t → Λ → ℤˣ,
    exp (-(β * (∑ s, Es (σ s) + ∑ s, tempE shift (σ s) (τ s) (σ (finRotate t s)))))

theorem layerZ_eq_trace (shift : D → Λ → Λ) (Es : (Λ → D → ℤˣ) → ℝ) (β : ℝ) {t : ℕ} (ht : 1 ≤ t) :
    layerZ shift Es β t = (transfer shift Es β ^ t).trace := by
  rw [← closed_walk_sum _ ht, layerZ, Finset.mul_sum]
  refine Finset.sum_congr rfl fun σ _ => ?_
  have hτ : ∑ τ : Fin t → Λ → ℤˣ, ∏ s, exp (-(β * tempE shift (σ s) (τ s) (σ (finRotate t s)))) =
      ∏ s, ∑ τ : Λ → ℤˣ, exp (-(β * tempE shift (σ s) τ (σ (finRotate t s)))) := by
    rw [Finset.prod_univ_sum]
    simp [Fintype.piFinset_univ]
  have hrot : ∏ s, exp (-(β * Es (σ (finRotate t s))) / 2) = ∏ s, exp (-(β * Es (σ s)) / 2) :=
    Equiv.prod_comp (finRotate t) (fun s => exp (-(β * Es (σ s)) / 2))
  simp only [transfer, Finset.prod_mul_distrib, Finset.prod_const, Finset.card_univ, Fintype.card_fin]
  rw [hrot, ← hτ, inv_pow, mul_assoc]
  congr 1
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun τ _ => ?_
  rw [← exp_sum, ← exp_sum, ← exp_add, ← exp_add]
  congr 1
  rw [← Finset.sum_div, Finset.sum_neg_distrib, Finset.sum_neg_distrib, ← Finset.mul_sum, ← Finset.mul_sum]
  ring

/-! ### Positivity -/

/-- The character of a set `A` of links. -/
def chi (A : Finset (Λ × D)) (σ : Λ → D → ℤˣ) : ℝ := ∏ l ∈ A, sgn (σ l.1 l.2)

/-- `e^{βu} = cosh β + u sinh β` for `u = ±1`. -/
theorem exp_sgn (β : ℝ) (u : ℤˣ) : exp (β * sgn u) = sinh β * sgn u + cosh β := by
  rcases sgn_cases u with h | h
  · rw [h, mul_one, mul_one, add_comm, cosh_add_sinh]
  · rw [h, mul_neg_one, mul_neg_one, neg_add_eq_sub, cosh_sub_sinh]

/-- The sign a gauge transformation `τ` contributes to `A`. -/
def psi (shift : D → Λ → Λ) (A : Finset (Λ × D)) (τ : Λ → ℤˣ) : ℤˣ :=
  ∏ l ∈ A, (τ (shift l.2 l.1) * τ l.1)

theorem psi_mul (shift : D → Λ → Λ) (A : Finset (Λ × D)) (τ τ' : Λ → ℤˣ) :
    psi shift A (τ * τ') = psi shift A τ * psi shift A τ' := by
  unfold psi
  rw [← Finset.prod_mul_distrib]
  refine Finset.prod_congr rfl fun l _ => ?_
  simp only [Pi.mul_apply]
  exact mul_mul_mul_comm _ _ _ _

/-- The sum of a sign character over a finite group is nonnegative. -/
theorem sum_psi_nonneg (shift : D → Λ → Λ) (A : Finset (Λ × D)) :
    0 ≤ ∑ τ : Λ → ℤˣ, sgn (psi shift A τ) := by
  by_cases h : ∀ τ, sgn (psi shift A τ) = 1
  · exact Finset.sum_nonneg fun τ _ => by rw [h τ]; exact zero_le_one
  · push Not at h
    obtain ⟨τ₀, h0⟩ := h
    have hm : sgn (psi shift A τ₀) = -1 := (sgn_cases _).resolve_left h0
    have hs : ∑ τ : Λ → ℤˣ, sgn (psi shift A (τ₀ * τ)) = ∑ τ : Λ → ℤˣ, sgn (psi shift A τ) :=
      Equiv.sum_comp (Equiv.mulLeft τ₀) (fun τ => sgn (psi shift A τ))
    simp only [psi_mul, sgn_mul, hm, ← Finset.mul_sum] at hs
    linarith

/-- The coefficient of the character `A` in the transfer matrix. -/
noncomputable def coeff (shift : D → Λ → Λ) (β : ℝ) (A : Finset (Λ × D)) : ℝ :=
  (layerCard Λ D)⁻¹ * exp (-(β * Fintype.card (Λ × D))) *
    (sinh β ^ A.card * cosh β ^ (Finset.univ \ A).card) * ∑ τ : Λ → ℤˣ, sgn (psi shift A τ)

theorem coeff_nonneg (shift : D → Λ → Λ) {β : ℝ} (hβ : 0 ≤ β) (A : Finset (Λ × D)) :
    0 ≤ coeff shift β A := by
  unfold coeff
  have h1 := (layerCard_pos (Λ := Λ) (D := D)).le
  have h2 := sinh_nonneg_iff.mpr hβ
  have h3 := (cosh_pos β).le
  have h4 := sum_psi_nonneg shift A
  positivity

theorem exp_tempE (shift : D → Λ → Λ) (β : ℝ) (σ σ' : Λ → D → ℤˣ) (τ : Λ → ℤˣ) :
    exp (-(β * tempE shift σ τ σ')) = exp (-(β * Fintype.card (Λ × D))) *
      ∑ A : Finset (Λ × D), (sinh β ^ A.card * cosh β ^ (Finset.univ \ A).card) *
        (chi A σ * chi A σ' * sgn (psi shift A τ)) := by
  have hE : -(β * tempE shift σ τ σ') = -(β * Fintype.card (Λ × D)) +
      ∑ l : Λ × D, β * sgn (σ l.1 l.2 * τ (shift l.2 l.1) * σ' l.1 l.2 * τ l.1) := by
    unfold tempE
    rw [← Fintype.sum_prod_type', Finset.sum_sub_distrib, ← Finset.mul_sum]
    simp only [Finset.sum_const, Finset.card_univ, nsmul_eq_mul, mul_one]
    ring
  rw [hE, exp_add, exp_sum]
  congr 1
  rw [Finset.prod_congr rfl fun l _ => exp_sgn β _, Finset.prod_add, Finset.powerset_univ]
  refine Finset.sum_congr rfl fun A _ => ?_
  have hc : ∏ l ∈ A, sgn (σ l.1 l.2 * τ (shift l.2 l.1) * σ' l.1 l.2 * τ l.1) =
      chi A σ * chi A σ' * sgn (psi shift A τ) := by
    unfold chi psi
    rw [sgn_prod, ← Finset.prod_mul_distrib, ← Finset.prod_mul_distrib]
    refine Finset.prod_congr rfl fun l _ => ?_
    simp only [sgn_mul]
    ring
  rw [Finset.prod_mul_distrib, Finset.prod_const, Finset.prod_const, hc]
  ring

theorem transfer_eq_sum (shift : D → Λ → Λ) (Es : (Λ → D → ℤˣ) → ℝ) (β : ℝ) :
    transfer shift Es β = ∑ A : Finset (Λ × D), coeff shift β A •
      Matrix.vecMulVec (fun σ => exp (-(β * Es σ) / 2) * chi A σ)
        (star fun σ => exp (-(β * Es σ) / 2) * chi A σ) := by
  ext σ σ'
  simp only [transfer, exp_tempE, Matrix.sum_apply, Matrix.smul_apply, Matrix.vecMulVec_apply,
    smul_eq_mul, star_trivial, coeff, Finset.mul_sum, Finset.sum_mul]
  rw [Finset.sum_comm]
  refine Finset.sum_congr rfl fun A _ => Finset.sum_congr rfl fun τ _ => ?_
  ring

theorem transfer_posSemidef (shift : D → Λ → Λ) (Es : (Λ → D → ℤˣ) → ℝ) {β : ℝ} (hβ : 0 ≤ β) :
    (transfer shift Es β).PosSemidef := by
  rw [transfer_eq_sum]
  exact Matrix.posSemidef_sum _ fun A _ =>
    (Matrix.posSemidef_vecMulVec_self_star _).smul (coeff_nonneg shift hβ A)

theorem transfer_trace_pos (shift : D → Λ → Λ) (Es : (Λ → D → ℤˣ) → ℝ) (β : ℝ) :
    0 < (transfer shift Es β).trace := by
  refine Finset.sum_pos (fun σ _ => ?_) Finset.univ_nonempty
  simp only [Matrix.diag_apply, transfer]
  have h1 := layerCard_pos (Λ := Λ) (D := D)
  have h2 : 0 < ∑ τ : Λ → ℤˣ, exp (-(β * tempE shift σ τ σ)) :=
    Finset.sum_pos (fun τ _ => exp_pos _) Finset.univ_nonempty
  positivity

/-- The sorted eigenvalues of the transfer matrix. -/
noncomputable def layerSpec (shift : D → Λ → Λ) (Es : (Λ → D → ℤˣ) → ℝ) {β : ℝ} (hβ : 0 ≤ β) : ℕ → ℝ :=
  specSeq (transfer_posSemidef shift Es hβ).1

theorem layerSpec_spectrum (shift : D → Λ → Λ) (Es : (Λ → D → ℤˣ) → ℝ) {β : ℝ} (hβ : 0 ≤ β) :
    Spectrum (layerSpec shift Es hβ) :=
  spectrum_specSeq (transfer_posSemidef shift Es hβ) (transfer_trace_pos shift Es β)

/-- **Trace formula.** -/
theorem layerZ_eq_partition (shift : D → Λ → Λ) (Es : (Λ → D → ℤˣ) → ℝ) {β : ℝ} (hβ : 0 ≤ β)
    {t : ℕ} (ht : 1 ≤ t) : layerZ shift Es β t = partition (layerSpec shift Es hβ) t := by
  rw [layerZ_eq_trace shift Es β ht, layerSpec, partition_specSeq _ ht]

end Layer
end YMFormal
