import YMFormal.Doubling

/-!
# The spectrum of a finite nonnegative transfer matrix

For a real positive semidefinite matrix `A` with `Tr A > 0`, the eigenvalues in decreasing order,
padded with zeros, form a `Spectrum`, and `Z(t) = Tr Aᵗ` is its partition function for `t ≥ 1`.
`specSeq hA 0` is the largest eigenvalue and `specSeq hA 1` the second, counted with multiplicity,
so `specSeq hA 1 / specSeq hA 0` is `λ₂/λ₁` of the handbook.
-/

namespace YMFormal

open Matrix

variable {n : Type*} [Fintype n] [DecidableEq n]

/-- `Tr Aᵗ = ∑ λᵢᵗ` for a real symmetric matrix. -/
theorem trace_pow_eq_sum_eigenvalues {A : Matrix n n ℝ} (hA : A.IsHermitian) (t : ℕ) :
    (A ^ t).trace = ∑ i, hA.eigenvalues i ^ t := by
  conv_lhs => rw [hA.spectral_theorem]
  rw [← map_pow, diagonal_pow, Unitary.conjStarAlgAut_apply, trace_mul_comm, ← mul_assoc]
  have hu : star (hA.eigenvectorUnitary : Matrix n n ℝ) * (hA.eigenvectorUnitary : Matrix n n ℝ) = 1 :=
    Unitary.coe_star_mul_self _
  rw [hu, one_mul, trace_diagonal]
  simp

/-- The eigenvalues in decreasing order, padded with zeros. -/
noncomputable def specSeq {A : Matrix n n ℝ} (hA : A.IsHermitian) : ℕ → ℝ := fun j =>
  if h : j < Fintype.card n then hA.eigenvalues₀ ⟨j, h⟩ else 0

theorem sum_eigenvalues_eq {A : Matrix n n ℝ} (hA : A.IsHermitian) (t : ℕ) :
    ∑ i, hA.eigenvalues i ^ t = ∑ k : Fin (Fintype.card n), hA.eigenvalues₀ k ^ t := by
  unfold IsHermitian.eigenvalues
  exact Equiv.sum_comp (Fintype.equivOfCardEq (Fintype.card_fin _)).symm (fun k => hA.eigenvalues₀ k ^ t)

theorem eigenvalues₀_nonneg {A : Matrix n n ℝ} (hA : A.PosSemidef) (k : Fin (Fintype.card n)) :
    0 ≤ hA.1.eigenvalues₀ k := by
  have := hA.eigenvalues_nonneg ((Fintype.equivOfCardEq (Fintype.card_fin _)) k)
  unfold IsHermitian.eigenvalues at this
  simpa using this

theorem partition_specSeq {A : Matrix n n ℝ} (hA : A.IsHermitian) {t : ℕ} (ht : 1 ≤ t) :
    partition (specSeq hA) t = (A ^ t).trace := by
  rw [trace_pow_eq_sum_eigenvalues hA, sum_eigenvalues_eq, partition]
  rw [tsum_eq_sum (s := Finset.range (Fintype.card n))]
  · rw [← Fin.sum_univ_eq_sum_range]
    refine Finset.sum_congr rfl fun k _ => ?_
    simp [specSeq, k.isLt]
  · intro j hj
    simp only [Finset.mem_range, not_lt] at hj
    simp [specSeq, Nat.not_lt.mpr hj, zero_pow (by omega : t ≠ 0)]

/-- A positive semidefinite matrix with positive trace has a `Spectrum`. -/
theorem spectrum_specSeq {A : Matrix n n ℝ} (hA : A.PosSemidef) (htr : 0 < A.trace) :
    Spectrum (specSeq hA.1) := by
  have hcard : 0 < Fintype.card n := by
    by_contra h0
    push Not at h0
    have : IsEmpty n := Fintype.card_eq_zero_iff.mp (by omega)
    simp [Matrix.trace] at htr
  have hanti := hA.1.eigenvalues₀_antitone
  have htop : ∀ k : Fin (Fintype.card n), hA.1.eigenvalues₀ k ≤ hA.1.eigenvalues₀ ⟨0, hcard⟩ :=
    fun k => hanti (Fin.le_iff_val_le_val.mpr (Nat.zero_le _))
  refine ⟨?_, ?_, ?_, ?_⟩
  · -- the largest eigenvalue is positive, since the trace is
    have hsum : A.trace = ∑ k : Fin (Fintype.card n), hA.1.eigenvalues₀ k := by
      have := trace_pow_eq_sum_eigenvalues hA.1 1
      simp only [pow_one] at this
      rw [this]
      simpa using sum_eigenvalues_eq hA.1 1
    by_contra hle
    push Not at hle
    have : ∑ k : Fin (Fintype.card n), hA.1.eigenvalues₀ k ≤ 0 :=
      Finset.sum_nonpos fun k _ => (htop k).trans (by simpa [specSeq, hcard] using hle)
    linarith
  · intro j
    unfold specSeq
    split_ifs with h
    · exact eigenvalues₀_nonneg hA _
    · exact le_refl 0
  · intro j
    have h0 : specSeq hA.1 0 = hA.1.eigenvalues₀ ⟨0, hcard⟩ := by simp [specSeq, hcard]
    rw [h0]
    unfold specSeq
    split_ifs with h
    · exact htop _
    · exact eigenvalues₀_nonneg hA _
  · refine summable_of_ne_finset_zero (s := Finset.range (Fintype.card n)) fun j hj => ?_
    simp only [Finset.mem_range, not_lt] at hj
    simp [specSeq, Nat.not_lt.mpr hj]

end YMFormal
