import Mathlib

/-!
# Casimir constant of `su(N)` (`ym:prop:casimir-sun`) and the slope `b_*`

Source: `documents/mono-notes/topics/yang-mills-lattice.tex`, Proposition `ym:prop:casimir-sun`,
Eq. `ym:slope` and the values after it.

Formalised:
- Eq. `ym:trace-AZB`: the trace of `Z ↦ A Z B` on `M_N(R)` is `tr A · tr B`.
- Eq. `ym:killing-gl`: `tr(ad X ∘ ad Y) = 2N tr(XY) - 2 tr X tr Y` on `M_N(R)`.
- The solving step for `C_𝔤`, and the values of `b_* = 11 C_𝔤 / (24π²)` for `SU(2)` and `SU(3)`.

Not formalised: that the real trace over `su(N)` equals the complex trace over `sl_N(ℂ)`
(complexification), and the reduction from `M_N(ℂ)` to `sl_N(ℂ)`.
-/

open Matrix

namespace YMFormal

variable {R : Type*} [CommRing R] {N : ℕ}

/-- The linear map `Z ↦ A Z B` on `N × N` matrices. -/
def sandwich (A B : Matrix (Fin N) (Fin N) R) : Matrix (Fin N) (Fin N) R →ₗ[R] Matrix (Fin N) (Fin N) R :=
  LinearMap.mulLeft R A ∘ₗ LinearMap.mulRight R B

/-- `ad X = [X, ·]` on `N × N` matrices. -/
def adjoint (X : Matrix (Fin N) (Fin N) R) : Matrix (Fin N) (Fin N) R →ₗ[R] Matrix (Fin N) (Fin N) R :=
  LinearMap.mulLeft R X - LinearMap.mulRight R X

/-- Eq. `ym:trace-AZB`. -/
theorem trace_sandwich (A B : Matrix (Fin N) (Fin N) R) :
    LinearMap.trace R _ (sandwich A B) = trace A * trace B := by
  sorry

/-- Eq. `ym:killing-gl`. -/
theorem trace_adjoint_comp (X Y : Matrix (Fin N) (Fin N) R) :
    LinearMap.trace R _ (adjoint X ∘ₗ adjoint Y) =
      2 * N * trace (X * Y) - 2 * trace X * trace Y := by
  sorry

/-- The last display of the proof of `ym:prop:casimir-sun`: if `K = 2N t` and
`K = -C (-(1/N) t)` for some `t ≠ 0`, then `C = 2N²`. -/
theorem casimir_of_killing {C t : ℝ} (hN : 0 < N) (ht : t ≠ 0)
    (h : 2 * N * t = -C * (-(1 / N) * t)) : C = 2 * N ^ 2 := by
  sorry

/-- `b_* = 11 C_𝔤 / (24π²)` with `C_𝔤 = 2N²`: `11/(3π²)` for `SU(2)`, `33/(4π²)` for `SU(3)`. -/
theorem slope_values :
    11 * (2 * 2 ^ 2 : ℝ) / (24 * Real.pi ^ 2) = 11 / (3 * Real.pi ^ 2) ∧
    11 * (2 * 3 ^ 2 : ℝ) / (24 * Real.pi ^ 2) = 33 / (4 * Real.pi ^ 2) := by
  sorry

end YMFormal
