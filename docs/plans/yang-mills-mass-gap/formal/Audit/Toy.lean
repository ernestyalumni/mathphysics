import Mathlib
import Audit.Export

/-! Development statements for the statement checker; not part of the benchmark. -/

namespace Toy

noncomputable def sq (x : ℝ) : ℝ := x ^ 2

theorem toy_true (a b : ℝ) (ha : 0 ≤ a) (hb : a ≤ b) : sq a ≤ sq b := by sorry
theorem toy_false (a b : ℝ) (hb : a ≤ b) : sq a ≤ sq b := by sorry
theorem toy_opt (a : ℝ) (s : ℝ := a ^ 2) (h : 0 ≤ a) : s = sq a := by sorry
theorem toy_sum (f : ℕ → ℝ) (h : ∀ j, 0 ≤ f j) (hs : Summable f) : 0 ≤ ∑' j, f j := by sorry
theorem toy_log (x : ℝ) (hx : 1 < x) : 0 < Real.log x ∧ Real.exp (Real.log x) = x := by sorry
theorem toy_let (n : ℕ) (hn : 1 ≤ n) : let m := (n : ℝ) / 2; m * 2 = n := by sorry

end Toy

#ym_export Toy.toy_true
#ym_export Toy.toy_false
#ym_export Toy.toy_opt
#ym_export Toy.toy_sum
#ym_export Toy.toy_log
#ym_export Toy.toy_let
