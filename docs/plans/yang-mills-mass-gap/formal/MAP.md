# Handbook ↔ Lean map

This file maps the Yang–Mills handbook to the Lean in this directory: `documents/mono-notes/topics/yang-mills-lattice.tex`, by equation and theorem label. The handbook does not name code, so the map lives here.

**Running the check.** All Lean is checked by `formal/check.sh --strict`, in Docker with no network. It allows no `sorry`, and only the axioms `propext`, `Classical.choice` and `Quot.sound`.

## Main results

| Handbook label | Lean (`YMFormal.…`) | What the Lean assumes beyond the handbook statement |
|---|---|---|
| `ym:purity-bounds` (Prop. `ym:prop:purity`) | `purity_bounds` | Spectrum `lam` summable, `0 ≤ lam j ≤ lam 0`, `0 < lam 0`, tail not ordered |
| `ym:gap-from-purity` (ratio form) | `ratio_pow_le_impurity` | none |
| `ym:gap-from-purity-log` | `exp_neg_gap_le`, `gap_from_purity` | `0 < lam 1` |
| `ym:defect-telescope` | `defect_telescope`, `defect_nonneg_and_last` | the four ratios `𝖯_k ≤ 1` as hypotheses |
| `ym:gap-from-defect` (ratio form) | `ratio_pow_le_defect` | `exp (-Δ) ≤ purity lam n` as input |
| `ym:gap-from-defect-log` | `gap_from_defect` | `0 < lam 1`, `0 < Δ` |
| rank one, after `ym:rank-one` | `rank_one_of_zero_defect` | — |
| `ym:doubling-cost`, `ym:cost-nonneg` | `cost`, `sideCost`; hypothesis `nonneg` | nonnegativity is a hypothesis; the handbook derives it from Prop. `ym:sc-prop:positivity` |
| Lemma `ym:lem:one-doubling` | hypothesis `estimate` | a hypothesis; the handbook derives it from `ym:purity-bounds` |
| Lemma `ym:lem:commuting` | `cost_swap` | — |
| `ym:defect-path`, `ym:cost-cube` | `pot_drop`, `pot_down_from_zero`, `pot_up_to_one`, `cost_le_defect_binary` | proved via the potential `φ(e) = 16 G(e)/2^{Σe}` rather than by choosing a path |
| `ym:cost-induction` | `cost_le_pow_three` | — |
| Theorem `ym:thm:recurrence` | `rectangle_recurrence` (exponents), `box_recurrence` (side lengths) | `nonneg`, `estimate` on boxes with positive sides |
| Prop. `ym:prop:conditional-gap` (ratio form) | `conditional_gap_ratio`, `conditional_gap_of_boxes` | the ratio bound `from_defect` at every scale as input |
| Prop. `ym:prop:conditional-gap`, log form | `conditional_gap` | `0 < δ_k` |
| `ym:log-z-bounds` | `log_partition_bound_wilson` (and the symmetric `log_partition_bound`) | measurable plaquette terms, probability measure |
| Prop. `ym:prop:strong-seed` | `seed_of_wilson_bound`, `seed_strong_coupling`, `conditional_gap_strong_coupling` | the plaquette bound on `log Z(x)` for every box |
| `ym:finite-correlation-bound` and the infinite-volume hypothesis | structure `LocalLimit`, theorem `LocalLimit.decay` | `finite_bound` and `tendsto` are fields, i.e. assumptions |
| Prop. `ym:prop:casimir-sun` | `trace_sandwich`, `trace_adjoint_comp`, `casimir_of_killing`, `slope_values` | matrix identities over a commutative ring; restriction to su(N) and complexification are not formalised |

## Not formalised

- **The trace formula** `ym:trace-formula`, and the identification of each `𝖯_k` as a purity. In the Lean these are the inputs `nonneg`, `estimate` and `from_defect`.
- **Positivity of the transfer operator,** Prop. `ym:sc-prop:positivity`.
- **Conjecture `ym:conj:defect`.**
- **The steps after `ym:physical-gap`.**

For the Z₂ gauge group all of the finite-volume steps are formalised, with no hypotheses: the trace formula in every direction (`z2Z_eq_partition`), positivity (`transfer_posSemidef`), the two doubling hypotheses (`doubling_cost_nonneg`, `doubling_cost_estimate`), the seed (`z2F_bounds`) and the conclusion (`z2_gap`). See `notes/YM-F02-z2-end-to-end.md`.
