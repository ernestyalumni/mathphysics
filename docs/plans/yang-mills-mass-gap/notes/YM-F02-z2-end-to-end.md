# YM-F02: Z₂ lattice gauge theory end to end in Lean

Author: Claude Opus 5.5, 2026-11-25 UTC. Branch `ym/YM-F02-claude-opus55-z2-end-to-end`. Request: Ernest, 2026-11-22, "go ahead with C and then B". Label: machine-checked; reviewed and accepted by Codex at `4a29e20` (`notes/YM-F02-codex-review-2026-11-25.md`, commit `6afbb49` on `ym/YM-F02-codex-z2`), no change to the proof requested.

## Result

`YMFormal.z2_gap`: for `0 ≤ β ≤ 1/(7480·192·n⁴)`, `n ≥ 1`, every direction `i` and every `k`, the transfer matrix in direction `i` on the transverse cube of side `2ᵏ n` satisfies `(λ₂/λ₁)ⁿ ≤ 1/2`. `z2_gap_log` is the form `log 2 / n ≤ -log(λ₂/λ₁)` when `λ₂ > 0`.

There is no hypothesis beyond the range of β: the trace formula in every direction, positivity, the recurrence and the seed are all proved. `formal/check.sh --strict` passes (107 theorems, no `sorry`, axioms propext, Classical.choice and Quot.sound; the count of 105 in commit `4a29e20` missed two with attributes).

**What it is not.** Z₂ is a toy gauge group, and the range of β is extremely strong coupling, where a gap has long been known. It is not a statement about SU(N), the continuum or the infinite-volume state. Its value is that every interface of the chain has been used once with nothing assumed.

## How it is proved

| Step | Where |
|---|---|
| `Z` defined from the lattice only: sites, links indexed by (site, direction), six planes, normalized counting measure, weight `e^{-β Σ(1-U_p)}` | `Z2.lean`: `z2Z`, `energy`, `plaq` |
| Sanity: `Z = 1` on the 1⁴ box, as in YM-T10 | `z2Z_one` |
| Slicing a configuration along any direction `i`, including sides of length 1 | `siteEquiv`, `configEquiv`, `energy_slice` |
| Sum over closed walks = trace of a power; trace formula with an arbitrary slice energy | `Layer.lean`: `closed_walk_sum`, `layerZ_eq_trace` |
| Positivity: `e^{βu} = cosh β + u sinh β`; the sum over temporal links is a character sum, which is 0 or the group order. The slice energy need not be gauge invariant | `transfer_posSemidef` |
| Sorted spectrum, `Tr Tᵗ = Σ λⱼᵗ` | `MatrixSpectrum.lean` |
| Doubling costs nonnegative; one-doubling estimate | `Doubling.lean` |
| Seed `Δ ≤ 192βn⁴` | `z2F_bounds`, `Seed.lean` |
| `(λ₂/λ₁)ᴺ ≤ e^{Δ(N)} - 1` on the transverse cube of side `N`, by the first doubled direction (YM-O05) | `z2_from_defect` |

The spatial side is `2ᵏ n`. The handbook's `ym:prop:conditional-gap` uses the last doubled direction and side `2ᵏ⁺¹ n`.

## Review outcome (Codex, 2026-11-25)

- **The lattice and `Z` are the intended ones,** including sides of length 1 and time extents 1 and 2. Independent enumeration in exact arithmetic agrees with the trace in all four directions on 1⁴, 2·1³ and 2²·1².
- **λ₂.** For this energy `T = PTP` with `P` the gauge average, so functions orthogonal to the gauge-invariant ones contribute only zeros. That is a paper observation, not formalised. With transverse side 1 the ratio is `tanh β`.
- **Not vacuous.** Positive β exist for every n, and entry 1 of the spectrum is a genuine eigenvalue.
- **Checker defect, fixed.** The strict check's axiom banner only tested for `sorryAx`. `check.sh --strict` now audits every theorem of the modules against the allowlist and fails if any is missing; a planted axiom and `native_decide` both fail it.
