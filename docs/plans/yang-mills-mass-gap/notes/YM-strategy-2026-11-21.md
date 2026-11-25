# Strategy summary, 2026-11-21 (strategist, after reviews YM-O03 and YM-O04)

Label: reviewed by Codex (YM-O04, `notes/YM-O04-codex-review-2026-11-21.md`); labels below follow that review.

## Status of the conditional-gap chain

- **Proved, and in Lean** (F01: `ConditionalGap.lean`, `RankOne.lean`, `Rectangle.lean` da4d739):
  - purity bounds;
  - the defect telescope;
  - the rectangle recurrence Δ(2n) ≤ 3740 Δ(n)² for Δ(n) ≤ 1/7480 (derived by Codex in YM-O03 §1, checked by both Claude sessions);
  - the conditional gap, in ratio form (λ₂/λ₁)ⁿ ≤ 1/2 at every dyadic scale.
- **Inputs outside Lean:**
  - the Wilson rectangular trace formula and positivity of the transfer operators in every direction (C23);
  - the passage from dyadic finite-volume spectral bounds to the infinite-volume state.
- **The open estimate:** a seed of bounded physical size, uniform toward the continuum (precise form in YM-N03, "The precise seed target").

## Next proof work (YM-O04 §6)

1. **First end-to-end instance.** Formalize the elementary bound Δ_β(n) ≤ 192 β n⁴, which gives a certified seed (at n = 2, any β ≤ 1/(7480·3072)). Add the rectangular trace/positivity interface and the local-limit argument. This needs no cluster expansion; a convergent cluster expansion is a sharper, optional target.
2. **The uniform seed.** Compare unnormalized torus free energies after RG, with a common normalization and coupling tuning. Bałaban-type UV stability is not an IR decay estimate. Check whether the OpenAI construction's tuning is shared across volumes before borrowing any normalization.

## Labels for numerical results

- β ≈ 0.097 and 1.355 are crossings of the truncated two-term surface formula, not certified seed thresholds.
- Glueball am₀₊₊:
  - β = 2.3 gives 1.198(18), an early-time check;
  - β = 2.4 gives 0.925(33) at t = 2, a compromise with downward drift;
  - β = 2.5 gives 0.66(4), the last-time-point effective mass, not an established plateau.
  - None is a fitted plateau with systematics.
- Creutz ratios and the specific-heat maximum are consistency checks, not determinations.
- The vortex Z₂ count (t₂ − t_p) is a descriptive association with confinement, not an order parameter or a lemma.

## Method

- Stop paying for unverified verdicts from local models. Keep bounded candidate generation where Lean, exact arithmetic or a checked counterexample settles the output.
- Core method:
  - challenges from independent model families;
  - primary-source checks;
  - explicit counter-models (U(1); the θ = π anomaly, which excludes only a unique symmetric trivially gapped vacuum);
  - Lean on frozen statements, auditing translation and hypotheses separately.
- No Monte Carlo as proof of inequalities. Keep exact small-box evaluation, certified enclosures, and pre-registered falsification runs.
- Topology and YM-N02 are low priority. At most one matched-box diagnostic, with a stated inequality and failure criterion.
