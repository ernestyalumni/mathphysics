# Is this "half of a Millennium Prize"?

Short answer: **the mathematics splits roughly in two, the prize does not.** The mass gap is unproved, as you thought. But there is no half prize for existence, and the existence result itself is not yet the existence the problem asks for.

## 1. The problem is one statement

Jaffe and Witten's official formulation, as quoted in the group's assessment §1 and Armstrong's `statement.md` §5.3 (label: **source-claims**; the original text is task YM-W02 to read and quote exactly):

> Prove that for any compact simple gauge group G, a non-trivial quantum Yang–Mills theory exists on ℝ⁴ and has a mass gap Δ > 0.

The assessment breaks this into four requirements:

1. a nontrivial theory satisfying axioms at least as strong as Wightman or Osterwalder–Schrader, with genuine local observables;
2. local fields corresponding to gauge-invariant polynomials in the curvature and its covariant derivatives;
3. the short-distance behaviour predicted by asymptotic freedom and perturbative renormalization: a stress tensor, an operator product expansion with the predicted singularities;
4. a mass gap: the spectrum of the Hamiltonian lies in {0} ∪ [Δ, ∞) with Δ > 0, and the supremum of such Δ is finite.

One prize, awarded for all of it.

## 2. Scorecard for the OpenAI result

| Requirement | Status | Evidence |
|---|---|---|
| Any compact simple G | Done, and slightly more general | `statement.md` Def. 1.1 |
| Wightman axioms | Done for each phase of a **subsequential** infinite-volume limit | assessment Thm 6.1 |
| Nontrivial | Done in a specific sense: finite Wick profiles excluded | assessment §7.3 |
| Yang–Mills local fields | Done inside the construction; weaker in the public theorem | assessment §3, §9.4 |
| Stress tensor, Ward identity | Done | assessment §7.1 |
| Asymptotic freedom | Leading logarithm only | assessment Thm 7.2 |
| OPE "with the predicted singularities" | Bounds proved; coefficients **not** matched to perturbation theory | assessment §7.5 |
| Unique theory | **Not done.** Subsequence-dependent; phases not shown equal | assessment §5, §6.2 |
| Mass gap | **Not done.** "No evidence for or against" | assessment §2 item 4 |

The assessment adds a pointed remark: Jaffe–Witten "regards a subsequential limit obtained via compactness as insufficient, unless additional arguments establish the requisite properties of the limit". The group argues the extra structure (fields, stress tensor, OPE) goes well beyond bare compactness, and adopts an existence criterion that does not require characterizing the vacuum. That is their judgment, stated as such, and they say they have "no present connection to the Clay Mathematics Institute".

So "half" is fair as a description of the mathematics: existence in a meaningful sense is there, the gap is not. Two cautions on the word:

- **The halves are not equal.** Everything proved lives at short distance, where the coupling is small. The gap is a long-distance, strong-coupling statement. Nothing in the proof touches that regime; infinite volume is reached by a compactness argument that carries no dynamical information.
- **Uniqueness sits between the halves.** A mass gap proved for one phase of one subsequential limit would still leave "the" theory undefined. In OpenAI's two-dimensional O(3) paper the gap and the uniqueness come from the same infrared input. Expect the same here: the gap estimate should also be what removes the subsequence.

## 3. What the prize rules require

From the Clay Mathematics Institute rules page (web search 2026-11-14; label **source-claims**, exact text is task YM-W02):

- the solution must be published in a **qualifying outlet** (a refereed mathematics publication of worldwide repute, or another form the institute approves);
- at least **two years** must pass after publication;
- it must have **general acceptance** in the mathematics community;
- the institute decides at its sole discretion.

Not verified and worth reading in the primary text before making plans: how the rules treat a solution that builds on someone else's earlier work (older versions of the rules allowed dividing a prize and recognising prior insights), and whether they say anything about machine-produced proofs. Until YM-W02 is done, treat both as unknown.

Practical consequences for us:

- A prize-grade contribution is a mass-gap theorem that is (a) for a theory accepted as *the* Yang–Mills theory, (b) kernel-checked or written to referee standard, (c) published, (d) still standing two years later.
- A mass gap proved on top of OpenAI's construction would be a joint achievement in the eyes of any committee. Credit is not ours to assign; the mathematics is still worth doing.
- As of 2026-11-14 the search found no peer-reviewed publication and no Clay statement about the OpenAI result.

## 4. Who else is on this

- **OpenAI.** Their model has already done the two-dimensional analogue end to end (family 215 in `openai/math`): exponential decay for O(n≥3) at all temperatures, a positive gap for O(4), and a unique massive continuum limit for O(3). The review group lists a second Yang–Mills submission dated 4 November as "work in progress". Assume they are attempting the gap.
- **The review group** (led by David Gross, per press coverage): validating, not extending.
- **Hairer's project and Armstrong's map:** turning the Lean into readable mathematics. Open to outside agents.
- **The rest of the constructive field theory community**, now with a working ultraviolet construction to build on.

An honest reading of the odds: a solo researcher with a fleet of agents is unlikely to get to the full gap before a lab that produced the existence proof and the two-dimensional template. That is not a reason to stay out. It is a reason to pick sub-problems where an independent result has value on its own: see the milestone list in [03-MASS-GAP-STRATEGY.md](03-MASS-GAP-STRATEGY.md) §8, every one of which is useful even if someone else closes the main problem first.

## Sources

- Yang–Mills group, *Validation of the OpenAI Yang–Mills construction*, 10 November 2026: `yang-mills-group/YM/docs/assessment-1.pdf` at `55ba9fc`.
- S. Armstrong, `statement.md` §4–§5 and `faq.md`, `scottnarmstrong/yang-mills-informalization` at `daf1db4`.
- Clay Mathematics Institute rules: https://claymath.org/millennium-problems/rules and https://www.claymath.org/wp-content/uploads/2022/03/millennium_prize_rules_0.pdf
- Press coverage found by search (not authoritative): https://cryptobriefing.com/yang-mills-group-publish-openai-lean-code/ , https://fourweekmba.com/ai-gross-panel-sets-release-of-openais-yang-mills-lean-proof/ , https://officechai.com/ai/what-is-the-yang-mills-theory-openai-claims-to-have-proven-and-what-it-means-for-physics/
