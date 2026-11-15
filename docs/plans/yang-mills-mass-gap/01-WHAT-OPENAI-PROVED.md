# What OpenAI proved, and how

Label for this whole file: **source-claims**, except lines marked otherwise. Written 2026-11-14 from the pinned sources listed in [README.md](README.md). The reader of record for wording is Armstrong's `statement.md`; the authority for truth is the Lean.

## 1. The theorem

Formal name: `YangMills.jointExistence` (in `Solution.lean`), proving `JointExistence` (in `ChallengeSpec.lean`). It is a projection of a stronger internal theorem `YangMills.constructed` (in `Construction.lean`).

Input: a closed subgroup G ⊂ U(m) whose real Lie algebra 𝔤 is simple and nonabelian. This covers SU(N), and also some disconnected groups.

Output: one family of Wilson lattice gauge theories on four-dimensional tori, with

- spacings a_N → 0 and physical side lengths L_p → ∞, both chosen by the construction;
- bare couplings β_{N,p} with β_{N,p} / log(1/a_N) → b_* > 0 uniformly in p;
- lattice "insertions" I_{N,p}(v, y): gauge-invariant functions of the links near site y, one for each field label v.

For every increasing sequence of volumes there is a subsequence and a continuum "hierarchy" such that:

| Clause | Content |
|---|---|
| P1 Lattice origin | Smeared lattice correlations of the insertions converge: first a_N → 0 at fixed volume along the **whole** spacing sequence, then volume → ∞ along the **subsequence**. Test functions vanish near coincident points. |
| P2 Phases and reconstruction | The limit S^mix is an average ∫ S^ξ dν(ξ) over a compact space X of "phases". One selected phase ξ_* carries a Wightman theory: Poincaré covariance, spectrum in the forward cone, vacuum unique up to scalars, tempered fields on a common domain, locality. |
| P3 Operator product expansion | One coefficient system gives finite-order expansions with remainder bounds r^{D+1−w−ε}, in the mixture and in every phase. |
| P4 Not free | A degree-4 field Θ has two-point function with r^8 log²(1/r) G_Θ(r) → 3 dim 𝔤 / π⁴ as r → 0, and G_Θ equals no finite sum of products of free mass-shell kernels near r = 0. |

The stronger theorem adds: the slope is the asymptotic-freedom value b_* = κ_G = 11 C_𝔤 / (24π²); the labels are gauge-invariant polynomials in curvature and covariant derivatives; Θ is the centred trace of a constructed conserved stress tensor; **every** phase reconstructs; correlations extend (non-uniquely) across coincident points; at short distance the rescaled theory converges to dim 𝔤 free Maxwell fields.

Normalization: ⟨X,Y⟩ = −(1/m) Re tr(XY), tr(ad X ad Y) = −C_𝔤 ⟨X,Y⟩, β = 1/g_0². For SU(N), C_𝔤 = 2N².

Checked by us on 2026-11-14 (label **reproduced**, script kept at `ledger/evidence/check_oneloop.py`):

- angular average over S³ of the comparison symbol h(k)|k|⁴ equals −(11/3) C_𝔤;
- ∫_{1/L<|k|≤1} |k|⁻⁴ d⁴k = 2π² log L;
- hence the per-step coefficient is κ_G log L with κ_G = 11 C_𝔤/(24π²);
- C_𝔤 = 8 for su(2) and 18 for su(3) in this normalization, so κ = 11N²/(12π²) and 1/g² ≃ (11N/24π²) log(1/a), the Gross–Wilczek–Politzer value.

## 2. What is not proved

These are stated as non-claims by the source README, the group's assessment and Armstrong's statement page.

1. **No mass gap.** Spectrum condition and unique vacuum only. The low-energy bounds are polynomial: spectral weight in (0,u] is O(u²), heat-kernel decay is O(s⁻²). Lean pointer: `YM/EveryPhase/Decay.lean:173–192`.
2. **No uniqueness of the infinite-volume limit.** Obtained by Arzelà–Ascoli. Different subsequences may give different limits.
3. **Phases not shown equal.** X is not shown to be a point. S^mix is not shown to cluster. No phase is shown to be a limit of lattice states.
4. **No universality.** The spacings, volumes, couplings and insertions are the construction's own. Nothing says another choice gives the same theory.
5. **Fixed order of limits.** No diagonal limit, no exchange, no rates.
6. **No measure on connections, no Wilson loops, no area law, no particles, no scattering.**
7. **Free fields only partly excluded.** Finite Wick profiles are excluded; continuous-mass generalized free fields are not.
8. **Perturbation theory not matched beyond leading order.** No two-loop term, no Λ parameter, OPE coefficients not compared.
9. The **public** theorem is weaker than the construction: OpenAI told the group that `jointExistence` "omitted some of the connections to the lattice construction". The group's own target `YMReview.Acceptance.combined_acceptance` closes that gap.

## 3. The proof in one page

The limit is taken in two stages that use completely different mathematics.

**Stage 1: ultraviolet, fixed volume (Parts I–III, about 70% of the Lean).** A Bałaban-type renormalization group runs from the fine lattice (level N) down to the torus at unit spacing (level 0). One block-spin step integrates the fine links against an exactly normalized gauge-covariant kernel. In the small-field region this is a Gaussian integral around a stationary background, organised by a cluster expansion. Large-field regions are kept and paid for by an energy lower bound. The one marginal term goes into the running coupling, which changes by κ_G log L_m per step up to a controlled error. The bare coupling is tuned so that the level-0 coupling is one fixed large number β̄ for every N and p. **The renormalization group stops at unit scale. Nothing is integrated out above it.**

Departures from Bałaban: a probabilistic block kernel with random tags, a spanning-forest gauge, block ratios that grow along the levels, a small-field threshold β^{−1/2} exp((log log β)²), an analytic family in a tuning parameter. The density format is closest to OpenAI's own O(4) paper (N2 §5).

The same steps transport observables and give, uniformly in the volume: chessboard bounds on every renormalized insertion, and convergence as a_N → 0.

**Stage 2: infrared, by compactness (Part IV, about 2%).** Reflection positivity (Osterwalder–Seiler) gives the chessboard inequality, hence local bounds uniform in the volume. A transfer-matrix representation plus a rank-free trace Hölder inequality gives equicontinuity. Arzelà–Ascoli gives a subsequence. No infrared analysis of the dynamics is done at all.

**Stage 3: structure of the limit (Parts V–VII, about 20%).**

- Ward identities from lattice integration by parts give a conserved stress tensor; comparison with free Maxwell theory fixes its normalization and integrates the infinitesimal rotations to O(4), without using any decay of correlations.
- The zero-energy subspace of the mixture's Osterwalder–Schrader space carries a commutative C*-algebra (commutativity uses O(4)); its characters are the phases.
- Each phase reconstructs to a Wightman theory. Vacuum uniqueness in every phase comes from the polynomial low-energy bound plus continuity in ξ and full support of ν, in place of a mass gap.
- The trace anomaly: β_k r_k⁴ Θ → −κ_G O_{F²} along scaling sequences, which with the free F² two-point function gives the 3 dim 𝔤/π⁴ law. The value of κ_G cancels out of this amplitude; it is the logarithm that proves non-freeness.
- The OPE comes from the same renormalization group.

**Stage 4: assembly (Part VIII).** Everything is shown to refer to one family and one limit table.

## 4. The interface that matters for the mass gap

The whole infrared stage consumes only three things from the ultraviolet stage (Armstrong README, "short, clean interface"):

1. a chessboard bound on every renormalized observable, uniform in the volume;
2. bounds on correlation kernels at separated points, uniform in the volume;
3. convergence of the kernels as the spacing tends to zero, at each fixed volume.

And the ultraviolet stage leaves behind, at unit scale, an explicit object: a lattice gauge theory on the unit lattice with density

    ρ_0(U) = exp( v|Λ| − β̄ 𝒦(U) − ℒ(U) − ℰ(U) ) 𝒢(U),

where 𝒦 is a normalized Wilson-like action, ℒ a coupling-independent loop term, ℰ a sum of small regular polymer activities (norm ≤ β̂^{−11/10}), 𝒢 a gas of large-field activities (norm ≤ exp(−P^{1/5}), P = exp((log log β̂)²)). Its total mass equals the original Wilson partition function exactly, at every level. **A mass gap for the continuum theory is a statement about this unit-scale lattice system, uniformly in the torus size.** See [03-MASS-GAP-STRATEGY.md](03-MASS-GAP-STRATEGY.md).

## 5. Where things are

| Need | Go to |
|---|---|
| Exact formal statement | `OpenAI-YM-20261028/ChallengeSpec.lean`, `PublicWilson.lean`, `PublicTheory.lean`, `Construction.lean` |
| Same statement in ordinary mathematics, clause by clause | `yang-mills-informalization/statement.md` §3 and §"Correspondence with the Lean source" |
| Whole proof in 10,000 words | `yang-mills-informalization/summary.md` |
| One unit in depth | `yang-mills-informalization/outline/units/<unit>.md` (32 units; I.1 … VIII.1) |
| Full prose proof, draft | `yang-mills-informalization/manuscript/main.pdf` (278 pages, AI-written, under review) |
| What the code's nicknames mean | `yang-mills-informalization/glossary.md` |
| Which declaration belongs to which unit | `yang-mills-informalization/data/declaration-units.tsv`, `key-statements.tsv` |
| Which ideas are old, which are new | `yang-mills-informalization/provenance.md` (§1, §4, §5, §6.3) |
| Independent validation | `YM/docs/assessment-1.pdf`; checks A1–A8 under `YM/submissions/openai-20261028/checks/` |
| How to rebuild | `YM/submissions/openai-20261028/docs/rerun.md`; Armstrong README "Compiled Lean files" |
| The mass-gap template in two dimensions | `openai/math/preprints/`: N1 `Exponential-decay-in-two-dimensional-classical-On-models-…`, N2 `Sharp-mass-bounds-for-the-two-dimensional-O4-model-…`, N3 `The-canonical-massive-continuum-limit-of-the-two-dimensional-O3-model-…` (LaTeX under each `build/`) |

## 6. Facts about the artifact

- Lean 4.32.0, Mathlib pinned in `lake-manifest.json`. A clean build took about 184 minutes on Apple Silicon.
- 250,319 project declarations; 100,168 are reached by the 17 audited results. Only the axioms `propext`, `Classical.choice`, `Quot.sound`.
- The group's Comparator runs were on macOS with `fake-landrun.sh`, so they "do not test isolation of potentially hostile build code". No second kernel implementation was used. Both are openings for us (tasks YM-A01, YM-A02).
- Names in the code are arbitrary (`Treasury`, `Jewel`, `Sheet`). Never infer meaning from a name. Use the glossary, then the definition.
- A second OpenAI submission dated 4 November 2026 is listed by the group as "work in progress". Its content is unknown to us. Task YM-W01 watches for it.
