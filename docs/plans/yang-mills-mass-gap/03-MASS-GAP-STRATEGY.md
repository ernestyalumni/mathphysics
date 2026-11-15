# Mass-gap research strategy

Author: strategist (Fable 5.2), 2026-11-14. Label for every mathematical assertion here that is not a quotation of a source: **conjectured**. This file is a plan, not a proof. It is revised only by the strategist role; workers propose changes through a `DISPUTE` or `FINDING` comment (see [06-AGENT-PROTOCOL.md](06-AGENT-PROTOCOL.md)).

## 1. What has to be shown

Fix G and the Wilson family of the construction. Two targets:

- **(GAP)** In the reconstructed theory of every phase ξ, the Hamiltonian satisfies spec(H_ξ) ⊂ {0} ∪ [Δ, ∞) with one Δ > 0.
- **(UNIQ)** The infinite-volume limit does not depend on the subsequence, and the phase space X is a single point.

Known relations (all **source-claims**, assessment §6.2): if S^mix clusters then all phases of that limit coincide; clustering does not give independence of the subsequence. A gap in every phase does not by itself make the phases equal.

What a gap means on the Euclidean side: for local gauge-invariant fields A, B, the truncated correlation ⟨A(0) B(x)⟩ − ⟨A⟩⟨B⟩ decays like e^{−Δ|x|}. What is proved today is s⁻² decay.

## 2. Where the gap has to come from

The ultraviolet construction ends at **unit scale** with a lattice gauge system whose density is

    ρ_0(U) = exp( v|Λ| − β̄ 𝒦(U) − ℒ(U) − ℰ(U) ) 𝒢(U)

on a torus of side L_p (see [01-WHAT-OPENAI-PROVED.md](01-WHAT-OPENAI-PROVED.md) §4). Three facts decide the strategy.

1. **β̄ is large and fixed.** The ultraviolet analysis needs weak coupling at every level including the last. So the unit-scale system is a weakly coupled four-dimensional lattice gauge theory with a Wilson-like action and small corrections.
2. **The partition function is preserved exactly.** ∫ρ_m equals the Wilson partition function at every level m. Ratios of partition functions on different tori are therefore the same whether computed on the fine lattice or at unit scale.
3. **Nothing has been integrated out above unit scale.** The physical mass in these units is expected to be of order exp(−β̄/(2κ_G)) up to powers of β̄ (physics expectation, not a theorem): exponentially small. The correlation length is exponentially large in β̄, and at that length the effective coupling is of order one.

So (GAP) is equivalent in spirit to the classic open problem: **a mass gap for four-dimensional lattice gauge theory at weak coupling, uniformly in the volume**, for a class of actions close to Wilson's. Strong coupling (small β) has been understood since Osterwalder–Seiler. The renormalization group here controls large β at short distances. No method controls the crossover, where there is no small parameter. Any proposal that does not say how it crosses that region is not a proposal.

## 3. The template: OpenAI's two-dimensional papers

Family 215 in `openai/math` does the whole job for sigma models. Label: **source-claims**; none of this is checked by us, and only part of the family is formalized.

| Paper | Result | Role in a gap proof |
|---|---|---|
| N1 (20 pp.) | Exponential decay of the spin two-point function for O(n≥3), every β, uniform over subgraphs and weakened bonds | The new infrared idea |
| N2 (102 pp.) | O(4): full transfer gap between c√β e^{−πβ} and C√β e^{−πβ}; unique periodic limit; positive full gap at every β | Infrared idea ⇒ full gap; ultraviolet RG ⇒ sharp scale |
| N3 (122 pp.) | O(3): canonical massive continuum limit, no subsequences, OS axioms, positive gap | The complete programme including uniqueness |

How N1 works, from its §1.3:

1. Rotate spins by a spatially varying angle; this gives exact identities for derivatives of the partition function.
2. Two different rotation axes produce a response along their commutator. The **fourth variation**, summed over a family of test functions, is a nonnegative square plus controlled local terms (N1 Prop. 2.2).
3. On an annulus of inner radius N with boundary layers pinned at relative angle t, test families at all spatial frequencies give a commutator field of size log N, so a quadratic response (log N)², against errors O(log N). Result: Z(π)/Z(0) ≥ 1 − O(1/log N), uniformly in which sites and bonds are present.
4. The signs of one spin component form a ferromagnetic Ising system. The probability that a sign cluster crosses the pinned annulus is exactly 1 − Z(π)/Z(0). Disjoint pinned annuli are independent, so many crossings are exponentially unlikely.

How N2 turns this into a full gap, from its abstract and §1.1: reflection positivity gives a positive transfer operator; "the finite-volume argument uses the complete trace, so multiplicities and invariant sectors enter automatically"; it "turns a small partition-function defect on one box into a spectral estimate uniform in larger circumferences" (N2 §2). A preliminary mixing estimate for all local observables (N2 §3) supplies the defect bound at any β. The exact renormalization map (N2 §§4–8) then carries a fixed-box test from a reference cutoff to all finer cutoffs (N2 §9).

**The Yang–Mills Lean contains the analogue of N2 §§4–8 and nothing else from this list.** The gap programme for Yang–Mills is the rest of the list.

## 4. Decomposition

| ID | Piece | Nature | Difficulty | Template |
|---|---|---|---|---|
| M0 | Precise statements: the class 𝒞(β̄) of unit-scale actions; conjectures L-twist, L-mix, L-defect below; the reduction theorem | Formulation | Moderate | N2 §1, §11 |
| M1 | Defect ⇒ gap: a small partition-function defect on one box gives a transfer gap uniform in larger boxes, for reflection-positive lattice gauge theory | Port | Moderate | N2 §2 |
| M2 | Transport: the defect on a fixed physical box is the same at every cutoff, so the gap bound is uniform in a_N and passes to each continuum phase | Port onto existing Lean | Hard, long | N2 §9, N3 §8 |
| M3 | **The infrared input**: prove the defect bound (or mixing, or the twist estimate that implies it) for the class 𝒞(β̄) | New mathematics | The open problem | N1; N2 §3 |
| M4 | Uniqueness: remove the subsequence and collapse the phases | Port | Hard | N2 Thm 1.1, N3 §§12–15 |
| M5 | Formalize in Lean on top of `OpenAI-YM-20261028`; write to referee standard | Engineering | Long | — |

M1, M2, M4 are conditional on M3 and can proceed in parallel with it. Together they give a theorem of the form "**if** the unit-scale system satisfies L-defect, **then** (GAP) and (UNIQ)". That conditional theorem is worth having even if we never prove M3: it tells everyone exactly what is missing, in checkable form.

### The three forms of the infrared conjecture (M0 drafts)

Let 𝒞(β̄) be the set of gauge-invariant, reflection-positive lattice systems on four-dimensional tori that the ultraviolet construction can output at level 0. Task YM-G01 extracts its exact definition from units III.9, III.10, III.12.

- **L-defect.** There are T_0, L_0 and δ < 1 such that for every system in 𝒞(β̄) on a box of spatial side L_0 the defect D(L_0, T_0) = Z(L_0³×T_0)² / Z(L_0³×2T_0) − 1 is at most δ. (Exact shape to be taken from N2 §2.)
- **L-mix.** Truncated correlations of all local gauge-invariant observables decay exponentially, with rate and constants uniform over 𝒞(β̄) and over the volume.
- **L-twist.** For a vortex container at scale N, the ratio of the partition function with a centre twist to the one without is at least 1 − ε(N) with ε(N) → 0, uniformly over 𝒞(β̄) and over deleting or weakening plaquettes.

Expected implications: L-twist ⇒ (with more work) L-mix ⇒ L-defect ⇒ gap. Each arrow is a task.

An elementary fact that shows why defects control gaps (derived here, label **conjectured** until YM-G02 checks it and compares with N2 §2). Let T be a positive trace-class operator with eigenvalues λ_0 ≥ λ_1 ≥ … ≥ 0 and Z(n) = Tr Tⁿ. Put t_i = (λ_i/λ_0)ⁿ and s = Σ_{i≥1} t_i. Then

    Z(n)² / Z(2n) = (1+s)² / (1 + Σ_{i≥1} t_i²) ≥ (1+s)² / (1+s) = 1 + s,

because t_i ≤ 1. So D(n) := Z(n)²/Z(2n) − 1 ≥ Σ_{i≥1} (λ_i/λ_0)ⁿ ≥ (λ_1/λ_0)ⁿ, and

    gap = log(λ_0/λ_1) ≥ (1/n) log(1/D(n)).

The weakness is visible at once: at fixed n the sum s grows with the spatial volume, because the number of low-lying states does. N2 §2.2 ("Increasing the spatial volume") is where the two-dimensional proof deals with that, using reflection positivity in the space direction as well. Porting that step to four dimensions is M1.

## 5. Routes to the infrared input (M3)

Ranked by my current judgment. Each route has a first task that is cheap and decisive.

### R1. A gauge version of N1 (primary)

Dictionary:

| N1 (2D spins) | 4D gauge theory |
|---|---|
| Spin σ_x ∈ S² | Link U_ℓ ∈ SU(2) |
| Edge twist A_e ∈ so(3); W(A) = twisted/untwisted partition function | Plaquette twist B_p ∈ su(2): U_p ↦ U_p exp(B_p) |
| Twist by angle π flips the embedded Ising sign | Twist by the centre element −1 = exp(π X), |X| = 1: a 't Hooft twist |
| Two rotation axes a, b with [a,b] = c | Two colour directions with [X_a, X_b] = X_c |
| Annulus between two pinned boundary layers | Vortex container: a thickened closed surface (Mack–Petkova) |
| Z(π)/Z(0) ≥ 1 − O(1/log N) | Vortex free energy of the container tends to zero |
| Embedded ferromagnetic Ising signs; crossing probability = 1 − Z(π)/Z(0) | Embedded Z₂ gauge theory of link signs; thin-vortex sheet through the container |
| Independent pinned annuli ⇒ exponential decay of ⟨σ_x·σ_y⟩ | Stacked containers linking a loop ⇒ area law for Wilson loops (Mack–Petkova, Tomboulis–Yaffe inequalities) |
| N2 §3: from charged observables to all local observables | From Wilson loops to plaquette–plaquette correlations: **no template in the Lean; must be built** |

Why this is plausible: both theories are marginal with logarithmic running, and N1's logarithm comes from summing fluctuation scales between the lattice and N, which is asymptotic freedom in a different language. In gauge theory a single-plane centre twist has zero classical cost ('t Hooft's twist-eating flat connections), so the twist ratio is decided by fluctuations alone. That is the regime where a fluctuation identity like N1 Prop. 2.2 has a chance.

Known obstacles, stated now so nobody rediscovers them late:

1. **Frustration.** Conditioning on |σ¹| gives a ferromagnetic Ising model. Conditioning SU(2) links on their SO(3) classes gives a Z₂ gauge theory whose plaquette couplings carry the sign of tr V_p, which cannot be gauged away where there are Z₂ monopoles. Ginibre and FKG inequalities then fail. Mack and Petkova proved their results for a modified model without monopoles. Stepping stone: the **positive-plaquette model** (constraint tr U_p > 0), where the embedded Z₂ theory is ferromagnetic. At large β̄ monopoles are rare (large-field regions cost e^{−cP²}), which may let them be treated as a dilute perturbation. Unproven.
2. **Planarity.** N1's crossing argument is planar. In four dimensions the objects are surfaces linking loops. The inequalities exist in the literature (Mack–Petkova 1979–80, Tomboulis–Yaffe 1985) but under hypotheses to be checked.
3. **Area law is not a gap.** Confinement of Wilson loops does not give decay of local observables. N2 needed a separate, long argument (§3) and it used two-dimensional tools.
4. **The class 𝒞(β̄).** N1 is uniform over weakened and deleted bonds, which is what absorbs the renormalization group's corrections. The gauge estimate needs the same robustness.

First tasks: YM-R01 (digest N1 with every identity written out), YM-G03 (derive the gauge fourth-variation identity symbolically and test on small lattices with exact arithmetic whether a sign-definite structure exists), YM-G04 (write the embedded Z₂ representation and locate the frustration exactly), YM-N02 (measure twist ratios on containers by Monte Carlo at several β and N).

### R2. Decimation sandwiches (Tomboulis 2007)

Tomboulis claimed confinement at all β for SU(2) by sandwiching the partition function and the twisted partition function between Migdal–Kadanoff decimations. Ito and Seiler pointed out a gap in the interpolation step. The exactly normalized block kernel of the OpenAI construction is a new tool of the same family. First task: YM-R06 (digest the claim and the objection, state precisely which inequality is missing, and ask whether exact partition-function preservation supplies it). Low cost, low prior.

### R3. Push the convergence region from the strong-coupling side

Osterwalder–Seiler cluster expansion; Bakry–Émery / Langevin arguments (Shen–Zhu–Zhu) give a gap for β below a constant depending on N; master-loop methods (Cao–Nissim–Sheffield) extend area-law regimes at large N. None reaches weak coupling, and β̄ is not ours to choose small. This route matters only if combined with a renormalization group continued into the infrared until the coupling is in the convergent region, which needs control of an intermediate regime with no small parameter. First task: YM-R07 (table of every rigorous gap or area-law result with its β range and group; label each range exactly). This is also the literature baseline for anything we write.

### R4. Finite-size criteria with rigorous computation

Results of Aizenman–Simon and Dobrushin–Shlosman type turn a finite-volume inequality into exponential decay. If such a criterion existed for 𝒞(β̄) on a box of feasible size, interval arithmetic could verify it. The box needed is of the order of the correlation length, which is exponentially large in β̄. Almost certainly infeasible directly; possibly feasible **after** renormalization to a scale where the coupling is moderate, if the effective action could be bounded rigorously there. First task: YM-G06 (feasibility estimate only: what box, what precision, what cost; no implementation).

### R5. Probabilistic and stochastic routes

Chatterjee's programme (strong mass gap plus centre symmetry implies confinement; the converse direction is what we need and is open); stochastic quantization of Yang–Mills in two and three dimensions (Chandra–Chevyrev–Hairer–Shen). Useful as a source of inequalities and as the language of the people who will referee. No first task beyond the literature table of R3.

### What I am not pursuing

Large-N and holographic arguments (not rigorous); numerical lattice spectroscopy as evidence for the gap (nobody doubts it); any approach that works in the continuum without a cutoff.

## 6. Experiments before proofs

Cheap computations that can kill or sharpen a route before anyone spends weeks on it. All under workstream B tooling.

| ID | Experiment | Decides |
|---|---|---|
| YM-N01 | SU(2) Wilson Monte Carlo, validated against published plaquette values | Baseline for everything else |
| YM-N02 | Twist-ratio Z_tw/Z on containers, several β, several scales | Whether L-twist has the 1 − O(1/log N) shape, a different shape, or fails at accessible sizes |
| YM-N03 | Defect D(L,T) in a model where the answer is known (2D O(3), then 3D Z₂ gauge) | Whether the defect criterion is numerically usable and how it scales with volume |
| YM-N04 | Z₂-monopole density and its correlation with twist response in SU(2), and the same in the positive-plaquette model | Whether frustration is a dilute perturbation at weak coupling |
| YM-G03 | Exact small-lattice evaluation of the gauge fourth variation | Whether the positive-square structure of N1 Prop. 2.2 survives |

## 7. Division of labour

| Tier | Does | Does not |
|---|---|---|
| T4 strategist (Fable-class) | Chooses routes, writes and retires packets, adjudicates disputes, writes the weekly synthesis, attempts the key new arguments | Bulk reading, bulk coding |
| T3 (Opus-class, or strongest Codex/Grok) | Mathematical digests of hard sources, derivations, proof attempts on bounded lemmas, red-team reviews | Decide strategy |
| T2 (Sonnet-class, Codex) | Rust units with tests, Monte Carlo, literature tables, reviews of T2 work, Lean statement transcription | Unsupervised proof claims |
| T1 (Haiku-class, local models) | Builds, reruns, extraction into fixed templates, watch tasks, formatting | Anything requiring a judgment about mathematics |

Every packet names its tier. A lower tier may attempt a higher-tier packet only as a *second opinion* whose output is labelled as such.

## 8. Milestones, each valuable on its own

| ID | Milestone | Standalone value | Depends on |
|---|---|---|---|
| G0 | Problem sheet: 𝒞(β̄) defined from the Lean, three conjectures stated, dictionary checked against N1–N3 | A precise open problem others can cite | YM-R01–R05, YM-G01 |
| G1 | Linux rebuild with real isolation, second-kernel check if available, spec red-team report | Fills two gaps the review group named | YM-O01, YM-A01–A03 |
| G2 | Defect ⇒ gap lemma for lattice gauge theory, written and stated in Lean | Reusable; first formal step beyond the submission | YM-G02 |
| G3 | Gauge fourth-variation identity, machine-checked on small lattices; verdict on sign structure | Decides R1 | YM-G03 |
| G4 | Twist-ratio numerics with error bars | Evidence for or against L-twist | YM-N01, YM-N02 |
| G5 | Conditional theorem: L-defect ⇒ (GAP) and (UNIQ), natural-language proof plus Lean statements | The reduction; publishable as such | G0, G2, M2, M4 |
| G6 | L-twist or area law for the positive-plaquette model at all β | A new theorem in lattice gauge theory | G3, YM-G04 |
| G7 | L-defect for 𝒞(β̄) | The mass gap | everything |

## 9. Kill and pivot rules

- If YM-G03 finds no sign-definite structure after two independent derivations by different model families, demote R1 below R2 and re-plan.
- If YM-N02 shows twist ratios that do **not** approach 1 with scale at fixed β in the scaling window, stop and find out why before any further proof work on R1; the physics expectation says they must.
- If the 4 November submission (or any later one) proves a piece of M1–M4, stop working on that piece the same day, read it, and move effort to what remains. YM-W01 exists for this.
- Any agent output claiming a proof of M3 or of a full gap is treated as **conjectured** and goes to red-team review by two different model families before the strategist reads it. The base rate for correct proofs of this statement is close to zero; act accordingly.
- No external announcement, upload or message about results without Ernest's explicit approval.

## 10. Reading list, in order

1. N1, all 20 pages. `openai/math/preprints/Exponential-decay-in-two-dimensional-classical-On-models-September-23-2026/` (LaTeX in `build/`).
2. N2 §1, §2, §3.10–3.12, §9, §11.
3. N3 §1, §8, §§12–15.
4. Armstrong `summary.md` Parts III–V; units III.9, III.10, III.12, IV.1, V.1, V.3.
5. Seiler, *Gauge Theories as a Problem of Constructive QFT and Statistical Mechanics* (1982): reflection positivity, transfer matrix, strong coupling.
6. Mack–Petkova (1979, 1980); Tomboulis–Yaffe (1985); Kovács–Tomboulis (1998–2002).
7. Tomboulis (2007) and Ito–Seiler's comment.
8. Chatterjee, *Yang–Mills for probabilists* and the confinement papers; Shen–Zhu–Zhu; Cao–Nissim–Sheffield.

Items 5–8 are cited from memory. YM-R07 must locate each one, confirm the reference, and record exactly which pages were read. No citation enters a deliverable until that is done.
