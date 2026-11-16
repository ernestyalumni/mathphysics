# Recommended working route: Wilson lattice, selective Lean, reusable Rust

Requested by Ernest, 16 November2026UTC. This is a research recommendation, not an assertion that the route closes the mass gap or approval to port the whole construction. No new Rust implementation or formal proof was produced by this note.

## Recommendation

Use Wilson lattice gauge theory as the working interface because the assessed existence construction starts there and the existing Rust exact-lattice tools already supply a test bed. Start with the precise finite-volume objects and the continuum bridge. Keep the choice of infrared estimate open: a proposed defect/twist mechanism must earn support from independent checks. Numerical confinement or a fixed-lattice gap is not the continuum theorem.

A target estimate must survive the cutoff and volume limits with physical normalization fixed. If one lattice time step represents physical time a, a dimensionless transfer gap m_lat corresponds to m_phys=m_lat/a. A uniform positive lower bound in physical units, together with the required convergence/domain/correlation control, is what matters. Pointwise positivity at each finite a and volume does not supply it. SU(2) is a useful development case; the full Clay target covers every compact simple gauge group.

## What to understand in Lean first

Read the public statement/definitions and bridges before implementation internals: PublicWilson.lean, PublicTheory.lean, ChallengeSpec.lean, then the selected-phase/EveryPhase decay and mixture endpoints. Existing R08 is the narrow reading packet. Establish the exact hypotheses, order of limits, selected phase versus phase average, and distinction between unique vacuum and unique extracted theory. Use Armstrong/Hairer exposition for navigation and verify the needed formal definitions in the pinned Lean source; AI informalization is not a second independent proof.

A full clean build/audit is a separate task with an explicit resource budget. Replaying the entire proof is not a prerequisite for each finite algebra experiment. Deep Lean reading is justified when a particular hypothesis/endpoint is needed by the proposed mass-gap argument.

## Rust companions to computational mathematics

| Mathematical operation | Existing owner / bounded next capability | Independent check |
|---|---|---|
| Lie-algebra normalization and Casimir constants | Rust/lie-groups, special_unitary | Definition-based trace calculation and separate convention check |
| Lattice cells, gauge words and twists | Rust/lattice-gauge, torus/cochain; later noncommutative group-word calculus | Gauge covariance, coboundary identities and deliberately invalid inputs |
| Partition functions and twist ratios | Existing Z2 enumeration; extend with a transfer-matrix or character-expansion representation | Solvable 2D formula and independently organized enumeration |
| Finite-volume spectral/defect identities | Existing spectral purity helper; genuine finite-group transfer operator is a next bounded extension | Complete finite spectra and independent direct correlation calculation |
| Exact coefficients and group integration | Existing qft-symbolics sphere unit; later bounded Haar/character identities | Analytic orthogonality identities and independent CAS oracle |

Translate mathematical inputs, assumptions, normalization and outputs into reusable APIs. Record the Lean theorem/definition locator as provenance in companion docs. Do not port tactic scripts into production algebra or infer correctness from names. Most compactness, measure-theoretic existence and limiting arguments are not symbolic computations with a direct Rust implementation.

Rust can discover counterexamples and emit checkable algebraic/combinatorial certificates. Lean may check a precisely specified certificate or formalize a proved bridge; that requires explicit connecting theorems. Rust tests, exact finite calculations, a checked certificate and a continuum mass-gap proof remain distinct evidence levels. CUDA follows a measured workload after the CPU mathematical contract is stable.

## First paired study/computation milestone

Read and explain Wilson measure → transfer operator → correlation decay → spectral gap, with source-located hypotheses and limit order. In parallel, extend the existing exact finite-group test bed to compute the actual transfer spectrum and compare it with a directly computed correlation function. Challenge the current defect criterion with known cases and counterexamples. Keep its connection to SU(2)/the 4D continuum explicit and incomplete until established. This milestone checks one useful bridge; it does not reproduce the full OpenAI proof.

## Primary sources checked

- Jaffe–Witten official text, Section4: https://www.claymath.org/wp-content/uploads/2022/06/yangmills.pdf
- Group assessment/status: https://github.com/yang-mills-group/YM
- Check scope and exclusions: https://raw.githubusercontent.com/yang-mills-group/YM/main/submissions/openai-20261028/docs/checks.md
- Formal theorem and entry points: https://github.com/yang-mills-group/OpenAI-YM-20261028

The group supports the October existence construction and explicitly excludes mass gap/infinite-volume uniqueness from the checked conclusions. These links support the scope statement; the choice to prioritize the lattice/selected-tool route is our recommendation.
