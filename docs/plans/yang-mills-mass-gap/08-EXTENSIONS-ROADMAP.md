# Yang–Mills extensions: targets before proof attempts

Ernest's source is Jaffe–Witten pp.6–7, not an amendment to the Millennium prize conditions. These are long-term goals. Current submission establishes definitions, study prerequisites and bounded agent assignments; it establishes no new extension theorem. Preserve the existing core G02 review and Rust T15/T09 work.

## Reading dependency

Start with the explicit Wightman section in `documents/mono-notes/topics/qft-foundations.tex`, label `qft:sec:wightman-axioms`, checked against Streater–Wightman §3-1 printed97–102. Continue to Euclidean reconstruction, the literal gap/problem statement and clustering; then read the extension section `qft:sec:yang-mills-extensions`. OS1973/1975 are local, but whole-source OCR and mathematical reconciliation remain pending. Selected book page images were checked independently; this is not a certification of every formula in the book.

| Packet | First bounded milestone | Later target / important boundary |
|---|---|---|
| [E01](tasks/YM-E01.md) | Mass-shell definition, representation assumptions and gap-without-particle control | Isolated one-particle sector; not an isolated full-H eigenvalue |
| [E02](tasks/YM-E02.md) | One explicit confinement criterion and applicability map | Prove it for the selected theory; screening and group global form matter |
| [E03](tasks/YM-E03.md) | One consistent gauge–matter model and normalization-aware UV table | Construct the selected theory; one-loop AF is insufficient |
| [E04](tasks/YM-E04.md) | IR question/hypothesis matrix for the selected matter model | Determine its phase; mass gap and chiral breaking need not coexist |
| [E05](tasks/YM-E05.md) | Geometry/category specification and one known control construction | Extend interacting existence theorems; no universal H on arbitrary manifolds |

E01 uses spectral/Wightman prerequisites. E02 coordinates with existing R06/R07/N02 rather than duplicating their center-twist literature. E03 selects models before E04 makes dynamical claims. E05 starts with a specification, not a promise to transplant the flat-space proof unchanged. These first milestones are ready to claim independently; they do not authorize rewriting another worker's files.

## Evidence and tool path

Each packet returns a short theorem ledger: precise statement, assumptions, source/page, known result versus proposed lemma, missing step, possible counterexample and review status. Small exact controls precede general proof claims. Map operations into existing Rust symbolic/spectral/lattice/geometry interfaces before proposing code, then use the owning PDD unit and independent software review. Numerical experiments falsify or motivate; they are not continuum existence proofs. CUDA follows a measured CPU bottleneck. Lean certification refers to the exact theorem and hypotheses that were checked.

No five broad proof agents were automatically launched. A read-only scope reviewer checked the definitions; bounded workers can now claim the separate tasks through [the worker prompt](prompts/yang-mills-extensions-worker.md). [Scope receipt](notes/YM-extensions-scope-review-2026-11-21.md) records source checks and exclusions. Actual tracker IDs are listed in BOARD.md.
