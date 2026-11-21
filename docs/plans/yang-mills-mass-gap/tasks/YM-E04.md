# YM-E04 — Infrared dynamics and chiral symmetry with matter

Status: ready for the first bounded definition/source milestone; no proof worker claimed and no theorem completed. Tier T3 for mathematical formulation; later extraction/exact checks may be delegated as bounded T1/T2 units. Independent substantive review required.

Depends on E03's model choice and coordinates with E02's confinement criterion. First milestone: an infrared question/hypothesis table. Specify exact global versus anomalous symmetries, explicit fermion masses, intended order parameters and limits; ask whether a gap/confinement/chiral breaking holds instead of requiring all at once.

If using a renormalized condensate, define its source/mass and volume limits; for example Σ=lim_(m→0+) lim_(V→∞) <barψψ>_ren is a selected-model target, not an automatic universal criterion. Under the hypotheses of Goldstone's theorem, spontaneous breaking of an exact continuous global symmetry produces massless modes. Distinguish the full physical spectrum from a proposed glue sector or pseudo-Goldstone masses with explicit breaking.

Acceptance: symmetry/anomaly/limit specification, Goldstone consistency check, source ledger and one bounded falsifiable next lemma. Numerical condensates or finite-volume spectra alone are not a four-dimensional continuum theorem. Prepare a reusable symbolic Ward-identity/spectral tool map only; no independent code implementation in this first packet.

Read [the extension roadmap](../08-EXTENSIONS-ROADMAP.md), [worker prompt](../prompts/yang-mills-extensions-worker.md), independent story `user_stories/story__yang-mills-extensions.md`, and the full canonical tracker thread before claiming. Work in your own ym/* worktree, announce manuscript scope, and preserve core G02 review. Follow STYLE.md: explicit algebra, locally repeated prerequisites and linked local equations. Separate theorem/proposal/exact/numerical/Lean evidence; Rust first, CUDA only after a measured workload. Implementation requires a later bounded ownership/PDD step, not authorization to rewrite existing tools.
