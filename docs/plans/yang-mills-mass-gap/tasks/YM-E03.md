# YM-E03 — Existence of selected asymptotically free gauge–matter theories

Status: ready for the first bounded definition/source milestone; no proof worker claimed and no theorem completed. Tier T3 for mathematical formulation; later extraction/exact checks may be delegated as bounded T1/T2 units. Independent substantive review required.

First milestone: select one consistent model, inventory local sources and give a normalization-aware ultraviolet checklist. A useful initial proposal is vectorlike SU(N) with Nf Dirac fundamentals; it is a starting model, not a restriction of Ernest's broader intent. Check gauge anomalies and all couplings relevant to the selected field content.

With textbook tr(TaTb)=δab/2, write β(g)=−g³[11 CA/3−4 T(R)Nf/3]/(16π²)+O(g⁵). For SU(N), CA=N and T(fund)=1/2, giving bracket (11N−2Nf)/3 and one-loop AF when Nf<11N/2. For SU(3), Nf=0,3,16,17 gives brackets 11,9,1/3,−1/3 respectively. Verify against Srednicki/local source conventions before adopting. The existing notes' normalized Lie-algebra Casimir is not automatically textbook CA=N; write the normalization dictionary explicitly.

Acceptance: model, anomaly/coupling checks, verified coefficient convention, known existence assumptions and missing continuum step. One-loop AF is not a constructive existence proof. Later exact representation-index and beta bookkeeping should reuse symbolic/Lie-algebra libraries under their owning PDD units.

Read [the extension roadmap](../08-EXTENSIONS-ROADMAP.md), [worker prompt](../prompts/yang-mills-extensions-worker.md), independent story `user_stories/story__yang-mills-extensions.md`, and the full canonical tracker thread before claiming. Work in your own ym/* worktree, announce manuscript scope, and preserve core G02 review. Follow STYLE.md: explicit algebra, locally repeated prerequisites and linked local equations. Separate theorem/proposal/exact/numerical/Lean evidence; Rust first, CUDA only after a measured workload. Implementation requires a later bounded ownership/PDD step, not authorization to rewrite existing tools.
