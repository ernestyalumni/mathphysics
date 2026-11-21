# YM-E01 — Isolated one-particle mass shell

Status: ready for the first bounded definition/source milestone; no proof worker claimed and no theorem completed. Tier T3 for mathematical formulation; later extraction/exact checks may be delegated as bounded T1/T2 units. Independent substantive review required.

Define H, P, M=sqrt(H²−P²), E_M and the physical Hilbert space. State a nonzero isolated positive mass projection together with the intended positive-energy Poincaré representation. Explain the moving-particle dispersion E(p)=sqrt(m²+|p|²) explicitly. Distinguish a lower energy gap, a mass-shell atom, isolation from the remaining mass spectrum, and asymptotic completeness.

First milestone: a source-grounded target/dependency sheet and a negative control. A generalized free scalar with positive continuous mass spectral density supported at m≥m0>0 and no atom is a useful control: a gap need not produce an isolated one-particle mass shell. State the spectral measure and assumptions rather than merely naming the example. Finite-box eigenvalues or exponential fits alone are insufficient. Reuse existing spectral/exact infrastructure in a future tool map; implement no chapter-specific solver.

Acceptance: precise spectral target, representation qualifications, checked control, theorem ledger and remaining constructive steps. Independent mathematical review required; particle existence is not yet proved.

Read [the extension roadmap](../08-EXTENSIONS-ROADMAP.md), [worker prompt](../prompts/yang-mills-extensions-worker.md), independent story `user_stories/story__yang-mills-extensions.md`, and the full canonical tracker thread before claiming. Work in your own ym/* worktree, announce manuscript scope, and preserve core G02 review. Follow STYLE.md: explicit algebra, locally repeated prerequisites and linked local equations. Separate theorem/proposal/exact/numerical/Lean evidence; Rust first, CUDA only after a measured workload. Implementation requires a later bounded ownership/PDD step, not authorization to rewrite existing tools.
