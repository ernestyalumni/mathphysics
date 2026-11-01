# Rust mathematical physics platform

Ernest approved the two user stories and implementation plan on 2026-10-31 and explicitly authorized subagents. The accepted meaning is preserved in [the program plan](plans/rust-physics-platform/README.md) and [original request](plans/rust-physics-platform/request.txt). Progress, concrete implementation and remaining gaps are in [EXECUTION.md](plans/rust-physics-platform/EXECUTION.md).

The program builds primarily Rust mathematical physics tools for the user's chosen agents, with explicit distinctions between computation, numerical evidence, exact identities, certificate replay and formal proof. It covers all 97 Srednicki chapters, their intermediate derivations/examples/prose/exercises, and a broader extensible CAS capability program. Language order is Rust, justified C++/CUDA C++, then Python; frontend code is TypeScript.

This record does not widen existing prompt ownership. Each implementation unit has its own bounded configuration and source prompt graph; original textbook/solution sources and independent tests remain separate. Full-book implementation, general CAS parity and automated PDD regeneration are not inferred from plan approval or metadata consistency.

## Unified notebook UX

The third story is [Unified mathematical physics notebook workspace](../user_stories/story__unified-notebook-workspace.md): one primary browser UI and double-click Linux/macOS apps, shared TypeScript notebook model and Rust services, a close Mathematica-style interface and an equally readable SNES pixel theme with sound controls. Both platforms connect available Claude Code, Codex, Grok Build, Hermes Agent and Pi through actual supported tool/MCP/CLI/API routes. The user explicitly requested this stream proceed alongside Srednicki tooling. [UX-DELIVERY.md](plans/rust-physics-platform/UX-DELIVERY.md) and UX01–UX10 coordinate the work; current frontend/desktop/harness readiness is not inferred from a story or reference download.

## Reusable research library requirement

The [reusable Rust research library story](../user_stories/story__reusable-rust-research-library.md) makes the architecture explicit: Srednicki is the first requirements/validation corpus, not a hierarchy of per-subsection production functions. Mathematical capabilities are implemented once behind reusable public Rust APIs, used by thin textbook/paper fixtures, agents and both independent UIs. Extending to other books and versioned arXiv research must not clone the engine. Production readiness and proof status require evidence; neither is inferred from a library name or passing source fixtures. See [architecture and dispatch policy](plans/rust-physics-platform/REUSABLE-LIBRARY.md).
