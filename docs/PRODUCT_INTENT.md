# Rust mathematical physics platform

Ernest approved the two user stories and implementation plan on 2026-10-31 and explicitly authorized subagents. The accepted meaning is preserved in [the program plan](plans/rust-physics-platform/README.md) and [original request](plans/rust-physics-platform/request.txt). Progress, concrete implementation and remaining gaps are in [EXECUTION.md](plans/rust-physics-platform/EXECUTION.md).

The program builds primarily Rust mathematical physics tools for the user's chosen agents, with explicit distinctions between computation, numerical evidence, exact identities, certificate replay and formal proof. It covers all 97 Srednicki chapters, their intermediate derivations/examples/prose/exercises, and a broader extensible CAS capability program. Language order is Rust, justified C++/CUDA C++, then Python; frontend code is TypeScript.

This record does not widen existing prompt ownership. Each implementation unit has its own bounded configuration and source prompt graph; original textbook/solution sources and independent tests remain separate. Full-book implementation, general CAS parity and automated PDD regeneration are not inferred from plan approval or metadata consistency.

## Unified notebook UX

The third story is [Unified mathematical physics notebook workspace](../user_stories/story__unified-notebook-workspace.md): one primary browser UI and double-click Linux/macOS apps, shared TypeScript notebook model and Rust services, a close Mathematica-style interface and an equally readable SNES pixel theme with sound controls. Both platforms connect available Claude Code, Codex, Grok Build, Hermes Agent and Pi through actual supported tool/MCP/CLI/API routes. The user explicitly requested this stream proceed alongside Srednicki tooling. [UX-DELIVERY.md](plans/rust-physics-platform/UX-DELIVERY.md) and UX01–UX10 coordinate the work; current frontend/desktop/harness readiness is not inferred from a story or reference download.

## Reusable research library requirement

The [reusable Rust research library story](../user_stories/story__reusable-rust-research-library.md) makes the architecture explicit: Srednicki is the first requirements/validation corpus, not a hierarchy of per-subsection production functions. Mathematical capabilities are implemented once behind reusable public Rust APIs, used by thin textbook/paper fixtures, agents and both independent UIs. Extending to other books and versioned arXiv research must not clone the engine. Production readiness and proof status require evidence; neither is inferred from a library name or passing source fixtures. See [architecture and dispatch policy](plans/rust-physics-platform/REUSABLE-LIBRARY.md).

<!-- pdd-intent-entry:yang-mills-construction-audit-rust-lattice-gauge-d3c6da77:start -->
## Yang-Mills construction audit, Rust lattice-gauge tooling and mass-gap research program

- Intent event: [`docs/intents/intent__yang-mills-construction-audit-rust-lattice-gauge-d3c6da77.md`](intents/intent__yang-mills-construction-audit-rust-lattice-gauge-d3c6da77.md)
- Change kind: `add`
- Supersedes: none
- Scope: `existing_pdd_change`
- Technology: `Rust`

> Given our repo: /media/propdev/Expansion/openclaw/.openclaw/workspace/repos/mathphysics and our attempts to make AI tooling for symbolic computation for high energy particle physics, namely quantum field theory, from a Rust first replacement for Sage Math, sympy, Cadabra2, Mathematica, consider these 2 repos: /media/propdev/Expansion/openclaw/.openclaw/workspace/repos/yang-mills-group/YM , /media/propdev/Expansion/openclaw/.openclaw/workspace/repos/yang-mills-group/OpenAI-YM-20261028 and these 2 accompanying repos: /media/propdev/Expansion/openclaw/.openclaw/workspace/repos/hairer/humanise-yang-mills and /media/propdev/Expansion/openclaw/.openclaw/workspace/repos/scottnarmstrong/yang-mills-informalization so it appears they have a 4D Yang Mills proof. Your job is learn as much as you can about what OpenAI did where apparently they constructed 4D pure Yang-Mills theory as a scaling limit of Wilson lattice gauge theory. As a user, I'd like to take on this effort to explore and prove or disprove parts or as a whole this proof and adopt as AI tools the techniques here and to also further developer our tools, again the Rust first replacement for Sage Math, Cadabra2, etc. for symbolic computation. Double check me if it's true that this is only half of a Millenium Prize, and the other half is showing mass gap. If this is true then as a user (do all work in mathphysics, with PDD), I want to be able to pursue research of mass gap and win the prize and resolve the questions about mass gap. Please as you are working, think strategically about .md files, instruction files, instructions, commands out to other less powerful LLMs tasks they could do, help me break down the tasks so less powerful LLMs can take up the work in parallel and be talking to each other, checking each other, collaborating with each other on the work; you as Fable 5.2 is to think strategically and plan out the research to be done by the other AI agents. Keep writing down .md files for progress so in case usage limit is reached, work can be handled off
<!-- pdd-intent-entry:yang-mills-construction-audit-rust-lattice-gauge-d3c6da77:end -->

<!-- pdd-intent-entry:yang-mills-program-amendment-local-network-agent-fb139368:start -->
## Yang-Mills program amendment: local-network agents, Rust-first discovery tools, Lean in Docker

- Intent event: [`docs/intents/intent__yang-mills-program-amendment-local-network-agent-fb139368.md`](intents/intent__yang-mills-program-amendment-local-network-agent-fb139368.md)
- Change kind: `clarify`
- Supersedes: none
- Scope: `existing_pdd_change`
- Technology: `Rust`

> Amendment by Ernest, 2026-11-15, to the Yang-Mills audit, tooling and mass-gap program in mathphysics. The review card in docs/plans/yang-mills-mass-gap/05-TOOLING-ROADMAP.md is approved with these changes. Agents MUST be allowed to talk to other local agents, including agents running on other machines on the same local network, because Ernest runs AI agents on this endeavor on that network. Agents MUST NOT post, publish, upload or message outside the local network without Ernest's explicit approval, with one standing exception: an agent that Ernest has given GitHub access MAY contribute to Hairer's public repository hairer/humanise-yang-mills on its own branch, following that repository's coordination rules. The symbolic computation and other AI tools MUST be Rust first, and CUDA C++ MUST be used only where a measured data-parallel workload benefits. The program MUST NOT rely completely on Lean: new Rust symbolic and exact-computation tools are the discovery and falsification loop, and Lean is the final certification. Lean MUST run in a Docker image, with the toolchain and compiled build stored on the Expansion drive. Worker tiers use the harnesses Ernest has: Pi, Hermes Agent, Codex, Claude Code, Grok Build, and local Qwen3.5 models when he is away. For example, a Qwen3.5 worker on a second machine on the local network claims packet YM-T05 through the shared tracker, pushes its branch to the canonical mathphysics checkout over the local network, and a Codex reviewer on the desktop accepts or rejects it. A written account of the construction and of the mass-gap programme, with derivations, propositions and theorems and a map to the Lean declarations, belongs in documents/mono-notes.
<!-- pdd-intent-entry:yang-mills-program-amendment-local-network-agent-fb139368:end -->


## Manually recorded research extension — 21 November 2026 UTC

Ernest asks our agents to investigate Jaffe–Witten's isolated-particle, confinement, gauge–matter, infrared-dynamics and four-manifold extensions, and to write the Wightman axioms explicitly from local sources. [Research intent](intents/intent__yang-mills-extensions-2026-11-21.md), [independent story](../user_stories/story__yang-mills-extensions.md), and [bounded roadmap](plans/yang-mills-mass-gap/08-EXTENSIONS-ROADMAP.md). These extend our research goals, not the Millennium prize requirements. This manual record is outside generated intent markers; no provider-backed apply or generated contract is claimed.
