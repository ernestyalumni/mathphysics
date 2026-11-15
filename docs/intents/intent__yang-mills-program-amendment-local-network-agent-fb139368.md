# Intent: Yang-Mills program amendment: local-network agents, Rust-first discovery tools, Lean in Docker

<!-- pdd-intent-id: yang-mills-program-amendment-local-network-agent-fb139368 -->
<!-- pdd-intent-sha256: fb13936825374f389dc07030fee7f44e760392e37321e4be87750c078499d7b9 -->

## Record

- Intent ID: `yang-mills-program-amendment-local-network-agent-fb139368`
- Kind: `clarify`
- Supersedes: none
- Approval ID: `yang-mills-program-amendment-local-network-agent-fb139368`
- Source kind: `inline`
- Source reference: `not applicable`
- Request SHA-256: `fb13936825374f389dc07030fee7f44e760392e37321e4be87750c078499d7b9`
- Project scope: `repository`
- Adoption scenario: `existing_pdd_change`

- Technology: `Rust`

## Original Request

> Amendment by Ernest, 2026-11-15, to the Yang-Mills audit, tooling and mass-gap program in mathphysics. The review card in docs/plans/yang-mills-mass-gap/05-TOOLING-ROADMAP.md is approved with these changes. Agents MUST be allowed to talk to other local agents, including agents running on other machines on the same local network, because Ernest runs AI agents on this endeavor on that network. Agents MUST NOT post, publish, upload or message outside the local network without Ernest's explicit approval, with one standing exception: an agent that Ernest has given GitHub access MAY contribute to Hairer's public repository hairer/humanise-yang-mills on its own branch, following that repository's coordination rules. The symbolic computation and other AI tools MUST be Rust first, and CUDA C++ MUST be used only where a measured data-parallel workload benefits. The program MUST NOT rely completely on Lean: new Rust symbolic and exact-computation tools are the discovery and falsification loop, and Lean is the final certification. Lean MUST run in a Docker image, with the toolchain and compiled build stored on the Expansion drive. Worker tiers use the harnesses Ernest has: Pi, Hermes Agent, Codex, Claude Code, Grok Build, and local Qwen3.5 models when he is away. For example, a Qwen3.5 worker on a second machine on the local network claims packet YM-T05 through the shared tracker, pushes its branch to the canonical mathphysics checkout over the local network, and a Codex reviewer on the desktop accepts or rejects it. A written account of the construction and of the mass-gap programme, with derivations, propositions and theorems and a map to the Lean declarations, belongs in documents/mono-notes.

## Must Stay Unchanged

- Agents MUST be allowed to talk to other local agents, including agents running on other machines on the same local network, because Ernest runs AI agents on this endeavor on that network.
- Agents MUST NOT post, publish, upload or message outside the local network without Ernest's explicit approval, with one standing exception: an agent that Ernest has given GitHub access MAY contribute to Hairer's public repository hairer/humanise-yang-mills on its own branch, following that repository's coordination rules.
- The symbolic computation and other AI tools MUST be Rust first, and CUDA C++ MUST be used only where a measured data-parallel workload benefits.
- The program MUST NOT rely completely on Lean: new Rust symbolic and exact-computation tools are the discovery and falsification loop, and Lean is the final certification.
- Lean MUST run in a Docker image, with the toolchain and compiled build stored on the Expansion drive.

## Examples

- Worker tiers use the harnesses Ernest has: Pi, Hermes Agent, Codex, Claude Code, Grok Build, and local Qwen3.5 models when he is away.
- For example, a Qwen3.5 worker on a second machine on the local network claims packet YM-T05 through the shared tracker, pushes its branch to the canonical mathphysics checkout over the local network, and a Codex reviewer on the desktop accepts or rejects it.

## Candidate Product Areas

- Strict bounded JSON discovery and exact existing Clifford trace bridge
- Offline external-cache workbench build and read-only local gateway launch
