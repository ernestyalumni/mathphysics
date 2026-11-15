# Source repositories

Five public repositories are the inputs of this program. They are **read-only**. Agents cite them as `<short name>@<pin>:<path>:<line>`.

| Short name | Repository | Pin (commit) | Size | What it is |
|---|---|---|---|---|
| `LEAN` | https://github.com/yang-mills-group/OpenAI-YM-20261028 | `3456e5e75cf8d4cd558955a5c7638fba2d68f464` | 144 MB | OpenAI's Lean 4 proof (28 Oct 2026 submission). 9,508 files, Lean 4.32.0. |
| `YM` | https://github.com/yang-mills-group/YM | `55ba9fc36293b3f4d52c6b1f813169a0be76ee98` | 163 MB | The review group's assessment (`docs/assessment-1.pdf`), examples, validation records. |
| `ARM` | https://github.com/scottnarmstrong/yang-mills-informalization | `daf1db452866dd773b17aa955d534e56196a4efc` | 94 MB | Armstrong's map: exact statement, summary, glossary, 32 unit outlines, data tables linking every declaration to a unit. |
| `HAIRER` | https://github.com/hairer/humanise-yang-mills | `16b070a9f5ebf4c83c210d8ccffba3be2c711c48` | 196 MB | Hairer's crowd-sourced monograph and its agent coordination rules. |
| `OAIMATH` | https://github.com/openai/math | `4eaee87548250b2682affda7d05d947c05554bae` | 3.7 GB | OpenAI's preprints. We need family 215 only: the two-dimensional O(n) papers that are the template for a mass-gap proof. |

## Does an agent need them cloned?

It depends on the task. Check the packet's "Inputs" line.

| Work | Needs |
|---|---|
| Rust tooling (`YM-T*`), numerics (`YM-N*`) | **None.** The oracles are in `ledger/evidence/` and in the packet. |
| Reading the mono-notes write-up, the plan, the strategy | **None.** |
| Micro-tasks (`bin/ym-micro`, planned: packet YM-T14) | **None on the worker.** The runner on the canonical host pastes the needed snippet into the prompt. |
| Statement, summary, glossary, unit outlines (`YM-R02`–`R05`, `R08`, `YM-A03`–`A05`) | `ARM` (94 MB) |
| Assessment (`YM-R08`, `YM-W*`) | `YM` (163 MB) |
| Reading Lean declarations | `LEAN` (144 MB of sources; no build needed to read) |
| Checking Lean (`YM-O01`, `A01`, `A02`) | `LEAN` plus the Docker image and compiled build: see "Lean in Docker" below. Canonical host only. |
| Two-dimensional template (`YM-R01`, `R06`) | `OAIMATH`, but only two directories: `preprints/Exponential-decay-in-two-dimensional-classical-On-models-September-23-2026/` and `preprints/Sharp-mass-bounds-for-the-two-dimensional-O4-model-September-23-2026/` |
| Contributing upstream | `HAIRER`, in a separate writable clone (see `06-AGENT-PROTOCOL.md` §11) |

## Where they are on the canonical host

Under `REPOS = <workspace>/repos/` on the Expansion drive:

    REPOS/yang-mills-group/OpenAI-YM-20261028      LEAN
    REPOS/yang-mills-group/YM                      YM
    REPOS/scottnarmstrong/yang-mills-informalization   ARM
    REPOS/hairer/humanise-yang-mills               HAIRER
    REPOS/openai/math                              OAIMATH
    REPOS/mathphysics                              this repository
    REPOS/yang-mills-work/                         writable: Lean working copy, downloads, logs (not a git repository)

## Getting them on another machine

    bin/ym-sources clone  <target-directory> [LEAN YM ARM HAIRER OAIMATH]   # default: LEAN YM ARM HAIRER
    bin/ym-sources verify <target-directory>                               # prints OK or MISMATCH per repository

`clone` fetches from GitHub over HTTPS without credentials and checks out the pin in detached-HEAD state. `OAIMATH` is 3.7 GB; clone it only if a packet needs it, or copy the two preprint directories from the canonical host over the local network.

Never `git pull` inside a pinned clone. A newer upstream revision is a new pin, and changing a pin is a strategist decision recorded in `HANDOFF.md`.

## Lean in Docker

Nothing Lean-related is installed on the host. Everything runs through `bin/ym-lean` (see the header of that script):

| Thing | Where |
|---|---|
| Image `ym-lean:4.32.0` | Built from `docker/lean/Dockerfile` (Debian slim, elan, Lean 4.32.0, the proof's pinned toolchain) |
| Writable working copy | `REPOS/yang-mills-work/OpenAI-YM-20261028-build` (a clone of `LEAN` at the pin) |
| Mathlib cache | Inside the working copy, fetched by `lake exe cache get` |
| Compiled project build | Release `lean-build-3456e5e` of `ARM` (2.2 GB download, 8.8 GB unpacked), SHA-256 checked |
| Logs | `REPOS/yang-mills-work/logs/` |

`bin/ym-lean offline <command>` runs with networking disabled. Use it for every check. `bin/ym-lean run` has network access and is only for fetching the Mathlib cache.

The existing images `claw-work/lean-agent:lean4.35.0-mathlib` and `monoclaw/lean4-agent:v4.35.0` carry Lean 4.35.0 and a different Mathlib. They cannot load this proof's compiled files. Use them for new, self-contained Lean work only.
