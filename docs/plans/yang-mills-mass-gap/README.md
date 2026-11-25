# Yang–Mills program: audit the construction, build the tools, attack the mass gap

Started 2026-11-14 by Ernest Yeung. Strategist: Fable 5.2 (Cyclonus). Gold version: `master` (protocol [06](06-AGENT-PROTOCOL.md) §4).
Original request, verbatim: [request.txt](request.txt). Amendment of 2026-11-15: [request-2026-11-15-amendment.txt](request-2026-11-15-amendment.txt). Both PDD intents are approved and recorded in `docs/PRODUCT_INTENT.md`; see [HANDOFF.md](HANDOFF.md) for what `pdd intent apply` did and did not do.

**If you are resuming this work, read [HANDOFF.md](HANDOFF.md) first.** It says what is done, what is running, and what to do next.

## The situation in six sentences

1. On 28 October 2026 OpenAI submitted a Lean 4 proof (about one million lines, 9,477 files) that a four-dimensional quantum Yang–Mills theory exists as a scaling limit of Wilson lattice gauge theory, for every compact matrix group with simple nonabelian Lie algebra.
2. A physicist-led review group rebuilt it, audited its axioms, ran Comparator, and wrote its own Lean checks; their assessment of 10 November says the result stands and satisfies Wightman axioms.
3. The proof does **not** show a mass gap, does **not** show the infinite-volume limit is unique, and the group says it gives "no evidence for or against a mass gap".
4. The Clay problem is one problem with one prize: existence **and** mass gap. Existence alone wins nothing, and the existence proved here has caveats the problem statement cares about. Details: [02-MILLENNIUM-PRIZE-STATUS.md](02-MILLENNIUM-PRIZE-STATUS.md).
5. OpenAI's own September–October 2026 papers on two-dimensional sigma models already contain the full template (ultraviolet renormalization group + an all-temperature decay theorem ⇒ massive, unique continuum limit); for Yang–Mills they have published only the ultraviolet half.
6. So the open piece is sharply located: an infrared estimate for four-dimensional lattice gauge theory at one fixed weak coupling, plus the transport of that estimate through the existing construction. Details: [03-MASS-GAP-STRATEGY.md](03-MASS-GAP-STRATEGY.md).

## Read order

| # | File | What it gives you |
|---|---|---|
| 0 | [HANDOFF.md](HANDOFF.md) | Current state and next actions. Always current. |
| 0a | [SOURCES.md](SOURCES.md) | The five source repositories, their pins, which task needs which, how to clone them, Lean in Docker. |
| 0b | `documents/mono-notes/topics/yang-mills-lattice.tex` (Section "Lattice Yang–Mills theory" of `master.pdf`) | The mathematics: definitions, derivations, the theorem, the mass-gap reduction with proofs, conjectures, and a map to Lean declarations. |
| 1 | [01-WHAT-OPENAI-PROVED.md](01-WHAT-OPENAI-PROVED.md) | The theorem, the proof architecture, the non-claims, where each piece lives in the four source repos. |
| 2 | [02-MILLENNIUM-PRIZE-STATUS.md](02-MILLENNIUM-PRIZE-STATUS.md) | Answer to "is this half the prize?" and what a prize-grade result needs. |
| 3 | [03-MASS-GAP-STRATEGY.md](03-MASS-GAP-STRATEGY.md) | The research plan: target statements, routes, ranking, kill criteria. |
| 4 | [04-AUDIT-PLAN.md](04-AUDIT-PLAN.md) | How we try to prove or disprove parts of the construction ourselves. |
| 5 | [05-TOOLING-ROADMAP.md](05-TOOLING-ROADMAP.md) | Rust-first tools this work needs and what each is validated against. |
| 6 | [06-AGENT-PROTOCOL.md](06-AGENT-PROTOCOL.md) | Roles, model tiers, claims, reviews, how agents talk and check each other. |
| 7 | [BOARD.md](BOARD.md) | The task queue: every packet, its tier, status, dependencies. |
| 8 | [tasks/](tasks/) | One self-contained packet per task. A worker needs only its packet plus the protocol. |
| 9 | [prompts/](prompts/) | Copy-paste prompts to launch a worker, a reviewer, or a coordinator. |
| 10 | [ledger/claims.md](ledger/claims.md) | Every claim we rely on, with its truth label and evidence. |

## Source repositories (all local, all read-only for us)

Paths are relative to `repos/` (the parent of this repository).

| Repo | Pin | Role |
|---|---|---|
| `yang-mills-group/OpenAI-YM-20261028` | `3456e5e` | The Lean source as submitted. |
| `yang-mills-group/YM` | `55ba9fc` | The review group's assessment (PDF), their extra Lean checks, evidence, rerun scripts. |
| `scottnarmstrong/yang-mills-informalization` | `daf1db4` | Armstrong's map: statement, 10k-word summary, 32 unit outlines, 278-page manuscript draft, glossary, provenance analysis, data tables. Best entry point for reading. |
| `hairer/humanise-yang-mills` | `16b070a9` | Hairer's crowd-sourced monograph project. Accepts agent contributions. Its coordination rules are the model for ours. |
| `openai/math` | `4eaee87` | OpenAI's 721 preprints. Family 215 holds the six sigma-model papers (N1–N6) that carry the mass-gap template. |

Never edit these checkouts. Cite them by repo, pin, path and line.

## Three workstreams

- **A. Audit.** Reproduce the build on Linux with real process isolation (the group's recorded runs were on macOS without it), red-team the formal statement for degenerate models, and check the informal accounts against the Lean. Output: findings we can send upstream.
- **B. Tools.** Extend the Rust mathematical-physics platform with what this subject needs: exact Lie-algebra data, the one-loop coefficient, free-field correlators, Haar integration and character expansions, a lattice Monte Carlo, twist free energies, a Lean dependency reader. Every tool is validated against a number or identity taken from the proof or the literature.
- **C. Mass gap.** Formulate the lattice infrared conjecture precisely, prove the conditional reduction (conjecture ⇒ gap and uniqueness in the continuum), and attack the conjecture by porting the idea of OpenAI's two-dimensional decay theorem to four-dimensional gauge theory.

## Truth standard

Every statement in this directory carries one of these labels, explicitly or by the file's header:

- **kernel-checked** — a Lean kernel accepted a proof of exactly this formal statement, and we name the declaration.
- **source-claims** — a pinned source says so; we have not checked it.
- **reproduced** — we recomputed a known result and it matched; evidence path given.
- **verified** — a named check ran, a second agent of a different model family reviewed it, evidence path given.
- **conjectured** — proposed, not proved.
- **inconclusive** — a check ran and could not decide.

"Proved" without a label is not allowed. A natural-language proof of a new theorem is **conjectured** until it is kernel-checked or accepted by a named human expert.

## Codex continuation, 2026-11-15

For the Mac Mini first milestone, original continuation request, selective stories and Srednicki reading task, see [07-CONTINUATION-ROADMAP.md](07-CONTINUATION-ROADMAP.md). Claude Code/Fable restart report: [CODEX-HANDOFF-2026-11-15.md](CODEX-HANDOFF-2026-11-15.md). Mac-side prompt: [prompts/mac-mini-bootstrap.md](prompts/mac-mini-bootstrap.md).

Daily public-source collection and versioned AI-review context: [watch/README.md](watch/README.md). Current SSH follow-up and commit attribution: [notes/YM-O02-hostname-trust-followup.md](notes/YM-O02-hostname-trust-followup.md).
