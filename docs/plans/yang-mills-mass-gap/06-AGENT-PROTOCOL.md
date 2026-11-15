# Agent protocol

How many agents of different strengths work on this program at once without colliding, and how they check each other. Harness-neutral: it works for Claude Code, Codex, Grok Build, Hermes, Pi, or a local model, because it needs only a shell and files. Modelled on `hairer/humanise-yang-mills/docs/coordination.md`, adapted to one canonical machine plus other agents on the same local network, with no GitHub in the loop.

Paths: `MP` = `/media/propdev/Expansion/openclaw/.openclaw/workspace/repos/mathphysics` (the canonical checkout). `PLAN` = `MP/docs/plans/yang-mills-mass-gap`.

## 1. Roles

| Role | Tier | Count | Job |
|---|---|---|---|
| Strategist | T4 | 1 | Owns `03-MASS-GAP-STRATEGY.md`, `BOARD.md` and the packets. Opens waves, adjudicates disputes, writes the weekly synthesis in `HANDOFF.md`. |
| Coordinator | T2–T3 | 1 per wave | Keeps the tracker tidy: checks claims have not expired, pairs each finished task with a reviewer, updates `BOARD.md` status. Makes no mathematical decisions. |
| Worker | per packet | many | Does exactly one packet. |
| Reviewer | ≥ packet tier | many | Checks exactly one deliverable. Never the author. Different model family from the author whenever one is available. |
| Red team | T3 | on demand | Tries to break a claimed result. Rewarded for findings, not for agreement. |
| Micro-worker | T0 | many | A small local model with no tools (planned; the runner `bin/ym-micro` is packet YM-T14 and does not exist yet). Receives one self-contained prompt from the runner, returns text in a fixed template. Its output is never evidence; it is raw material that a T2+ reviewer accepts or discards. |

Tiers are defined in `03-MASS-GAP-STRATEGY.md` §7 and mapped to Ernest's harnesses in §12 below. A packet's tier is a floor.

## 2. Where things live

| Thing | Home | Versioned |
|---|---|---|
| Task definition | `PLAN/tasks/<ID>.md` | yes |
| Queue and status summary | `PLAN/BOARD.md` | yes |
| Live conversation: claims, progress, questions, reviews | Local tracker, one work item per task (`PLAN/bin/ym`) | no (lives in `MP/.pdd/work_items/`) |
| Deliverable notes and digests | `PLAN/notes/<ID>-<slug>.md` | yes |
| Evidence: scripts, logs, receipts | `PLAN/ledger/evidence/<ID>/` | yes if small; large outputs stay outside git with a receipt inside |
| Claims we rely on | `PLAN/ledger/claims.md` | yes |
| Findings against the sources | `PLAN/ledger/findings/` | yes |
| Disputes | `PLAN/ledger/disputes/` | yes |
| Code | the owning crate under `MP/Rust/…`, through PDD | yes |
| Mathematics for keeps | LaTeX under `MP/documents/` once a result is reviewed | yes |

The tracker is the conversation. Files are the memory. If something matters after the task closes, it must be in a file.

## 3. The tracker

`PLAN/bin/ym` wraps Ernest's local PDD work items (`pdd work`), always against the canonical checkout so that every agent sees the same thread regardless of which worktree it runs in.

    ym list                          # all items
    ym show 7                        # item 7 with its comments
    ym say 7 <session> "<text>"      # add a comment as <session>
    ym new "<title>" <packet-file>   # strategist/coordinator only, canonical host only

From another machine on the local network the same commands work once `YM_HOST` is set (§10). Writers are serialized with a lock, so several agents may comment at the same moment.

Session names: `<harness>-<model>-<4 random chars>`, for example `codex-gpt56-k3pd`. One session, one packet.

Every comment starts with one keyword in capitals:

| Keyword | Meaning |
|---|---|
| `CLAIM until <UTC time ≤ 24h ahead>` | I am working on this. Include a two-sentence plan. |
| `RENEW until <UTC>` | Still working. Only before expiry, only by the holder. |
| `PROGRESS` | A fact others can use now. Not a diary. |
| `QUESTION @<role or session>` | I need an answer to continue. State what you will assume if nobody answers. |
| `ANSWER` | Reply to a question; quote it. |
| `FINDING` | Something that contradicts a source, the packet, or another agent's result. Include evidence path. |
| `SUBMIT` | Deliverable is at these paths, on this branch, at this commit. Lists what was and was not verified. |
| `REVIEW accept` / `REVIEW changes` / `REVIEW reject` | Reviewer's verdict with numbered reasons. |
| `DISPUTE` | Author and reviewer disagree after one round. Summarise both positions. Goes to the strategist. |
| `RELEASE` | I stop. Hand-off state in two lines. |

Claim rules:

1. Read the whole thread first. If there is an unexpired `CLAIM` by someone else, pick another task.
2. Post `CLAIM`. Then **read the thread again**. If an earlier `CLAIM` appeared, post `RELEASE` and leave. Earliest comment wins.
3. An expired claim is free. Read the previous holder's notes before taking it.
4. Never work on a packet without a posted claim. A branch name is not a claim.

## 4. Workspace

Each session works in its own git worktree and branch, so parallel agents never share a working directory:

    cd $MP
    git worktree add ../mathphysics-wt/<ID>-<session> -b ym/<ID>-<session> feat/yang-mills-mass-gap-program

An agent on another machine clones the canonical checkout instead and pushes its `ym/...` branch back (§10). Commit to that branch only. Never commit to `master`. Never push anywhere unless the packet says so. The coordinator or strategist merges accepted branches into `feat/yang-mills-mass-gap-program`. Remove the worktree after merge.

Read-only inputs: the four Yang–Mills repositories and `openai/math` under `repos/`. Do not write in them, do not run build tools inside them unless the packet is a build packet (YM-O01, YM-A01, YM-A02).

## 5. What a deliverable must contain

Every `SUBMIT` points to a note `PLAN/notes/<ID>-<slug>.md` with these sections, in this order:

1. **Result.** Three to ten lines. What is now known that was not.
2. **Labels.** For each statement in the result, one truth label from `README.md` §"Truth standard".
3. **Evidence.** Exact commands run and where their output is. Source citations as `repo@pin:path:line` or paper, section, page, equation.
4. **Not done / not checked.** Explicit. "Nothing" is almost never true.
5. **Questions for the strategist.** Optional.

Rules that reviewers enforce:

- A citation means the author opened the source at that place. Remembered references are marked `(from memory, unverified)` and cannot support a result.
- A name in the Lean code is never evidence of what a declaration means. Unfold the definition.
- A number is given with how it was computed. Two independent methods for anything another task will depend on.
- A failed or inconclusive attempt is a valid deliverable. Say so plainly. Do not dress it up.
- Do not fill a gap in an argument with a plausible sentence. Write `GAP:` and what is missing.

## 6. Review

- Every packet is reviewed before the board says `done`. Reviewer ≠ author, and a different model family when available.
- The reviewer reruns the evidence commands, opens at least three of the cited source locations chosen by themselves, and tries one thing the author did not.
- Verdict vocabulary: `accept`, `changes` (numbered, each one checkable), `reject` (the approach is wrong; say why).
- One round of changes, then either `accept` or `DISPUTE`. No silent rewriting of the author's work.
- Reviewers do not approve mathematics they did not follow. "I could not check step 4" is a correct review comment.
- For any packet in workstream C that claims a new lemma: two reviewers, and one of them is a red-team pass whose instruction is to find the error.

## 7. Talking to each other

- Cross-task questions go on the **other** task's thread as `QUESTION @<session>`, with a link back.
- If two tasks need the same definition or convention, the first to need it writes it into `PLAN/notes/conventions.md` and posts `PROGRESS` on both threads. Later tasks use it or open a `DISPUTE`.
- Two agents given the same packet on purpose (independent derivation) must not read each other's thread until both have posted `SUBMIT`. The packet says when this applies.
- Disagreement is information. A reviewer who cannot reproduce a result posts `FINDING`, not a quiet fix.

## 8. Stopping and handing off

Stop when the claim expires, the budget is spent, or you are blocked. Before stopping:

1. Commit what you have on your branch.
2. Write or update the note with a **Resume here** section: what is done, the next concrete step, what you would try next and why.
3. Post `RELEASE` with the branch name and commit.

The strategist updates `HANDOFF.md` at the end of every session. If you are the last one working and the strategist is not available, append a dated paragraph to `HANDOFF.md` §"Log" yourself.

## 9. Hard limits

- Nothing leaves the **local network** without Ernest's explicit approval: no pushes to outside remotes, issues, emails, uploads, posts. Talking to other agents on this machine or on the local network is allowed and expected (§10). One standing exception, approved 2026-11-15: contributions to `hairer/humanise-yang-mills` under §11.
- Local network addresses, host names, user names and credentials never appear in this repository, which is public. They live in the untracked file `MP/.pdd/ym-lan.env`.
- No credentials in prompts, notes or commits.
- No edits to the pinned source repositories.
- No claim of a proof of the mass gap, of confinement, or of any M3 statement outside the red-team process of `03-MASS-GAP-STRATEGY.md` §9.
- No new dependency, toolchain or large download unless the packet authorises it.
- `Data/Private/` paths never appear in this repository.

## 10. Agents on the local network

Ernest runs agents on several machines of one local network. All of them work on one tracker and one canonical repository.

**Canonical host.** The machine that holds `MP`. It holds the tracker (`MP/.pdd/work_items/`), the pinned sources, the Lean image and the compiled build.

**Setup on another machine (once).** Ernest provides SSH access to the canonical host. Then:

    git clone ssh://<user>@<canonical-host>/<absolute path of MP> mathphysics
    cd mathphysics && git checkout feat/yang-mills-mass-gap-program
    mkdir -p .pdd && cat > .pdd/ym-lan.env <<'EOT'
    YM_HOST=<user>@<canonical-host>
    YM_REMOTE_YM=<absolute path of MP>/docs/plans/yang-mills-mass-gap/bin/ym
    EOT
    docs/plans/yang-mills-mass-gap/bin/ym list          # must print the same list as on the canonical host

`.pdd/` is untracked, so the address stays out of git. Sources, if the packet needs them: `bin/ym-sources clone <dir>` (from GitHub) or copy from the canonical host.

**Working.** Same protocol as §3–§8 with two differences:

1. Branch from the remote: `git checkout -b ym/<ID>-<session> origin/feat/yang-mills-mass-gap-program`.
2. `SUBMIT` requires the branch to be on the canonical host: `git push origin ym/<ID>-<session>`. Push only `ym/*` branches. Never push `master` or `feat/yang-mills-mass-gap-program`; git on the canonical host refuses a push to its checked-out branch, and that refusal is intended.

**Talking.** Agent-to-agent conversation is the tracker thread: `QUESTION @<session>` and `ANSWER`, on the thread of the task the question is about. It works identically across machines. No side channels: a decision that is not on a thread or in a file did not happen.

**Model endpoints.** A llama.cpp server binds `127.0.0.1` by default. To let other machines use it, Ernest changes `server.host` in the profile; agents do not change bind addresses or firewall settings themselves. Endpoint URLs go in the worker's own untracked configuration.

**What must stay on the canonical host.** `ym new`, merges into the program branch, Lean builds, edits of `BOARD.md` and `HANDOFF.md` by the coordinator or strategist.

## 11. Contributing to Hairer's repository

Approved by Ernest on 2026-11-15 for an agent that he has given GitHub access. It is the only outward channel.

1. A contribution starts as one of our reviewed deliverables (status `done` on the board, reviewed by a different model family). Nothing unreviewed goes upstream.
2. Work in a **separate writable clone** (`REPOS/yang-mills-work/humanise-yang-mills`), never in the pinned `HAIRER` checkout. One branch per contribution.
3. Follow their rules, not ours: `HAIRER/docs/coordination.md` (task issue, self-claim, expiry, pull request). Read `HAIRER/AGENTS.md` and `CONTRIBUTING.md` at their current upstream revision first.
4. Every upstream comment and pull request says that it was written by an AI agent working for Ernest, names the model, and cites the Lean pin.
5. No credentials in prompts or files. If the agent's GitHub login is not working, it stops and says so. It does not run `gh auth login`.
6. Record each upstream action (URL, date, session) in `ledger/upstream.md`.

State on 2026-11-15: the `gh` token on the canonical host is invalid, so no agent here can post yet. No contribution has been made.

## 12. Harnesses and tiers on this setup

| Tier | Typical model | Harness on this machine | Non-interactive form | Use for |
|---|---|---|---|---|
| T4 | Fable-class | Claude Code | interactive | Strategy, packets, disputes, the key arguments |
| T3 | Opus-class; strongest Codex or Grok model | Claude Code, Codex, Grok Build | `claude -p`, `codex exec` | Digests of hard sources, derivations, red team |
| T2 | Sonnet-class; Codex | Claude Code, Codex, Pi, Hermes Agent | `claude -p`, `codex exec`, `pi -p`, `hermes -z` | Rust units with tests, Monte Carlo, reviews, tables |
| T1 | Haiku-class; Qwythos or Qwen 9B **with tools** through Pi or Hermes | Pi, Hermes Agent against the local llama.cpp server | `pi -p --model <local alias>` | Builds, reruns, filling fixed templates, one Rust test at a time |
| T0 | Qwen3.5 9B **without tools** | `bin/ym-micro` (plain HTTP to the local server; **planned**, packet YM-T14) | `ym-micro run …` | Micro-tasks: one snippet in, one template out |

Facts about the local models that shape the packets (from `claw-work/Deployments/Scripts/llama-cpp-server`): context 32,768 tokens for the Qwen3.5 profile, one server slot (requests are serial), about 40 tokens per second on the RTX 3060. Therefore:

- A T1 or T0 task must fit in about 6,000 tokens of input and 1,500 of output. The packet or the runner supplies the material; the model does not browse the repository.
- A T1 or T0 task must have an acceptance check that a script can run (a test passes, a file has the required sections, quoted identifiers occur in the quoted snippet).
- Local-model output is labelled `unreviewed` and cannot support a claim in `ledger/claims.md` until a T2+ reviewer accepts it.
- One local server means local workers run one at a time. Parallelism comes from the cloud harnesses and from other machines.

Launch recipes: `prompts/launch.md`.
