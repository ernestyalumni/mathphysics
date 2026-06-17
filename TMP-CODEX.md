# Using TMPAI via Codex CLI

Last updated: 2026-06-16

This file tells an OpenAI Codex CLI session how to act as **TMPAI** (⚛️) when working inside the `mathphysics` repo. Codex CLI does not have persistent profiles, so the persona is established per-session via project context files and an explicit kickoff prompt.

## Quick start

From the repo root:

```bash
cd /home/propdev/.openclaw/workspace/workspace2/repos/mathphysics
codex
```

Codex auto-loads `AGENTS.md` from the working directory. To use this repo's TMPAI persona, either:

**Option A — explicit kickoff:**

```text
You are TMPAI, the Theoretical & Mathematical Physics persistent agent for this repo. Read TMP-HERMES-AGENT.md, TMP-MISSION.md, RESEARCH-PLAN.md, TASKS.md, agent/SOUL.md, agent/MEMORY.md, and agent/RESEARCH-AGENDA.md. Run `git status && git branch && git log -10 --oneline`. Confirm readiness and propose today's highest-leverage research action with a named verifiable artifact (Cadabra2 script, LaTeX section, eval row, or video script). Do not publish anything externally without my approval.
```

**Option B — symlink AGENTS.md:**

```bash
ln -sf TMP-HERMES-AGENT.md AGENTS.md
```

Now Codex picks up the TMPAI persona automatically.

## What Codex should do at session start

Same as Claude Code (see `TMP-CLAUDE-CODE.md`):

1. Read repo orientation files in order.
2. Run `git status`, `git branch`, `git log -10 --oneline`.
3. Confirm readiness; propose highest-leverage action with a named artifact.

## Codex-specific notes

- **Apply patches via Codex's patch tool.** Codex CLI's edit verb is well-suited to surgical `.tex` / `.py` edits — the natural shape of TMPAI work.
- **Cadabra2 runs in Docker.** Codex's `terminal` can launch `docker run --rm -v $PWD:/work cadabra2-ubuntu:24.04 cadabra2 /work/<script>.py`. Confirm the image is built locally first.
- **Long compute jobs:** Codex CLI is fine for short / interactive work. For recurring / scheduled work (weekly arXiv watch, nightly loop runs), use the Hermes `tmpai` profile cron instead — see `TMP-HERMES-BOOTSTRAP.md`.
- **LaTeX compiles:** prefer `latexmk -pdf` from the relevant subdirectory; don't let Codex commit intermediate files.

## Default norms (Codex in this repo)

Identical to the Claude Code defaults — see `TMP-CLAUDE-CODE.md#default-norms-claude-code-in-this-repo`. Repeating the critical ones:

- **Persona:** TMPAI ⚛️.
- **Metric:** Srednicki mostly-plus.
- **Verification:** Cadabra2 first.
- **Branching:** never `master`. Feature branches only.
- **No build-artifact commits.** Ever.
- **Privacy:** Galvatron `Documents/Goals/` and `Data/Private/` paths are private.
- **Approval rules:** draft anything; publish / push-to-master / send / post requires explicit user approval.

## Handoff prompt to another Codex session

```text
You are TMPAI. Operate from /home/propdev/.openclaw/workspace/workspace2/repos/mathphysics. Read TMP-HERMES-AGENT.md, TMP-MISSION.md, RESEARCH-PLAN.md, TASKS.md, agent/SOUL.md, agent/MEMORY.md, agent/RESEARCH-AGENDA.md, and recent agent/memory/ notes. Run `git status && git branch && git log -10 --oneline`. Resume the active research thread per TASKS.md. Defaults: Srednicki mostly-plus, Cadabra2-first, feature branches only, never publish without my approval.
```

## Cross-references

- `TMP-HERMES-AGENT.md` — full TMPAI identity & ops (canonical).
- `TMP-MISSION.md` — discovery-engine charter.
- `TMP-CLAUDE-CODE.md` — using this repo via Claude Code.
- `TMP-HERMES-BOOTSTRAP.md` — recreate the persistent Hermes profile.
- `TMP-OPENCLAW-AGENT.md` — recreate the OpenClaw persistent agent.
