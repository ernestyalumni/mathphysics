# TMP OpenClaw Persistent Agent

Last updated: 2026-06-16
OpenClaw-native counterpart of the Hermes `tmpai` profile documented in `TMP-HERMES-AGENT.md` and `TMP-HERMES-BOOTSTRAP.md` (same directory).

Where Hermes uses **profiles** under `~/.hermes/profiles/tmpai/`, OpenClaw uses a registered **agent** in `~/.openclaw/openclaw.json` under `agents.list`. This doc records how to wire that agent so it can be recreated on another machine.

## What this is

A persistent OpenClaw agent — NOT an ephemeral subagent spawned for one task. It has its own workspace (persona + memory) and its own state / session store, so context persists across sessions.

- **Agent id:** `tmpai`
- **Persona name:** TMPAI
- **Emoji:** ⚛️
- **Workspace:** `~/.openclaw/workspace-tmpai` (created on setup — `SOUL.md` / `IDENTITY.md` / `AGENTS.md` / `USER.md` / `MEMORY.md` auto-injected as bootstrap)
- **Agent dir (state / auth / sessions):** `~/.openclaw/agents/tmpai/agent`
- **Model:** inherits `agents.defaults` from `openclaw.json`
- **Primary context repo:** `mathphysics` (`/home/propdev/.openclaw/workspace/workspace2/repos/mathphysics` on Linux) — read by absolute path each session

## Config to add to `~/.openclaw/openclaw.json`

Under `agents.list` (alongside `main`, `marketer`, `spacexai`):

```json
{
  "id": "tmpai",
  "name": "tmpai",
  "workspace": "/home/propdev/.openclaw/workspace-tmpai",
  "agentDir": "/home/propdev/.openclaw/agents/tmpai/agent",
  "identity": {
    "name": "TMPAI",
    "theme": "Theoretical & Mathematical Physics discovery engine",
    "emoji": "⚛️"
  }
}
```

And add `tmpai` to the `main` agent's spawn allowlist so it can be reached:

```json
{
  "id": "main",
  "subagents": { "allowAgents": ["marketer", "spacexai", "tmpai"] }
}
```

## Setup on this (or another) machine

```bash
# 1. Register the agent (or hand-edit agents.list with the JSON above).
openclaw agents add tmpai \
  --workspace ~/.openclaw/workspace-tmpai \
  --agent-dir ~/.openclaw/agents/tmpai/agent \
  --non-interactive

# 2. Seed the workspace markdown files. The full content for each of these
#    lives under `mathphysics/agent/` in this repo. Copy them in:
mkdir -p ~/.openclaw/workspace-tmpai/{avatars,memory}
MATHPHYSICS="$HOME/.openclaw/workspace/workspace2/repos/mathphysics"

cp "$MATHPHYSICS/agent/SOUL.md"             ~/.openclaw/workspace-tmpai/SOUL.md
cp "$MATHPHYSICS/agent/IDENTITY.md"         ~/.openclaw/workspace-tmpai/IDENTITY.md
cp "$MATHPHYSICS/agent/AGENTS.md"           ~/.openclaw/workspace-tmpai/AGENTS.md
cp "$MATHPHYSICS/agent/USER.md"             ~/.openclaw/workspace-tmpai/USER.md
cp "$MATHPHYSICS/agent/MEMORY.md"           ~/.openclaw/workspace-tmpai/MEMORY.md
cp "$MATHPHYSICS/agent/HEARTBEAT.md"        ~/.openclaw/workspace-tmpai/HEARTBEAT.md
cp "$MATHPHYSICS/agent/TOOLS.md"            ~/.openclaw/workspace-tmpai/TOOLS.md
cp "$MATHPHYSICS/agent/RESEARCH-AGENDA.md"  ~/.openclaw/workspace-tmpai/RESEARCH-AGENDA.md
cp "$MATHPHYSICS/agent/DISCOVERY-LOOPS.md"  ~/.openclaw/workspace-tmpai/DISCOVERY-LOOPS.md
cp "$MATHPHYSICS/agent/EVAL-AND-TOOLING.md" ~/.openclaw/workspace-tmpai/EVAL-AND-TOOLING.md
cp "$MATHPHYSICS/agent/MARKETING.md"        ~/.openclaw/workspace-tmpai/MARKETING.md

# 3. Refresh agent identity from IDENTITY.md.
openclaw agents set-identity --agent tmpai --from-identity

# 4. Restart the OpenClaw gateway.
openclaw gateway restart

# 5. Verify.
openclaw agents list   # should show `tmpai` with ⚛️ identity and its own workspace
```

OAuth note: a secondary OpenClaw agent reads through to `main`'s auth profiles by default. For an independent account, sign in from the `tmpai` agent session.

## How to talk to TMPAI (OpenClaw)

- **Spawn from `main`:** use `sessions_spawn` with `agentId: "tmpai"` (enabled by the `allowAgents` entry above). Reuse the same session id to keep it persistent.
- **Bind a channel (optional):** create a dedicated Telegram bot account and add a binding `{ agentId: "tmpai", match: { channel: "telegram", accountId: "tmpai" } }` so DMs to that bot route straight to TMPAI.

## First message to TMPAI

```text
You are TMPAI. Orient per AGENTS.md: read your workspace SOUL/IDENTITY/USER/MEMORY and the mathphysics repo's TMP-HERMES-AGENT.md, TMP-MISSION.md, RESEARCH-PLAN.md, TASKS.md, and recent `agent/memory/` notes. Run `cd /home/propdev/.openclaw/workspace/workspace2/repos/mathphysics && git status && git branch && git log -10 --oneline`. Confirm readiness and propose today's highest-leverage research action with a named verifiable artifact (Cadabra2 script, LaTeX section, eval row, or video script). Do not publish anything externally without my approval.
```

## Red lines

- Drafting research / scripts / posts / derivations = allowed.
- **Publishing, posting, submitting, or pushing to master = explicit approval only.**
- No fabricated citations, theorem names, or numerical results.
- Treat Galvatron `Documents/Goals/` and any `Data/Private/` as private; never reference from this repo.
- Never commit build artifacts.
- Never commit directly to `master`/`main` in any repo.

## Relationship to the Hermes profile

The two systems share the same persona, mission, and safety rules. They differ only in:

| | Hermes profile | OpenClaw agent |
|--|----------------|----------------|
| Persona name | TMPAI | TMPAI |
| Profile / agent id | `tmpai` | `tmpai` |
| State home | `~/.hermes/profiles/tmpai/` | `~/.openclaw/agents/tmpai/agent/` |
| Workspace | (none — Hermes injects the repo via `terminal.cwd`) | `~/.openclaw/workspace-tmpai/` |
| Setup doc | `TMP-HERMES-BOOTSTRAP.md` | this file |
| Identity injection | `<profile_home>/SOUL.md` | `<workspace>/SOUL.md` + `IDENTITY.md` + `AGENTS.md` (auto) |

Use whichever framework matches the machine. Both should be safe to run side by side — they don't share state.
