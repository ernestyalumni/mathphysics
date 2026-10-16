# AGENTS.md - mathphysics repo (TMPAI workspace)

This file is the entry-point for AI agents working in this repo. It's auto-loaded by Hermes (when `terminal.cwd` is here), Claude Code, Codex, and other AGENTS.md-aware tools.

This repo is the persistent context for **TMPAI** (⚛️) — the Theoretical & Mathematical Physics discovery agent. The repo itself holds research artifacts (LaTeX, Cadabra2, Jupyter notebooks); the `agent/` subdirectory holds the persona, memory, and research-engine bookkeeping.

Repo root: `/home/propdev/.openclaw/workspace/workspace2/repos/mathphysics`
Persona: TMPAI ⚛️
Hermes profile: `tmpai` (see `TMP-HERMES-BOOTSTRAP.md`)
OpenClaw agent: `tmpai` (see `TMP-OPENCLAW-AGENT.md`)
Claude Code: see `TMP-CLAUDE-CODE.md`
Codex CLI: see `TMP-CODEX.md`

## Every Session (orientation)

Before doing anything else:

1. Read the repo-root identity & mission docs:
   - `TMP-HERMES-AGENT.md` — TMPAI identity & ops (canonical)
   - `TMP-MISSION.md` — discovery-engine charter
   - `RESEARCH-PLAN.md` — active paper plan
   - `TASKS.md` — immediate priorities
2. Read the persona files in `agent/`:
   - `agent/SOUL.md`, `agent/IDENTITY.md`, `agent/USER.md`, `agent/MEMORY.md`
3. Read the research-engine docs in `agent/`:
   - `agent/RESEARCH-AGENDA.md` — multi-quarter agenda
   - `agent/DISCOVERY-LOOPS.md` — registered AI loops
   - `agent/EVAL-AND-TOOLING.md` — registered evals
   - `agent/MARKETING.md` — content calendar & drafts
4. Read recent `agent/memory/YYYY-MM-DD.md` notes (today + yesterday) if present.
5. Run `git status && git branch && git log -10 --oneline`.

Don't dump large `.tex` bodies into chat. Use `read_file` with `offset` / `limit`; cite by `path:line`.

## Standard Workflows

- **Daily research sprint** — review active branch + `TASKS.md`, give 1–3 highest-leverage actions, each pinned to a verifiable artifact (Cadabra2 script, LaTeX section, eval row, video draft). End with one concrete action started or staged.
- **Derive + verify** — show derivation step-by-step; produce a Cadabra2 script that verifies; update the relevant LaTeX section; commit on a `feat/` branch.
- **Discovery-loop run** — look up the loop in `agent/DISCOVERY-LOOPS.md`, run proposer → critic → verifier → logger, persist result to `agent/discoveries/YYYY-MM-DD-<slug>.md`.
- **Eval build** — `evals/<name>/{README,problems.jsonl,verifier.py,harness.py}` at repo root; 10 problems; symbolic ground truth; pinned baseline model.
- **Short-form draft** — `agent/marketing/drafts/YYYY-MM-DD-<slug>.md` with 60-90s script + visual beats + closing hook; mark `STATUS: DRAFT — awaiting approval`.

## OCR / Document Parsing for Math & Physics Sources

Whenever a task involves parsing, OCR'ing, or extracting text/equations from a document (PDF, `.djvu`, `.ps` / PostScript, scanned papers, textbooks) — or improving/debugging/extending that capability — use the Monoclaw `ocr-compare` pipeline as the single source of truth:

`/media/propdev/Expansion/openclaw/.openclaw/workspace/repos/Monoclaw/Python/ocr-compare/`

Before doing anything, read both files in that directory; they override assumptions:

- `README.md` — workflow and commands.
- `AGENTS.md` — hard assumptions, conventions, load-bearing pins.

Pipeline summary: Marker is the authoritative backbone with full page coverage; Nougat is a per-equation cross-check with cleaner arXiv-style LaTeX but can drop/repeat pages, so it is untrusted by default. `combine.py` reconciles by `\tag{N}` and writes `reconciled/equations.json`; `resolve.py` is the optional vision-judge tier (Grok, or Claude Code via manifest).

Hard rules:

- This box has an NVIDIA GPU. Scripts install the CUDA torch build only; do not add CPU fallbacks.
- Machine-local paths/GPU settings live in `scripts/env.sh` / `scripts/_paths.sh` (`OCR_STORAGE`, `CUDA_VISIBLE_DEVICES=0`). Never hardcode paths; change locations via env vars only.
- Never commit venvs, weights, `out/`, or OCR-generated bulk outputs.
- Nougat pins are load-bearing: `transformers==4.38.2`, `albumentations==1.3.1`, `pypdfium2==4.30.0`. Do not bump without re-validating.
- Big PDFs/books OOM under `marker_single`; use `scripts/run_ocr_large.sh`, not `scripts/run_ocr.sh`.

Input normalization: the pipeline is PDF-in only. For non-PDF inputs, convert to PDF first, then run the pipeline:

- `.djvu` → `ddjvu -format=pdf in.djvu out.pdf` (or `djvups | ps2pdf`).
- `.ps` → `ps2pdf in.ps out.pdf`.
- images → `img2pdf` or ImageMagick `convert`.

If the converter is missing, install it (`apt`) before proceeding or report the blocker; never feed non-PDF input directly into the pipeline. The venvs and weights are expected to already exist at `OCR_STORAGE`; do not rerun `scripts/setup_envs.sh` unless they are missing or broken.

On parse completion, inspect `reconciled/equations.json` and report page coverage, Nougat dropped-page ranges, conflict count, and output locations (`marker` markdown, Nougat `.mmd`, `reconciled/merged.md`). For books, run `scripts/triage_conflicts.py` and separate cosmetic auto-resolutions from genuine symbol disagreements before escalating to a vision judge.

## Persistent Artifacts

Create / update under `agent/` (persona-local) or repo root (public-facing) as needed:

Under `agent/`:
- `agent/RESEARCH-AGENDA.md` — multi-quarter agenda, ranked threads.
- `agent/DISCOVERY-LOOPS.md` — registered loops with prompts, results, runs.
- `agent/EVAL-AND-TOOLING.md` — registered evals, harnesses, baselines.
- `agent/MARKETING.md` — content calendar, drafts, published index.
- `agent/discoveries/YYYY-MM-DD-<slug>.md` — one file per loop-discovered or verified result.
- `agent/marketing/drafts/YYYY-MM-DD-<slug>.md` — short-form / long-form drafts.
- `agent/memory/YYYY-MM-DD.md` — daily raw notes.
- `agent/MEMORY.md` — curated long-term wisdom (main session only).
- `agent/HEARTBEAT.md` — optional periodic checklist.

At repo root:
- `evals/<name>/` — eval directories (public, when promoted out of drafts).
- `RESEARCH-PLAN.md` — current active paper plan.
- `TASKS.md` — current immediate priorities.

Don't pre-create empty directories or files until they're needed.

## Memory Discipline

- Daily raw notes → `agent/memory/YYYY-MM-DD.md` (create `agent/memory/` on first use).
- Curated long-term wisdom → `agent/MEMORY.md` (main session only; do NOT load in shared contexts).
- Every discovery loop run is logged whether it passed or failed — failure logs are themselves valuable data.
- Periodically distill daily notes into `agent/MEMORY.md`.
- **Write it down.** Mental notes don't survive session restarts; files do.

## Git Rules

- **Never commit / push to `master`** — Ernest merges manually. Work on feature branches (`feat/`, `fix/`, `chore/`, `experiment/`).
- **Never commit build artifacts.** That includes:
  - LaTeX: `.aux`, `.log`, `.toc`, `.out`, `.synctex.gz`, `.fls`, `.fdb_latexmk`, intermediate PDFs.
  - C++ / CMake: `build/`, `_build/`, `cmake-build-*/`, `CMakeCache.txt`, `CMakeFiles/`, `Makefile`, `*.o`, `*.a`, `*.so`, compiled binaries.
  - Generated outputs: CSVs, PNGs, logs produced by running scripts.
  - Python: `__pycache__/`, `*.pyc`, `.ipynb_checkpoints/`.
- Don't compensate with `.gitignore` — just don't stage them in the first place.
- Cadabra2 script + LaTeX update for the same result should commit together when they form one logical change.

## Red Lines

- Drafting derivations, scripts, posts, evals = allowed.
- **Publishing, posting, submitting, or pushing to master = explicit approval only.**
- No fabricated citations, theorem names / attributions, or numerical results.
- Treat Galvatron `Documents/Goals/` and any `Data/Private/` as private — never reference from this repo.
- Re-derivations are valuable; do NOT claim originality on known results.
- A "discovery loop" without a verifier step is not a loop, just a brainstorm — do not log it under `agent/discoveries/`.

Document significant decisions here (or in `agent/MEMORY.md`) so future runs have continuity.
