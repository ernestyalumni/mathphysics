# DELEGATE.md — Tasks for Claude Code / Codex / Other Coding Agents

**Last updated:** 2026-04-15
**Context:** Preparing for Wolfram conference talk (April 17) + ongoing MHV amplitudes research.
**Branch convention:** Create a feature branch from `master`. Never push directly to `master`.

Read `amplitudes/TALK.md`, `amplitudes/TASKS.md`, and `amplitudes/DISCOVERY-PLAN.md` for full context.

---

## How to Pick Up a Task

1. Read this file top-to-bottom
2. Pick a task marked 🟡 AVAILABLE
3. Create a branch: `git checkout -b feat/<task-slug>`
4. Do the work
5. Commit with descriptive messages
6. Push the branch — Ernest merges manually
7. Update this file: change 🟡 to ✅ and add date

---

## CODING TASK 1: Verify All 8 Scripts Pass (P0 — Conference Critical)
**Status:** 🟡 AVAILABLE
**Difficulty:** Easy
**Time estimate:** 30 min
**Branch:** `feat/verify-scripts`

Run each script in `amplitudes/scripts/` (01 through 08) and capture output.
These are standalone Python scripts requiring `numpy` and `sympy` (no Cadabra2 needed).

```bash
cd amplitudes/scripts
for f in 0*.py; do
    echo "=== $f ==="
    python3 "$f" 2>&1
    echo ""
done
```

**Deliverables:**
- Fix any script that fails (document what you changed and why)
- Create `amplitudes/talk/script_outputs/` directory
- Save each script's stdout to `amplitudes/talk/script_outputs/01_output.txt`, etc.
- Create `amplitudes/talk/script_outputs/SUMMARY.md` with pass/fail status per script

**Notes:**
- Scripts 01–07 are self-contained verification (color algebra, spinor math, etc.)
- Script 08 reproduces SMGA paper Table 1 — compare output against SMGA Eq.(17)–(20)
- If a script has import errors, check if it needs `sys.path` adjustment for the `MHVamplitudes` package

---

## CODING TASK 2: Build Slide Deck (P0 — Conference Critical)
**Status:** 🟡 AVAILABLE
**Difficulty:** Medium
**Time estimate:** 2–3 hours
**Branch:** `feat/slide-deck`

Create a presentation following the outline in `amplitudes/TALK.md`.

**Option A: LaTeX Beamer** (preferred for academic audience)
- Create `amplitudes/talk/slides.tex` using beamer
- Use the existing `00-preamble.tex` macros for physics notation
- Compile to `amplitudes/talk/slides.pdf`

**Option B: reveal.js / HTML** (preferred for VC audience)
- Create `amplitudes/talk/slides.html`
- Single self-contained HTML file with embedded CSS

**Required slides (minimum):**

1. Title: "OpenClaw: Agent-Driven Exploration of MHV Scattering Amplitudes"
2. The problem: researchers drown in notation, literature, routine computation
3. What is OpenClaw (bullet points — see `amplitudes/TASKS.md` TASK-004 spec)
4. Agent workflow: read paper → derive equations → write verification script → document
5. Key equation: Parke-Taylor formula (use LaTeX)
   `A_n(1^-,2^-;3^+...n^+) = \langle 12\rangle^4 / (\langle 12\rangle\langle 23\rangle\cdots\langle n1\rangle)`
6. Key equation: SMGA Eq.(16) stripped amplitude formula
   `A_{1...n}|_{R_1} = (1/2^{n-2}) \prod_{m=2}^{n-1} (\text{sg}_{m,m+1} + \text{sg}_{1,2...m})`
7. What we built: repo structure (show file tree of amplitudes/)
8. Dixon–Srednicki bridge: notation translation example
9. Verification: 8/8 scripts passing (show captured output from TASK 1)
10. SMGA reproduction: comparison table (from script 08 output)
11. Berends-Giele recursion: standard BG vs SMGA BG (diagram or side-by-side)
12. What works / what doesn't (honest assessment)
13. Future directions: NMHV, gravitons, loop-level, Cadabra2 package
14. Links: github.com/openclaw/openclaw, docs.openclaw.ai

**For the VC/general version** (separate file `amplitudes/talk/pitch-slides.html` or `.tex`):
- Add: "Why this matters" slide (AI doing real science, not just chatting)
- Add: Market opportunity slide (every research lab needs this)
- Add: Team/traction slide
- Remove: heavy equations (keep Parke-Taylor as the "wow" equation, skip SMGA details)
- Add: "What OpenAI did vs what we did" comparison

**Deliverables:**
- `amplitudes/talk/slides.tex` + `amplitudes/talk/slides.pdf` (academic)
- `amplitudes/talk/pitch-slides.tex` or `.html` (VC/general) — optional but valuable

---

## CODING TASK 3: OpenClaw Overview for Academics (P0)
**Status:** 🟡 AVAILABLE
**Difficulty:** Easy
**Time estimate:** 30 min
**Branch:** `feat/openclaw-overview`

Write `amplitudes/talk/openclaw-overview.md` — a 1-page summary of OpenClaw for physics/math academics.

**Must cover:**
- What it does: persistent AI assistant with shell access, file system, web, memory, sub-agents
- How it differs from ChatGPT (no persistence, no tools, no agents), Copilot (code-only), Wolfram Alpha (no autonomy)
- Key features: reads papers (LaTeX source, not just PDFs), writes derivations, runs verification scripts, maintains reading lists and task boards, spawns sub-agents for parallel work
- Open-source, self-hosted on your own hardware (Mac mini, Linux box, VPS)
- Links: https://docs.openclaw.ai, https://github.com/openclaw/openclaw, https://discord.com/invite/clawd

**Tone:** Accessible to a physicist who has used ChatGPT but doesn't know what an "agent framework" is.

---

## CODING TASK 4: Notation Rosetta Stone (P1)
**Status:** 🟡 AVAILABLE
**Difficulty:** Medium (requires reading source LaTeX)
**Time estimate:** 1–2 hours
**Branch:** `feat/notation-rosetta`

Create `amplitudes/talk/notation-rosetta.md` — a clean comparison table of conventions.

**Read these source files first:**
- Dixon: look for `amplitudes/09-dixon-srednicki-bridge.tex` existing content
- Srednicki: `amplitudes/Srednicki/individual_chapters/ch50_*.tex` and `ch60_*.tex`
- SMGA paper: the docstring in `amplitudes/MHVamplitudes/src/mhvamplitudes/amplitudes/berends_giele.py` has an excellent convention comparison — use that

**Table columns:** Srednicki | Dixon (1310.5353) | SMGA (2602.12176) | Elvang-Huang (1308.1697)

**Rows to include:**
- Metric signature
- Spinor notation (λ vs λ̃, dotted vs undotted, which is "left" vs "right")
- Angle bracket ⟨ij⟩ definition
- Square bracket [ij] definition
- Trace normalization: Tr(TᵃTᵇ) = ½δ^{ab} vs δ^{ab}
- Structure constant convention: [Tᵃ,Tᵇ] = if^{abc}Tᶜ vs i√2 f^{abc}Tᶜ
- Momentum decomposition: p^{αα̇} = λ^α λ̃^{α̇} (which index is which)
- Polarization vector conventions
- Color-ordered amplitude notation

**Deliverable:** `amplitudes/talk/notation-rosetta.md` (markdown table)

---

## CODING TASK 5: Extend SMGA Reproduction to n=7,8,9,10 (P3 — Post-Conference)
**Status:** 🟡 AVAILABLE
**Difficulty:** Medium
**Time estimate:** 2–3 hours
**Branch:** `feat/smga-higher-n`

Extend `amplitudes/scripts/08_smga_stripped_amplitudes.py` to verify SMGA Eq.(16) for higher multiplicity.

**Steps:**
1. Read `amplitudes/MHVamplitudes/src/mhvamplitudes/amplitudes/berends_giele.py` to understand the recursion
2. Read `amplitudes/MHVamplitudes/src/mhvamplitudes/kinematics/` for the spinor-helicity infrastructure
3. Extend script 08 to compute A_{1234567}, A_{12345678}, A_{123456789}, A_{1234567890}
4. For each n, pick several random kinematic configurations in region R₁
5. Compare: (a) direct formula Eq.(16) vs (b) Berends-Giele recursion output
6. Report agreement/disagreement and precision

**Deliverables:**
- Updated `amplitudes/scripts/08_smga_stripped_amplitudes.py` (or new `09_smga_higher_n.py`)
- Print clear comparison table for each n
- Note: if n≥9 is slow, document timings

---

## CODING TASK 6: pytest Suite Verification + CI Prep (P2)
**Status:** 🟡 AVAILABLE
**Difficulty:** Medium
**Time estimate:** 1–2 hours
**Branch:** `feat/pytest-verify`

Run the full pytest suite inside the Cadabra2 Docker container and fix any failures.

**Steps:**
1. Build Docker image if needed (see `repos/Monoclaw/Deployments/DockerBuilds/Physics/Cadabra2/`)
2. Run container with mathphysics mounted:
   ```
   MATHPHYSICS_DIR=/path/to/mathphysics ./run.sh cli
   ```
3. Inside container:
   ```
   cd /mathphysics/amplitudes/MHVamplitudes
   pip install -e .
   pytest -s ./tests/
   ```
4. Fix any failures
5. Document results in `amplitudes/MHVamplitudes/TESTING.md` (append a "Last Run" section)

**Bonus:** Create `.github/workflows/test.yml` for GitHub Actions CI (run pytest on push/PR using the Cadabra2 Docker image).

---

## CODING TASK 7: Berends-Giele Explainer Section (P1)
**Status:** 🟡 AVAILABLE
**Difficulty:** Hard (requires physics understanding)
**Time estimate:** 3–4 hours
**Branch:** `feat/bg-explainer`

Add a new section to `amplitudes/05-bcfw.tex` (or create `amplitudes/10-berends-giele.tex`) explaining the Berends-Giele recursion.

**Source material:**
- The docstring in `amplitudes/MHVamplitudes/src/mhvamplitudes/amplitudes/berends_giele.py` — this is excellent and should be the basis
- SMGA paper Appendix A (equations A7, A17)
- TASKS.md TASK-008 spec

**Must include:**
1. Standard BG recursion (off-shell currents, 3/4-gluon vertices from Lagrangian)
2. How SMGA modifies it for (2,2) signature + half-collinear limit
3. Why propagators cancel (piecewise-constant result)
4. The vertex functions V and V̄ (SMGA eq. 7) — what they mean physically
5. Side-by-side: "Standard BG step" vs "SMGA BG step"
6. Connection to BCFW (both are recursions, different inputs)

**If creating a new .tex file**, update `amplitudes/99-master.tex` to include it.

**Deliverable:** New/updated .tex file with the BG section + update to 99-master.tex

---

## CODING TASK 8: LaTeX Compile Clean-Up (P0)
**Status:** 🟡 AVAILABLE
**Difficulty:** Easy
**Time estimate:** 30 min
**Branch:** `feat/latex-cleanup`

Clean compile all LaTeX documents and fix any errors.

```bash
cd amplitudes

# Remove stale build artifacts
rm -f *.aux *.log *.toc *.out *.fls *.fdb_latexmk

# Compile master document (twice for ToC)
pdflatex 99-master.tex && pdflatex 99-master.tex

# Compile Srednicki companion
cd Srednicki
rm -f *.aux *.log *.toc *.out
pdflatex Srednicki_Companion_Full.tex && pdflatex Srednicki_Companion_Full.tex
```

**Deliverables:**
- Clean `99-master.pdf` and `Srednicki_Companion_Full.pdf` (don't commit PDFs — just confirm they build)
- Fix any LaTeX errors (missing packages, undefined references, etc.)
- Document page counts and any warnings in a commit message

---

## CODING TASK 9: Interactive Spinor-Helicity Demo (P2 — Nice to Have)
**Status:** 🟡 AVAILABLE
**Difficulty:** Medium-Hard
**Time estimate:** 3–4 hours
**Branch:** `feat/interactive-demo`

Create a single-file HTML page (`amplitudes/talk/demo.html`) that lets someone enter momenta and see spinor-helicity calculations live.

**Features:**
- Input: 4-momenta p₁,...,pₙ (or generate random massless momenta)
- Compute: spinor brackets ⟨ij⟩, [ij], Mandelstam s_{ij}
- Compute: Parke-Taylor MHV amplitude for chosen helicity config
- Compute: SMGA stripped amplitude via Eq.(16) for single-minus config
- Display: formatted equations + numerical results
- All in-browser JavaScript, no server needed

**Use:** numpy-equivalent JS library or just raw array math (2x2 matrices, nothing heavy).

---

## Notes for Agents

- **Python version:** 3.10+ (scripts use modern syntax)
- **Dependencies for scripts:** `numpy`, `sympy` (both pip-installable)
- **Cadabra2 scripts** (in `amplitudes/Srednicki/`) require the Docker container — see `amplitudes/MHVamplitudes/TESTING.md`
- **LaTeX:** standard `pdflatex` with common packages (amsmath, amssymb, hyperref, etc.)
- **Don't commit:** `.aux`, `.log`, `.toc`, `.out`, `.fls`, `.fdb_latexmk`, `.pdf`, `__pycache__/`, `.pyc`
- **Do commit:** `.py`, `.tex`, `.md`, `.html`, `.css`, `.js`

---

## Priority Summary

| # | Task | Priority | Time | Good For |
|---|------|----------|------|----------|
| 1 | Verify scripts | P0 | 30m | Codex |
| 2 | Slide deck | P0 | 2-3h | Claude Code |
| 3 | OpenClaw overview | P0 | 30m | Claude Code or Codex |
| 4 | Notation Rosetta | P1 | 1-2h | Claude Code |
| 5 | SMGA higher-n | P3 | 2-3h | Codex |
| 6 | pytest + CI | P2 | 1-2h | Codex |
| 7 | BG explainer LaTeX | P1 | 3-4h | Claude Code |
| 8 | LaTeX compile | P0 | 30m | Codex |
| 9 | Interactive demo | P2 | 3-4h | Codex |
