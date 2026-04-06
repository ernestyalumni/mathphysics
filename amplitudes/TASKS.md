# TASKS.md — Agent Task Board for Talk Prep + Research

**Talk:** "OpenClaw: Agent-Driven Exploration of MHV Scattering Amplitudes"  
**Deadline:** April 17, 2026 (conference starts)  
**Branch:** `feat/wolfram-talk-prep`  
**See also:** `TALK.md` (talk outline), `DISCOVERY-PLAN.md` (long-term research plan)

---

## Priority Levels

- **P0** — Must be done before the talk (April 17)
- **P1** — Should be done, significantly strengthens the talk
- **P2** — Nice to have, do if time permits
- **P3** — Post-talk / long-term research

---

## P0: Talk Essentials

### TASK-001: Slide Deck / Presentation Materials
- **What:** Create slides (PDF or HTML) for the talk
- **Format:** ~15–20 slides for 45-min, ~8–10 for 15-min
- **Content:** Follow structure in `TALK.md`
- **Must include:**
  - What is OpenClaw (architecture diagram or screenshot)
  - Agent workflow diagram (read paper → derive → verify → document)
  - Key equations: Parke-Taylor, SMGA Eq.(16), spinor-helicity basics
  - Screenshots of: repo structure, a verification script output, bridge document PDF
  - Honest "what works / what doesn't" slide
- **Output:** `amplitudes/talk/slides.pdf` or `amplitudes/talk/slides.html`
- **Status:** Not started
- **Depends on:** TASK-002, TASK-003

### TASK-002: Verify All Existing Scripts Pass
- **What:** Run all 8 scripts in `amplitudes/scripts/` and confirm they pass
- **Why:** Can't present verification results if scripts are broken
- **Steps:**
  1. `cd amplitudes/scripts && python3 01_color_fierz_identity.py`
  2. Repeat for 02–08
  3. Fix any failures
  4. Document output (capture stdout for slides)
- **Output:** All scripts green; captured output in `amplitudes/talk/script_outputs/`
- **Status:** Last run: unknown (commit message says "all 8 passing" but re-verify)

### TASK-003: Compile LaTeX Notes PDF
- **What:** Clean compile of `99-master.tex` and `Srednicki/Srednicki_Companion_Full.tex`
- **Why:** Need clean PDFs as talk artifacts and handouts
- **Steps:**
  1. Remove LaTeX intermediary files (`.aux`, `.log`, `.toc`, `.out`, `.fls`, `.fdb_latexmk`)
  2. `cd amplitudes && pdflatex 99-master.tex && pdflatex 99-master.tex` (twice for ToC)
  3. Same for Srednicki companion
  4. Verify no errors, check page count
- **Output:** Clean `99-master.pdf`, `Srednicki_Companion_Full.pdf`
- **Status:** PDFs exist but may be stale

### TASK-004: OpenClaw Overview for Academics
- **What:** Write a 1-page summary of what OpenClaw is, aimed at physics/math academics
- **Why:** Most attendees won't know what an "AI agent framework" is
- **Must cover:**
  - What it does (persistent AI assistant with tools, memory, sub-agents)
  - How it differs from ChatGPT / Copilot / Wolfram Alpha
  - Open-source, self-hosted
  - Link to docs: https://docs.openclaw.ai and https://github.com/openclaw/openclaw
- **Output:** `amplitudes/talk/openclaw-overview.md`
- **Status:** Not started

---

## P1: Strengthens the Talk

### TASK-005: Run SMGA Stripped Amplitude Reproduction
- **What:** Run `scripts/08_smga_stripped_amplitudes.py` and verify it reproduces SMGA Table 1 values
- **Why:** "We reproduced their results" is a strong statement for the talk
- **Steps:**
  1. Run the script, capture output
  2. Compare against SMGA paper Table 1 (stripped amplitudes A_{123} through A_{123456})
  3. If mismatches, debug and fix
  4. Document: which values match, to what precision
- **Output:** Verified reproduction with comparison table
- **Status:** Script exists, needs verification run
- **Depends on:** TASK-002

### TASK-006: Notation Rosetta Stone
- **What:** Create a clean comparison table: Srednicki ↔ Dixon ↔ SMGA ↔ Elvang-Huang notation
- **Why:** Notation translation is a key selling point of agent-assisted research
- **Must include:**
  - Spinor conventions (which is λ, which is λ̃, dotted vs undotted)
  - Metric signature
  - Trace normalization (Tr(TᵃTᵇ) = ½δ vs δ)
  - Bracket notation
  - Momentum conservation sign conventions
- **Output:** Table in `09-dixon-srednicki-bridge.tex` §1 (may already partially exist) + standalone `amplitudes/talk/notation-rosetta.md`
- **Status:** Partially done in bridge document

### TASK-007: Finish Dixon–Srednicki Bridge Sections 5–6
- **What:** Complete the bridge document (sections on polarization vectors and soft gluon limit)
- **Why:** Demonstrates depth of agent-driven derivation work
- **See:** `TASK-dixon-srednicki-bridge.md` for full spec
- **Status:** Sections 1–4 drafted; 5–6 need work

### TASK-008: Add Berends-Giele Recursion Section
- **What:** Add a section to the bridge doc or Ch.05 on Berends-Giele recursion
- **Why:** This is the backbone of the SMGA proof; understanding it is essential
- **Must include:**
  - Derivation of SMGA Eq.(A7) from first principles
  - Connection to BCFW
  - Python implementation (may already be in `08_smga_stripped_amplitudes.py`)
- **Output:** New section in bridge doc or `05-bcfw.tex`
- **Status:** Not started

---

## P2: Nice to Have

### TASK-009: Interactive Demo Page
- **What:** Simple HTML page showing spinor-helicity calculations live (like the Cosmos GNC HTML demo)
- **Why:** Impressive for a live demo during the talk
- **Could show:** Enter momenta → compute spinor brackets → compute Parke-Taylor → compare with Berends-Giele
- **Output:** `amplitudes/talk/demo.html`
- **Status:** Not started

### TASK-010: Cadabra2 Computation for Srednicki Ch.48 + Ch.50
- **What:** Complete the Cadabra2 computations for the spinor-helicity core chapters
- **Why:** Fills the gap between Ch.36 (done) and Ch.60 (done)
- **Output:** `amplitudes/Srednicki/individual_chapters/ch48_*.py`, `ch50_*.py`
- **Status:** Not started (Ch.34, 35, 36, 60 done)

### TASK-011: Comparison: Agent Workflow vs Wolfram Language
- **What:** Prepare a slide or talking point comparing OpenClaw agent workflow with Mathematica/Wolfram Language
- **Why:** Audience is Wolfram-adjacent; they'll ask
- **Key differences:**
  - Wolfram: symbolic CAS, deterministic, user drives every step
  - OpenClaw: autonomous agent, reads papers, writes code, maintains state
  - Complementary, not competing (agent could use Wolfram Language as a tool)
- **Output:** Talking points in `TALK.md` or `amplitudes/talk/wolfram-comparison.md`
- **Status:** Not started

---

## P3: Post-Talk / Long-Term

### TASK-012: Reproduce SMGA Results for n=7,8,9,10
- **What:** Extend `08_smga_stripped_amplitudes.py` to higher multiplicity
- **Why:** Validates SMGA Eq.(16) beyond what they showed in the paper
- **Status:** Depends on TASK-005

### TASK-013: Explore Multi-Minus Amplitudes (NMHV)
- **What:** Investigate whether similar simplifications exist for NMHV amplitudes in kinematic limits
- **See:** `DISCOVERY-PLAN.md` Phase 3
- **Status:** Not started

### TASK-014: Explore Graviton Amplitudes
- **What:** SMGA mentions direct generalization to gravity; investigate
- **Status:** Not started

### TASK-015: Write Up Results for arXiv
- **What:** If discoveries are made, prepare a paper
- **Status:** Blocked on TASK-012–014

---

## For Other OpenClaw Agents

If you're picking up work from this task board:

1. **Read first:** `TALK.md`, `DISCOVERY-PLAN.md`, `STUDY-GUIDE.md`, `READING-LIST.md`
2. **Check branch:** Work on `feat/wolfram-talk-prep` (or create a sub-branch)
3. **Never push to master** — Ernest merges manually
4. **LaTeX builds:** `cd amplitudes && pdflatex 99-master.tex` (run twice for ToC)
5. **Scripts:** All in `amplitudes/scripts/`, standalone Python, require numpy + sympy
6. **Cadabra2:** Scripts in `amplitudes/Srednicki/` need Docker container `cadabra2-ubuntu:24.04`
7. **Key source files:**
   - SMGA paper: `workspace2/Data/Public/papers/physics/arXiv-2602.12176v2/SMGA.tex`
   - Dixon review: `workspace2/Data/Public/papers/physics/arXiv-1310.5353v1/ModAmpIntro.tex`
   - Srednicki LaTeX: `workspace2/Data/Public/books/Physics/Srednicki-QuantumFieldTheory/Srednicki_QFT.tex`
8. **Don't commit:** build artifacts, `.aux`, `.log`, `.toc`, `.out`, `.fls`, `.fdb_latexmk`, compiled binaries
9. **Update this file** when you complete a task (change Status to ✅ Done with date)

---

Last updated: 2026-04-06
