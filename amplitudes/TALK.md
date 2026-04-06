# Talk: OpenClaw — Agent-Driven Exploration of MHV Scattering Amplitudes

**Conference:** Bridging Wolfram Models and Mathematical Physics  
**Venue:** Austin, TX (Capital Factory)  
**Dates:** April 17–19, 2026  
**Speaker:** Ernest Yeung  
**Format:** 45-min talk + Q&A, or 15-min talk + extended Q&A  

## Abstract (draft)

We present OpenClaw, an open-source autonomous AI agent framework, and its application to
exploratory research in MHV (maximally helicity violating) scattering amplitudes. Using persistent
AI agents equipped with symbolic computation tools (Cadabra2, SymPy) and access to primary
literature (Srednicki QFT, Dixon 2013, Guevara et al. 2025), we demonstrate a workflow where
agents independently read papers, derive equations step-by-step, write verification scripts, and
build bridging documents between different conventions in the literature. We discuss what works,
what doesn't, and what this means for computational mathematical physics research.

## Talk Structure (suggested)

### Part 1: OpenClaw — What It Is (~20 min for 45-min format / ~8 min for 15-min)

1. **The problem**: Researchers drown in literature, notation translation, and routine computation
2. **What OpenClaw is**: Open-source autonomous AI agent platform
   - Persistent memory across sessions (MEMORY.md, daily logs)
   - Tool use: shell, web, file system, code execution
   - Sub-agent spawning for parallel work
   - Cron jobs for background/periodic tasks
   - Multi-channel: CLI, web, Discord, Telegram, WhatsApp
3. **The agent-driven research workflow**:
   - Agent reads LaTeX source of papers (not just PDFs)
   - Agent writes derivations in LaTeX
   - Agent writes Python/Cadabra2 verification scripts
   - Agent maintains reading lists, task boards, study guides
   - Human steers; agent executes
4. **Live demo or screenshots**: show the repo structure, a script running, the bridge document

### Part 2: MHV Scattering Amplitudes — The Physics (~20 min / ~5 min)

1. **Quick MHV primer** (for non-experts in the room):
   - Spinor-helicity formalism: massless momenta → 2-component spinors λ, λ̃
   - Parke-Taylor formula: A_n(1⁻,2⁻;3⁺...n⁺) = ⟨12⟩⁴ / (⟨12⟩⟨23⟩...⟨n1⟩)
   - Why it's beautiful: n-point tree amplitude in one line
2. **The SMGA discovery** (arXiv:2602.12176):
   - Single-minus amplitudes: long thought to vanish at tree level
   - GPT-5.2 Pro conjectured closed-form in half-collinear regime
   - What "AI discovery in physics" looks like in practice
3. **What we've built so far**:
   - Dixon–Srednicki bridge document (notation translation, full derivations)
   - 8 verification scripts (color algebra, spinor construction, MHV amplitudes, soft limits, SMGA stripped amplitudes)
   - Cadabra2 computations for Srednicki Ch.34–36, Ch.60
   - Modular LaTeX notes (~100 pages)
   - Study guide and reading list for future agents
4. **What's next**: reproduce SMGA numerically, explore extensions

### Part 3: Q&A / Discussion (~5 min / ~remaining time)

- How does this compare to Wolfram Language/Mathematica for symbolic physics?
- Connections to Wolfram models: discrete/computational approaches to physics
- Can agents handle genuinely creative mathematical reasoning?
- Open-source vs proprietary AI for research

## Key Points to Emphasize

- **OpenClaw is open-source** — anyone in the room can use it
- **Agents as research infrastructure**, not one-shot chatbots
- **Notation translation** is a huge unsolved problem in physics — agents can help
- **Reproducibility**: every derivation has a verification script
- **Honest about limitations**: agents make mistakes, need human oversight, can't (yet) do creative leaps
- **Connection to conference theme**: computational/discrete methods for mathematical physics

## Artifacts to Prepare

See `TASKS.md` for the full task list.

---

Created: 2026-04-06
