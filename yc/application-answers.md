# YC Application — Draft Answers

Paste-ready drafts for the standard YC form questions. `[FILL]` marks facts only Ernest
can supply — do not submit with any `[FILL]` remaining. Keep answers in your own voice;
these are 90% drafts, not scripts. YC reads for clarity and honesty, not polish.

---

## Company name

TMP Labs *(working name — alternatives: Amplitude Labs, Verifier Physics. Pick one and
stop thinking about it; YC does not care about names.)*

## Describe what your company does (50 characters)

AI agents that discover and verify new physics.

## What is your company going to make? Please describe your product and what it does or will do.

We build a discovery engine for theoretical and mathematical physics: fleets of AI agents
that read the literature, conjecture new results, and — the part that matters — **verify
every claim with symbolic and numerical computation** (Cadabra2, SymPy, FORM) before it's
allowed to be called a result. LLMs are already good enough to conjecture real physics: in
Feb 2026, a GPT-5.2-conjectured closed form for single-minus gluon amplitudes was published
by Strominger's group (arXiv:2602.12176). The bottleneck isn't ideas — it's trustworthy
verification at scale. That's what we automate.

Two products come out of the same engine: (1) the research engine itself, which does the
work of a team of graduate students — our milestones are arXiv papers it materially
produced; (2) a physics-hallucination eval suite and verifier API, where every test item is
machine-gradable by our symbolic backends — sold to frontier AI labs, who currently pay
heavily for exactly this kind of hard, un-gameable reasoning eval.

## Why did you pick this idea to work on? Do you have domain expertise in this area? How do you know people need what you're making?

I hold a graduate degree in Theoretical and Mathematical Physics from LMU Munich [FILL:
exact degree/program wording], and I've spent [FILL: N] years building GPU/CUDA numerical
and simulation software in industry [FILL: 1-line strongest industry credential]. This is
the field I was trained for, now attackable with tools that didn't exist two years ago.

Demand evidence: frontier labs (OpenAI, Google DeepMind, Anthropic) are publicly racing to
demonstrate AI-driven science — the 2602.12176 paper *is* OpenAI demonstrating demand — and
labs pay for expert-constructed hard evals (FrontierMath-style) precisely because they can't
generate trustworthy ones internally. Meanwhile theory groups run on grad-student labor
that costs ~$80k/yr/person and is chronically scarce.

## What's new about what you're making? What substitutes do people resort to because it doesn't exist (or they don't know about it)?

New: verifier-in-the-loop discovery. Everyone else demos an LLM writing plausible physics;
we wire the LLM to symbolic engines that mechanically check every sign, index, and limit,
so wrong conjectures die in seconds instead of surviving into papers. Substitutes today:
grad students (scarce, slow), Mathematica in expert hands (not autonomous), and raw LLMs
(hallucinate silently — wrong signs and index placements that read as correct).

## Who are your competitors? What do you understand about your business that they don't?

Frontier labs' internal science teams (OpenAI's work with Strominger et al.), Google
DeepMind (AlphaProof/AlphaEvolve lineage), academic groups using LLMs ad hoc, and
math-adjacent startups (Harmonic, Axiom-style formal-math efforts). What we understand:
formal proof (Lean) is the wrong substrate for theoretical physics — physicists work in
computer algebra, and the verification layer that matters is symbolic computation with
correct conventions (metric signatures, spinor conventions), which is exactly the
unglamorous expertise we have and labs keep getting silently wrong. Also: labs will pay
for evals *from outside* precisely because internal evals leak into training.

## How do or will you make money? How much could you make?

Near-term: (1) eval-suite licensing + private benchmark runs for frontier labs — deals in
this market run $[FILL: research current FrontierMath/Epoch-style deal sizes]00k+/yr per
lab; (2) verifier API + research-agent platform for physics and physics-adjacent R&D
groups (national labs, quant firms, aerospace) on seats + compute. Long-term: the lab that
owns verified AI-driven discovery in physics owns a piece of every downstream result. The
honest version: this starts as a high-margin evals/tooling business with a research lab
attached, and the research lab is the moat and the marketing.

## How far along are you?

- Working spinor-helicity + amplitude computation stack (Python package with tests, 8
  verification scripts: color/Fierz identities, Parke-Taylor, soft limits, SMGA stripped
  amplitudes) in a public repo [FILL: link].
- Reproduction of arXiv:2602.12176's results with independent verification of their
  Eq.(16) to n = 10, beyond the paper's published checks [status: in progress this week —
  update before submitting].
- Symbolic verification layer over Cadabra2/SymPy with known-answer tests.
- Agent-orchestration layer (multi-agent task board, autonomous sessions) already running
  this repo day-to-day.
- Related shipping proof: built claw-dj, an autonomous DJ agent, at the H Company
  hackathon [FILL: result/placement], solo, in [FILL: timeframe].

## How long have each of you been working on this? Full-time?

Working on the physics-agent stack since [FILL: month] 2026, alongside contract/interview
pipeline; going full-time on it through the batch. [Adjust to reality — do not overstate.]

## Who writes code, or does other technical work on your product?

Ernest Yeung — 100%. All research, code, and agent infrastructure. [If applying solo,
answer the solo-founder question head-on: you've shipped an entire working stack alone in
months; you're open to a cofounder with [FILL: complementary skill] and have [FILL: any
actual candidate/network, else delete].]

## What is your company going to do next? (or: What's your plan for the batch?)

Batch goal: ship the first arXiv paper with a genuinely new, machine-verified result in
gluon amplitudes (regions beyond R₁ / NMHV extensions — concrete candidate results already
mapped), and close 2 paid eval-suite pilots with AI labs. Demo-day story: "our agents found
and verified something humans hadn't, and labs pay us to test their models against it."

## Equity / legal

[FILL: not yet incorporated / incorporated as X. If not incorporated: say so plainly;
YC handles incorporation at funding. 100% Ernest.]

---

## Founder video

1 minute, unlisted YouTube. Script: `founder-video-script.md`. YC's guidance: no slides,
no production value, just the founder talking. Face the light, phone at eye level.

## Demo video (optional but we include it)

Script: `demo-video-script.md`, recorded off the TASK-009 dashboard.
