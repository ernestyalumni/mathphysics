# Pitch Deck — 10 Slides

YC's application doesn't require a deck, but you want one ready for the interview and for
any investor conversation this triggers. Build as HTML deck per the Monoclaw short-form
pipeline (or Keynote — whatever's fastest). One idea per slide, big type, no paragraphs.

---

## 1. Title
**TMP Labs** — AI agents that discover and verify new physics.
Ernest Yeung · ernestyalumni@gmail.com · [repo link]

## 2. The moment (why now)
Feb 2026: an LLM conjectured a new closed-form gluon amplitude formula → published with
Andrew Strominger (arXiv:2602.12176). *Speaker note: this is the single most important
slide — an outside, maximally credible proof that AI can source real theoretical physics.*

## 3. The problem
Ideas are no longer the bottleneck. **Trust is.** LLMs produce physics that looks right
and is silently wrong: signs, indices, tensor structure, invalid limits. Today's
verification layer is grad students — scarce, slow, ~$80k/yr each — or nothing.

## 4. The product
A discovery engine with a verifier in the loop: agents read → conjecture → **symbolic
engines (Cadabra2/SymPy/FORM) mechanically check** → only verified claims survive.
Diagram: conjecture → verify → publish loop, with the FAIL branch prominent.

## 5. Demo
Screenshot of the TASK-009 dashboard: the Strominger-paper formula verified to n=10 live;
the flipped-sign version caught in seconds. *Speaker note: "this system cannot be
bullshitted" — say it out loud.*

## 6. Business
Two revenue lines off one engine:
- **Evals & private benchmarks for frontier labs** — physics-hallucination suite where
  every item is machine-gradable; labs already pay six figures for hard reasoning evals
  that can't leak into training.
- **Verifier API / research agents** for physics-adjacent R&D (national labs, aerospace,
  quant). *The research lab is the moat and the marketing; evals are the margin.*

## 7. Milestones (traction slide, kept honest)
- Working amplitude + verification stack, public, tested. ✅
- SMGA reproduction + Eq.(16) verified beyond the paper's published range. [status]
- Next: first arXiv paper with a new machine-verified result (regions beyond R₁, NMHV).
- Long arc: agent fleets running Nobel/Fields-scale problem programs.

## 8. Why me
LMU Munich grad work in Theoretical & Mathematical Physics + years of industry GPU/CUDA
numerical software + demonstrated agent-infrastructure builder (this stack solo in months;
claw-dj autonomous DJ built at H Company hackathon). Rare combo: can do the physics AND
build the fleet.

## 9. Competition
Frontier labs' internal science pushes (they validate the market, and they'll buy evals
externally *because* internal ones leak), formal-math startups (Lean is the wrong
substrate for physics — physicists live in computer algebra), ad-hoc academic LLM use (no
verification discipline).

## 10. The ask / vision
YC batch goals: 1 arXiv paper the engine materially produced + 2 paid lab pilots.
Vision: the institution where AI does physics that humans then get to stand on —
discovery at industrial scale.

---

*Design: follow the dataviz/artifact-design skills if built as HTML; dark, sparse, one
hero visual per slide (helicity diagrams on 2–4, dashboard on 5).*
