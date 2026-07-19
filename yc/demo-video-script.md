# Demo Video — ~2-Minute Script

Prereq: TASK-009 dashboard running (`tasks/TASK-009-demo-dashboard.md`). Screen recording
with voiceover (QuickTime/OBS — see the Monoclaw short-form pipeline doc for OBS setup).
Everything shown must be computed live on screen — no pre-baked numbers narrated as live.

---

**[0:00 — repo/terminal visible]**

"This is a discovery engine for theoretical physics. I'll show you the loop that matters:
an AI conjectures a formula, and a symbolic engine decides — mechanically — whether it's
true."

**[0:15 — dashboard, formula panel]**

"In February, a formula for single-minus gluon amplitudes was conjectured by an LLM and
published by Strominger's group at Harvard. Here it is. The paper verified it in a
half-collinear regime. My engine is going to check it independently, right now, further
than the paper did."

**[0:35 — start the run; verification rows streaming]**

"Each row is a random point in the kinematic regime. On every point, the engine computes
the amplitude two completely different ways — the conjectured closed form, and direct
Berends-Giele recursion from first principles — and compares. Three particles… up through
ten. Relative errors at ten to the minus ten. The formula holds."

**[1:05 — the money shot: wrong-conjecture toggle]**

"Now the important part. Here's the same formula with one sign flipped — the kind of error
language models make constantly, and the kind that reads as correct to a human skimming a
paper. Watch."

*[toggle → red FAIL within seconds]*

"Caught instantly. This engine cannot be bullshitted. That's the point: LLM ideas are
cheap now — machine-checked truth is the product."

**[1:35 — close, repo/tests visible]**

"Everything here is reproducible — one command, open code, every claim backed by a script.
Next milestones: new results in regimes nobody has mapped, published on arXiv, and this
same verifier sold as an eval suite to AI labs. I'm Ernest — this is TMP Labs."

**[~1:55 end]**

---

## Recording notes

- Rehearse the toggle beat once; the FAIL must land within ~5s on screen or trim the gap.
- Keep the mouse still while talking; move only when directing attention.
- If n=10 is slow live, it's fine to show 3…8 live and say "it runs to ten" only if it
  actually does — never narrate beyond what the recording shows.
