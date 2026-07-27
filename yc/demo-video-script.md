# Product Demo — 90-Second Public Script

**Company name:** Explore The Universe (finalized 2026-07-25).

**Prerequisite met 2026-07-25:** TASK-009 shipped `amplitudes/scripts/09_verification_demo.py`
(see [`DEMO.md`](../amplitudes/DEMO.md)). Command:

```bash
uv run --project amplitudes/MHVamplitudes python amplitudes/scripts/09_verification_demo.py
```

Runtime is well under a second; narration below assumes screen-recording the real output, not
a mockup. Record at a normal, unhurried speaking pace — the terminal output alone takes under
a second, so there's no need to rush the narration to fit; let the numbers sit on screen.

## 0–3 seconds — Title card

On-screen text only, no narration yet: **"Explore The Universe"** with the one-line tagline
"AI agents that verify theoretical physics." Cut to terminal.

## 3–13 seconds — Problem

Show a correct amplitude expression next to a subtly wrong one.

“AI can generate physics that looks right while being wrong — a sign, a convention, a
normalization factor. Explore The Universe separates proposing an expression from verifying
it.”

## 13–28 seconds — Scope

Show the closed-form implementation and specialized recurrence.

“This example uses a known all-`n` formula from arXiv:2602.12176 and a separately implemented
recurrence described in the same paper. This is a reproducibility test, not a new physics
claim.”

## 28–58 seconds — Run

Run the one-command demo and show `n = 3..10`, seeds, trial counts, and residuals.

“Across seeded kinematics through ten particles, both code paths agree exactly — residual zero
every time. The focused unit and integration suite has 43 passing tests, and the full suite is
64 passing with no known failures as of today.”

## 58–75 seconds — Deliberate failure

Show the broken-candidate section of the same run.

“This second candidate has one deliberate bug — a wrong normalization exponent, built from the
same verified building blocks. The verifier catches it immediately and shows red, not green. A
verifier that can't say no isn't a verifier.”

## 75–89 seconds — Evidence

Open the generated report (`amplitudes/results/demo-latest.md`).

“Every run records its git commit, branch, backend versions, seeds, and result, so a
researcher can audit exactly what was checked and reproduce it.”

## 89–93 seconds — Close

“Today this reproduces known work. Next, the same loop tests bounded new conjectures and keeps
failures as data.”

Cut back to the title card: **"Explore The Universe"**.

## Rules

- Actual terminal output only — this script was written against a real run on 2026-07-25;
  re-verify output matches before recording.
- Do not say “independent proof,” “from first principles,” “beyond the paper,” or “discovered.”
- Do not show private application links, notifications, credentials, or employer material.
- If the demo behavior differs, update the narration before recording.
