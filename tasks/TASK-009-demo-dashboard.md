# TASK-009 — Live demo artifact: verification dashboard

Status: TODO
Priority: P0 for the YC application (deadline 2026-07-27) — see `../yc/README.md`
Depends on: TASK-003 (its JSON output is the data source); degraded demo possible from TASK-001 alone
Estimated size: one session (3–5 h)
Branch suggestion: `feat/task-009-demo`

## Goal

A single self-contained demo Ernest can screen-record for the YC application (and reuse for
SpaceX-wedge portfolio + shorts): show the loop *LLM-conjectured formula → symbolic/numeric
engine verifies it live → verdict with error bars*. Concretely: SMGA Eq.(16) being verified
across n and random kinematic points in real time.

## Requirements

1. **One command to run.** `uv run python amplitudes/dashboard/serve.py` (there is an
   existing `amplitudes/dashboard/` — inspect and extend it rather than starting over).
2. **Left panel:** the conjectured formula, rendered (MathJax/KaTeX inline is fine).
3. **Right panel:** streaming verification runs — n, sampled point, |Eq.(16) − BG| relative
   error, running max — with a green/red pass/fail verdict per n.
4. **A "wrong conjecture" toggle:** flip one sign in the formula and show the verifier
   catching it immediately. This is the money shot — it proves the system can't be
   bullshitted, which is the entire thesis.
5. Terminal fallback: `--headless` mode printing the same run as a clean table, in case
   the recording is done in a terminal instead of a browser.

## Non-goals

No auth, no deployment, no styling polish beyond legible. It must be honest: real
computations running live, no pre-baked numbers presented as live.

## Definition of done

- Fresh checkout of the branch → documented one command → dashboard runs and completes a
  full n = 3…10 verification pass in under ~2 minutes on the M5 MacBook.
- `amplitudes/dashboard/DEMO-SCRIPT.md` written: 90-second walkthrough of what to click
  and say while recording (coordinate with `../yc/demo-video-script.md`).
