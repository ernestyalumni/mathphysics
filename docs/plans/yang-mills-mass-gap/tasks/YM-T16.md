# YM-T16 — A weight edit that makes the local 9B referee catch more planted errors (claw-work)

Tier T2 · Workstream B (tooling for the program itself) · Depends on: YM-T14 (`bin/ym-micro`, the proof units and mutants) · Review: T2, different family · Budget: two GPU nights

## Why
The local 9B models are useful only as flagging referees (`notes/local-models.md` §§4a, 6). They miss errors that need bookkeeping over several steps (probe row h). Ernest's request, 2026-11-17: manipulate the weights of a local model (Qwen3.8 9B Q5) directly to improve its accuracy on this program's work. The code lives in the claw-work repository, `Rust/weight-space/` (Rust host code, CUDA C++ for every GPU computation, llama.cpp for inference); this packet owns the Yang–Mills side: the cases, the reading of the results, and what the program does with them.

## Inputs
claw-work `Rust/weight-space/docs/PHASE2.md` (design and pre-registered protocol) and `scripts/referee-experiment.sh`; `bin/ym-micro export-cases`; the units in the YM-T14 state directory.

## Do
1. `bin/ym-micro export-cases --state <YM-T14 state> --out <dir outside git>`: discovery and held-out referee case files, split by proof (rules in the docstring). Record the two SHA-256 values.
2. Run claw-work's `scripts/referee-experiment.sh` on the 3060 when no other server holds it.
3. Report in `notes/YM-T16-weight-space-referee.md`: the gate statistics first (do missed planted errors go with shorter reasoning?), then held-out mutants located and false alarms for control, edited and prompt-only arms, per case as well as pooled.
4. If the edited model is better than both controls on held-out mutants without more false alarms, propose adding it as a fourth referee in the union-of-flags rule of `notes/local-models.md` §6; otherwise say plainly that the edit did not help.

## Acceptance
Discovery and held-out proofs are disjoint (the split file shows it). Every number in the note comes from files written by the runs, with their paths and hashes. Negative results are a valid outcome.

## Do not
Use held-out outputs to choose anything. Treat the edited model's verdicts as evidence: like every local model, its INCORRECT flags are cheap pointers for a stronger reviewer, and its CORRECT verdicts carry no weight.
