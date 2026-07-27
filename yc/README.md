# Public YC Companion — Explore The Universe Fall 2026

**Drafted:** 2026-07-19
**Official on-time deadline:** 2026-07-27 at 8:00 PM PT
**Internal submit target:** 2026-07-25

This directory contains public-safe product language and rehearsal material for Explore The Universe's YC application. Authenticated application links, equity, employment, user/revenue facts, and other private founder information do not belong in this public research repo.

## Product thesis

Explore The Universe is an AI research system for theoretical physics that generates candidate calculations, checks them with independent symbolic or numerical verifiers, and records why each result passed, failed, or remained inconclusive.

The first wedge is scattering amplitudes. The current prototype implements two separate code paths from arXiv:2602.12176v2—a closed-form expression and the paper's specialized recurrence—and compares them through `n = 10`. This is a reproducibility milestone, not an extension beyond the paper's all-`n` result and not an original discovery.

## Files

| File | Purpose | Status |
|---|---|---|
| `application-answers.md` | Public-safe answer language; private facts omitted | DRAFT |
| `founder-video-script.md` | Bullet rehearsal guide; do not recite | DRAFT |
| `demo-video-script.md` | Product-demo shot list | Ready (TASK-009 done 2026-07-25) |
| `pitch-deck.md` | Ten-slide outline | DRAFT, lower priority |

## Technical critical path

1. ~~TASK-000 — reproducible developer/test entry point.~~ Done 2026-07-25.
2. ~~TASK-009 — passing and deliberately failing demo cases plus evidence ledger.~~ Done
   2026-07-25 — see `amplitudes/DEMO.md`.
3. TASK-002/TASK-003 — persistent reproduction reports with run metadata.
4. Record videos only after the behavior and numbers are frozen.

## Truth rules

- The 2026 paper reports that GPT-5.2 Pro conjectured its key formula and that another model proved it; the authors checked it using recurrence and consistency relations.
- Our current code reproduces known work using separately implemented code paths derived from that paper.
- Do not say our `n = 10` checks go beyond an all-`n` result.
- Do not claim users, revenue, pilots, market prices, or original physics without evidence.
- Keep the long-term scientific ambition distinct from the near-term product proof.
- Publishing, uploading, or submitting requires Ernest's approval.
