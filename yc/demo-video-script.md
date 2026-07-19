# Product Demo — 90-Second Public Script

**Prerequisite:** TASK-009 produces one command, a deliberate failure case, and an evidence report.

## 0–10 seconds — Problem

Show a correct amplitude expression next to a one-sign perturbation.

“AI can generate physics that looks right while being wrong by one sign. TMP separates proposing an expression from verifying it.”

## 10–25 seconds — Scope

Show the closed-form implementation and specialized recurrence.

“This example uses a known all-`n` formula from arXiv:2602.12176 and a separate code implementation of the recurrence described in the same paper. This is a reproducibility test, not a new physics claim.”

## 25–55 seconds — Run

Run the one-command demo and show `n = 3..10`, seeds or point counts, and residuals.

“Across seeded kinematics through ten particles, both code paths agree. The focused unit and integration suite currently has 43 passing tests.”

## 55–72 seconds — Deliberate failure

Run the perturbed expression.

“When I flip one sign, the verifier fails and records the residual. A useful verifier has to show red results, not just green checks.”

## 72–86 seconds — Evidence

Open the generated report.

“Every run records its code version, backend, seeds, scope, and result so a researcher can audit what was actually checked.”

## 86–90 seconds — Close

“Today this reproduces known work. Next, the same loop tests bounded new conjectures and keeps failures as data.”

## Rules

- Actual terminal output only.
- Do not say “independent proof,” “from first principles,” “beyond the paper,” or “discovered.”
- Do not show private application links, notifications, credentials, or employer material.
- If the demo behavior differs, update the narration before recording.
