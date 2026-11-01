# PDD adoption evidence — 2026-10-31

## Scope and completed setup

Initial PDD ownership is one unit: `prompts/rotations_python.prompt` maps to
`rotations.py`. The nested `.pddrc` and `architecture.json` define that mapping.
All ten approved constraints and both examples appear verbatim in the prompt.
The prompt also preserves the public API, return types, numerical tolerances,
singular branch choices, validation boundaries and no-input-mutation behavior.
The repository `PDD.md`, agent entry point and local README explain future work.
No production source or package dependency metadata was changed.

The exact Claude Code request was recovered from its local PDD planning record.
Replanning returned the same ID:
`rotation-library-conventions-as-pdd-source-22448aa3`.
The user's takeover/go-ahead approves that proposed meaning and scope.
`pdd intent apply --characterized` preserved the complete request in
`docs/PRODUCT_INTENT.md` and the immutable `docs/intents/intent__*.md` event.

## Executed evidence

- Existing root uv environment, CPython 3.13:
  `PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=Python/rotations .venv/bin/python -m unittest discover -s Python/rotations/tests -v`
  passed **32 tests**: the prior 26 plus six independent convention checks.
  The complete run includes rigid-body and localhost HTTP integration tests.
- New tests cover explicit ZYX matrix coefficients and extraction order,
  quarter-turn quaternion components, canonical quaternion signs across branches,
  identity axis and canonical ranges, CK components/Goldstein transpose,
  nonfinite CK rejection, no projection or input mutation, and execution in an
  isolated Python process with site packages disabled.
- `pdd preprocess prompts/rotations_python.prompt --output /tmp/mathphysics-rotations-expanded.prompt --pdd-tags`
  passed without model inference.
- `pdd sync rotations --dry-run --json` resolved only `rotations.py`, the correct
  prompt and `tests/test_rotations.py`. It reported no conflicts or failures,
  but **UNBASELINED**, with `ok=false`; this is path-resolution evidence only.
- `git diff --check` passed. Production rotation, dynamics and server source
  remain byte-for-byte unchanged relative to HEAD.

## Automation and remaining limits

The automatic architecture workflow completed its first PRD-analysis step using
its configured Anthropic agent and entered step 2. Codex stopped only that
invocation and its descendants, then finished the bounded single-module setup
directly. The workflow was interrupted, not successfully applied; its runtime
intent record may still say `applying`. Do not interpret that stale status as a
running background job or mark it `applied` manually. The saved first-step
workflow state reports $0.9980925; incomplete later-step usage is not included
in that figure.

The supported `pdd baseline` command was attempted with Codex as reviewer. It
refused with `Input/Output Error: base and head must contain .pdd/repository-id`.
No baseline was fabricated. Repository identity initialization and the required
Git history must be handled through the fork's documented workflow before
claiming a current baseline.

The acceptance story is agent-authored from the approved independent request.
Its generated contract, semantic story validation, actual synchronization and
fresh regeneration have **not** completed. No generated contract was hand-edited.
Passing characterization tests establishes the adopted implementation's behavior;
it does not prove that the prompt regenerates an equivalent implementation.
Keep subsequent adoption work confined to this unit, preserve the full request,
and inspect current CLI help before resuming the interrupted intent workflow.

Changes are on `feat/pdd-rotations-adoption`, uncommitted and unpushed.
The pre-existing mono-notes build report remains unrelated and untouched.
