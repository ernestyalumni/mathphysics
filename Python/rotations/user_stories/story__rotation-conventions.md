<!-- pdd-story-prompts: prompts/rotations_python.prompt -->
# Rotation conventions

## Story
As a researcher comparing Goldstein and engineering attitude conventions, I want regeneration to preserve the approved physical rotation actions, angle ordering, canonical representatives and rejection boundaries, so that a convention error cannot silently change my calculations.

## Source
The user approved continuing Claude Code's proposed adoption on 2026-10-31.
The full independent request is preserved in
[the intent event](../docs/intents/intent__rotation-library-conventions-as-pdd-source-22448aa3.md).

## Acceptance evidence
`tests/test_rotations.py` and `tests/test_rotation_conventions.py` check physical
quarter turns, explicit ZYX matrix entries, singular choices, canonical ranges,
Hamilton/CK signs, invalid inputs, no projection and standard-library execution.
The containing Python suite checks existing dynamics and HTTP integration.

The generated story contract and semantic model validation remain pending if
the provider-backed adoption workflow fails. This story is agent-authored from
the approved request; it is not a claim that a model generated a contract.
