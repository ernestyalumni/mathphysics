# Prompt-driven development in mathphysics

Registered PDD scopes use nested configuration and architecture files:

- `Python/rotations`: rotation conversion module.
- `Rust/qft-symbolics`: fourteen bounded arithmetic/operator/Lorentz/CLI units; manual adoption, full regeneration unverified. See its PDD.md.
- `Rust/physics-contracts`: typed metadata contracts and JSON validation CLI.
- `Rust/physics-tools`: bounded JSON bridge to existing exact QFT operations.
- `TypeScript/physics-workbench`: six bounded browser TS/CSS units, including PDF.js rendering; manually implemented, automated regeneration unverified. See its PDD.md.

Run commands from the applicable unit directory. Initial agent implementation and path/preprocessing checks do not imply a successful full PDD regeneration; inspect each unit's evidence.
LaTeX notes, research notebooks, validation oracles, unmapped native library units,
`rigid_body.py` and `simulator_server.py` remain outside the prompt graph.
Do not infer ownership from the presence of a repository-level PDD guide.

Ernest describes desired behavior in ordinary language. Agents preserve that
meaning, run `pdd intent plan` before mutation, and apply the exact approved
intent after meaning approval. Approval from the current conversation persists;
do not ask him to choose prompt filenames, dev units, story slugs or commands.
For future behavior changes, edit the prompt before regenerating managed code.
Keep independent acceptance tests and numerical oracles outside generated output.

Before PDD work, read these canonical Monoclaw playbooks (sibling repository):

- `../Monoclaw/docs/pdd/PDD_INTENT_FIRST_WORKFLOW.md`
- `../Monoclaw/docs/pdd/PDD_NATURAL_LANGUAGE_AGENT_PLAYBOOK.md`
- `../Monoclaw/docs/pdd/PDD_WITH_ANY_AGENT_HARNESS.md`
- `../Monoclaw/docs/pdd/PDD_PROMPT_GRAPH_AND_ARTIFACT_BUNDLES.md`

Use Ernest's `../PromptDrivenDevelopment/pdd/.venv/bin/pdd`; installed help
is authoritative. Run commands from the applicable nested scope, with
`PDD_LOCAL_ONLY=1 PDD_SKIP_UPDATE_CHECK=1`. Local-only disables GitHub access;
model calls can still use configured external providers. Do not change provider
or authentication configuration to evade a failed run.

From the repository root, run the independent and integration tests with:

```sh
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=Python/rotations .venv/bin/python -m unittest discover -s Python/rotations/tests -v
```

The existing uv environment supplies SciPy for dynamics verification; the
rotation module itself uses only Python's standard library. HTTP tests need a
localhost test socket. Native/CAS suites are additional gates when changing
shared mathematical behavior, not evidence supplied by the Python suite.

Current adoption evidence and remaining verification: [PDD_ADOPTION.md](Python/rotations/docs/PDD_ADOPTION.md).

Adoption preserves the existing implementation. Configuration and passing
characterization tests alone do not prove that a fresh PDD regeneration passes.
Record provider failures and synchronization status honestly; never mark a
failed intent run as applied by editing its state. Never hand-edit a generated
story contract.
