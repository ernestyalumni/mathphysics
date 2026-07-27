# TESTING.md - Running Tests for MHVamplitudes

**Last verified:** 2026-07-25, Python 3.12.12, `uv` 0.9.8, fresh checkout (`.venv` and
`uv.lock` deleted before each run below).

## Focused MHV command (no Docker required)

From this directory:

```bash
uv run pytest -q tests/unit_tests/amplitudes tests/integration_tests/test_smga_reproduction.py
```

Result as of 2026-07-25: `43 passed`, no warnings. `numpy` and `pytest` are declared in
`pyproject.toml` (`dependencies` / `[dependency-groups] dev`), so `uv` resolves them
automatically — no `--with` flags, no `--no-project`, no manual environment needed.

Add `-m "not slow"` to skip the `n = 8, 9, 10` cases (exponential-time recursion), or
`-m slow` to run only those.

## Full-suite command (also no Docker required, currently)

```bash
uv run pytest -q tests/
```

Result as of 2026-07-25: `64 passed`. The spinor/gamma-matrix tests
(`tests/unit_tests/spinors/`) do not currently import `cadabra2`, so the full suite does not
require the Cadabra2 Docker image today.

## Cadabra2 (only needed for future Cadabra2-backed verification work)

`src/mhvamplitudes/spinors/cadabra2/declarations.py` imports `cadabra2`, but no test currently
exercises that module. If/when Cadabra2-backed verifier tests are added, use the Docker
workflow below.

Prerequisite: Cadabra2 Docker image built and available (see Monoclaw repo for build
instructions: `repos/Monoclaw/Deployments/DockerBuilds/Physics/Cadabra2/`).

1. Start the Cadabra2 Docker container with the mathphysics repo mounted. This uses an
   environment variable to avoid hardcoding paths, ensuring portability across machines or
   OpenClaw instances.

   Assuming the Docker run script is at
   `~/.openclaw/workspace/repos/Monoclaw/Deployments/DockerBuilds/Physics/Cadabra2/run.sh`:

   ```
   MATHPHYSICS_DIR=/path/to/your/mathphysics/repo ./run.sh cli
   ```

   Example (on MS-7885 desktop):
   ```
   MATHPHYSICS_DIR=/home/propdev/.openclaw/workspace/workspace2/repos/mathphysics ./run.sh cli
   ```

2. Inside the Docker container:

   ```
   cd /mathphysics/amplitudes/MHVamplitudes/
   pytest -s ./tests/
   ```

   The `-s` flag captures full output, useful for debugging Cadabra2/SymPy interactions.

## Notes

- If adding more tests, ensure they cover key mathematical objects (e.g., spinor brackets,
  Schouten identities) to align with the discovery plan.
- For CI/CD: Consider GitHub Actions workflow to automate this on push/PR.

This doc can be expanded as the project grows.
