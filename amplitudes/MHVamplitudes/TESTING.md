# TESTING.md - Running Tests for MHVamplitudes

## Prerequisites

- Cadabra2 Docker image built and available (see Monoclaw repo for build instructions: `repos/Monoclaw/Deployments/DockerBuilds/Physics/Cadabra2/`).

- mathphysics repo cloned locally.

## Running Tests

1. Start the Cadabra2 Docker container with the mathphysics repo mounted. This uses an environment variable to avoid hardcoding paths, ensuring portability across machines or OpenClaw instances.

   Assuming the Docker run script is at `~/.openclaw/workspace/repos/Monoclaw/Deployments/DockerBuilds/Physics/Cadabra2/run.sh`:

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

## Notes

- The `-s` flag on pytest captures full output, useful for debugging Cadabra2/SymPy interactions.

- If adding more tests, ensure they cover key mathematical objects (e.g., spinor brackets, Schouten identities) to align with the discovery plan.

- For CI/CD: Consider GitHub Actions workflow to automate this on push/PR.

This doc can be expanded as the project grows.