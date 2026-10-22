#!/usr/bin/env bash
set -euo pipefail
repo_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
viewer_dir="$repo_dir/JavaScript/rigid-body-viewer"
if [[ "${1:-}" == "--no-build" ]]; then
  shift
else
  if [[ ! -d "$viewer_dir/node_modules" ]]; then
    npm --prefix "$viewer_dir" ci --no-audit --no-fund
  fi
  npm --prefix "$viewer_dir" run build
fi
if [[ -x "$repo_dir/.venv/bin/python" ]] && "$repo_dir/.venv/bin/python" -c 'import numpy' 2>/dev/null; then
  exec "$repo_dir/.venv/bin/python" "$repo_dir/Python/rotations/simulator_server.py" "$@"
else
  exec uv run --project "$repo_dir/Python/rotations" "$repo_dir/Python/rotations/simulator_server.py" "$@"
fi
