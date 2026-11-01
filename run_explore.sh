#!/usr/bin/env bash
# Independent Explore Physics application; notebook is not loaded.
set -euo pipefail
repo_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
export WORKBENCH_APP=explore
export WORKBENCH_PORT=${EXPLORE_PORT:-8877}
export WORKBENCH_CACHE_DIR=${EXPLORE_CACHE_DIR:-"${XDG_CACHE_HOME:-$HOME/.cache}/mathphysics-explore"}
exec "$repo_dir/run_workbench.sh" "$@"
