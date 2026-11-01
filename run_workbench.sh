#!/usr/bin/env bash
# Build outside the repository, then serve approved roots read-only on loopback.
set -euo pipefail
repo_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd -P)
frontend_dir=${WORKBENCH_FRONTEND_DIR:-"$repo_dir/TypeScript/physics-workbench"}
cache_dir=${WORKBENCH_CACHE_DIR:-"${XDG_CACHE_HOME:-$HOME/.cache}/mathphysics-workbench"}
books_dir=${WORKBENCH_BOOK_ROOT:-"$(dirname -- "$(dirname -- "$repo_dir")")/Data/Public/books/Physics/Srednicki-QuantumFieldTheory"}
port=${WORKBENCH_PORT:-8876}
build_only=false
while (($#)); do
  case "$1" in
    --port) [[ $# -ge 2 ]] || { echo 'Missing --port value' >&2; exit 2; }; port=$2; shift 2 ;;
    --build-only) build_only=true; shift ;;
    *) echo "Usage: $0 [--port N] [--build-only]" >&2; exit 2 ;;
  esac
done
[[ $port =~ ^[0-9]+$ ]] && ((${#port} <= 5 && 10#$port <= 65535)) || { echo 'Invalid port' >&2; exit 2; }
[[ -f "$frontend_dir/package-lock.json" ]] || { echo "Frontend missing: $frontend_dir (integrate TypeScript/physics-workbench first)" >&2; exit 1; }
npm_offline=${WORKBENCH_NPM_OFFLINE:-1}
[[ $npm_offline == 0 || $npm_offline == 1 ]] || { echo 'WORKBENCH_NPM_OFFLINE must be 0 or 1' >&2; exit 2; }
root_args=(--root "mathphysics=$repo_dir")
if [[ -d "$books_dir" ]]; then
  books_dir=$(cd -- "$books_dir" && pwd -P)
  root_args+=(--root "books=$books_dir")
else
  echo "Optional Srednicki corpus unavailable: $books_dir; serving mathphysics root alone (set WORKBENCH_BOOK_ROOT to override)." >&2
fi
mkdir -p -- "$cache_dir"
cache_dir=$(cd -- "$cache_dir" && pwd -P)
case "$cache_dir/" in "$repo_dir/"*) echo 'Build cache must be outside repository' >&2; exit 1 ;; esac
build_dir=$(mktemp -d "$cache_dir/build.XXXXXXXX")
mkdir -- "$build_dir/frontend"
# Unique snapshots keep builds and Vite configuration caches out of source.
(tar -C "$frontend_dir" --exclude='./node_modules' --exclude='./dist' --exclude='./.git' --exclude='./test-results' --exclude='./playwright-report' -cf - .) | tar -C "$build_dir/frontend" -xf -
if [[ -n ${WORKBENCH_NODE_MODULES:-} ]]; then
  [[ -d "$WORKBENCH_NODE_MODULES" ]] || { echo 'WORKBENCH_NODE_MODULES must be a directory' >&2; exit 1; }
  deps_dir=$(cd -- "$WORKBENCH_NODE_MODULES" && pwd -P)
  ln -s -- "$deps_dir" "$build_dir/frontend/node_modules"
else
  # Network is permitted only by explicit user opt-in; installs stay in cache.
  npm_args=(ci --ignore-scripts --no-audit --no-fund)
  if [[ $npm_offline == 1 ]]; then
    npm_args+=(--offline)
  else
    echo 'WORKBENCH_NPM_OFFLINE=0: online npm dependency installation explicitly enabled.' >&2
  fi
  if ! (cd -- "$build_dir/frontend" && npm "${npm_args[@]}"); then
    echo 'Dependency install failed. For an explicit online bootstrap, rerun with WORKBENCH_NPM_OFFLINE=0, or set WORKBENCH_NODE_MODULES to matching installed dependencies.' >&2
    exit 1
  fi
fi
(cd -- "$build_dir/frontend" && WORKBENCH_DIST="$build_dir/dist" npm run build)
export CARGO_TARGET_DIR="$cache_dir/rust-target"
cargo build --offline --locked --manifest-path "$repo_dir/Rust/physics-gateway/Cargo.toml"
if [[ ${WORKBENCH_APP:-notebook} == explore ]]; then
  cp -- "$build_dir/dist/explore.html" "$build_dir/dist/index.html"
fi
printf 'Frontend: %s\nGateway: %s\n' "$build_dir/dist" "$CARGO_TARGET_DIR/debug/physics-gateway"
if $build_only; then exit 0; fi
printf 'Application: http://127.0.0.1:%s\nRead-only root: mathphysics=%s\n' "$port" "$repo_dir"
if ((${#root_args[@]} > 2)); then printf 'Read-only root: books=%s\n' "$books_dir"; fi
exec "$CARGO_TARGET_DIR/debug/physics-gateway" --static-dir "$build_dir/dist" "${root_args[@]}" --port "$port"
