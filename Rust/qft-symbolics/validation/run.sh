#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
repo="$(cd ../.. && pwd)"
cache="${XDG_CACHE_HOME:-$HOME/.cache}/mathphysics-rust"
export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$cache/target}"
export PYTHONDONTWRITEBYTECODE=1
run="${MATHPHYSICS_RUN_DIR:-$cache/runs/$(date -u +%Y%m%dT%H%M%SZ)}"
mkdir -p "$run"
python="${PYTHON:-$repo/amplitudes/MHVamplitudes/.venv/bin/python}"
image="${CADABRA_IMAGE:-cadabra2-ubuntu:24.04}"
{
  date -u
  git rev-parse HEAD
  git branch --show-current
  git status --short
  rustc --version
  cargo --version
  "$python" --version
  docker image inspect "$image" --format '{{.Id}}'
  sha256sum Cargo.toml Cargo.lock corpus.json
  find src tests validation docs -type f \( -name '*.rs' -o -name '*.py' -o -name '*.sh' -o -name '*.tex' \) -print0 |
    sort -z | xargs -0 sha256sum
} > "$run/provenance.txt"
cargo fmt --check 2>&1 | tee "$run/format.log"
cargo clippy --offline --all-targets -- -D warnings 2>&1 | tee "$run/clippy.log"
cargo test --offline -- --nocapture 2>&1 | tee "$run/tests-debug.log"
cargo test --offline --release 2>&1 | tee "$run/tests-release.log"
cargo build --offline 2>&1 | tee "$run/build.log"
"$python" validation/corpus_check.py | tee "$run/corpus.json"
binary="$CARGO_TARGET_DIR/debug/qft-symbolics"
"$binary" demo | tee "$run/demo.json"
"$binary" tensor-cases | tee "$run/tensor-cases.json"
"$binary" schouten | tee "$run/schouten.json"
if "$binary" schouten --wrong-sign > "$run/negative-control.json"; then
  echo "ERROR: deliberately wrong candidate passed" >&2
  exit 1
else
  code=$?
  test "$code" = 1
fi
"$python" validation/python_parity.py "$binary" | tee "$run/python-parity.json"
docker run --rm --network none \
  -v "$PWD/validation:/validation:ro" "$image" \
  python3 /validation/cadabra_oracle.py | tee "$run/cadabra.json"
docker run --rm --network none \
  -v "$PWD/validation:/validation:ro" -v "$run:/evidence:ro" "$image" \
  python3 /validation/tensor_oracle.py /evidence/tensor-cases.json | tee "$run/cadabra-tensors.json"
docker run --rm --network none \
  -v "$PWD/validation:/validation:ro" "$image" \
  python3 /validation/dimension_oracle.py | tee "$run/cadabra-dimension.json"
printf '\nValidation complete. Evidence outside repository: %s\n' "$run"
