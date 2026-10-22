#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
module="$PWD"
repo="$(cd ../.. && pwd)"
python="${ROTATIONS_PYTHON:-$repo/.venv/bin/python}"
run="${ROTATIONS_RUN_DIR:-$(mktemp -d /tmp/mathphysics-goldstein-XXXXXXXX)}"
mkdir -p "$run"
export PYTHONDONTWRITEBYTECODE=1
export PYTHONPATH="$module${PYTHONPATH:+:$PYTHONPATH}"
"$python" -m unittest discover -s tests -v > "$run/python-tests.log" 2>&1
cargo build --offline --manifest-path "$repo/Rust/lie-groups/Cargo.toml" --example goldstein_fixtures --example rigid_body_fixtures --target-dir "$run/rust-target" > "$run/rust-build.log" 2>&1
rust_binary="$run/rust-target/debug/examples/goldstein_fixtures"
"$rust_binary" > "$run/rust-numerical.jsonl"
"$rust_binary" --symbolic > "$run/rust-symbolic.jsonl"
rust_dynamics="$run/rust-target/debug/examples/rigid_body_fixtures"
"$rust_dynamics" > "$run/rust-dynamics.jsonl"
cmake -S "$repo/CppMath/LieGroups" -B "$run/cpp" -DCMAKE_BUILD_TYPE=Release > "$run/cpp-configure.log" 2>&1
cmake --build "$run/cpp" --target goldstein_rotation_interchange rigid_body_interchange -j2 > "$run/cpp-build.log" 2>&1
cpp_binary="$run/cpp/goldstein_rotation_interchange"
"$cpp_binary" > "$run/cpp-numerical.jsonl"
cpp_dynamics="$run/cpp/rigid_body_interchange"
"$cpp_dynamics" > "$run/cpp-dynamics.jsonl"
"$python" validation/numeric_parity.py "$run/rust-numerical.jsonl" "$run/cpp-numerical.jsonl" | tee "$run/numeric-parity.json"
"$python" validation/dynamics_oracle.py "$run/rust-dynamics.jsonl" "$run/cpp-dynamics.jsonl" | tee "$run/dynamics.json"
"$python" validation/sympy_oracle.py "$run/rust-symbolic.jsonl" > "$run/sympy.jsonl" 2> "$run/sympy-progress.log"
if test "${WITH_CADABRA:-1}" = 1; then
  image="${CADABRA_IMAGE:-cadabra2-ubuntu:24.04}"
  docker image inspect "$image" --format '{{.Id}}' > "$run/cadabra-image.txt"
  docker run --rm --pull=never --network none -v "$module/validation:/validation:ro" "$image" python3 /validation/cadabra_oracle.py | tee "$run/cadabra.json"
else
  printf '%s\n' '{"status":"SKIPPED","reason":"WITH_CADABRA=0"}' > "$run/cadabra.json"
fi
"$python" validation/provenance.py "$repo" "$run" "$rust_binary" "$cpp_binary" "$rust_dynamics" "$cpp_dynamics" > "$run/provenance.json"
printf 'PASS. Evidence: %s\n' "$run"
