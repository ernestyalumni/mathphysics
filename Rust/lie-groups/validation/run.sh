#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
repo="$(cd ../.. && pwd)"
cpp="$repo/CppMath/LieGroups"
cache="${XDG_CACHE_HOME:-$HOME/.cache}/mathphysics-lie-groups"
export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-$cache/target}"
export PYTHONDONTWRITEBYTECODE=1
run="${MATHPHYSICS_LIE_RUN_DIR:-$cache/runs/$(date -u +%Y%m%dT%H%M%SZ)}"
sage_image="${SAGE_IMAGE:-sagemath/sagemath:latest}"
cadabra_image="${CADABRA_IMAGE:-cadabra2-ubuntu:24.04}"
sympy="${SYMPY_ROOT:-$repo/../sympy/sympy}"
mkdir -p "$run"
{
  date -u
  git rev-parse HEAD
  git branch --show-current
  git status --short
  rustc --version
  cargo --version
  c++ --version
  cmake --version
  docker image inspect "$sage_image" --format '{{.Id}}'
  docker image inspect "$cadabra_image" --format '{{.Id}}'
  git -C "$sympy" rev-parse HEAD
  sha256sum Cargo.toml Cargo.lock "$cpp/CMakeLists.txt"
  find src tests examples validation docs "$cpp/include" "$cpp/tests" ../qft-symbolics/src ../qft-symbolics/tests -type f \
    \( -name '*.rs' -o -name '*.cpp' -o -name '*.hpp' -o -name '*.py' -o -name '*.sh' -o -name '*.tex' \) \
    -print0 | sort -z | xargs -0 sha256sum
} > "$run/provenance.txt"
cargo fmt --check 2>&1 | tee "$run/format.log"
cargo clippy --offline --all-targets -- -D warnings 2>&1 | tee "$run/clippy.log"
cargo test --offline -- --nocapture 2>&1 | tee "$run/tests-debug.log"
cargo test --offline --release 2>&1 | tee "$run/tests-release.log"
cargo run --offline --example frankel 2>&1 | tee "$run/example.log"
cargo build --offline
binary="$CARGO_TARGET_DIR/debug/lie-groups"
"$binary" verify | tee "$run/verify.json"
if "$binary" verify --wrong-sign > "$run/negative-control.json"; then exit 1; else test "$?" = 1; fi
if "$binary" unknown-command > "$run/unsupported.json"; then exit 1; else test "$?" = 2; fi
for mode in debug release sanitize; do
  build="$cache/cpp-$mode"
  args=(-DCMAKE_BUILD_TYPE=Debug)
  if test "$mode" = release; then args=(-DCMAKE_BUILD_TYPE=Release); fi
  if test "$mode" = sanitize; then args+=(-DLIE_GROUPS_SANITIZE=ON); fi
  cmake -S "$cpp" -B "$build" "${args[@]}" 2>&1 | tee "$run/cpp-$mode-configure.log"
  cmake --build "$build" -j2 2>&1 | tee "$run/cpp-$mode-build.log"
  ctest --test-dir "$build" -V 2>&1 | tee "$run/cpp-$mode-tests.log"
done
"$binary" symbolic > "$run/rust-symbolic.json"
"$cache/cpp-release/lie_groups_interchange" symbolic > "$run/cpp-symbolic.json"
"$binary" fixtures > "$run/rust-fixtures.json"
"$cache/cpp-release/lie_groups_interchange" fixtures > "$run/cpp-fixtures.json"
python3 validation/numeric_parity.py "$run/rust-fixtures.json" "$run/cpp-fixtures.json" | tee "$run/numeric-parity.json"
docker run --rm --pull=never --network none \
  -v "$PWD/validation:/validation:ro" -v "$run:/evidence:ro" -v "$sympy:/sympy:ro" \
  "$sage_image" sage -python /validation/sympy_oracle.py /evidence/rust-symbolic.json /evidence/cpp-symbolic.json | tee "$run/sympy.json"
docker run --rm --pull=never --network none \
  -v "$PWD/validation:/validation:ro" -v "$run:/evidence:ro" \
  "$sage_image" sage -python /validation/sage_oracle.py /evidence/rust-symbolic.json /evidence/cpp-symbolic.json | tee "$run/sage.json"
docker run --rm --pull=never --network none \
  -v "$PWD/validation:/validation:ro" "$cadabra_image" \
  python3 /validation/cadabra_oracle.py | tee "$run/cadabra.json"
mkdir -p "$run/latex"
for pass in 1 2; do
  pdflatex -interaction=nonstopmode -halt-on-error -output-directory="$run/latex" \
    docs/frankel-lie-groups.tex > "$run/latex/compile-$pass.stdout" 2>&1
done
if grep -E "undefined references|Citation .* undefined|Overfull" "$run/latex/compile-2.stdout"; then
  echo "LaTeX reference/layout regression" >&2
  exit 1
fi
# Run the old dependent kernel's complete existing validation too.
MATHPHYSICS_RUN_DIR="$run/qft-regression" bash ../qft-symbolics/validation/run.sh > "$run/qft-regression-console.log" 2>&1
printf '\nPASS. Evidence: %s\n' "$run"
