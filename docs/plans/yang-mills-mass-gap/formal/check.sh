#!/usr/bin/env bash
# Compile the YMFormal modules against Mathlib, in Docker, with no network.
# Usage: formal/check.sh [MODULE ...]     (default: all modules, in dependency order)
# Exit status is nonzero if any module fails. A `sorry` is reported as a warning by Lean;
# `--strict` makes any `sorry` an error and prints the axioms of the main theorems (Axioms.lean).
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
image=${YM_FORMAL_IMAGE:-claw-work/lean-agent:lean4.35.0-mathlib}
strict=0
if [ "${1:-}" = "--strict" ]; then strict=1; shift; fi
modules=("$@")
if [ ${#modules[@]} -eq 0 ]; then
  modules=(YMFormal.Purity YMFormal.Defect YMFormal.ConditionalGap YMFormal.Casimir YMFormal.RankOne YMFormal.Rectangle YMFormal.Seed)
fi
docker run --rm --network none -u "$(id -u):$(id -g)" -v "$here":/work -w /work \
  --entrypoint bash "$image" -c '
set -uo pipefail
strict=$1; shift
mkdir -p build
export LEAN_PATH=/work/build:$LEAN_PATH
status=0
for module in "$@"; do
  source=${module//./\/}.lean
  mkdir -p "build/$(dirname "$source")"
  output=$(lean -R /work -o "build/${module//./\/}.olean" "$source" 2>&1); code=$?
  printf "%s\n" "$output"
  if [ $code -ne 0 ]; then echo "FAIL $module"; status=1; continue; fi
  if [ "$strict" = 1 ] && grep -q "declaration uses .sorry." <<<"$output"; then
    echo "FAIL $module (sorry)"; status=1; continue
  fi
  echo "OK   $module"
done
if [ "$strict" = 1 ] && [ $status -eq 0 ] && [ -f Axioms.lean ]; then
  axioms=$(lean Axioms.lean 2>&1); printf "%s\n" "$axioms"
  if grep -q "sorryAx" <<<"$axioms" || ! grep -q "depends on axioms" <<<"$axioms"; then
    echo "FAIL axioms"; status=1
  else
    echo "OK   axioms (propext, Classical.choice, Quot.sound only)"
  fi
fi
exit $status' check "$strict" "${modules[@]}"
