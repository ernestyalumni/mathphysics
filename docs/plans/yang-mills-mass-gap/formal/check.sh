#!/usr/bin/env bash
# Compile the YMFormal modules against Mathlib, in Docker, with no network.
# Usage: formal/check.sh [MODULE ...]     (default: all modules, in dependency order)
# Exit status is nonzero if any module fails. A `sorry` is reported as a warning by Lean;
# `--strict` makes any `sorry` an error and audits every theorem of the modules: each must be reported by
# `#print axioms`, and only propext, Classical.choice and Quot.sound are allowed.
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
image=${YM_FORMAL_IMAGE:-claw-work/lean-agent:lean4.35.0-mathlib}
strict=0
if [ "${1:-}" = "--strict" ]; then strict=1; shift; fi
modules=("$@")
if [ ${#modules[@]} -eq 0 ]; then
  modules=(YMFormal.Purity YMFormal.Defect YMFormal.ConditionalGap YMFormal.Casimir YMFormal.RankOne YMFormal.Rectangle YMFormal.Seed YMFormal.Doubling YMFormal.MatrixSpectrum YMFormal.Layer YMFormal.Z2)
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
if [ "$strict" = 1 ] && [ $status -eq 0 ]; then
  # Every theorem and lemma of the checked modules, attributes included; each must be reported, and
  # every axiom reported must be on the allowlist.
  audit=$(mktemp -d)
  for module in "$@"; do
    echo "import $module"
  done > "$audit/All.lean"
  names=$(for module in "$@"; do
    sed -nE "s/^(@\[[^]]*\] *)?(private +|protected +)?(theorem|lemma) +([^ :({]+).*/\4/p" "${module//./\/}.lean"
  done)
  for name in $names; do echo "#print axioms YMFormal.$name"; done >> "$audit/All.lean"
  axioms=$(lean "$audit/All.lean" 2>&1); code=$?
  wanted=$(wc -w <<<"$names")
  reported=$(grep -cE "depends on axioms|does not depend on any axioms" <<<"$axioms")
  other=$(tr "\n" " " <<<"$axioms" | grep -oE "depends on axioms: \[[^]]*\]" | sed -E "s/.*\[(.*)\]/\1/" | tr "," "\n" | sed "s/ //g" | grep -vxE "propext|Classical\.choice|Quot\.sound|" | sort -u)
  if [ $code -ne 0 ] || [ "$wanted" -eq 0 ] || [ "$reported" -ne "$wanted" ] || [ -n "$other" ]; then
    printf "%s\n" "$axioms" | grep -E "error" | head
    echo "FAIL axioms (lean exit $code; $reported of $wanted reported; other axioms: ${other:-none})"; status=1
  else
    echo "OK   axioms: $reported theorems, only propext, Classical.choice, Quot.sound"
  fi
fi
exit $status' check "$strict" "${modules[@]}"
