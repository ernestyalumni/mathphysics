#!/bin/bash
# Build the letter master and the 1200pt by 700pt screen edition.
# Intermediates are removed. PDFs stay next to the source.
set -euo pipefail
cd "$(dirname "$0")"
pass() {
  local job="$1"
  shift
  pdflatex -interaction=nonstopmode -halt-on-error -jobname="$job" "$@"
  pdflatex -interaction=nonstopmode -halt-on-error -jobname="$job" "$@"
}
pass master master.tex
pass screen screen.tex
rm -f master.aux master.log master.out master.toc \
      screen.aux screen.log screen.out screen.toc \
      *.fls *.fdb_latexmk *.synctex.gz
