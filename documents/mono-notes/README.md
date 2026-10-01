# Mono notes

Thermodynamics is the first part. The source is shared. `master.tex` is US Letter, single column. `screen.tex` defines `\screenedition` and inputs `master.tex`, which selects a \(1200\,\mathrm{pt}\times 700\,\mathrm{pt}\) page and two columns. That is the same split as `Propulsion/documents/notes/preamble.tex`.

```sh
bash build.sh
```

The script runs pdfLaTeX twice for each edition and deletes the auxiliary files. It leaves `master.pdf` and `screen.pdf`.

`legacy/` is the 2008 Kittel--Kroemer manuscript and the 2015 `Propulsion/LaTeXandpdfs/thermo.tex`, split by topic and included in one progression. There is no Thermodynamics (Revisited) part. Correction environments mark a later check that disagrees with a displayed step. The archival files were not edited. New derivations are in `topics/`. Read `topics/90-resolutions.tex` for the arguments.

Build tools: `pdflatex` and the packages named in `preamble.tex`.
