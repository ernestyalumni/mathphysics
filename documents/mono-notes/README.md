# Mathematical Physics — collected notes and worked problems

This is the single mathematical-physics master. The established
`documents/mono-notes/` directory now covers mathematics and physics together,
expanding the previous thermodynamics-only master. `master.tex` assembles
subject files; `screen.tex` uses the same content with a screen layout.

## Read and build

- `master.pdf`: US Letter, single column, with wide pages for tensor-network diagrams.
- `screen.pdf`: 1200 × 700 TeX points, two columns; the tensor-network chapter uses the full page width.

From this directory:

```sh
bash build.sh
```

To publish outside the source tree:

```sh
bash build.sh --output /absolute/path/to/artifacts
```

Requires Python 3, pdfLaTeX, and the packages in `preamble.tex`. The builder
runs three passes per edition with shell escape disabled in a temporary
working directory. It publishes both PDFs only after both compile, together
with `build-report.json` containing source hashes, PDF hashes, page counts,
reference/citation warnings and overfull-box measurements. Intermediates stay
outside the repository. Generated PDFs and build reports should not be staged.

## Subject organization

The master interleaves material from the two former discipline manuscripts:

1. Linear algebra, algebra, sets, probability, and integral transforms.
2. Categories, computational algebraic geometry, topology, and graphs.
3. Manifolds, Lie groups, forms, connections, holonomy, and Morse theory.
4. Classical and geometric mechanics, rigid-body flight, and waves.
5. Thermodynamics and statistical mechanics.
6. Quantum mechanics, tensor networks, conformal and quantum field theory,
   the Standard Model, and general relativity.

These are collected working notes, with uneven coverage and unfinished
passages inherited from the manuscripts. A successful LaTeX build does not
certify the mathematics. Existing long displays and diagrams still need
layout review; the build report makes remaining overflow explicit.

`topics/rigid-body-foundations.tex` adds the bridge from rotations and
momentum balances to a constant-mass six-DOF model, with a Holm--Schmah--Stoica
reading sequence. It states frame and physical assumptions locally. The
Holm page-reconciled parse is beside the PDF, in
`Holm-Schmah-Stoica-GeometricMechanicsAndSymmetry-2009/ocr-compare/reconciled/resolved.md`.
These notes do not copy the book.

## Source and provenance

The historical `LaTeX_and_pdfs/AGAT_dump.tex` and `DGDT_dump.tex` remain
unchanged. Their mathematical bodies are included as editable subject files
in `legacy/mathematics/`, split at original part headings. They retain original
prose, exercises, and unfinished calculations. Source hashes and import
transformations are recorded in `legacy/mathematics/provenance.json`.

The import removes archival title pages, abstracts, donation text, comments,
and embedded column-layout commands. It prefixes labels and bibliography
keys by source; repairs blank paragraphs inside displays/TikZ and a delimiter
alignment error; fixes the RecordTypesProduct and RightMPSoperation reference typos; and gives
repeated equation labels distinct targets. It does not resolve mathematical
errors or complete unfinished proofs. Both archives and thermodynamics use
one bibliography in `backmatter.tex`.

`topics/thermodynamics.tex` preserves the existing thermodynamic progression,
including its Kittel--Kroemer and Propulsion excerpts and the correction
arguments in `topics/90-resolutions.tex`.

`import_legacy.py` reproduces the initial import from the two archives. It
**overwrites** the imported topic files, so do not rerun it after editing those
files without comparing and preserving the edits first. New work belongs in
shared topic files, with no recurring need to regenerate the archives.

Other manuscripts, including `EM_dump.tex`, `Analysis.tex`, the Apostol
solutions, Frankel problem files, and Srednicki/amplitudes notes, remain
available in their established locations. They have not yet been merged into
this master; incorporate passages by subject after checking duplication and
references. Do not silently replace personal exercise work with generated
solutions.

## Extending the master

Read `STYLE.md` before adding derivations. Add a topic under `topics/`, input it
at the appropriate subject location in `master.tex`, and add its citations to
`backmatter.tex`. Both editions then receive the same mathematics. Keep
exercises and solutions inline, state assumptions before a derivation, show
the algebra, and preserve historical content with explicit corrections.
