# Handbook reading guide: finite-torus experiments

## Result

Both editions now have a shared entry section, **What can a finite torus tell us about Yang–Mills theory?** (`ym:sec:finite-guide`), after the QFT foundations and before the full lattice chapter. Start on PDF page **540** (portrait) or **291** (screen). The front status page links to it. The section derives the binary plaquette sum, flux projection, enumerator identity, first-axis purity bound and a one-slice spectrum. It compares these with scalar and winding-channel measurements. The detailed chapter distinguishes a bound on physical length from a leading-log estimate and accounts for repeated link incidences when moving twists. The general kernel definition now requires symmetry and the trace-class argument invokes positivity.

## Labels

Exposition of known lattice methods and previously reviewed finite results, not a new continuum theorem. The cohomology and flux check is for the sign group. The stated finite-box formal theorem is for Z₂; general compact-group positivity is a paper argument. The pair interpretation of the scalar fit is a hypothesis, not an eigenstate identification.

## Evidence

- T10 `e7f58af`, `YM-T10-z2-exact.md`: enumerator conventions and flux values; O06 `b71a8e9`, `YM-O06-codex-review-2026-11-25.md`: independent sector-sum check, threshold correction, twelve replayed fits and compact-group positivity review.
- F02 `4a29e20`, independent review `6afbb49`, `YM-F02-codex-review-2026-11-25.md`: first-axis geometry, β range, ratio and its zero case. No mass eigenvalue inferred from a plaquette average.
- Independently evaluated the stored untwisted 2⁴ histogram: its counts sum to 2⁴⁵, twice the mean weight is 96, and at β=1.3926034652×10⁻⁶ the absolute discrepancy between −log Z and 1/7480 is below 2×10⁻¹⁵.
- From `documents/mono-notes`, `bash build.sh --output "$OUT"`: 630/338 pages, no LaTeX reference/citation warnings; 47/3 overfull boxes across the collected volume. Rechecked the source and PDF hashes against the build report and the section destinations with `pdfinfo -dests`; checked the rendered section text in each edition. Only the two final reading PDFs, not build intermediates or the report, accompany the sources.

## Not done / not checked

No new simulation, Lean proof, continuum-limit argument, or claim that small twist cost implies a seed. The existing overfull boxes in the collected notes were not a layout-cleanup target.
