# PROGRESS.md — Srednicki Cadabra2 Chapter Status

## Completed Chapters

| Chapter | Topic | Script | Export/PDF | Notes |
|---------|-------|--------|------------|-------|
| Ch.34 | Left/Right Weyl Spinors | ✅ `ch34_left_right_spinors.py` | ✅ | Base chapter — metric (+,-,-,-), ε^{12}=+1 established |
| Ch.35 | σ-matrix Algebra | ✅ `ch35_sigma_algebra.py` | ✅ | σ^μ completeness, trace identities |
| Ch.36 | Weyl Lagrangian | ✅ `ch36_weyl_lagrangian.py` | ✅ | EOM, kinetic term, mass term |
| Ch.60 | Spinor Helicity / MHV | ✅ `ch60_spinor_helicity.py` | — | Twistors, angle/square brackets, polarization vectors, Fierz; no export script yet |

## Completed Chapters (cont.)

| Chapter | Topic | Script | Export/PDF | Notes |
|---------|-------|--------|------------|-------|
| Ch.37 | Canonical Quantization of Spinor Fields I | ✅ `ch37_canonical_quantization.py` | ✅ | CARs, mode expansion, basis spinors u^s v^s, normal-ordered H, spin-statistics |
| Ch.38 | Spinor Technology | ✅ `ch38_spinor_technology.py` | ✅ | σ completeness, 2-trace, 4-trace (+2iε sign), index gymnastics, Fierz |

## Renamed/Corrected Chapters (2026-03-30)

The following files were previously mislabeled by a prior AI. File renames, internal
reference fixes, and AGENTS.md roadmap corrections were applied on this date.

| Old Name | New Name | Correct Chapter | Content |
|----------|----------|-----------------|---------|
| `ch48_massless_spinor_helicity.py` | `ch50_massless_spinor_helicity.py` | Ch.50 | Massless particles & spinor-helicity (Parke-Taylor) |
| `ch48_spinor_helicity.py` | `ch50_spinor_helicity_core.py` | Ch.50 | Core spinor-helicity formalism |
| `ch48_spinor_helicity.tex/.pdf` | `ch50_spinor_helicity.tex/.pdf` | Ch.50 | LaTeX notes (title + framed box references corrected to Ch.50) |
| `ch48_export_latex.py` | `ch50_export_latex.py` | Ch.50 | Export script (all ch48 refs corrected to ch50) |
| `ch76_nonabelian_gauge.py` | `ch69_nonabelian_gauge.py` | Ch.69 | Nonabelian gauge theory (NOT Ch.76 = Anomaly in QED) |
| `ch77_group_representations.py` | `ch70_group_representations.py` | Ch.70 | Group representations (NOT Ch.77 = Fujikawa anomaly) |
| `ch80_feynman_rules_nonabelian.py` | `ch72_feynman_rules_nonabelian.py` | Ch.72 | Feynman rules for nonabelian gauge theory (NOT Ch.80 = Large-N limit) |

**Note:** The real Ch.80 (Large-N limit / planar diagrams / double-line notation) is still
missing. See NEXT_STEPS.md for priority ordering.

**Actual content of mislabeled chapters:**
- Srednicki Ch.48 = Spin Sums (Yukawa trace methods, spin-averaged amplitudes) — NO script yet
- Srednicki Ch.76 = The Anomaly in QED (axial current, Ward identities) — NO script yet
- Srednicki Ch.77 = Anomalies and the Path Integral (Fujikawa method) — NO script yet
- Srednicki Ch.80 = The Large-N Limit (planar diagrams, double-line notation) — NO script yet

## Docker Verification Status (2026-03-30)

| Script | Docker Result | Notes |
|--------|--------------|-------|
| `ch50_massless_spinor_helicity.py` | ✅ EXIT 0 | All sections pass |
| `ch50_spinor_helicity_core.py` | ✅ EXIT 0 | All sections pass |
| `ch69_nonabelian_gauge.py` | ✅ EXIT 0 | All sections pass |
| `ch70_group_representations.py` | ✅ EXIT 0 | All sections pass |
| `ch72_feynman_rules_nonabelian.py` | ✅ EXIT 0 | All sections pass |
| `ch60_spinor_helicity.py` | ⚠️ §60.E fails | §60.A–D pass; §60.E Fierz identity preexisting bug (completeness relation fails for hardcoded 4-component spinors); see below |

**ch60 bug details (2026-03-30):**
- Fixed crash: `pol_vector_plus()` called with 9 args instead of 4 (line 603)
- Fixed physics bug: `pol_vector_minus()` used wrong spinor; corrected to eps_-^mu = [k|sigma^mu|q> / (sqrt2 <kq>) in (+,-,-,-) convention
- Fixed transversality check: metric sign in k_lower_d was wrong (need (+,-,-,-) lowering)
- Fixed normalization assertions: eps_+ . eps_+^* = -1 (was "should be 0"), eps_+ . eps_-^* = 0 (was "should be -1")
- §60.E preexisting bug: Fierz identity fails because `massless_spinors` 4-component spinors do not
  satisfy the completeness relation (-k/ = |k]<k| + |k><k|); requires rewriting §60.E. Out of scope.

## Not Started

| Chapter | Topic | Priority | Notes |
|---------|-------|----------|-------|
| Ch.37 export | LaTeX export for Ch.37 | high | needs `ch37_export_latex.py` |
| Ch.38 export | LaTeX export for Ch.38 | high | needs `ch38_export_latex.py` |
| Ch.33 | Lorentz group representations | high | missing start of spinor chain; prerequisite for Ch.34 |
| Ch.50 export | PDF recompile (title fix) | high | `ch50_spinor_helicity.tex` title fixed; needs pdflatex rerun |
| Ch.60 export | LaTeX export for Ch.60 | medium | `ch60_spinor_helicity.py` exists, export script missing |
| Ch.80 | Large-N limit / color decomposition | medium | real Ch.80 content; no script exists |
| BCFW recursion | On-shell recursion | low | after Ch.50 + Ch.60 are solid |

## Last Worked On

**2026-03-30** — Phase 1 relabeling complete. Verified Srednicki TOC via .mmd parse:
- Ch.48 = Spin Sums; Ch.50 = Massless Particles & Spinor-Helicity
- Ch.69 = Nonabelian Gauge Theory; Ch.70 = Group Representations; Ch.72 = Feynman Rules
- Ch.76 = Anomaly in QED; Ch.77 = Fujikawa method; Ch.80 = Large-N limit
Renamed 7 files (4 ch48->ch50, 1 ch76->ch69, 1 ch77->ch70, 1 ch80->ch72).
Fixed all internal chapter refs in .tex and .py files. Updated AGENTS.md roadmap.
Branch: feat/srednicki-relabel-and-verify.

**2026-03-27** — Verified ch37 + ch38 scripts run clean in `cadabra2-ubuntu:24.04` Docker.
Fixed bugs: ch37 `Ex()` free-index-sum crash (mode expansion), ch38 unavailable
`DiagonalMetric`/`EpsilonTensor` API calls, ch38 4-trace e sign (correct: +2i for
Tr[sigma sigma-bar sigma sigma-bar] with mostly-plus metric). All 9 numerical checks pass.
Wrote `ch37_export_latex.py` and `ch38_export_latex.py`; produced 6-page PDFs for both
chapters. Committed to `feat/srednicki-ch37-ch38-breakdown`, pushed.

## Notes / Caveats

- `ch36_weyl_lagrangian.aux/.out/.toc` are LaTeX build artifacts in the working
  tree — not committed, safe to delete or ignore
- `MHV_research_survey.md` is a research notes file — not a chapter script
- `pol_fix_patch.py` is a one-off utility patch — not part of the chapter sequence
- Ch.60 was done out of order (before Ch.37/38) as a preview of MHV machinery
- `ch38_lsz_feynman_rules.py` appears to be an alternate/early draft — not the canonical ch38 script
