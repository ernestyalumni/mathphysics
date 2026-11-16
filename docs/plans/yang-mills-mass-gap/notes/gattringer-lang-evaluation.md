# Gattringer–Lang, *Quantum Chromodynamics on the Lattice* (2010): evaluation for this program

Author: strategist (Claude Opus 5.5), 2026-11-16. Requested by Ernest. Status: evaluation from the PDF's text layer; OCR (Marker + Nougat) queued behind another job on the GPU.

Source: Ernest's copy, `Data/Public/books/Physics/[Christof_Gattringer,_Christian_B._Lang]_Quantum_C(BookFi.org).pdf`, 351 PDF pages. Parse workspace: `Data/Public/books/Physics/GattringerLang-QCDOnTheLattice-2010/` (text layer, per-page text, page map, contents, 978 display-equation tags, `STATUS.md`). Page numbers below are **printed** pages, located in the text layer by the strategist (PDF page = printed page + 8 to 12, varying; use `text-layer/page_map.json`).

## Verdict

**Useful, but not as a guide to OpenAI's proof.** It is a physicist's textbook on lattice QCD as a numerical method. It contains no constructive field theory: no Bałaban renormalization group, no Wightman axioms, no uniform estimates, no proof of a mass gap. It assumes a discrete spectrum where that is convenient (p. 5) and sends the reader to Reed–Simon for the continuous spectrum (p. 3).

What it is good for, in order of value to us:

| Use | Where | Serves |
|---|---|---|
| **Pointers to the rigorous literature on exactly our open hypothesis.** It cites Osterwalder–Seiler 1978 and Seiler's LNP 159 for the Hilbert space of lattice gauge theory, and Lüscher 1977 for a self-adjoint, strictly positive transfer matrix | p. 21, refs. [12]–[15] of ch. 1 | The hypothesis "T ≥ 0 in every direction" left unproved in mono-notes Prop. `ym:prop:defect`; packet YM-G02. See the study-route note §3 |
| Conventions: Wilson action (β/3) Σ Re tr(1 − U_p), β = 6/g² for SU(3), general SU(N) at its eq. (3.93) | pp. 44–45 | Same form as the Lean's Wilson energy; agrees with the coupling dictionary β = 2m/g² now in the notes (claim C22) |
| Haar measure, invariance, the SU(3) integrals (3.24)–(3.29), integrating out a shared link (Fig. 3.1) | §3.1, pp. 44–49 | YM-T09 (matrix integration), YM-T15 test fixtures |
| Gauge invariance, gauge fixing, Wilson and Polyakov loops; glueball operators from loops | §§3.2–3.3, pp. 49–58 | Observables for YM-N01/N02 |
| Strong-coupling expansion of the Wilson loop (area law at small β), with the remark that Wilson's derivation was proved rigorously by Osterwalder–Seiler | §3.4.1, pp. 59–62 | Route R3 (strong-coupling side) in `03-MASS-GAP-STRATEGY.md` §5 |
| Monte Carlo: Metropolis, heat bath, over-relaxation, error analysis | ch. 4, pp. 73–101 | YM-N01 (SU(2) baseline), YM-N02 (twist ratios): the reference for implementing T06–T07 |
| Extracting masses from correlators: effective-mass curves, fits, and the spectral density as a Laplace transform C(t) = ∫ρ(E)e^{−tE} dE | §6.3, eq. (6.64), pp. 143–150 | The numerical counterpart of Prop. `ym:prop:gap-decay` and Lemma `ym:lem:layer`; how a lattice "glueball mass" is measured |
| Real-space renormalization group: blocking, RG flow of couplings, saddle-point analysis | §9.3, pp. 227–236 | Physical intuition for OpenAI's block-spin steps (unit III.1 in Armstrong's map); not their method |
| Centre (Z₃) symmetry, Polyakov loop, deconfinement | §12.1, pp. 301–311 | Centre symmetry is the symmetry behind our 't Hooft twists (route R1) |

Not useful for us: chapters 5, 7, 8, 10, 11 (fermions, chiral symmetry, dynamical quarks, hadron structure) and §§12.2–12.3 (chemical potential). Pure Yang–Mills has no quarks.

## Relation to OpenAI's proof

The book and the Lean start from the same object, Wilson's lattice measure, and then go opposite ways. The book takes it as a computational device and computes. The Lean proves that a scaling limit of it exists. For reading the Lean, the book supplies vocabulary and sign conventions (links, plaquettes, Haar measure, β), not arguments. For the mass gap, the book's chapter 6 is how physicists *measure* the gap numerically; it says nothing on how to *prove* one.

## Comparison with the other lattice books on disk

Topic counts in the text layers (DeGrand–DeTar 2006 and Smit 2002 PDFs; Montvay–Münster is DjVu and was not counted):

| Topic | Gattringer–Lang | DeGrand–DeTar | Smit |
|---|---:|---:|---:|
| transfer matrix | 2 | 9 | 2 |
| "mass gap" | 0 | 0 | 3 |
| area law | 0 | 2 | 4 |
| vortices | 0 | 0 | 4 |
| Haar measure | 25 | 4 | 2 |
| heat bath / over-relaxation | 26 / 16 | 16 / 7 | 0 / 0 |
| Polyakov loop | 39 | 17 | 0 |
| finite volume / finite size | 5 | 39 | 22 |

Reading of the table: Gattringer–Lang is the best of the three for Haar integration and Monte Carlo; Smit is the most conceptual about confinement and the gap; DeGrand–DeTar is strongest on finite-volume analysis. None replaces the rigorous sources in the study-route note.

## A new item for the watch list

While checking the Osterwalder–Seiler citation, a search found arXiv:2606.19362 (Faizal, Shabir, 9 June 2026, hep-lat/hep-th), "Reflection-Positive Construction of a Four-Dimensional SU(N) Yang-Mills Theory with Mass Gap and Confinement", with four companion papers in Int. J. Geom. Meth. Mod. Phys. (2026). The abstract claims a uniform transfer-operator gap from strong-coupling character expansions, carried to the continuum by a multiscale analysis, plus universality and agreement with the weak-coupling continuum. **Label: source-claims, unreviewed.** Only the abstract was read. Clay still listed the problem as unsolved on 2026-11-16 (Codex's check, `07-CONTINUATION-ROADMAP.md`). The crossing from strong to weak coupling is the step with no known small parameter (`03-MASS-GAP-STRATEGY.md` §2), so that is where a reader should look first. Packet YM-W04.

## Not done

- OCR not run yet (GPU busy with another session's Marker job); no equation-level reconciliation; no figure extraction.
- Equation numbers and page locators above come from the PDF text layer, which drops some symbols; formulas were not transcribed into our notes.
- Lüscher 1977, Osterwalder–Seiler 1978 and Seiler's LNP 159 were **not read**; their content is known here only from titles, abstracts and how Gattringer–Lang and later papers cite them.
