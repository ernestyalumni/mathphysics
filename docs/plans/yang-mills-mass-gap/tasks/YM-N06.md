# YM-N06 — Spectral lower bound on the short-side cost

Tier T2 · Workstream B/C · Depends on: YM-F03 (short-side recurrence), YM-N05 (glueball and torelon tools) · Review: T2

## Why

The seed of the short-side recurrence is κ_β(n) ≤ 1/80. The direct measurement of κ_β(4) by integrating plaquette differences (YM-F03) has errors of 0.06 to 0.11 after ten GPU pairs. That error is set by the volume and not by κ, so more pairs do not resolve the threshold at n = 4, and n = 8 is out of reach.

## What

κ(n) is a function of the spectrum of the transfer operator in the short direction on the spatial torus (2n)³. Measured energies of the lightest states therefore bound κ(n) from below at a precision of about 10⁻³. Measure the torelon and 0++ glueball energies on (2n)³ × 16 for n = 4, 6, 8 at couplings on both sides of the expected threshold and report the lower bound.

## Acceptance

- The pre-registration in `notes/YM-N06-spectral-seed.md` is committed before any stream is analyzed.
- Energies come from `glueball_fit` with the windows fixed in that note; every printed window is reported.
- The note states the bound as a necessary condition for the seed and never as evidence that the seed holds.

## Rules

- Monte Carlo certifies nothing.
- States left out of the spectral sum only raise κ; an estimate that adds states by a dispersion assumption is labeled as an estimate.
