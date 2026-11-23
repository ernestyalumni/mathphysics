# YM-N05: the 0++ glueball against box size at β = 2.5, with variational fits

Label: reviewed by Codex (YM-O05, `notes/YM-O05-codex-review-2026-11-22.md` §4–5: all eight fit outputs reproduced byte for byte); corrections applied 2026-11-23. The torelon section was added after that review and is unreviewed. It answers YM-O04 §3: fit windows chosen before the new runs, binned-jackknife covariance, correlated fits, volume varied. Code: `Rust/lattice-gauge/examples/glueball_cpu.rs` and `glueball_fit.rs`; the GPU streams come from `glueball_gpu`. Raw streams and fit outputs are kept privately (about 770 MB), not in the repository.

## Method

- **Operators.** Zero-momentum spatial plaquettes after 0, 4, 8 and 16 APE steps (α = 0.5).
- **Correlators.** The 4×4 matrix C(t), vacuum-subtracted.
- **Variational projection.** Solve C(1)v = λC(0)v on the full sample and keep the top vector, held fixed.
- **Fits.** One-state cosh fits to v^T C(t) v with the jackknife covariance (20 bins per replica, replicas pooled), over windows of 2 to 5 time points.
- **What the errors include.** In the jackknife refits both v and the inverse covariance W are held fixed. The errors do not propagate reselection of the basis, re-estimation of W, the choice of window, or excited-state bias. W comes from a noisy inverse, so χ²/dof is a fit-quality indicator, not a calibrated p-value. An independent re-analysis (YM-O05 §4) found no material shift at 20⁴ from recomputing v per deletion or from equal-size blocks.
- **Reporting rule.**
  - Fixed before the new runs: primary window [2,5], check window [3,6].
  - Acceptance rule (a window counts only if χ²/dof < 2): added after seeing the old 16⁴ stream, which is part of the combined 16⁴ result. This is a pilot-informed reporting rule, not a blind validation.
  - The reported value is the primary window if it passes, otherwise the check window. Otherwise there is no result.

## Results

| β | L³ × T | Measurements (replicas) | [2,5] (χ²/dof) | [3,6] (χ²/dof) | Reported | Teper, hep-th/9812187 Table 10 |
|---|---|---|---|---|---|---|
| 2.5 | 6³×16 | 80 000 (2) | 0.969(22) (1.9) | 0.84(6) (4.2) | 0.969(22), see below | |
| 2.5 | 8³×16 | 120 000 (3) | 0.714(7) (31) | 0.630(12) (11) | none | |
| 2.5 | 10³×16 | 120 000 (3) | 0.624(7) (31) | 0.565(9) (7.0) | none | |
| 2.5 | 12³×16 | 120 000 (3) | 0.677(7) (11) | 0.636(12) (3.8) | none | |
| 2.5 | 14³×16 | 90 000 (3) | 0.726(8) (3.6) | 0.687(18) (1.2) | 0.687(18) | |
| 2.5 | 16³×16 | 100 000 (5) | 0.725(9) (4.4) | 0.691(17) (1.8) | 0.691(17) | 0.660(22) on 16⁴ |
| 2.5 | 20³×20 | 52 000 (5) | 0.730(11) (0.43) | 0.711(22) (0.02) | 0.730(11) | 0.723(23) on 20⁴ |
| 2.4 | 16³×16 | 20 000 (1) | 0.927(30) (0.01) | 0.936(61) (0.6) | 0.927(30) | 0.936(23) on 16⁴ |

Earlier single streams with the same fit and 40 bins: at β = 2.3 on 10³×16, [1,5] gives 1.207(14) (20 bins: 1.208(14); Teper 1.216(23) on 10⁴); at β = 2.4 on 14³×16 (a different ensemble from the 16⁴ row), [2,5] gives 0.936(28) (20 bins: 0.936(29); Teper 0.959(20) on 14⁴).

- **Matched-volume comparisons.** At 20⁴ and at β = 2.4 on 16⁴ the primary window passes and agrees with Teper within one error. At 16⁴ (β = 2.5) only the check window passes; it lies 1.1 combined errors above Teper's. These check the measurement at matched volumes; they do not validate the one-state fit model.
- **No result for L = 8, 10, 12.**
  - Both designated windows fail. The periodic (cosh) effective masses keep falling with t (at L = 10: 0.635, 0.582, 0.529, 0.485 at t = 2–5). Other windows were not promoted after the fact.
  - Hypothesis: a centre-neutral pair of oppositely wound spatial flux tubes (torelons) mixes with the scalar channel and lies below the glueball in small boxes. The truncated long-string form for non-interacting strings, 2aE ≈ 2(a²σL − π/(3L)) with Teper's a√σ = 0.1834(26), gives 0.28, 0.46, 0.63 and 0.79 at L = 8, 10, 12, 14. Failure of two one-state fits does not identify the lowest level, and short strings, interactions and finite T are not included.
  - The single-torelon energies were then measured directly (next section). Pair operators, which would test the mixing itself, were not measured.
- **L = 6.** The primary window passes narrowly (χ²/dof 1.9), but the check window fails and the cosh effective masses fall to 0.43 by t = 5. Passing one window is not an identification of the state; it is not a glueball mass.

## Direct test: torelon energies (2026-11-22)

`examples/torelon_scan.rs` writes, from the same configurations, the plaquette streams above and the zero-momentum single-torelon operator (centre-charged, so it does not mix with the plaquettes): a₀ of the Polyakov loop winding once around a spatial direction, averaged over the three directions and the transverse positions, after the same smearing. The runs are new and independent: β = 2.5, T = 16 (20 for L = 20); 160 000, 228 000, 150 000, 105 000, 66 000 and 6 000 measurements at L = 8, 10, 12, 14, 16, 20. The same fit and reporting rule were used.

| L | E_T, window | 2E_T (no interaction) | a√σ from E_T = σL − π/(3L) | plaquette channel, late window (not pre-registered) |
|---|---|---|---|---|
| 8 | 0.1925(9), [2,5] | 0.385 | 0.201 | [5,8]: 0.470(24) (χ²/dof 1.6), still falling |
| 10 | 0.2368(13), [3,6] | 0.474 | 0.185 | [5,8]: 0.484(17) (0.18) |
| 12 | 0.3123(22), [3,6] | 0.625 | 0.183 | [4,7]: 0.591(18) (0.30) |
| 14 | 0.4023(27), [2,5] | 0.805 | 0.185 | [4,7]: 0.624(26) (0.45) |
| 16 | 0.4851(42), [2,5] | 0.970 | 0.186 | [2,5]: 0.740(10) (0.72), primary window passes |
| 20 | 0.660(17), [2,5] | 1.32 | 0.189 | [2,5]: 0.78(4) (1.5) |

- **Consistent with the hypothesis, not a test of mixing.** Twice the measured single-torelon energy lies below the glueball (≈ 0.72–0.74) for L ≤ 12 and above it from L = 14.
  - At L = 10 the late-window plaquette mass equals 2E_T within errors. The window was chosen after seeing the data.
  - At L = 8 it is still decreasing at t = 8.
  - At L = 12 it lies 1.6–1.9 errors below 2E_T. That would be consistent with an attractive torelon–torelon interaction, which was not measured.
  - At L = 14 the late window is noisy and not decisive.
- **String tension.** For L ≥ 10 the truncated long-string form gives a√σ = 0.183–0.189. Teper's Table 10 gives 0.1834(26) on 16⁴ and 0.1881(28) on 20⁴; ours are 0.186 and 0.189 at L = 16 and 20. L = 8 lies off the form.
- **16⁴ again.** On this independent 16⁴ ensemble the primary window passes and gives 0.740(10). The ensemble of the table above gave 0.725(9) on [2,5] (rejected, χ²/dof 4.4) and 0.691(17) on [3,6]. The 16⁴ value therefore depends on the window at the level of 0.03–0.05.
  - The values at 14⁴, 16⁴ and 20⁴ on [2,5] (0.718–0.740) show no volume dependence at that resolution.
  - Teper's 16⁴ value 0.660(22) is lower than this ensemble's.

## Relation to the seed

No defect was measured in this scan, so it gives no statement about where a seed holds.

What the chain does give is a necessary condition (YM-O05 §5):
- In the defect telescope at side n, take time as the first doubled direction. Its cost p carries weight 8, so 8p ≤ Δ(n).
- The purity bound gives e^{−n m_box} ≤ e^p − 1, where m_box is the true first excitation of the transfer operator on the spatial torus of side n.
- So a seed Δ(n) ≤ 1/7480 requires n · m_box ≥ −log(e^{1/59840} − 1) = 10.9994.

A seed box must therefore be at least about eleven of its own correlation lengths wide.
- A certified upper bound on m_box below 10.9994/n would exclude that n. Effective masses from finite-T fits are not such a certificate.
- Passing the condition bounds neither the multiplicities nor the defect.
- The scan does not test the conditional-gap chain, which is proved.

## Status of the earlier labels (YM-O04 §3)

- The β = 2.5 value 0.66(4) quoted earlier was a last-time-point effective mass. It is replaced by 0.691(17) on 16⁴ ([3,6]) and 0.730(11) on 20⁴ ([2,5]).
- The β = 2.4 value 0.925(33) at t = 2 is replaced by 0.927(30) ([2,5], 16⁴).
- Not done:
  - two-state fits;
  - varying t₀ in the variational step;
  - an independent check of the GPU update against the CPU update at these volumes. That check was done at smaller volumes in YM-N01: GPU and CPU plaquettes agree.
