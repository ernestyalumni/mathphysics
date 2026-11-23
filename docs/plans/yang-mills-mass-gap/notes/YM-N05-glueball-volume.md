# YM-N05: the 0++ glueball against box size at β = 2.5, with variational fits

Label: unreviewed (claude-opus55-strat, 2026-11-22). It answers YM-O04 §3: fit windows fixed in advance, binned-jackknife covariance, correlated fits, all windows reported, volume varied. Code: `Rust/lattice-gauge/examples/glueball_cpu.rs` and `glueball_fit.rs`; the GPU streams come from `glueball_gpu`. Raw streams and fit outputs are kept privately (about 770 MB), not in the repository.

## Method

- **Operators.** Zero-momentum spatial plaquettes after 0, 4, 8 and 16 APE steps (α = 0.5).
- **Correlators.** The 4×4 matrix C(t), vacuum-subtracted.
- **Variational projection.** Solve C(1)v = λC(0)v on the full sample and keep the top vector, held fixed.
- **Fits.** One-state cosh fits to v^T C(t) v with the jackknife covariance (20 bins per replica, replicas pooled).
- **Reporting rule.**
  - Fixed before the new runs: primary window [2,5], check window [3,6].
  - Acceptance rule (a window counts only if χ²/dof < 2): added after seeing the old 16⁴ stream, before any new data were analysed.
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

Earlier single streams with the same fit: at β = 2.3 on 10³×16, [1,5] gives 1.207(14) (Teper 1.216(23) on 10⁴); at β = 2.4 on 14³×16, [2,5] gives 0.936(28) (Teper 0.959(20) on 14⁴).

- **Agreement where the box is large.** At 20⁴ and at β = 2.4 on 16⁴ the primary window passes and agrees with Teper within one error. At 16⁴ (β = 2.5) only the later window passes; its value lies 1.1 combined errors above Teper's.
- **No result for L = 8, 10, 12.**
  - Every window fails. The effective masses keep falling with t (at L = 10: 0.635, 0.582, 0.529, 0.485 at t = 2–5).
  - Expected cause: a pair of flux tubes wound around a spatial direction (a torelon pair) lies below the glueball in small boxes. In the Nambu–Goto estimate E(L) ≈ σL − π/(3L), with a√σ ≈ 0.185 at β = 2.5 (measured below), the pair energy is about 0.28, 0.46, 0.63 and 0.79 at L = 8, 10, 12, 14. It crosses the glueball between L = 12 and 14, which is where acceptable fits begin.
  - The torelon energies were then measured directly (next section).
- **L = 6.** The primary window passes narrowly (χ²/dof 1.9), but the later window fails and the effective masses fall to 0.43 by t = 5. Read it as a mixture, not a glueball mass.

## Direct test: torelon energies (2026-11-22)

`examples/torelon_scan.rs` writes, from the same configurations, the plaquette streams above and the zero-momentum torelon operator: a₀ of the Polyakov loop winding once around a spatial direction, averaged over the three directions and the transverse positions, after the same smearing. The runs are new and independent: β = 2.5, T = 16 (20 for L = 20); 160 000, 228 000, 150 000, 105 000, 66 000 and 6 000 measurements at L = 8, 10, 12, 14, 16, 20. The same fit and reporting rule were used.

| L | E_T, window | 2E_T | a√σ from E_T = σL − π/(3L) | plaquette channel, late window (not pre-registered) |
|---|---|---|---|---|
| 8 | 0.1925(9), [2,5] | 0.385 | 0.201 | [5,8]: 0.470(24) (χ²/dof 1.6), still falling |
| 10 | 0.2368(13), [3,6] | 0.474 | 0.185 | [5,8]: 0.484(17) (0.18) |
| 12 | 0.3123(22), [3,6] | 0.625 | 0.183 | [4,7]: 0.591(18) (0.30) |
| 14 | 0.4023(27), [2,5] | 0.805 | 0.185 | [4,7]: 0.624(26) (0.45) |
| 16 | 0.4851(42), [2,5] | 0.970 | 0.186 | [2,5]: 0.740(10) (0.72), primary window passes |
| 20 | 0.660(17), [2,5] | 1.32 | 0.189 | [2,5]: 0.78(4) (1.5) |

- **The prediction holds.** 2E_T lies below the glueball (≈ 0.72–0.74) for L ≤ 12 and above it from L = 14.
  - At L = 10 the late-time plaquette mass equals 2E_T within errors.
  - At L = 8 it is still decreasing at t = 8.
  - At L = 12 it lies 1.6–1.9 errors below 2E_T. That would be consistent with an attractive torelon–torelon interaction, which was not measured.
  - At L = 14 the late window is noisy and not decisive.
- **String tension.** For L ≥ 10 the leading Nambu–Goto form gives a√σ = 0.183–0.189, so a√σ ≈ 0.185 at β = 2.5. This replaces the value quoted from memory; L = 8 lies off the form.
- **16⁴ again.** On this independent 16⁴ ensemble the primary window passes and gives 0.740(10). The ensemble of the table above gave 0.725(9) on [2,5] (rejected, χ²/dof 4.4) and 0.691(17) on [3,6]. The 16⁴ value therefore depends on the window at the level of 0.03–0.05.
  - The values at 14⁴, 16⁴ and 20⁴ on [2,5] (0.718–0.740) show no volume dependence at that resolution.
  - Teper's 16⁴ value 0.660(22) is lower than this ensemble's.

## What this says about box size

- At β = 2.5 the lowest 0++ level in the box is the infinite-volume glueball only once the torelon pair lies above it: 2E_T(L) > m, between L = 12 and 14, which is L√σ ≈ 2.4. The pair energy, not the glueball, sets the scale.
  - Below that, the lowest level the plaquette operator sees is a finite-volume state, and it is lower than the glueball. The finite-box "gap" is not monotone in L: the t = 2 value dips from 0.98 (L = 6) to 0.62 (L = 10) and returns to 0.73 (L ≥ 14).
- For the seed, this is a lower bound on the scale, not an estimate of it. Δ(n) ≤ 1/7480 asks for finite-size effects in log Z at the 10⁻⁴ level. Merely getting the right lowest level is a much weaker requirement.
- It is consistent with the program's statement that the seed must be met at a physical size of a few flux-tube widths, uniformly in β.
- It does not test the conditional-gap chain, which is proved.

## Status of the earlier labels (YM-O04 §3)

- The β = 2.5 value 0.66(4) quoted earlier was a last-time-point effective mass. It is replaced by 0.691(17) on 16⁴ ([3,6]) and 0.730(11) on 20⁴ ([2,5]).
- The β = 2.4 value 0.925(33) at t = 2 is replaced by 0.927(30) ([2,5], 16⁴).
- Not done:
  - two-state fits;
  - varying t₀ in the variational step;
  - an independent check of the GPU update against the CPU update at these volumes. That check was done at smaller volumes in YM-N01: GPU and CPU plaquettes agree.
