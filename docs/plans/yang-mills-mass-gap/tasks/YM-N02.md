# YM-N02 — Twist ratios on vortex containers

Tier T3 design, T2 implementation · Workstream B/C · Depends on: YM-N01, YM-G04 · Blocks: milestone G4 · Review: T3

## Goal
Measure R(N, β) = Z_twisted / Z for a centre twist confined to a container of linear scale N, for β in and around the scaling window and N as large as the hardware allows, and see what function of N it follows.

## Do
1. Design note first (reviewed before any long run): exact container geometry (start with the standard 't Hooft twist on an asymmetric torus, then a Mack–Petkova style container with fixed boundary), the estimator (chain of intermediate ensembles flipping one twisted plaquette stack at a time, or an equivalent published method; locate and cite the SU(2) twist free-energy papers of de Forcrand and von Smekal and of Kovács and Tomboulis), error budget.
2. Validate on cases with exact answers from YM-T10 (2D; Z₂).
3. Run. Fit R against 1 − c/log N, 1 − c N^{−p}, and exp(−c e^{−σN²}) style forms; report which are excluded.
4. Repeat with randomly weakened plaquettes (b_p ∈ [0, β]) to probe the uniformity that the renormalization-group corrections would require.

## Acceptance
Validation passes before production. Every number has an error bar and an autocorrelation estimate. The conclusion section says what the data exclude, not what they suggest.
