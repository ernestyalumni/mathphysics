# YM-N01 — SU(2) Wilson Monte Carlo baseline (Rust)

Tier T2 · Workstream B · Depends on: YM-T05; PDD intent approved · Crate `Rust/lattice-gauge` · Review: T2

Heat-bath (Kennedy–Pendleton or Creutz) plus over-relaxation for SU(2) as unit quaternions, d = 4, Wilson action S = β Σ_p (1 − ½ tr U_p). Seeded, reproducible, multi-threaded by checkerboard. Observables: plaquette average, Wilson loops, Creutz ratios, with jackknife errors and integrated autocorrelation times.

Validation: plaquette averages at β = 2.2, 2.3, 2.4, 2.5 on 12⁴ or larger against published values (locate a table in the literature, cite it, and state the tolerance); leading strong-coupling behaviour ⟨½ tr U_p⟩ ≈ β/4 at small β; weak-coupling 1 − 3/(4β) at large β. Negative control: a deliberately broken update (no accept step where one is needed) must fail the plaquette test.

Receipt: `PLAN/ledger/evidence/YM-N01/receipt.md` with seeds, parameters, revision, wall time, table of results against references.
