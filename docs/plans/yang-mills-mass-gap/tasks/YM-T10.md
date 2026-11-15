# YM-T10 — Exact small-lattice evaluation (Rust)

Tier T2–T3 · Workstream B · Depends on: YM-T05; PDD intent approved · Crate `Rust/lattice-gauge` · Review: T3

## Behaviour wanted
- Full enumeration of Z₂ (then Z_N) gauge theory on small lattices in d = 2, 3, 4: partition functions as exact polynomials in e^{±β}, with and without centre twists, with per-plaquette couplings.
- SU(2) in d = 2 by the character expansion, exact in terms of ratios of modified Bessel functions (symbolic coefficients, numeric evaluation with error bound).

## Fixtures
- 2D Z₂ gauge theory on a torus: closed form (plaquettes are independent up to the global constraint); derive it in the test comment and assert equality of polynomials.
- 3D Z₂ gauge theory on a 2×2×2 torus against a brute-force oracle written independently in Python.
- 2D SU(2): Wilson loop = (I₂(β)/I₁(β))^{area} in the infinite-volume normalization (state the β convention and check it against the convention table before asserting).
- Twist: in 2D, Z_tw/Z in closed form.

Serves YM-G03 phase B, YM-G04, YM-N03.
