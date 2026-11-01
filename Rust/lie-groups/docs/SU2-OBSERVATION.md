# Reusable SU(2) observation API — LIB01/LIB02

`lie_groups::su2_observation::observe(axis, t, slice)` is a direct Rust API.
It returns typed `Observation`, `Chart`, `Slice`, `Checks`, and `Error` values.
Complex numbers are `[real, imaginary]` pairs; matrices are row-major pairs.
There is no serde, HTTP, frontend, IO, or chapter-specific dependency in the
module. The existing lie-groups dependency on qft-symbolics is unchanged.

The gateway owns strict JSON parsing and response serialization only. Its JSON
keys, errors and their validation order, numerical status, and
`formal_proof:false` are preserved. Main/body guards and HTTP tests are unchanged.

## Convention and domain

The observation convention is `U=a0 I+i a.sigma`, with
`a0=cos(t/2)` and `a=-axis*sin(t/2)`. Existing `su2::quaternion_matrix`
uses `U=q0 I-i q.sigma`: reuse it with `q=[a0,-a1,-a2,-a3]`.
That helper supports raw coordinates without normalizing or asserting group
membership. Tests independently expand Pauli matrices and evolve a +x spinor;
an additional test compares the sign conversion against the existing algebra.

Do not replace this path with `Su2Algebra::exp(axis*t)`: it uses the length of
its argument as phase and divides by that length, changing tolerated norm drift
and the legacy tiny-angle behavior. Do not force the strict `Su2` constructor:
its default tolerance is 32 eps, whereas this API accepts
`abs(hypot(axis)-1)<=1e-12`. The supplied axis is returned unchanged. Accepted
slightly nonunit axes produce approximate matrices with actual residuals;
these are numerical observations, not exact membership certificates.

All inputs must be finite, `abs(t)<=1e6`, and the axis must pass the norm
check. Slice values outside [-1,1] are valid empty slices. Stable stereography
uses `-axis*cot(t/4)`. `Chart::Pole` includes both exact f64 multiples of 4pi
and numerical off-chart results (any coordinate cannot be represented finitely).
Near-pole values are never clipped. Signed U/-U remain distinct even when their
Bloch actions agree. Bloch values act on initial +x, not the SU(2) S3 itself.

## Evidence and reproduction

From the repository root, with external build products:

```sh
CARGO_TARGET_DIR=/tmp/mathphysics-foundation-su2/target cargo test --offline --locked --manifest-path Rust/lie-groups/Cargo.toml
CARGO_TARGET_DIR=/tmp/mathphysics-foundation-su2/target cargo test --offline --locked --manifest-path Rust/physics-gateway/Cargo.toml
CARGO_TARGET_DIR=/tmp/mathphysics-foundation-su2/target cargo run --offline --locked --manifest-path Rust/lie-groups/Cargo.toml --example observation
```

Library tests cover independent Pauli/Bloch oracles, convention conversion,
0/pi/2pi/4pi, negative angles, finite near-pole versus off-chart/subnormal
values, axis drift on both sides of the tolerance, and slice boundaries.
The eight existing gateway HTTP tests pass independently and remain untouched.
Both affected crates pass scoped Clippy with warnings denied; changed Rust
sources pass rustfmt checking.

A disposable external crate at
`/tmp/mathphysics-foundation-su2/external-client` depends only on lie-groups
and successfully calls the API with fresh axis `[3,-5,7]/sqrt(83)`, checking
Bloch against independently written Rodrigues action. The maintained source
example is `Rust/lie-groups/examples/observation.rs`; the /tmp client is evidence,
not the canonical maintained example. Logs/build products are under
`/tmp/mathphysics-foundation-su2`.

The read-only PDD plan used the exact user request, intent ID
`implement-lib01-library-boundary-fix-and-lib02-e-66c144db`.
Its deterministic extraction missed constraints; the full request governs the
prompt-first implementation. No intent apply, provider-backed generation,
automated synchronization, commits, or broader LIB01/LIB02 program completion
is claimed. Only gateway/lie-groups lockfiles were refreshed offline for the
current shared dependency tree; parent may perform final lock synchronization.
