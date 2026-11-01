# Numerical SU(2) laboratory endpoint

`POST /api/su2` accepts exactly `{"axis":[0,0,1],"t":1.5707963267948966,"slice":0}`.
The axis must have Euclidean norm within 1e-12 of one; it is returned unchanged,
never normalized. Numbers must be finite, |t| <= 1e6, and t is a dimensionless
rotation angle in radians. Missing/extra fields, wrong shapes/types, zero or
nonunit axes and out-of-range t return HTTP 400 with the gateway JSON error.
Only POST is supported (405 otherwise). Existing Host/Origin restrictions and
1 MiB declared/streamed request-body limits apply before computation.

The response uses convention `pauli-plus-i-v1`:

```
U = a0 I + i(a1 sigma1 + a2 sigma2 + a3 sigma3)
a0 = cos(t/2), (a1,a2,a3) = -axis sin(t/2)
alpha = [a0,a3], beta = [a2,a1]
matrix = [[[a0,a3],[a2,a1]],[[-a2,a1],[a0,-a3]]]
```

All complex numbers are `[real,imaginary]` pairs. `matrix` contains two rows,
each with two such pairs. The response also contains the original `t`, `axis`,
`coordinates:[a0,a1,a2,a3]`, `bloch:[x,y,z]`, `stereographic`, `slice`, `checks`,
`evidence_kind:"numerical"` and `formal_proof:false`. These are f64 observations,
never exact algebra or formal proof. Sign is retained: 0 and 4pi give I while
2pi gives -I within roundoff. Existing lie-groups quaternion matrices use
minus-i vector coordinates; this standalone module avoids silently reusing
that convention and adds no dependencies.

`bloch` is the action U rho U-dagger on the initial +x density matrix, on S2.
It is not all of SU(2), whose unit coordinates lie on S3. A positive z-axis
rotation initially moves +x toward +y; U and -U have the same Bloch action.

The stereographic map is `(a1,a2,a3)/(1-a0)`, excluding the north pole I.
`{status:"finite",point:[x,y,z]}` uses the stable expression `-axis*cot(t/4)`
to avoid subtracting nearly equal numbers near I. Exact f64 multiples of
4pi return `{status:"pole",point:null}`. This classification also covers
unrepresentable chart coordinates at extreme subnormal t; `pole` therefore
means numerical off-chart as well as the excluded mathematical pole. No
finite component is silently clipped or emitted as null. No near-pole visual
cutoff is applied; ordinary near-pole points can have very large coordinates.
At tiny t, rounded a0 may equal one even though the stable chart is finite.
The chart and rounded coordinate readout then agree only within f64 limits.

`slice` describes the whole a0=c intersection of the unit S3, independently
of whether the selected state lies on it. Its `coordinate` is input c:

- |c|<1: `{status:"sphere",radius:sqrt((1-c)*(1+c))}`.
- |c|=1: `{status:"point",radius:0}`.
- |c|>1: `{status:"empty",radius:null}` (valid input, including +/-1.2).

`checks.unitarity_error` is the maximum complex entry magnitude of U-dagger U-I;
`checks.determinant_error` is |det(U)-1|. Both use the returned matrix. Axis
acceptance tolerance is distinct from roundoff; tolerated axis discrepancies
are not hidden by normalizing them. Acceptance tests use 2e-11 absolute matrix
comparison tolerance and 3e-12 residual bounds for their unit-axis fixtures.

## Implementation and evidence

Approved story: [mathematical physics exploration](../../../user_stories/story__mathematical-physics-exploration.md),
scenarios 2–5. Exact approved request and read-only scoped plan:
[su2-request.txt](su2-request.txt), [su2-intent-plan.json](su2-intent-plan.json).
Plan ID `implement-rust-su-2-lab-for-approved-mathphysics-da41ac10`.
The planner recognized the gateway but omitted most compactly written constraints;
the complete request and prompt govern. Its generic scope/example suggestions
are resolved by the user's explicit exclusive scope and approved story.

Prompt updated before manual implementation, as explicitly authorized. The
existing gateway prompt owns main.rs and this module. PDD configuration and
architecture are unchanged within the exclusive write boundary. No intent apply,
provider generation, regenerated baseline, synchronization or generated story
contract is claimed. No frontend work or commits performed by this task.

Verification on 2026-10-31, from Rust/physics-gateway:

```sh
cargo fmt --check
CARGO_TARGET_DIR=/tmp/mathphysics-explore-rust/target cargo clippy --offline --locked --all-targets -- -D warnings
CARGO_TARGET_DIR=/tmp/mathphysics-explore-rust/target cargo test --offline --locked
```

All passed. Eight real-server HTTP tests: five existing gateway regressions and
three SU(2) tests with many cases, including explicit Pauli construction,
independent complex products and density-matrix Bloch oracle; 0/pi/2pi/4pi,
negative and maximum t; noncommuting x/y product signs; mixed axes; chart poles,
near-pole finite values; slice boundaries and outside; malformed/nonfinite/nonunit
inputs; methods, Host/Origin and declared/chunked body limits. Each fixture owns
and stops its own ephemeral server. No frontend integration result is claimed.

Evidence: `/tmp/mathphysics-explore-rust/{test.log,clippy.log,fmt-check.log}`;
plan raw stdout/stderr and exact request are in the same directory. Binary:
`/tmp/mathphysics-explore-rust/target/debug/physics-gateway`.
