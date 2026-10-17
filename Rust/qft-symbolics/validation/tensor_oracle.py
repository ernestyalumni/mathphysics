"""Independent real Cadabra tensor transformations; compare to Rust's case verdicts."""
import json
import sys
import cadabra2 as c
from cadabra2 import Ex

c.__cdbkernel__ = c.create_scope()
c.Indices(Ex(r"{a,b,c,d,m,n}"), Ex(r"L, position=fixed"))
c.Integer(Ex(r"{a,b,c,d,m,n}"), Ex("0..3"))
c.Metric(Ex(r"g_{a b}"))
c.InverseMetric(Ex(r"g^{a b}"))
c.KroneckerDelta(Ex(r"\delta^{a}_{b}"))
c.Symmetric(Ex(r"S_{a b}"))
c.AntiSymmetric(Ex(r"K_{a b}"))
c.Symmetric(Ex(r"S^{a b}"))
c.AntiSymmetric(Ex(r"K^{a b}"))
c.AntiSymmetric(Ex(r"K_{a}^{b}"))
c.AntiSymmetric(Ex(r"K^{a}_{b}"))

cases = [
    ("symmetric_swap", r"S_{a b}-S_{b a}", True),
    ("antisymmetric_swap", r"K_{a b}+K_{b a}", True),
    ("antisymmetric_wrong_sign", r"K_{a b}-K_{b a}", False),
    ("symmetric_antisymmetric_contraction", r"K_{a b} S^{a b}", True),
    ("antisymmetric_identical_vectors", r"K_{a b} V^{a} V^{b}", True),
    ("antisymmetric_distinct_vectors", r"K_{a b} V^{a} W^{b}", False),
    ("metric_lowering", r"g_{a b} V^{b}-V_{a}", True),
    ("inverse_metric", r"g_{a b} g^{b c}-\delta^{c}_{a}", True),
    ("metric_trace_four", r"g_{a b} g^{a b}-4", True),
    ("metric_trace_wrong_two", r"g_{a b} g^{a b}-2", False),
    ("dummy_variance_orientation", r"A^{a} B_{a}-A_{b} B^{b}", True),
    ("dummy_rename_exchange_sign", r"K_{a b} V^{a} W^{b}+W^{c} K_{c d} V^{d}", True),
]
with open(sys.argv[1]) as f:
    rust = json.load(f)
assert rust["status"] == "PASS", rust
rust_rows = {r["id"]: r for r in rust["checks"]}
assert len(rust_rows) == len(cases) == len(rust["checks"])
rows = []
for name, source, expected in cases:
    ex = Ex(source)
    for _ in range(3):
        c.distribute(ex)
        c.eliminate_metric(ex)
        # Cadabra retains a mixed metric as g_a^b. Explicitly identify this
        # inverse-metric composition with the identity before taking traces.
        c.substitute(ex, Ex(r"g_{a}^{b} -> \delta^{b}_{a}"))
        c.substitute(ex, Ex(r"g^{a}_{b} -> \delta^{a}_{b}"))
        c.eliminate_kronecker(ex)
        c.sort_product(ex)
        c.canonicalise(ex)
        c.rename_dummies(ex)
        c.collect_terms(ex)
    zero = ex == 0
    assert zero == expected, (name, str(ex), expected)
    assert rust_rows[name]["residual_zero"] == zero, (name, str(ex), rust_rows[name])
    assert rust_rows[name]["expected_zero"] == expected
    rows.append({"id": name, "residual": str(ex), "zero": zero})
print(json.dumps({"status": "PASS", "backend": "cadabra2", "cases": rows}, indent=2))
