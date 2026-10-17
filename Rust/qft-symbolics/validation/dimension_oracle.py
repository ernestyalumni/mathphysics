"""Formal-d contracted gamma identities, using real Cadabra transformations."""
import json
import cadabra2 as c
from cadabra2 import Ex

c.__cdbkernel__ = c.create_scope()
c.Indices(Ex(r"{a,b,c,m,n}"), Ex(r"L, position=fixed"))
c.Integer(Ex(r"{a,b,c,m,n}"), Ex("1..d"))
c.Metric(Ex(r"h_{a b}"))
c.InverseMetric(Ex(r"h^{a b}"))
c.KroneckerDelta(Ex(r"\delta^{a}_{b}"))
c.GammaMatrix(Ex(r"\Gamma_{#}"), Ex("metric=h"))
sources = [
    ("47.20", r"-\Gamma^{m}\Gamma^{a}\Gamma^{b}\Gamma_{m}+4 h^{a b}+(d-4)\Gamma^{a}\Gamma^{b}", True),
    ("47.21", r"-\Gamma^{m}\Gamma^{a}\Gamma^{b}\Gamma^{c}\Gamma_{m}-2\Gamma^{c}\Gamma^{b}\Gamma^{a}-(d-4)\Gamma^{a}\Gamma^{b}\Gamma^{c}", True),
    ("47.20_wrong_sign", r"-\Gamma^{m}\Gamma^{a}\Gamma^{b}\Gamma_{m}-4 h^{a b}+(d-4)\Gamma^{a}\Gamma^{b}", False),
]
rows = []
for name, source, expected_zero in sources:
    expr = Ex(source)
    for _ in range(12):
        c.distribute(expr)
        c.sort_product(expr)
        c.join_gamma(expr)
        c.distribute(expr)
        c.eliminate_metric(expr)
        c.substitute(expr, Ex(r"h_{a}^{b} -> \delta^{b}_{a}"))
        c.substitute(expr, Ex(r"h^{a}_{b} -> \delta^{a}_{b}"))
        c.eliminate_kronecker(expr)
        c.sort_product(expr)
        c.canonicalise(expr)
        c.rename_dummies(expr)
        c.collect_terms(expr)
    assert (expr == 0) == expected_zero, (name, str(expr))
    rows.append({"identity": name, "residual": str(expr), "zero": expr == 0})
print(json.dumps({"status": "PASS", "convention": "Cadabra Clifford h=-eta; minus sign converts h-lowered sandwich to eta-lowered sandwich", "checks": rows}))
