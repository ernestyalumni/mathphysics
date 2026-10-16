"""Independent live Cadabra checks. Run in the existing image, network disabled.
No reimplementation of symbolic transformations in Python: Cadabra performs them.
The Clifford form h equals -eta for Srednicki.
"""
import json
import cadabra2 as c
from cadabra2 import Ex

def spinor_residual(sign):
    expr = Ex(
        r"(a d-b c)(e h-f g)"
        + ("+" if sign == 1 else "-") + r"(a f-b e)(g d-h c)"
        + r"+(a h-b g)(c f-d e)"
    )
    c.distribute(expr)
    c.sort_product(expr)
    c.collect_terms(expr)
    return expr

c.__cdbkernel__ = c.create_scope()
good = spinor_residual(1)
bad = spinor_residual(-1)
assert good == 0, str(good)
assert bad != 0, "negative control unexpectedly passed"
c.Indices(Ex(r"{a,b}"), Ex(r"position=fixed"))
c.Symmetric(Ex(r"h_{a b}"))
c.GammaMatrix(Ex(r"\Gamma_{#}"), Ex(r"metric=h"))
anti = Ex(r"\Gamma_{a}\Gamma_{b}+\Gamma_{b}\Gamma_{a}-2 h_{a b}")
c.join_gamma(anti)
c.distribute(anti)
c.canonicalise(anti)
c.collect_terms(anti)
assert anti == 0, str(anti)
print(json.dumps({
    "backend": "cadabra2",
    "status": "PASS",
    "schouten_residual": str(good),
    "wrong_sign_residual": str(bad),
    "clifford_anticommutator_residual": str(anti),
    "module": c.__file__,
}))
