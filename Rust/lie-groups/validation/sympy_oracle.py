"""Independent SymPy from the user's checkout, not a vendored/reimplemented CAS."""
import json
import sys
from pathlib import Path
sys.path.insert(0, "/sympy")  # read-only local checkout mounted by run.sh
import sympy as s
from sympy.physics.matrices import msigma
from sympy.algebras.quaternion import Quaternion

checks = 0
def require(ok, message):
    global checks
    if not ok:
        raise AssertionError(message)
    checks += 1
symbols = {name: s.Symbol(name, real=True) for name in
           ("v0", "v1", "v2", "t0", "t1", "t2", "s", "c", "theta")}
symbols["I"] = s.I
def expr(x):
    return s.sympify(x, locals=symbols)  # ONLY trusted, fixed-symbol producer output
def mat(rows):
    return s.Matrix([[expr(x) for x in row] for row in rows])
def equal(a, b, message):
    require(all(s.simplify(x) == 0 for x in a-b), message)

E = [s.Matrix(3, 3, lambda a, b: -s.LeviCivita(i, a, b)) for i in range(3)]
sigma = [msigma(i+1) for i in range(3)]
T = [-s.I*p/2 for p in sigma]
F = [-s.I*p for p in sigma]
theta = symbols["theta"]
circle = {symbols["s"]: s.sin(theta), symbols["c"]: s.cos(theta)}
for path in sys.argv[1:]:
    data = json.loads(Path(path).read_text())
    require(data["status"] == "COMPUTED", path)
    require(data["epsilon3"] == [int(s.LeviCivita(i,j,k)) for i in range(3) for j in range(3) for k in range(3)], "epsilon")
    for name, expected in (("E",E),("sigma",sigma),("T",T),("F",F)):
        require(len(data[name])==3, name+" count")
        for actual, wanted in zip(data[name], expected):
            equal(mat(actual), wanted, name)
    for i in range(3):
        for j in range(3):
            target = s.zeros(3)
            for k in range(3):
                target += s.LeviCivita(i,j,k)*E[k]
            equal(mat(data["E"][i])*mat(data["E"][j])-mat(data["E"][j])*mat(data["E"][i]), target, "commutator")
    v = s.Matrix([symbols[f"v{i}"] for i in range(3)])
    w = s.Matrix(s.symbols("w0:3", real=True))
    equal(mat(data["hat"])*w, v.cross(w), "cross product")
    equal(mat(data["r_exp_z"]).subs(circle), (theta*E[2]).exp(), "SO3 matrix exponential")
    for key, basis, scale in (("u_exp_full_z",F,"1"),("u_exp_half_z",T,"1/2")):
        actual = mat(data[key]).subs({
            symbols["s"]: s.sin(expr(scale)*theta),
            symbols["c"]: s.cos(expr(scale)*theta)})
        equal(actual, (theta*basis[2]).exp(), key)
    require(expr(data["full_argument_scale"]) == 1 and expr(data["half_argument_scale"]) == s.Rational(1,2), "argument normalization")
    t = [symbols[f"t{i}"] for i in range(3)]
    squared = sum(x*x for x in t)
    q = [(1-squared)/(1+squared)] + [2*x/(1+squared) for x in t]
    expected_u = q[0]*s.eye(2) - s.I*sum((q[i+1]*sigma[i] for i in range(3)), s.zeros(2))
    u = mat(data["chart_u"]["numerator"])/expr(data["chart_u"]["denominator"])
    r = mat(data["chart_r"]["numerator"])/expr(data["chart_r"]["denominator"])
    equal(u, expected_u, "arbitrary symbolic SU2 chart")
    equal(u.H*u, s.eye(2), "SU2 unitarity")
    require(s.cancel(u.det()-1) == 0, "SU2 determinant")
    equal(r, Quaternion(*q).to_rotation_matrix(), "SymPy quaternion SO3 conversion")
    equal(r.T*r, s.eye(3), "SO3 orthogonality")
    require(s.cancel(r.det()-1) == 0, "SO3 determinant")
    require(E[0]*E[1]-E[1]*E[0] != -E[2], "wrong-sign control")
    require(msigma(3).det() == -1, "Pauli is not SU2 control")
print(json.dumps({"status":"PASS", "backend":"sympy", "version":s.__version__,
                  "source":s.__file__, "producer_records":len(sys.argv)-1, "checks":checks}))
