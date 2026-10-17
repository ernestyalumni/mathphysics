"""Real SageMath exact checks. Do not rely solely on a named group's validator."""
import json
import sys
from pathlib import Path
from sage.all import (QQ, QQbar, SR, I, var, matrix, identity_matrix, diagonal_matrix,
    sage_eval, SO, SU, LieAlgebra, FiniteRankFreeModule, EuclideanSpace)
from sage.version import version

checks = 0
def require(ok, message):
    global checks
    if not bool(ok):
        raise AssertionError(message)
    checks += 1
names = {n: var(n, domain="real") for n in
         ("v0","v1","v2","t0","t1","t2","s","c","theta")}
names["I"] = I
def expr(x):
    return sage_eval(x, locals=names)  # ONLY trusted producer expressions, not model input
def mat(rows):
    return matrix(SR, [[expr(x) for x in row] for row in rows])
def equal(a,b,message):
    require(all(x.simplify_full() == 0 for x in (a-b).list()), message)
# Actual alternating form, not just a Python permutation implementation.
module = FiniteRankFreeModule(QQ, 3, name="V")
module.basis("e")
epsilon = module.alternating_form(3, name="epsilon")
epsilon[0,1,2] = 1
eps_values = [int(epsilon[i,j,k]) for i in range(3) for j in range(3) for k in range(3)]
# The metric volume form agrees in a positively oriented Cartesian orthonormal frame.
space = EuclideanSpace(3, start_index=0)
volume = space.metric().volume_form()
require(volume[0,1,2] == 1, "oriented Euclidean volume form")
E = [matrix(QQ, [[0,0,0],[0,0,-1],[0,1,0]]),
     matrix(QQ, [[0,0,1],[0,0,0],[-1,0,0]]),
     matrix(QQ, [[0,-1,0],[1,0,0],[0,0,0]])]
sigma = [matrix(SR,[[0,1],[1,0]]),matrix(SR,[[0,-I],[I,0]]),matrix(SR,[[1,0],[0,-1]])]
T = [-I*p/2 for p in sigma]
F = [-I*p for p in sigma]
# This explicit compact rational model is NOT silently Sage's split Cartan A1.
L = LieAlgebra(QQ, {("e1","e2"):{"e3":1},("e2","e3"):{"e1":1},
                   ("e3","e1"):{"e2":1}}, names=("e1","e2","e3"))
x,y,z = L.gens()
require(x.bracket(y)==z and y.bracket(z)==x and z.bracket(x)==y, "abstract compact structure constants")
require(L.killing_form_matrix()==-2*identity_matrix(QQ,3), "negative-definite compact Killing form")
reflection=diagonal_matrix(QQ,[-1,1,1])
reflection_accepted = False
try:
    SO(3,QQ)(reflection)
    reflection_accepted = True
except (TypeError,ValueError):
    pass
require(reflection.det()==-1, "reflection must not be our SO3 element")

for path in sys.argv[1:]:
    data=json.loads(Path(path).read_text())
    require(data["status"]=="COMPUTED",path)
    require(data["epsilon3"]==eps_values,"epsilon form components")
    for key, expected in (("E",E),("sigma",sigma),("T",T),("F",F)):
        require(len(data[key])==3,key+" count")
        for a,b in zip(data[key],expected): equal(mat(a),b,key)
    for i in range(3):
        for j in range(3):
            target=sum((epsilon[i,j,k]*E[k] for k in range(3)),matrix(QQ,3,3))
            equal(mat(data["E"][i])*mat(data["E"][j])-mat(data["E"][j])*mat(data["E"][i]),target,"bracket")
    v=matrix(SR,3,1,[names[f"v{i}"] for i in range(3)])
    expected_hat=names["v0"]*E[0]+names["v1"]*E[1]+names["v2"]*E[2]
    equal(mat(data["hat"]),expected_hat,"abstract hat coordinates")
    s,c=names["s"],names["c"]
    equal(mat(data["r_exp_z"]),matrix(SR,[[c,-s,0],[s,c,0],[0,0,1]]),"Rodrigues expression")
    equal(mat(data["u_exp_full_z"]),matrix(SR,[[c-I*s,0],[0,c+I*s]]),"Frankel diagonal exponential")
    require(expr(data["full_argument_scale"])==1 and expr(data["half_argument_scale"])==QQ(1)/2,"half-angle metadata")
    t=[names[f"t{i}"] for i in range(3)]
    squared=sum(a*a for a in t)
    q=[(1-squared)/(1+squared)]+[2*a/(1+squared) for a in t]
    expected=q[0]*identity_matrix(SR,2)-I*sum((q[i+1]*sigma[i] for i in range(3)),matrix(SR,2,2))
    u=mat(data["chart_u"]["numerator"])/expr(data["chart_u"]["denominator"])
    r=mat(data["chart_r"]["numerator"])/expr(data["chart_r"]["denominator"])
    equal(u,expected,"symbolic SU2 chart")
    equal(u.conjugate_transpose()*u,identity_matrix(SR,2),"exact unitarity")
    require((u.det()-1).simplify_full()==0,"unit determinant")
    equal(r.transpose()*r,identity_matrix(SR,3),"exact orthogonality")
    require((r.det()-1).simplify_full()==0,"proper rotation")
    # Independent conjugation representation, including sign and factor of two.
    for j in range(3):
        column=sum((r[i,j]*sigma[i] for i in range(3)),matrix(SR,2,2))
        equal(u*sigma[j]*u.conjugate_transpose(),column,"Pauli conjugation covering")
    substitution={t[0]:0,t[1]:0,t[2]:QQ(1)/2}
    u0=matrix(QQbar, u.subs(substitution));r0=matrix(QQ,r.subs(substitution))
    SU(2,QQbar)(u0)
    SO(3,QQ)(r0)
    require(u0.det()==1 and r0.det()==1,"native group positive witness with explicit determinant gate")
    rejected=False
    try: SU(2,QQbar)(matrix(QQbar,sigma[2]))
    except (TypeError,ValueError): rejected=True
    require(rejected,"native SU2 Pauli negative witness")
print(json.dumps({"status":"PASS","backend":"sage","version":version,
    "producer_records":len(sys.argv)-1,"checks":checks,
    "diagnostic_SO_QQ_reflection_accepted":reflection_accepted,
    "warning":"Native SO label was not used as determinant-one ground truth; determinants independently checked."}))
