"""Cadabra-style indexed Euclidean structure constants and epsilon contractions."""
import json
import cadabra2 as c
from cadabra2 import Ex
c.__cdbkernel__=c.create_scope()
c.Indices(Ex(r"{a,b,c,d,e,i,j,k,l,m}"),Ex("euclidean, position=free"))
c.Integer(Ex(r"{a,b,c,d,e,i,j,k,l,m}"),Ex("0..2"))
c.KroneckerDelta(Ex(r"\delta{#}"))
c.EpsilonTensor(Ex(r"\epsilon_{a b c}"),Ex(r"delta=\delta"))
c.SelfNonCommuting(Ex(r"E_{a}"))
# E_i E_j - E_j E_i = epsilon_ijk E_k is an algebra declaration, not
# something Cadabra can infer merely from the string 'E'. The following
# establishes the Jacobi consistency of those declared structure constants.
cases=[
 ("epsilon_contraction",r"\epsilon_{a b c}\epsilon_{a d e}-\delta_{b d}\delta_{c e}+\delta_{b e}\delta_{c d}",True),
 ("epsilon_trace",r"\epsilon_{a b c}\epsilon_{a b d}-2\delta_{c d}",True),
 ("structure_constant_jacobi",r"\epsilon_{i j m}\epsilon_{m k l}+\epsilon_{j k m}\epsilon_{m i l}+\epsilon_{k i m}\epsilon_{m j l}",True),
 ("wrong_trace_sign",r"\epsilon_{a b c}\epsilon_{a b d}+2\delta_{c d}",False),
]
rows=[]
for name,source,expected in cases:
    ex=Ex(source)
    c.epsilon_to_delta(ex)
    c.expand_delta(ex)
    for _ in range(3):
        c.distribute(ex)
        c.eliminate_kronecker(ex)
        c.sort_product(ex)
        c.canonicalise(ex)
        c.rename_dummies(ex)
        c.collect_terms(ex)
    if (ex==0)!=expected: raise AssertionError((name,str(ex)))
    rows.append({"id":name,"residual":str(ex),"zero":ex==0})
print(json.dumps({"status":"PASS","backend":"cadabra2","cases":rows,
    "scope":"3D Euclidean epsilon/structure constants; not an automatic group-element factory"}))
