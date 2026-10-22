"""Indexed identities supporting Rodrigues; no group membership inferred from names."""
import json
import sys
import cadabra2 as c
from cadabra2 import Ex
c.__cdbkernel__=c.create_scope()
c.Indices(Ex(r'{a,b,c,d,e,i,j,k,l,m}'),Ex('euclidean, position=free'))
c.Integer(Ex(r'{a,b,c,d,e,i,j,k,l,m}'),Ex('0..2'))
c.KroneckerDelta(Ex(r'\delta{#}'))
c.EpsilonTensor(Ex(r'\epsilon_{i j k}'),Ex(r'delta=\delta'))
cases=[
 ('cross_cross',r'\epsilon_{i j k}\epsilon_{k l m}n_j n_l x_m-n_i n_j x_j+n_j n_j x_i',True),
 ('hat_cube',r'\epsilon_{i a j}\epsilon_{j b k}\epsilon_{k c l}n_a n_b n_c+n_a n_a\epsilon_{i b l}n_b',True),
 ('axis_fixed',r'\epsilon_{i j k}n_j n_k',True),
 ('wrong_triple_sign',r'\epsilon_{i j k}\epsilon_{k l m}n_j n_l x_m+n_i n_j x_j-n_j n_j x_i',False),
]
rows=[]
for name,source,expected in cases:
    ex=Ex(source)
    c.epsilon_to_delta(ex)
    c.expand_delta(ex)
    for _ in range(4):
        c.distribute(ex)
        c.eliminate_kronecker(ex)
        c.sort_product(ex)
        c.canonicalise(ex)
        c.rename_dummies(ex)
        c.collect_terms(ex)
    if (ex==0)!=expected: raise AssertionError((name,str(ex)))
    rows.append({'id':name,'zero':ex==0,'residual':str(ex)})
print(json.dumps({'status':'PASS','backend':'cadabra2','python':sys.version,'module':c.__file__,'scope':'formal 3D Euclidean epsilon contractions, arbitrary commuting n and x','cases':rows}))
