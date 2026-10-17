"""Rust/C++ f64 comparison and independent matrix Taylor reference, no skipped cases."""
import json
import math
import sys
from pathlib import Path

def require(ok, message):
    if not ok:
        raise AssertionError(message)
def multiply(a,b):
    return [[sum(a[i][k]*b[k][j] for k in range(len(b))) for j in range(len(b[0]))] for i in range(len(a))]
def exponential(a):
    n=len(a)
    term=[[float(i==j) for j in range(n)] for i in range(n)]
    result=[row[:] for row in term]
    # All fixtures have small generator norms; 80 terms makes truncation negligible
    # compared with f64 arithmetic. This algorithm does not use sin/cos/Rodrigues.
    for k in range(1,81):
        term=[[x/k for x in row] for row in multiply(term,a)]
        result=[[result[i][j]+term[i][j] for j in range(n)] for i in range(n)]
    return result
def qmatrix(q):
    a,b,c,d=q
    return [[a-1j*d,-c-1j*b],[c-1j*b,a+1j*d]]
def maximum(a,b):
    return max(abs(x-y) for ra,rb in zip(a,b) for x,y in zip(ra,rb))
rust=json.loads(Path(sys.argv[1]).read_text())
cpp=json.loads(Path(sys.argv[2]).read_text())
require(rust["status"]==cpp["status"]=="COMPUTED", "producer failure")
require(len(rust["cases"])==len(cpp["cases"])==25, "fixture count")
errors=[]
for i,(a,b) in enumerate(zip(rust["cases"],cpp["cases"])):
    require(a["id"]==b["id"]==i, "fixture identity")
    require(a["v"]==b["v"], "input parity")
    x,y,z=a["v"]
    e=[[0.,-z,y],[z,0.,-x],[-y,x,0.]]
    u=[[complex(0,-z/2),complex(-y/2,-x/2)],
       [complex(y/2,-x/2),complex(0,z/2)]]
    expected_r=exponential(e);expected_u=exponential(u)
    errors.append(maximum(a["rotation"],b["rotation"]))
    errors.extend(abs(x-y) for x,y in zip(a["quaternion"],b["quaternion"]))
    for actual in (a,b):
        require(len(actual["rotation"])==3 and all(len(row)==3 for row in actual["rotation"]), "rotation dimensions")
        require(len(actual["quaternion"])==4, "quaternion dimensions")
        errors.append(maximum(actual["rotation"],expected_r))
        errors.append(maximum(qmatrix(actual["quaternion"]),expected_u))
require(all(math.isfinite(x) and x<=1e-12 for x in errors), "parity/reference error")
print(json.dumps({"status":"PASS","cases":25,"skipped":0,"tolerance":1e-12,
                  "max_absolute_error":max(errors),"reference":"independent 80-term matrix exponential"}))
