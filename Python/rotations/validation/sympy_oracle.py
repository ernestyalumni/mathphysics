"""Exact explicit-matrix oracle, independent of numerical implementation."""
import json
import sys
import sympy as s
from sympy.algebras.quaternion import Quaternion

rows=[]
def check(name,residual):
    values=list(residual) if isinstance(residual,s.MatrixBase) else [residual]
    values=[s.trigsimp(s.expand(x)) for x in values]
    if any(x!=0 for x in values): raise AssertionError((name,values))
    rows.append(name)
    print(name, file=sys.stderr, flush=True)

def rx(t): return s.Matrix([[1,0,0],[0,s.cos(t),-s.sin(t)],[0,s.sin(t),s.cos(t)]])
def ry(t): return s.Matrix([[s.cos(t),0,s.sin(t)],[0,1,0],[-s.sin(t),0,s.cos(t)]])
def rz(t): return s.Matrix([[s.cos(t),-s.sin(t),0],[s.sin(t),s.cos(t),0],[0,0,1]])
def hat(v): return s.Matrix([[0,-v[2],v[1]],[v[2],0,-v[0]],[-v[1],v[0],0]])

# Exact three-half-angle polynomial quotient, rather than heuristic trig
# simplification. Constraints are independently imposed for all three angles.
cs=s.symbols('c0 s0 c1 s1 c2 s2',real=True)
c0,s0,c1,s1,c2,s2=cs
basis=s.groebner([c0*c0+s0*s0-1,c1*c1+s1*s1-1,c2*c2+s2*s2-1],*cs,domain=s.QQ)
def quotient_check(name,residual):
    values=list(residual) if isinstance(residual,s.MatrixBase) else [residual]
    if any(basis.reduce(s.expand(value))[1]!=0 for value in values):
        raise AssertionError(name)
    rows.append(name)
    print(name,file=sys.stderr,flush=True)
def erx(c,t): return s.Matrix([[1,0,0],[0,c,-t],[0,t,c]])
def ery(c,t): return s.Matrix([[c,0,t],[0,1,0],[-t,0,c]])
def erz(c,t): return s.Matrix([[c,-t,0],[t,c,0],[0,0,1]])
full=[(c0*c0-s0*s0,2*c0*s0),(c1*c1-s1*s1,2*c1*s1),(c2*c2-s2*s2,2*c2*s2)]
qzs=[Quaternion(c0,0,0,s0),Quaternion(c1,0,0,s1),Quaternion(c2,0,0,s2)]
products={}
for name,middle,last,qmiddle,qlast in [
 ('zxz',erx,erz,Quaternion(c1,s1,0,0),Quaternion(c2,0,0,s2)),
 ('zyz',ery,erz,Quaternion(c1,0,s1,0),Quaternion(c2,0,0,s2)),
 ('zyx',ery,erx,Quaternion(c1,0,s1,0),Quaternion(c2,s2,0,0))]:
    product=erz(*full[0])*middle(*full[1])*last(*full[2])
    products[name]=product
    q=qzs[0]*qmiddle*qlast
    w,x,y,z=q.a,q.b,q.c,q.d
    qm=s.Matrix([[w*w+x*x-y*y-z*z,2*(x*y-w*z),2*(x*z+w*y)],[2*(x*y+w*z),w*w-x*x+y*y-z*z,2*(y*z-w*x)],[2*(x*z-w*y),2*(y*z+w*x),w*w-x*x-y*y+z*z]])
    quotient_check(name+'_quaternion_norm',w*w+x*x+y*y+z*z-1)
    quotient_check(name+'_quaternion_product',product-qm)
    quotient_check(name+'_orthogonality',product.T*product-s.eye(3))
    quotient_check(name+'_determinant',product.det()-1)
# The x/y dictionary changes full first/last angles by -/+ pi/2.
quotient_check('zxz_zyz_shift',products['zxz']-erz(full[0][1],-full[0][0])*ery(*full[1])*erz(-full[2][1],full[2][0]))
quotient_check('Goldstein_half_angle_ZXZ',(c1*(c0*c2-s0*s2))**2-(1+s.trace(products['zxz']))/4)
phi,theta,psi=s.symbols('phi theta psi',real=True)
# Generic unit quaternion: all unit quaternions, with no rational-chart
# denominator blow-up. Real symbols and unit-norm ideal are explicit.
w,x,y,z=s.symbols('w x y z',real=True)
v=s.Matrix(s.symbols('v1:4',real=True))
qvars=(w,x,y,z,*v)
qbasis=s.groebner([w*w+x*x+y*y+z*z-1],*qvars,domain=s.QQ)
def ring_check(name,residual,ideal):
    values=list(residual) if isinstance(residual,s.MatrixBase) else [residual]
    for value in values:
        numerator,denominator=s.fraction(s.cancel(value))
        # All rational denominators in these checks are either native q.norm^2
        # (equal to 1 under the ideal) or h*(1+h^2), with real h != 0.
        for part in (s.re(numerator),s.im(numerator)):
            if ideal.reduce(s.expand(part))[1]!=0: raise AssertionError((name,part))
    rows.append(name)
    print(name,file=sys.stderr,flush=True)
def qcheck(name,residual): ring_check(name,residual,qbasis)
q=Quaternion(w,x,y,z)
R=s.Matrix([[w*w+x*x-y*y-z*z,2*(x*y-w*z),2*(x*z+w*y)],
            [2*(x*y+w*z),w*w-x*x+y*y-z*z,2*(y*z-w*x)],
            [2*(x*z-w*y),2*(y*z+w*x),w*w-x*x-y*y+z*z]])
qcheck('native_SymPy_quaternion',R-q.to_rotation_matrix())
qcheck('quaternion_orthogonality',R.T*R-s.eye(3))
qcheck('quaternion_determinant',R.det()-1)
qcheck('hat_conjugation',R*hat(v)*R.T-hat(R*v))
alpha,beta=w+s.I*z,y+s.I*x;gamma,delta=-s.conjugate(beta),s.conjugate(alpha)
K=s.Matrix([[alpha,beta],[gamma,delta]])
U=s.Matrix([[w-s.I*z,-y-s.I*x],[y-s.I*x,w+s.I*z]])
qcheck('CK_U_dagger',K-U.conjugate().T)
qcheck('CK_unitarity',K.conjugate().T*K-s.eye(2))
# Goldstein's printed polynomial in the four CK entries.
AG=s.Matrix([[(alpha**2-gamma**2+delta**2-beta**2)/2,s.I*(gamma**2-alpha**2+delta**2-beta**2)/2,gamma*delta-alpha*beta],
 [s.I*(alpha**2+gamma**2-beta**2-delta**2)/2,(alpha**2+gamma**2+beta**2+delta**2)/2,-s.I*(alpha*beta+gamma*delta)],
 [beta*delta-alpha*gamma,s.I*(alpha*gamma+beta*delta),alpha*delta+beta*gamma]])
qcheck('CK_Goldstein_equals_R_transpose',AG-R.T)
t=s.trace(R)
qcheck('SO3_Cayley_Hamilton',R**3-t*R**2+t*R-s.eye(3))
l=s.symbols('lambda',real=True)
# CharacteristicPolynomial canonicalizes generator assumptions; substitute
# its actual generator before collecting the four coefficients.
characteristic=R.charpoly()
spectral_residual=characteristic.as_expr().subs(characteristic.gen,l)-(l-1)*(l*l-(t-1)*l+1)
qcheck('SO3_characteristic_spectrum',spectral_residual.expand().coeff(l,0))
for degree in range(1,4):
    qcheck('SO3_spectrum_coefficient_'+str(degree),spectral_residual.expand().coeff(l,degree))
# Axis-angle relations in a second explicit real polynomial quotient.
n=s.Matrix(s.symbols('n1:4',real=True));N=hat(n)
qw,qs=s.symbols('qw qs',real=True)
h=s.symbols('h',real=True,nonzero=True)
translation=s.Matrix(s.symbols('t1:4',real=True))
point=s.Matrix(s.symbols('p1:4',real=True))
axisbasis=s.groebner([n.dot(n)-1,qw*qw+qs*qs-1],*n,qw,qs,h,*translation,*point,domain=s.QQ)
def acheck(name,residual): ring_check(name,residual,axisbasis)
acheck('hat_cube',N**3+N)
sinphi,cosphi=2*qw*qs,qw*qw-qs*qs
Rod=s.eye(3)+sinphi*N+(1-cosphi)*N**2
Q=Quaternion(qw,*(qs*ni for ni in n))
acheck('Rodrigues_quaternion',Rod-Q.to_rotation_matrix())
acheck('axis_fixed',Rod*n-n)
acheck('axis_trace',s.trace(Rod)-1-2*cosphi)
acheck('Goldstein_clockwise_sign',Rod.T-(s.eye(3)-sinphi*N+(1-cosphi)*N**2))
# Strong Chasles, h=tan(phi/2) != 0; phi=pi follows by the cot=0
# limit, covered numerically. Pure translation handled separately in notes.
Rh=s.eye(3)+2*h/(1+h*h)*N+2*h*h/(1+h*h)*N**2
distance=n.dot(translation);perp=translation-distance*n
center=(perp+N*perp/h)/2
acheck('Chasles_axis_point',(s.eye(3)-Rh)*center-perp)
acheck('Chasles_screw_action',translation+Rh*point-(center+Rh*(point-center)+distance*n))
# General matrix CH, beyond the SO3-specialized polynomial.
m=s.Matrix(3,3,s.symbols('m0:9'));poly=m.charpoly(l)
check('general_3x3_Cayley_Hamilton',sum((coefficient*m**(3-i) for i,coefficient in enumerate(poly.all_coeffs())),s.zeros(3)))
wrong=rz(s.pi/2)*hat(s.Matrix([1,0,0]))*rz(s.pi/2)-hat(rz(s.pi/2)*s.Matrix([1,0,0]))
if wrong==s.zeros(3): raise AssertionError('negative control failed')
# Reversing a quarter-turn or omitting CK conjugacy minus must fail.
if rz(s.pi/2)==rz(-s.pi/2): raise AssertionError('sign control failed')
wrongK=s.Matrix([[s.sqrt(2)/2,s.sqrt(2)/2],[s.sqrt(2)/2,s.sqrt(2)/2]])
if wrongK.conjugate().T*wrongK==s.eye(2): raise AssertionError('CK control failed')
print(json.dumps({'status':'PASS','backend':'sympy','version':s.__version__,'checks':rows,'negative_controls':['missing_conjugation_inverse','clockwise_sign','missing_CK_minus']}))
# Optional Rust-owned symbolic-ring exports: producer is trusted fixed syntax,
# never untrusted expressions. Domains and arities are checked explicitly.
if len(sys.argv)>1:
    import pathlib
    exports=[json.loads(line) for line in pathlib.Path(sys.argv[1]).read_text().splitlines() if line.strip()]
    expected_keys={(name,argument) for name in ('zxz','zyz','zyx') for argument in ('full','half')}
    if {(row['convention'],row['argument']) for row in exports}!=expected_keys:
        raise AssertionError('missing symbolic conventions')
    local={str(symbol):symbol for symbol in cs};local['I']=s.I
    def parse(value):
        expression=s.sympify(value,locals=local)
        if expression.free_symbols-set(cs): raise AssertionError('unexpected symbols')
        return expression
    for row in exports:
        name=row['convention'];argument=row['argument']
        triples=full if argument=='half' else [(c0,s0),(c1,s1),(c2,s2)]
        rotations={'x':erx,'y':ery,'z':erz}
        expected=rotations[name[0]](*triples[0])*rotations[name[1]](*triples[1])*rotations[name[2]](*triples[2])
        got=s.Matrix([[parse(value) for value in line] for line in row['matrix']])
        quotient_check('rust_'+row['name']+'_matrix',got-expected)
        if 'quaternion' in row:
            factors=[]
            for axis,(co,si) in zip(name,[(c0,s0),(c1,s1),(c2,s2)]):
                factors.append(Quaternion(co,si if axis=='x' else 0,si if axis=='y' else 0,si if axis=='z' else 0))
            q=factors[0]*factors[1]*factors[2]
            quotient_check('rust_'+row['name']+'_quaternion',s.Matrix([parse(v) for v in row['quaternion']])-s.Matrix([q.a,q.b,q.c,q.d]))
            ck=s.Matrix([[q.a+s.I*q.d,q.c+s.I*q.b],[-q.c+s.I*q.b,q.a-s.I*q.d]])
            residual=s.Matrix([[parse(value) for value in line] for line in row['cayley_klein']])-ck
            quotient_check('rust_'+row['name']+'_CK_real',s.re(residual))
            quotient_check('rust_'+row['name']+'_CK_imag',s.im(residual))
    print(json.dumps({'status':'PASS','backend':'sympy','scope':'Rust independent MultiTrig ring exports','checks':rows[-15:]}))
