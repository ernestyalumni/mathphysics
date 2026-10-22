"""RH rotations of column vectors; R maps body components to inertial components.

Euler inputs are always product order: ZXZ/ZYZ (phi, theta, psi), ZYX
(heading, pitch, roll). Goldstein's coordinate matrix is A_G = R.T.
Numerical imports reject reflections, nonfinite data and nonunit quaternions.
"""
from __future__ import annotations
import math
from dataclasses import dataclass

Matrix = tuple[tuple[float, ...], ...]
I3 = ((1., 0., 0.), (0., 1., 0.), (0., 0., 1.))

def _shape(a):
    if not a or not a[0] or any(len(row)!=len(a[0]) for row in a):
        raise ValueError('nonempty rectangular matrix required')
    _finite(x for row in a for x in row)
    return len(a),len(a[0])

def transpose(a):
    _shape(a)
    return tuple(zip(*a))

def multiply(a, b):
    _,cols=_shape(a)
    rows,_=_shape(b)
    if cols!=rows: raise ValueError('incompatible matrix dimensions')
    return tuple(tuple(sum(x*y for x, y in zip(row, col)) for col in zip(*b)) for row in a)

def apply(a, x):
    _,cols=_shape(a)
    if cols!=len(x): raise ValueError('incompatible matrix/vector dimensions')
    _finite(x)
    return tuple(sum(ai*xi for ai, xi in zip(row, x)) for row in a)

def det(a):
    if _shape(a)!=(3,3): raise ValueError('3x3 determinant required')
    return (a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
            -a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
            +a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]))

def _finite(values):
    if not all(math.isfinite(x) for x in values):
        raise ValueError('finite real values required')

def _positive_tolerance(tolerance):
    if not math.isfinite(tolerance) or tolerance <= 0 or tolerance > 1e-4:
        raise ValueError('tolerance must be finite and in (0, 1e-4]')

def validate_matrix(a, tolerance=1e-10):
    _positive_tolerance(tolerance)
    if len(a) != 3 or any(len(row) != 3 for row in a):
        raise ValueError('3x3 matrix required')
    _finite(x for row in a for x in row)
    gram = multiply(transpose(a), a)
    if max(abs(gram[i][j]-I3[i][j]) for i in range(3) for j in range(3)) > tolerance or abs(det(a)-1) > tolerance:
        raise ValueError('proper orthogonal matrix required')
    return tuple(tuple(float(x) for x in row) for row in a)

def hat(x):
    if len(x) != 3:
        raise ValueError('3-vector required')
    _finite(x)
    a,b,c=x
    return ((0.,-c,b),(c,0.,-a),(-b,a,0.))

def elementary(axis, phi):
    _finite((phi,))
    c,s=math.cos(phi),math.sin(phi)
    if axis == 'x': return ((1.,0.,0.),(0.,c,-s),(0.,s,c))
    if axis == 'y': return ((c,0.,s),(0.,1.,0.),(-s,0.,c))
    if axis == 'z': return ((c,-s,0.),(s,c,0.),(0.,0.,1.))
    raise ValueError('axis must be x, y or z')

def from_euler(convention, first, middle, last):
    """Intrinsic moving-axis product R_axis1(first) R_axis2(middle) R_axis3(last)."""
    if convention not in ('zxz','zyz','zyx'):
        raise ValueError('supported conventions: zxz, zyz, zyx')
    return multiply(multiply(elementary(convention[0],first),elementary(convention[1],middle)),elementary(convention[2],last))

@dataclass(frozen=True)
class EulerAngles:
    first: float
    middle: float
    last: float
    singular: bool

    @property
    def values(self): return (self.first,self.middle,self.last)

def to_euler(a, convention, singular_tolerance=1e-12):
    """Canonical middle range [0,pi] (proper Euler), [-pi/2,pi/2] (ZYX).

    At lock choose last=0; the singular flag reports this choice. Within the
    chosen singular tolerance reconstruction differs by O(tolerance).
    """
    a=validate_matrix(a)
    _positive_tolerance(singular_tolerance)
    if convention in ('zxz','zyz'):
        st=math.hypot(a[0][2],a[1][2])
        theta=math.atan2(st,a[2][2])
        if st > singular_tolerance:
            if convention=='zxz':
                return EulerAngles(math.atan2(a[0][2],-a[1][2]),theta,math.atan2(a[2][0],a[2][1]),False)
            return EulerAngles(math.atan2(a[1][2],a[0][2]),theta,math.atan2(a[2][1],-a[2][0]),False)
        phi=math.atan2(a[1][0],a[0][0])
        if convention=='zyz' and a[2][2]<0:
            phi=math.atan2(-a[1][0],-a[0][0])
        return EulerAngles(phi,theta,0.,True)
    if convention=='zyx':
        ct=math.hypot(a[0][0],a[1][0])
        theta=math.atan2(-a[2][0],ct)
        if ct > singular_tolerance:
            return EulerAngles(math.atan2(a[1][0],a[0][0]),theta,math.atan2(a[2][1],a[2][2]),False)
        return EulerAngles(math.atan2(-a[0][1],a[1][1]),theta,0.,True)
    raise ValueError('supported conventions: zxz, zyz, zyx')

def validate_quaternion(q, tolerance=1e-10):
    _positive_tolerance(tolerance)
    if len(q)!=4: raise ValueError('scalar-first Hamilton quaternion required')
    _finite(q)
    if abs(sum(x*x for x in q)-1)>tolerance: raise ValueError('unit quaternion required')
    return tuple(float(x) for x in q)

def quaternion_product(q,p):
    w,x,y,z=validate_quaternion(q)
    v,a,b,c=validate_quaternion(p)
    return validate_quaternion((w*v-x*a-y*b-z*c,w*a+x*v+y*c-z*b,w*b-x*c+y*v+z*a,w*c+x*b-y*a+z*v))

def from_quaternion(q):
    w,x,y,z=validate_quaternion(q)
    return validate_matrix(((1-2*(y*y+z*z),2*(x*y-w*z),2*(x*z+w*y)),
            (2*(x*y+w*z),1-2*(x*x+z*z),2*(y*z-w*x)),
            (2*(x*z-w*y),2*(y*z+w*x),1-2*(x*x+y*y))))

def to_quaternion(a):
    """Largest-diagonal extraction, canonical w>=0; q/-q represent the same R."""
    a=validate_matrix(a)
    tr=sum(a[i][i] for i in range(3))
    if tr>0:
        s=2*math.sqrt(max(0.,tr+1))
        q=(s/4,(a[2][1]-a[1][2])/s,(a[0][2]-a[2][0])/s,(a[1][0]-a[0][1])/s)
    else:
        i=max(range(3),key=lambda k:a[k][k]);j=(i+1)%3;k=(i+2)%3
        s=2*math.sqrt(max(0.,1+a[i][i]-a[j][j]-a[k][k]))
        v=[0.,0.,0.];v[i]=s/4;v[j]=(a[j][i]+a[i][j])/s;v[k]=(a[k][i]+a[i][k])/s
        q=((a[k][j]-a[j][k])/s,*v)
    norm=math.hypot(*q)
    q=tuple(x/norm for x in q)
    return tuple(-x for x in q) if q[0]<0 else q

def from_axis_angle(axis, phi):
    if len(axis)!=3: raise ValueError('3-vector required')
    _finite((*axis,phi))
    n=math.hypot(*axis)
    if abs(n-1)>1e-10: raise ValueError('unit axis required')
    # 2 sin(phi/2)^2 avoids cancellation in 1-cos(phi) near zero.
    h=hat(axis);h2=multiply(h,h);s=math.sin(phi);d=2*math.sin(phi/2)**2
    return validate_matrix(tuple(tuple(I3[i][j]+s*h[i][j]+d*h2[i][j] for j in range(3)) for i in range(3)))

def to_axis_angle(a):
    """Canonical phi in [0,pi]; identity uses arbitrary axis (1,0,0).

    Quaternion extraction avoids division by sin(phi) near pi.
    """
    q=to_quaternion(a);s=math.hypot(*q[1:])
    if s==0: return ((1.,0.,0.),0.)
    return (tuple(x/s for x in q[1:]),2*math.atan2(s,q[0]))

def quaternion_to_ck(q):
    """Goldstein K=U(q)^dagger, alpha=w+iz, beta=y+ix."""
    w,x,y,z=validate_quaternion(q)
    alpha,beta=complex(w,z),complex(y,x)
    return (alpha,beta,-beta.conjugate(),alpha.conjugate())

def ck_to_quaternion(parameters, tolerance=1e-10):
    _positive_tolerance(tolerance)
    if len(parameters)!=4: raise ValueError('four complex parameters required')
    alpha,beta,gamma,delta=parameters
    _finite(y for x in parameters for y in (x.real,x.imag))
    if abs(gamma+beta.conjugate())>tolerance or abs(delta-alpha.conjugate())>tolerance:
        raise ValueError('CK conjugacy constraints required')
    return validate_quaternion((alpha.real,beta.imag,beta.real,alpha.imag),tolerance)

def ck_to_goldstein(parameters):
    return transpose(from_quaternion(ck_to_quaternion(parameters)))

def cayley_hamilton_residual(a):
    """Actual CH polynomial R^3-tR^2+tR-I for R in SO(3), t=tr(R)."""
    a=validate_matrix(a);t=sum(a[i][i] for i in range(3));a2=multiply(a,a);a3=multiply(a2,a)
    return tuple(tuple(a3[i][j]-t*a2[i][j]+t*a[i][j]-I3[i][j] for j in range(3)) for i in range(3))
