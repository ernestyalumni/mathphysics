"""Ideal COM-based 6-DOF rigid body, with body-to-inertial Hamilton attitude.

State is (r_I, v_I, q, omega_B). Wrenches are constant over one step;
body force is evaluated at the half-step attitude. Units must be consistent.
The implicit midpoint momentum/Cayley update is second order, preserves free
rotational energy and the complete inertial spin vector to solve roundoff,
and has finite time-step phase error. No aerodynamic/contact/gravity model
is implicitly supplied: gravity can be an explicit inertial force m*g.
"""
from __future__ import annotations
from dataclasses import dataclass, field
import math
import numpy as np
import rotations

ZERO=(0.,0.,0.)
IDENTITY_Q=(1.,0.,0.,0.)

class ConvergenceError(RuntimeError):
    """The nonlinear step failed its residual/finite-value contract."""

def _vector(value,size,name):
    a=np.asarray(value,dtype=float)
    if a.shape!=(size,) or not np.all(np.isfinite(a)):
        raise ValueError(f'{name} must have {size} finite real components')
    return tuple(float(x) for x in a)

def _qmul(q,p):
    # Hamilton product on arbitrary quaternions (the pure angular-velocity
    # quaternion is not a unit quaternion).
    w,x,y,z=q;v,a,b,c=p
    return np.array((w*v-x*a-y*b-z*c,w*a+x*v+y*c-z*b,
                     w*b-x*c+y*v+z*a,w*c+x*b-y*a+z*v))

def _hat(v):
    return np.asarray(rotations.hat(v))

def _finite_result(value,name):
    if not np.all(np.isfinite(value)):
        raise ConvergenceError(f'nonfinite {name}')
    return value

@dataclass(frozen=True)
class RigidBodyModel:
    mass: float
    inertia_body: tuple[tuple[float,...],...]
    _j: np.ndarray=field(init=False,repr=False,compare=False)
    _inverse: np.ndarray=field(init=False,repr=False,compare=False)

    def __post_init__(self):
        if not math.isfinite(self.mass) or self.mass<=0:
            raise ValueError('positive finite mass required')
        j=np.asarray(self.inertia_body,dtype=float)
        if j.shape!=(3,3) or not np.all(np.isfinite(j)):
            raise ValueError('finite 3x3 body inertia required')
        scale=float(np.max(np.abs(j)))
        if scale==0 or np.max(np.abs(j-j.T))>1e-12*scale:
            raise ValueError('symmetric body inertia required')
        # Roundoff-sized asymmetry is imported as its symmetric part.
        j=j/2+j.T/2
        try: values=np.linalg.eigvalsh(j)
        except np.linalg.LinAlgError as exc: raise ValueError("inertia eigensolve failed") from exc
        if not np.all(np.isfinite(values)): raise ValueError("finite principal moments required")
        if not np.all(values>0): raise ValueError('positive definite inertia required')
        if values[-1]>values[0]+values[1]+1e-12*scale:
            raise ValueError('inertia principal moments must satisfy triangle inequalities')
        try: inverse=np.linalg.inv(j)
        except np.linalg.LinAlgError as exc: raise ValueError('numerically invertible inertia required') from exc
        if not np.all(np.isfinite(inverse)):
            raise ValueError('inertia inverse must be finite at float64 precision')
        j.setflags(write=False);inverse.setflags(write=False)
        object.__setattr__(self,'mass',float(self.mass))
        object.__setattr__(self,'inertia_body',tuple(tuple(float(x) for x in row) for row in j))
        object.__setattr__(self,'_j',j)
        object.__setattr__(self,'_inverse',inverse)

    @classmethod
    def from_box(cls,mass,dimensions):
        a,b,c=_vector(dimensions,3,'box dimensions')
        if min(a,b,c)<=0: raise ValueError('positive box dimensions required')
        return cls(mass,((mass*(b*b+c*c)/12,0.,0.),
                         (0.,mass*(a*a+c*c)/12,0.),
                         (0.,0.,mass*(a*a+b*b)/12)))

@dataclass(frozen=True)
class RigidBodyState:
    position: tuple[float,...]=ZERO
    velocity: tuple[float,...]=ZERO
    quaternion: tuple[float,...]=IDENTITY_Q
    omega_body: tuple[float,...]=ZERO

    def __post_init__(self):
        for name in ('position','velocity','omega_body'):
            object.__setattr__(self,name,_vector(getattr(self,name),3,name))
        q=_vector(self.quaternion,4,'quaternion')
        squared=sum(x*x for x in q)
        if not math.isfinite(squared) or abs(squared-1)>1e-10:
            raise ValueError('unit scalar-first Hamilton quaternion required')
        # Accepted roundoff is normalized; q/-q sign is never canonicalized.
        norm=math.sqrt(squared)
        object.__setattr__(self,'quaternion',tuple(x/norm for x in q))

    def as_vector(self):
        return (*self.position,*self.velocity,*self.quaternion,*self.omega_body)

    @classmethod
    def from_vector(cls,value):
        value=_vector(value,13,'state')
        return cls(value[:3],value[3:6],value[6:10],value[10:])

@dataclass(frozen=True)
class StepInfo:
    iterations: int
    residual: float
    residual_bound: float

def rhs(model,state,force_body=ZERO,torque_body=ZERO,force_inertial=ZERO):
    """13-component continuous RHS, independent of the discrete integrator."""
    f=np.asarray(_vector(force_body,3,'body force'))
    tau=np.asarray(_vector(torque_body,3,'body torque'))
    fi=np.asarray(_vector(force_inertial,3,'inertial force'))
    omega=np.asarray(state.omega_body);q=np.asarray(state.quaternion)
    acceleration=(np.asarray(rotations.from_quaternion(q))@f+fi)/model.mass
    qdot=.5*_qmul(q,(0.,*omega))
    omegadot=model._inverse@(tau-np.cross(omega,model._j@omega))
    result=np.concatenate((state.velocity,acceleration,qdot,omegadot))
    return tuple(float(x) for x in _finite_result(result,'RHS'))

def step_with_info(model,state,dt,force_body=ZERO,torque_body=ZERO,
                   force_inertial=ZERO,*,tolerance=1e-12,max_iterations=20):
    """One implicit midpoint momentum / Cayley attitude / midpoint-force step.

    Newton uses the analytic Jacobian and stops only if its residual infinity
    norm <= tolerance*(1+||h0||inf+dt*||tau||inf). The 1 supplies an explicit
    absolute tolerance in the chosen momentum units. No time-step cap is
    imposed; a large or ill-conditioned step can raise ConvergenceError.
    """
    if not math.isfinite(dt) or dt<=0: raise ValueError('positive finite dt required')
    if not math.isfinite(tolerance) or tolerance<=0 or tolerance>1e-4:
        raise ValueError('tolerance must be finite and in (0,1e-4]')
    if not isinstance(max_iterations,int) or isinstance(max_iterations,bool) or not 1<=max_iterations<=100:
        raise ValueError('max_iterations must be an integer in [1,100]')
    f=np.asarray(_vector(force_body,3,'body force'))
    tau=np.asarray(_vector(torque_body,3,'body torque'))
    fi=np.asarray(_vector(force_inertial,3,'inertial force'))
    omega0=np.asarray(state.omega_body);h0=model._j@omega0
    with np.errstate(over='ignore',invalid='ignore',divide='ignore'):
        h1=h0+dt*(tau-np.cross(omega0,h0))
        bound=tolerance*(1+np.linalg.norm(h0,ord=np.inf)+dt*np.linalg.norm(tau,ord=np.inf))
        if not math.isfinite(bound): raise ConvergenceError('nonfinite residual bound')
        for iteration in range(max_iterations+1):
            _finite_result(h1,'momentum iterate')
            hm=(h0+h1)/2;wm=model._inverse@hm
            residual=h1-h0-dt*tau+dt*np.cross(wm,hm)
            _finite_result(residual,'Newton residual')
            residual_norm=float(np.linalg.norm(residual,ord=np.inf))
            if residual_norm<=bound: break
            if iteration==max_iterations:
                raise ConvergenceError(f'Newton did not converge: residual {residual_norm:g} > {bound:g}')
            jacobian=np.eye(3)+dt/2*(_hat(wm)-_hat(hm)@model._inverse)
            _finite_result(jacobian,'Newton Jacobian')
            try: delta=np.linalg.solve(jacobian,residual)
            except np.linalg.LinAlgError as exc:
                raise ConvergenceError('singular Newton Jacobian') from exc
            h1=h1-delta
        # Cayley increment: tan(angle/2)=dt*|omega_mid|/2.
        increment=np.array((1.,*(dt*wm/2)))
        _finite_result(increment,'attitude increment')
        incnorm=math.hypot(*increment)
        increment/=incnorm
        half_scalar=math.sqrt((1+increment[0])/2)
        half=np.array((half_scalar,*(increment[1:]/(2*half_scalar))))
        q0=np.asarray(state.quaternion)
        qmid=_qmul(q0,half)
        qmid/=math.hypot(*qmid)
        acceleration=(np.asarray(rotations.from_quaternion(qmid))@f+fi)/model.mass
        v0=np.asarray(state.velocity);v1=v0+dt*acceleration
        r1=np.asarray(state.position)+dt*(v0+v1)/2
        q1=_qmul(q0,increment)
        q1/=math.hypot(*q1)
        omega1=model._inverse@h1
        output=np.concatenate((r1,v1,q1,omega1))
        _finite_result(output,'updated state')
    return (RigidBodyState(tuple(r1),tuple(v1),tuple(q1),tuple(omega1)),
            StepInfo(iteration,residual_norm,float(bound)))

def step(model,state,dt,force_body=ZERO,torque_body=ZERO,force_inertial=ZERO,
         *,tolerance=1e-12,max_iterations=20):
    return step_with_info(model,state,dt,force_body,torque_body,force_inertial,
                          tolerance=tolerance,max_iterations=max_iterations)[0]

def diagnostics(model,state):
    omega=np.asarray(state.omega_body);v=np.asarray(state.velocity)
    p=model.mass*v
    h=np.asarray(rotations.from_quaternion(state.quaternion))@(model._j@omega)
    rotational=float(omega@(model._j@omega)/2)
    translational=float(model.mass*(v@v)/2)
    _finite_result(np.concatenate((p,h,[rotational,translational,rotational+translational])), 'diagnostics')
    orbital=h+np.cross(state.position,p)
    _finite_result(orbital,'angular momentum about origin')
    return {'energy':rotational+translational,'rotational_energy':rotational,
            'translational_energy':translational,
            'angular_momentum_inertial':h.tolist(),
            'linear_momentum_inertial':p.tolist(),
            'total_angular_momentum_about_origin':orbital.tolist(),
            'quaternion_norm':math.hypot(*state.quaternion)}
