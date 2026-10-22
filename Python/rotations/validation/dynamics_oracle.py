"""Independent DOP853/exact solutions and native-fixture comparison for 6-DOF.

No discrete integrator internals are reused by the reference ODE or closed
solutions. Exact symmetric-top attitude includes both precession and spin.
"""
from __future__ import annotations
import json
import math
import pathlib
import sys
import numpy as np
from scipy.integrate import solve_ivp
import scipy
import rigid_body as rb

def hamilton(q,p):
    w,x,y,z=q;v,a,b,c=p
    return np.array([w*v-x*a-y*b-z*c,w*a+x*v+y*c-z*b,
                     w*b-x*c+y*v+z*a,w*c+x*b-y*a+z*v])

def quaternion_axis(axis,phi):
    return np.array([math.cos(phi/2),*(math.sin(phi/2)*np.asarray(axis))])

def rotation_matrix(q):
    # Independent homogeneous quaternion matrix, evaluated on near-unit
    # DOP853 state; no step-time re-normalization or shared conversion code.
    w,x,y,z=q
    return np.array([[w*w+x*x-y*y-z*z,2*(x*y-w*z),2*(x*z+w*y)],
                     [2*(x*y+w*z),w*w-x*x+y*y-z*z,2*(y*z-w*x)],
                     [2*(x*z-w*y),2*(y*z+w*x),w*w-x*x-y*y+z*z]])

def reference(model,state,time,force_body=(0,0,0),torque_body=(0,0,0),force_inertial=(0,0,0),*,samples=None):
    j=np.asarray(model.inertia_body)
    f,tau,fi=map(np.asarray,(force_body,torque_body,force_inertial))
    def ode(_time,y):
        position,velocity,q,omega=y[:3],y[3:6],y[6:10],y[10:]
        angular=np.linalg.solve(j,tau-np.cross(omega,j@omega))
        return np.concatenate((velocity,(rotation_matrix(q)@f+fi)/model.mass,
                               hamilton(q,np.array([0.,*omega]))/2,angular))
    result=solve_ivp(ode,(0.,time),state.as_vector(),method='DOP853',rtol=2e-13,atol=2e-15,t_eval=samples)
    if not result.success or not np.all(np.isfinite(result.y)):
        raise AssertionError(result.message)
    return result.y if samples is not None else result.y[:,-1]

def attitude_error(q,p):
    # Sign-invariant quaternion chordal distance. For small error this is
    # half the physical angle; avoids acos(trace) cancellation near identity.
    q=np.asarray(q);p=np.asarray(p)
    return float(min(np.linalg.norm(q-p),np.linalg.norm(q+p)))

def exact_free_spherical(state,time):
    omega=np.asarray(state.omega_body);speed=np.linalg.norm(omega)
    inc=quaternion_axis(omega/speed,speed*time) if speed else np.array([1.,0.,0.,0.])
    return hamilton(state.quaternion,inc)

def exact_free_symmetric_top(model,state,time):
    j=np.asarray(model.inertia_body)
    if np.max(np.abs(j-np.diag(np.diag(j))))>1e-14 or abs(j[0,0]-j[1,1])>1e-14:
        raise ValueError('reference expects principal inertia diag(I,I,C)')
    transverse,axial=j[0,0],j[2,2]
    omega=np.asarray(state.omega_body)
    momentum=rotation_matrix(state.quaternion)@(j@omega)
    magnitude=np.linalg.norm(momentum)
    precession=quaternion_axis(momentum/magnitude,magnitude*time/transverse) if magnitude else np.array([1.,0.,0.,0.])
    rate=(axial-transverse)*omega[2]/transverse
    spin=quaternion_axis((0,0,1),-rate*time)
    q=hamilton(hamilton(precession,state.quaternion),spin)
    co,si=math.cos(rate*time),math.sin(rate*time)
    wb=np.array([co*omega[0]-si*omega[1],si*omega[0]+co*omega[1],omega[2]])
    return q,wb

def integrate(model,state,dt,steps,**wrench):
    for _ in range(steps): state=rb.step(model,state,dt,**wrench)
    return state

def fixtures():
    q0=(1.,0.,0.,0.);p=(.1,-.2,.3);v=(.3,.1,-.2)
    sphere=((2.,0.,0.),(0.,2.,0.),(0.,0.,2.))
    asym=((1.,0.,0.),(0.,2.,0.),(0.,0.,2.5))
    off=((2.,.2,.1),(.2,3.,.3),(.1,.3,4.))
    specs=[
        ('free_sphere',sphere,(.2,-.4,.7),.01,100,(0,0,0),(0,0,0),q0),
        ('free_asymmetric',asym,(.2,1.,.3),.005,400,(0,0,0),(0,0,0),q0),
        ('free_offdiag',off,(.3,-.4,.8),.01,100,(0,0,0),(0,0,0),q0),
        ('principal',asym,(0,0,.7),.01,100,(0,0,0),(0,0,0),q0),
        ('forced_translation',sphere,(0,0,0),.01,100,(2.,-1.,.5),(0,0,0),q0),
        ('forced_torque',sphere,(0,0,0),.01,100,(0,0,0),(.1,-.2,.3),q0),
        ('forced_coupled',off,(.3,-.4,.8),.002,100,(.2,-.1,.3),(.03,-.02,.01),q0),
        ('nonidentity_offdiag',off,(.3,-.4,.8),.01,100,(0,0,0),(0,0,0),(math.cos(.3),math.sin(.3),0,0)),
    ]
    return [(name,rb.RigidBodyModel(2.,j),rb.RigidBodyState(p,v,q,w),dt,count,{'force_body':force,'torque_body':torque})
            for name,j,w,dt,count,force,torque,q in specs]

def run_reference_checks():
    rows=[]
    for name,model,initial,dt,count,wrench in fixtures():
        final=integrate(model,initial,dt,count,**wrench)
        ref=reference(model,initial,dt*count,**wrench)
        errors={'position':float(np.linalg.norm(np.asarray(final.position)-ref[:3])),
                'velocity':float(np.linalg.norm(np.asarray(final.velocity)-ref[3:6])),
                'attitude':attitude_error(final.quaternion,ref[6:10]),
                'omega':float(np.linalg.norm(np.asarray(final.omega_body)-ref[10:]))}
        if max(errors.values())>1e-4: raise AssertionError((name,errors))
        rows.append({'name':name,'dt':dt,'steps':count,'DOP853_errors':errors})
    # Nontrivial symmetric top and offdiagonal inertia: second-order full state,
    # not just an invariant scalar. Decreasing dt should divide errors by ~4.
    model=rb.RigidBodyModel(2.,((2.,0,0),(0,2.,0),(0,0,3.)))
    initial=rb.RigidBodyState(quaternion=(math.cos(.2),0,math.sin(.2),0),omega_body=(.4,-.7,1.2))
    expected_q,expected_w=exact_free_symmetric_top(model,initial,2.)
    dop=reference(model,initial,2.)
    if attitude_error(expected_q,dop[6:10])>2e-11 or np.linalg.norm(expected_w-dop[10:])>2e-11:
        raise AssertionError('closed symmetric-top solution disagrees with independent ODE')
    errors=[]
    for dt in (.04,.02,.01):
        state=integrate(model,initial,dt,round(2/dt))
        errors.append(attitude_error(state.quaternion,expected_q))
    ratios=[errors[i]/errors[i+1] for i in range(2)]
    if not all(3.8<ratio<4.2 for ratio in ratios): raise AssertionError(('symmetric-top order',errors,ratios))
    rows.append({'name':'symmetric_top_full_attitude_order','dt':[.04,.02,.01],
                 'errors':errors,'ratios':ratios,'closed_form_DOP853_attitude_error':attitude_error(expected_q,dop[6:10])})
    return rows

def native_parity(paths):
    results=[]
    for path in paths:
        records=[json.loads(line) for line in pathlib.Path(path).read_text().splitlines() if line.strip()]
        if not records: raise AssertionError('empty native dynamics producer')
        maximum=0.
        for row in records:
            model=rb.RigidBodyModel(row['mass'],row['inertia'])
            value=row['initial']
            initial=rb.RigidBodyState(value['r'],value['v'],value['q'],value['omega'])
            wrench=row['wrench']
            final=integrate(model,initial,row['dt'],row['steps'],force_body=wrench['force'],torque_body=wrench['torque'],force_inertial=wrench.get('force_inertial',(0,0,0)))
            output=row['final']
            error=max(np.max(np.abs(np.asarray(final.position)-output['r'])),
                      np.max(np.abs(np.asarray(final.velocity)-output['v'])),
                      np.max(np.abs(np.asarray(final.omega_body)-output['omega'])),
                      attitude_error(final.quaternion,output['q']))
            diagnostic=rb.diagnostics(model,final)
            error=max(error,abs(float(row['energy'])-diagnostic['energy']),
                      float(np.max(np.abs(np.asarray(row['H_I'])-diagnostic['angular_momentum_inertial']))))
            if not math.isfinite(error) or error>2e-10: raise AssertionError((row['name'],error))
            maximum=max(maximum,float(error))
        results.append({'producer':str(path),'cases':len(records),'max_error':maximum})
    return results

def symbolic_midpoint_checks():
    import sympy as sym
    w=sym.Matrix(sym.symbols('w0:3',real=True))
    h=sym.Matrix(sym.symbols('h0:3',real=True))
    tau=sym.Matrix(sym.symbols('tau0:3',real=True))
    dt=sym.symbols('dt',real=True)
    omega=sym.Matrix([[0,-w[2],w[1]],[w[2],0,-w[0]],[-w[1],w[0],0]])
    minus=sym.eye(3)-dt*omega/2
    plus=sym.eye(3)+dt*omega/2
    # h1=plus^-1 minus h0; Rinc=minus^-1 plus. The two
    # factors commute as polynomials in the same skew matrix.
    Rinc=minus.inv()*plus
    hstep=plus.inv()*minus
    residual=Rinc*hstep-sym.eye(3)
    if any(sym.cancel(value)!=0 for value in residual):
        raise AssertionError('Cayley inertial momentum identity')
    energy=sym.expand(w.dot(-dt*w.cross(h)))
    work=sym.expand(w.dot(dt*(tau-w.cross(h)))-dt*w.dot(tau))
    if energy!=0 or work!=0: raise AssertionError('midpoint energy/work identity')
    return ['Cayley_inverse_momentum_update','free_midpoint_energy','midpoint_torque_work']

if __name__=='__main__':
    print(json.dumps({'status':'PASS','reference':'independent scipy.integrate.solve_ivp DOP853 + closed full-attitude solutions',
                      'numpy':np.__version__,'scipy':scipy.__version__,
                      'checks':run_reference_checks(),'symbolic_midpoint_checks':symbolic_midpoint_checks(),'native_parity':native_parity(sys.argv[1:])}))
