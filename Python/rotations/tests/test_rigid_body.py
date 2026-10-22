import math
import unittest
import numpy as np
import rigid_body as rb
import rotations
from validation.dynamics_oracle import (attitude_error,exact_free_spherical,
    exact_free_symmetric_top,reference,integrate)

class RigidBody(unittest.TestCase):
    def setUp(self):
        self.model=rb.RigidBodyModel(2.,((2.,.2,.1),(.2,3.,.3),(.1,.3,4.)))
        self.initial=rb.RigidBodyState((.1,-.2,.3),(.3,.1,-.2),
            (math.cos(.3),math.sin(.3),0,0),(.3,-.4,.8))

    def test_model_box_and_immutability(self):
        model=rb.RigidBodyModel.from_box(3.,(2.,3.,4.))
        np.testing.assert_allclose(np.diag(model.inertia_body),(6.25,5.,3.25))
        with self.assertRaises(ValueError): model._j[0,0]=9
        # Equality of two principal moments is allowed, as is a planar
        # realizability boundary (largest moment equals sum of others).
        rb.RigidBodyModel(1.,((1,0,0),(0,1,0),(0,0,2)))

    def test_invalid_models_states_controls(self):
        for mass in (0.,-1.,float('nan'),float('inf')):
            with self.assertRaises(ValueError): rb.RigidBodyModel(mass,np.eye(3))
        for j in (np.zeros((3,3)),np.diag((-1,2,2)),np.diag((1,1,3)),
                  ((1,.3,0),(0,1,0),(0,0,1)),np.eye(2),np.full((3,3),float('nan'))):
            with self.assertRaises(ValueError): rb.RigidBodyModel(1.,j)
        for q in ((0,0,0,0),(2,0,0,0),(1,2,3),(float('nan'),0,0,1)):
            with self.assertRaises(ValueError): rb.RigidBodyState(quaternion=q)
        for dt in (0.,-1.,float('nan'),float('inf')):
            with self.assertRaises(ValueError): rb.step(self.model,self.initial,dt)
        for tol in (0.,-1.,float('nan'),float('inf'),1.):
            with self.assertRaises(ValueError): rb.step(self.model,self.initial,.01,tolerance=tol)
        for iterations in (0,-1,101,1.1,True):
            with self.assertRaises(ValueError): rb.step(self.model,self.initial,.01,max_iterations=iterations)
        with self.assertRaises(ValueError): rb.step(self.model,self.initial,.01,force_body=(1,2))
        with self.assertRaises(ValueError): rb.step(self.model,self.initial,.01,torque_body=(0,float('inf'),0))

    def test_rhs_and_state_dictionary(self):
        state=self.initial
        self.assertEqual(rb.RigidBodyState.from_vector(state.as_vector()),state)
        derivative=np.asarray(rb.rhs(self.model,state,force_body=(.2,-.1,.3),torque_body=(.03,-.02,.01)))
        np.testing.assert_allclose(derivative[:3],state.velocity)
        # qdot is tangent to S3.
        self.assertLess(abs(np.dot(derivative[6:10],state.quaternion)),1e-15)
        np.testing.assert_allclose(derivative[10:],np.linalg.solve(self.model.inertia_body,
            np.array((.03,-.02,.01))-np.cross(state.omega_body,np.asarray(self.model.inertia_body)@state.omega_body)))

    def test_free_energy_full_inertial_momentum_and_linear_motion(self):
        before=rb.diagnostics(self.model,self.initial)
        state=self.initial
        max_residual=0.
        for _ in range(1000):
            state,info=rb.step_with_info(self.model,state,.01,tolerance=2e-14)
            self.assertLessEqual(info.residual,info.residual_bound)
            max_residual=max(max_residual,info.residual)
        after=rb.diagnostics(self.model,state)
        self.assertLess(abs(after['rotational_energy']-before['rotational_energy']),2e-11)
        np.testing.assert_allclose(after['angular_momentum_inertial'],before['angular_momentum_inertial'],atol=2e-11,rtol=0)
        np.testing.assert_allclose(state.velocity,self.initial.velocity,atol=0,rtol=0)
        np.testing.assert_allclose(state.position,np.asarray(self.initial.position)+10*np.asarray(self.initial.velocity),atol=2e-12,rtol=0)
        self.assertAlmostEqual(after['quaternion_norm'],1.,places=14)

    def test_spherical_and_principal_axis_full_attitude(self):
        for model,omega in [(rb.RigidBodyModel(1.,2*np.eye(3)),(.2,-.4,.7)),
                            (rb.RigidBodyModel(1.,np.diag((1.,2.,2.5))),(0.,0.,.7))]:
            initial=rb.RigidBodyState(quaternion=self.initial.quaternion,omega_body=omega)
            expected=exact_free_spherical(initial,2.)
            errors=[]
            for dt in (.04,.02,.01):
                state=integrate(model,initial,dt,round(2/dt))
                errors.append(attitude_error(state.quaternion,expected))
                np.testing.assert_allclose(state.omega_body,omega,atol=1e-13)
            self.assertGreater(errors[-1],1e-8) # invariant conservation does not remove phase error
            self.assertTrue(all(3.9<errors[i]/errors[i+1]<4.1 for i in range(2)))

    def test_symmetric_top_full_attitude_and_omega(self):
        model=rb.RigidBodyModel(1.,np.diag((2.,2.,3.)))
        initial=rb.RigidBodyState(quaternion=self.initial.quaternion,omega_body=(.4,-.7,1.2))
        q,w=exact_free_symmetric_top(model,initial,2.)
        dop=reference(model,initial,2.)
        self.assertLess(attitude_error(q,dop[6:10]),2e-11)
        np.testing.assert_allclose(w,dop[10:],atol=2e-11,rtol=0)
        previous=None
        for dt in (.04,.02,.01):
            state=integrate(model,initial,dt,round(2/dt))
            error=attitude_error(state.quaternion,q)
            if previous is not None: self.assertTrue(3.8<previous/error<4.2)
            previous=error

    def test_intermediate_axis_flip_and_DOP853(self):
        model=rb.RigidBodyModel(1.,np.diag((1.,2.,2.5)))
        initial=rb.RigidBodyState(omega_body=(.01,1.,.01))
        state=initial;minimum=1.
        for _ in range(4800):
            state=rb.step(model,state,.005,tolerance=2e-14)
            minimum=min(minimum,state.omega_body[1])
        self.assertLess(minimum,-.9)
        dop=reference(model,initial,24.)
        self.assertLess(attitude_error(state.quaternion,dop[6:10]),1e-4)
        np.testing.assert_allclose(state.omega_body,dop[10:],atol=1e-4,rtol=0)

    def test_constant_force_and_torque(self):
        model=rb.RigidBodyModel(2.,2*np.eye(3));initial=rb.RigidBodyState()
        state=integrate(model,initial,.01,100,force_body=(2,-1,.5))
        np.testing.assert_allclose(state.velocity,(1.,-.5,.25),atol=1e-13)
        np.testing.assert_allclose(state.position,(.5,-.25,.125),atol=1e-13)
        state=integrate(model,initial,.01,100,torque_body=(.1,-.2,.3))
        np.testing.assert_allclose(state.omega_body,(.05,-.1,.15),atol=2e-14)
        # Work-energy identity of the midpoint angular momentum step.
        nextstate=rb.step(self.model,self.initial,.02,torque_body=(.2,-.1,.4),tolerance=2e-14)
        de=rb.diagnostics(self.model,nextstate)['rotational_energy']-rb.diagnostics(self.model,self.initial)['rotational_energy']
        work=.02*np.dot((np.asarray(self.initial.omega_body)+nextstate.omega_body)/2,(.2,-.1,.4))
        self.assertAlmostEqual(de,work,places=12)

    def test_midpoint_rotating_body_force_and_offdiagonal_DOP853(self):
        f=(.2,-.1,.3);tau=(.03,-.02,.01)
        expected=reference(self.model,self.initial,1.,f,tau)
        errors=[]
        for dt in (.04,.02,.01):
            state=integrate(self.model,self.initial,dt,round(1/dt),force_body=f,torque_body=tau)
            errors.append(np.linalg.norm(np.asarray(state.as_vector())-expected))
        self.assertTrue(all(3.8<errors[i]/errors[i+1]<4.2 for i in range(2)))
        # Inertial force is a separate explicit input.
        state=integrate(self.model,rb.RigidBodyState(),.01,100,force_inertial=(0,0,-19.62))
        np.testing.assert_allclose(state.position,(0,0,-4.905),atol=2e-13)

    def test_quaternion_sign_continuity_and_newton_failure(self):
        negative=rb.RigidBodyState(quaternion=(-1.,0,0,0),omega_body=(0,0,.2))
        state=rb.step(self.model,negative,.01)
        self.assertLess(state.quaternion[0],0.)
        with self.assertRaises(rb.ConvergenceError):
            rb.step(self.model,self.initial,.5,tolerance=1e-16,max_iterations=1)
        with self.assertRaises(rb.ConvergenceError):
            rb.step(self.model,self.initial,1e308,torque_body=(1e308,0,0))

if __name__=='__main__': unittest.main()
