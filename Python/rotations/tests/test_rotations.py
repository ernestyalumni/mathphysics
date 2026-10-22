import math
import random
import unittest
import rotations as r

class Rotations(unittest.TestCase):
    def close(self,a,b,tol=5e-12):
        self.assertLessEqual(max(abs(x-y) for ra,rb in zip(a,b) for x,y in zip(ra,rb)),tol)

    def test_independent_quarter_turns(self):
        # Direct physical actions fix handedness and row/column convention.
        self.close(r.elementary('z',math.pi/2),((0,-1,0),(1,0,0),(0,0,1)))
        self.close(r.elementary('x',math.pi/2),((1,0,0),(0,0,-1),(0,1,0)))
        self.close(r.elementary('y',math.pi/2),((0,0,1),(0,1,0),(-1,0,0)))
        self.assertGreater(abs(r.apply(r.elementary('z',math.pi/2),(1,0,0))[1]-(-1)),1.)

    def test_noncommuting_and_zxz_zyz_dictionary(self):
        a,b,c=.43,.79,-.27
        self.close(r.from_euler('zxz',a,b,c),r.from_euler('zyz',a-math.pi/2,b,c+math.pi/2))
        x,y=r.elementary('x',a),r.elementary('y',b)
        self.assertGreater(max(abs(v-w) for ra,rb in zip(r.multiply(x,y),r.multiply(y,x)) for v,w in zip(ra,rb)),.1)

    def test_branch_stress_reconstruction(self):
        rng=random.Random(476)
        cases=0
        for convention in ('zxz','zyz','zyx'):
            mids=[0.,1e-14,1e-10,.7,math.pi-1e-10,math.pi-1e-14,math.pi]
            if convention=='zyx': mids=[-math.pi/2,-math.pi/2+1e-14,-math.pi/2+1e-10,0.,math.pi/2-1e-10,math.pi/2-1e-14,math.pi/2]
            triples=[(.4,m,-.9) for m in mids]
            triples += [tuple(rng.uniform(-8,8) for _ in range(3)) for _ in range(1000)]
            for angles in triples:
                matrix=r.from_euler(convention,*angles)
                recovered=r.to_euler(matrix,convention)
                self.close(matrix,r.from_euler(convention,*recovered.values))
                if convention=='zyx': self.assertLessEqual(abs(recovered.middle),math.pi/2)
                else: self.assertTrue(0<=recovered.middle<=math.pi)
                q=r.to_quaternion(matrix)
                self.close(matrix,r.from_quaternion(q))
                self.close(matrix,r.from_quaternion(tuple(-x for x in q)))
                self.close(r.transpose(matrix),r.ck_to_goldstein(r.quaternion_to_ck(q)))
                axis,phi=r.to_axis_angle(matrix)
                self.close(matrix,r.from_axis_angle(axis,phi))
                self.close(r.cayley_hamilton_residual(matrix),((0,0,0),)*3)
                cases+=1
        self.assertEqual(cases,3021)

    def test_axis_angle_independent_vector_action_and_spectrum(self):
        rng=random.Random(462)
        for phi in (0.,1e-14,1e-9,.31,math.pi-1e-12,math.pi,math.pi+1e-12,2*math.pi):
            for _ in range(30):
                n=[rng.uniform(-1,1) for _ in range(3)];norm=math.hypot(*n);n=[x/norm for x in n]
                a=r.from_axis_angle(n,phi)
                vector=[rng.uniform(-2,2) for _ in range(3)]
                dot=sum(x*y for x,y in zip(n,vector))
                cross=(n[1]*vector[2]-n[2]*vector[1],n[2]*vector[0]-n[0]*vector[2],n[0]*vector[1]-n[1]*vector[0])
                reference=tuple(math.cos(phi)*vector[i]+(1-math.cos(phi))*dot*n[i]+math.sin(phi)*cross[i] for i in range(3))
                self.close((r.apply(a,vector),),(reference,))
                self.close((r.apply(a,n),),(n,))
                self.assertAlmostEqual(sum(a[i][i] for i in range(3)),1+2*math.cos(phi),places=12)
                recovered,angle=r.to_axis_angle(a)
                self.close(a,r.from_axis_angle(recovered,angle))
                self.assertTrue(0<=angle<=math.pi)

    def test_quaternion_composition(self):
        a=r.from_euler('zxz',.4,-.8,.2);b=r.from_euler('zyx',-.5,.3,-.9)
        q=r.quaternion_product(r.to_quaternion(a),r.to_quaternion(b))
        self.close(r.from_quaternion(q),r.multiply(a,b))

    def test_conjugation_physical_action(self):
        b=r.from_euler('zyx',.4,-.8,.9);a=r.from_euler('zxz',.3,.7,.1);x=(.3,-.6,.4)
        new=r.multiply(r.multiply(b,a),r.transpose(b))
        self.close((r.apply(new,r.apply(b,x)),),(r.apply(b,r.apply(a,x)),))
        self.close(r.multiply(r.multiply(b,r.hat(x)),r.transpose(b)),r.hat(r.apply(b,x)))
        # Missing inverse must be detected.
        wrong=r.multiply(r.multiply(b,a),b)
        self.assertGreater(max(abs(x-y) for ra,rb in zip(new,wrong) for x,y in zip(ra,rb)),.1)

    def test_chasles_screw_action(self):
        axis=(2/3,-1/3,2/3);translation=(.3,-.8,.6);point=(-.9,.2,1.1)
        distance=sum(a*b for a,b in zip(axis,translation))
        perp=tuple(translation[i]-distance*axis[i] for i in range(3))
        cross=r.apply(r.hat(axis),perp)
        for phi in (.01,.3,math.pi,2*math.pi-.01):
            matrix=r.from_axis_angle(axis,phi)
            center=tuple((perp[i]+cross[i]/math.tan(phi/2))/2 for i in range(3))
            lhs=tuple(a+b for a,b in zip(translation,r.apply(matrix,point)))
            rotated=r.apply(matrix,tuple(point[i]-center[i] for i in range(3)))
            rhs=tuple(center[i]+rotated[i]+distance*axis[i] for i in range(3))
            self.close((lhs,),(rhs,))

    def test_singular_metadata(self):
        for convention,middle in [('zxz',0),('zxz',math.pi),('zyz',0),('zyz',math.pi),('zyx',math.pi/2),('zyx',-math.pi/2)]:
            out=r.to_euler(r.from_euler(convention,.4,middle,-.3),convention)
            self.assertTrue(out.singular)
            self.assertEqual(out.last,0.)

    def test_generic_shape_rejection(self):
        for bad in ((),((),),((1,2),(3,))):
            with self.assertRaises(ValueError): r.transpose(bad)
            with self.assertRaises(ValueError): r.multiply(bad,r.I3)
        with self.assertRaises(ValueError): r.multiply(((1,2),),((3,4),))
        with self.assertRaises(ValueError): r.apply(r.I3,(1,2))
        with self.assertRaises(ValueError): r.apply(r.I3,(1,2,3,4))
        with self.assertRaises(ValueError): r.det(((1,2),(3,4)))
        self.assertEqual(r.multiply(((1,2),),((3,),(4,))),((11,),))
        self.assertEqual(r.apply(((1,2),),(3,4)),(11,))

    def test_invalid_imports(self):
        bad_matrices=[((-1,0,0),(0,1,0),(0,0,1)),((1,1,0),(0,1,0),(0,0,1)),((1,0),(0,1))]
        for value in (float('nan'),float('inf'),-float('inf')):
            bad_matrices.append(((value,0,0),(0,1,0),(0,0,1)))
            for constructor,args in [(r.elementary,('z',value)),(r.from_quaternion,((value,0,0,1),)),(r.from_axis_angle,((1,0,0),value))]:
                with self.assertRaises(ValueError): constructor(*args)
        for a in bad_matrices:
            with self.assertRaises(ValueError): r.validate_matrix(a)
        for q in ((0,0,0,0),(2,0,0,0),(1,0,0)):
            with self.assertRaises(ValueError): r.from_quaternion(q)
        for tol in (0,-1,float('nan'),float('inf'),1):
            with self.assertRaises(ValueError): r.validate_matrix(r.I3,tol)
        with self.assertRaises(ValueError): r.from_axis_angle((0,0,0),.3)
        with self.assertRaises(ValueError): r.from_axis_angle((1.00000000004,0,0),math.pi)
        with self.assertRaises(ValueError): r.from_quaternion((0,1.00000000004,0,0))
        with self.assertRaises(ValueError): r.from_euler('xyz',1,2,3)
        with self.assertRaises(ValueError): r.ck_to_quaternion((1,1,1,1))
        ck=r.quaternion_to_ck(r.to_quaternion(r.from_euler('zxz',.3,.4,.5)))
        with self.assertRaises(ValueError): r.ck_to_quaternion((ck[0],ck[1],-ck[2],ck[3]))

if __name__=='__main__': unittest.main()
