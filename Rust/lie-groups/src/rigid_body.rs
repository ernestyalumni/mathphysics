//! COM-centered numerical 6DOF dynamics, f64. r/v are inertial; omega and
//! inertia are body components. Scalar-first Hamilton q maps body to inertial.
//! No gravity, damping or external wrench is implicit. A step holds its wrench
//! fixed. Accepted quaternion roundoff is normalized without changing lift sign.
use crate::{default_tolerance, su2, vector, Error, Matrix, Su2};

pub type Vec3 = [f64; 3];
pub type Mat3 = Matrix<f64, 3, 3>;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct State {
    pub position: Vec3,
    pub velocity: Vec3,
    pub quaternion: [f64; 4],
    pub omega_body: Vec3,
}
impl Default for State {
    fn default() -> Self {
        Self {
            position: [0.; 3],
            velocity: [0.; 3],
            quaternion: [1., 0., 0., 0.],
            omega_body: [0.; 3],
        }
    }
}
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Wrench {
    pub force_body: Vec3,
    pub torque_body: Vec3,
    pub force_inertial: Vec3,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Derivative {
    pub position: Vec3,
    pub velocity: Vec3,
    pub quaternion: [f64; 4],
    pub omega_body: Vec3,
}
#[derive(Clone, Copy, Debug)]
pub struct StepOptions {
    pub tolerance: f64,
    pub max_iterations: usize,
}
impl Default for StepOptions {
    fn default() -> Self {
        Self {
            tolerance: 1e-12,
            max_iterations: 20,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct StepResult {
    pub state: State,
    pub iterations: usize,
    pub residual: f64,
}
fn finite(v: &Vec3) -> Result<(), Error> {
    if v.iter().all(|x| x.is_finite()) {
        Ok(())
    } else {
        Err(Error::NonFinite)
    }
}
fn inf_norm(v: &Vec3) -> f64 {
    v.iter().fold(0., |a, b| a.max(b.abs()))
}
fn add(a: Vec3, b: Vec3) -> Vec3 {
    std::array::from_fn(|i| a[i] + b[i])
}
fn sub(a: Vec3, b: Vec3) -> Vec3 {
    std::array::from_fn(|i| a[i] - b[i])
}
fn scale(a: Vec3, x: f64) -> Vec3 {
    a.map(|v| v * x)
}
fn normalized(mut q: [f64; 4]) -> Result<[f64; 4], Error> {
    let n = q[0].hypot(q[1]).hypot(q[2]).hypot(q[3]);
    if !n.is_finite() || n == 0. {
        return Err(Error::NonFinite);
    }
    for x in &mut q {
        *x /= n;
    }
    Ok(q)
}
fn validated(mut s: State) -> Result<State, Error> {
    finite(&s.position)?;
    finite(&s.velocity)?;
    finite(&s.omega_body)?;
    Su2::try_from_unit_quaternion(s.quaternion, 1e-10)?;
    s.quaternion = normalized(s.quaternion)?;
    Ok(s)
}
fn validate_wrench(w: Wrench) -> Result<(), Error> {
    finite(&w.force_body)?;
    finite(&w.torque_body)?;
    finite(&w.force_inertial)
}
fn rotation(q: [f64; 4]) -> Result<crate::So3<f64>, Error> {
    Su2::try_from_unit_quaternion(q, default_tolerance())?.to_so3()
}
fn hat(v: Vec3) -> Mat3 {
    crate::So3Algebra::new(v).matrix()
}
// Pivoted nonsymmetric 3x3 Newton solve; J itself uses Cholesky below.
fn solve3(mut a: Mat3, mut b: Vec3) -> Result<Vec3, Error> {
    for k in 0..3 {
        let mut p = k;
        for i in k + 1..3 {
            if a[(i, k)].abs() > a[(p, k)].abs() {
                p = i;
            }
        }
        if !a[(p, k)].is_finite() || a[(p, k)] == 0. {
            return Err(Error::Domain("singular Newton Jacobian"));
        }
        if p != k {
            for j in k..3 {
                let v = a[(p, j)];
                a[(p, j)] = a[(k, j)];
                a[(k, j)] = v;
            }
            b.swap(p, k);
        }
        for i in k + 1..3 {
            let f = a[(i, k)] / a[(k, k)];
            for j in k + 1..3 {
                a[(i, j)] -= f * a[(k, j)];
            }
            b[i] -= f * b[k];
        }
    }
    let mut x = [0.; 3];
    for i in (0..3).rev() {
        let mut v = b[i];
        for j in i + 1..3 {
            v -= a[(i, j)] * x[j];
        }
        x[i] = v / a[(i, i)];
    }
    finite(&x)?;
    Ok(x)
}
fn largest_symmetric_eigenvalue(a: &Mat3) -> f64 {
    let p1 = a[(0, 1)].powi(2) + a[(0, 2)].powi(2) + a[(1, 2)].powi(2);
    if p1 == 0. {
        return a[(0, 0)].max(a[(1, 1)]).max(a[(2, 2)]);
    }
    let q = a.trace() / 3.;
    let mut b = a - &Mat3::identity().scaled(q);
    let p = ((b[(0, 0)].powi(2) + b[(1, 1)].powi(2) + b[(2, 2)].powi(2) + 2. * p1) / 6.).sqrt();
    if p == 0. {
        return q;
    }
    b = b.scaled(1. / p);
    let angle = (b.determinant() / 2.).clamp(-1., 1.).acos() / 3.;
    q + 2. * p * angle.cos()
}

#[derive(Clone, Debug)]
pub struct Model {
    mass: f64,
    inertia: Mat3,
    scale: f64,
    chol: Mat3,
    inverse: Mat3,
}
impl Model {
    /// Full COM inertia in body axes, SPD and satisfying principal-moment
    /// triangle inequalities. Symmetry roundoff <=1e-12 relative is averaged.
    pub fn new(mass: f64, mut inertia: Mat3) -> Result<Self, Error> {
        if !mass.is_finite() || mass <= 0. {
            return Err(Error::Domain("mass must be positive finite"));
        }
        let mut scale: f64 = 0.;
        for x in inertia.rows().iter().flatten() {
            if !x.is_finite() {
                return Err(Error::NonFinite);
            }
            scale = scale.max(x.abs());
        }
        if scale == 0. {
            return Err(Error::Domain("inertia must be SPD"));
        }
        for i in 0..3 {
            for j in i + 1..3 {
                if (inertia[(i, j)] / scale - inertia[(j, i)] / scale).abs() > 1e-12 {
                    return Err(Error::Domain("inertia must be symmetric"));
                }
                let mean = inertia[(i, j)] / 2. + inertia[(j, i)] / 2.;
                inertia[(i, j)] = mean;
                inertia[(j, i)] = mean;
            }
        }
        let a = Mat3::from_rows(std::array::from_fn(|i| {
            std::array::from_fn(|j| inertia[(i, j)] / scale)
        }));
        let mut chol = Mat3::zero();
        for i in 0..3 {
            for j in 0..=i {
                let mut x = a[(i, j)];
                for k in 0..j {
                    x -= chol[(i, k)] * chol[(j, k)];
                }
                if i == j {
                    if !x.is_finite() || x <= 0. {
                        return Err(Error::Domain("inertia must be SPD"));
                    }
                    chol[(i, j)] = x.sqrt();
                } else {
                    chol[(i, j)] = x / chol[(j, j)];
                }
            }
        }
        if 2. * largest_symmetric_eigenvalue(&a) > a.trace() + 1e-12 {
            return Err(Error::Domain(
                "inertia violates physical triangle condition",
            ));
        }
        let mut model = Self {
            mass,
            inertia,
            scale,
            chol,
            inverse: Mat3::zero(),
        };
        for j in 0..3 {
            let mut e = [0.; 3];
            e[j] = 1.;
            let col = model.solve_inertia(e)?;
            for (i, value) in col.into_iter().enumerate() {
                model.inverse[(i, j)] = value;
            }
        }
        Ok(model)
    }
    pub fn mass(&self) -> f64 {
        self.mass
    }
    pub fn inertia(&self) -> &Mat3 {
        &self.inertia
    }
    pub fn solve_inertia(&self, b: Vec3) -> Result<Vec3, Error> {
        finite(&b)?;
        let mut y = [0.; 3];
        let mut x = [0.; 3];
        for i in 0..3 {
            let mut v = b[i] / self.scale;
            for (j, value) in y.iter().enumerate().take(i) {
                v -= self.chol[(i, j)] * value;
            }
            y[i] = v / self.chol[(i, i)];
        }
        for i in (0..3).rev() {
            let mut v = y[i];
            for (j, value) in x.iter().enumerate().skip(i + 1) {
                v -= self.chol[(j, i)] * value;
            }
            x[i] = v / self.chol[(i, i)];
        }
        finite(&x)?;
        Ok(x)
    }
    pub fn momentum_body(&self, state: State) -> Result<Vec3, Error> {
        let s = validated(state)?;
        let h = self.inertia.mul_vector(&s.omega_body);
        finite(&h)?;
        Ok(h)
    }
    pub fn momentum_inertial(&self, state: State) -> Result<Vec3, Error> {
        let s = validated(state)?;
        rotation(s.quaternion)?.apply(&self.momentum_body(s)?)
    }
    pub fn kinetic_energy(&self, state: State) -> Result<f64, Error> {
        let s = validated(state)?;
        let e = 0.5
            * (self.mass * vector::dot(&s.velocity, &s.velocity)
                + vector::dot(&s.omega_body, &self.momentum_body(s)?));
        if !e.is_finite() {
            Err(Error::NonFinite)
        } else {
            Ok(e)
        }
    }
    pub fn rhs(&self, state: State, wrench: Wrench) -> Result<Derivative, Error> {
        let s = validated(state)?;
        validate_wrench(wrench)?;
        let force = add(
            rotation(s.quaternion)?.apply(&wrench.force_body)?,
            wrench.force_inertial,
        );
        let acceleration = force.map(|v| v / self.mass);
        finite(&acceleration)?;
        let mut dq = su2::hamilton(
            &s.quaternion,
            &[0., s.omega_body[0], s.omega_body[1], s.omega_body[2]],
        );
        for x in &mut dq {
            *x *= 0.5;
        }
        if dq.iter().any(|x| !x.is_finite()) {
            return Err(Error::NonFinite);
        }
        let angular = self.solve_inertia(sub(
            wrench.torque_body,
            vector::cross(&s.omega_body, &self.momentum_body(s)?),
        ))?;
        Ok(Derivative {
            position: s.velocity,
            velocity: acceleration,
            quaternion: dq,
            omega_body: angular,
        })
    }
    /// Implicit midpoint angular momentum + Cayley attitude increment + midpoint
    /// translation. Finite dt>0; no hidden cap. Failed Newton/overflow is an Error.
    pub fn step(
        &self,
        state: State,
        dt: f64,
        wrench: Wrench,
        options: StepOptions,
    ) -> Result<StepResult, Error> {
        let s = validated(state)?;
        validate_wrench(wrench)?;
        if !dt.is_finite() || dt <= 0. {
            return Err(Error::Domain("dt must be positive finite"));
        }
        if !options.tolerance.is_finite()
            || options.tolerance <= 0.
            || options.tolerance > 1e-4
            || options.max_iterations == 0
            || options.max_iterations > 100
        {
            return Err(Error::Domain("invalid Newton options"));
        }
        let h0 = self.momentum_body(s)?;
        let mut h1 = add(
            h0,
            scale(
                sub(wrench.torque_body, vector::cross(&s.omega_body, &h0)),
                dt,
            ),
        );
        let threshold =
            options.tolerance * (1. + inf_norm(&h0) + dt * inf_norm(&wrench.torque_body));
        if !threshold.is_finite() {
            return Err(Error::NonFinite);
        }
        let mut iteration = 0;
        let (wm, residual) = loop {
            finite(&h1)?;
            let hm = scale(add(h0, h1), 0.5);
            let wm = self.solve_inertia(hm)?;
            let f = sub(
                sub(h1, h0),
                scale(sub(wrench.torque_body, vector::cross(&wm, &hm)), dt),
            );
            finite(&f)?;
            let residual = inf_norm(&f);
            if residual <= threshold {
                break (wm, residual);
            }
            if iteration >= options.max_iterations {
                return Err(Error::Evaluation(
                    "implicit midpoint Newton did not converge".into(),
                ));
            }
            let jac = &Mat3::identity() + &(&hat(wm) - &(&hat(hm) * &self.inverse)).scaled(dt / 2.);
            h1 = sub(h1, solve3(jac, f)?);
            iteration += 1;
        };
        let increment = normalized([1., dt * wm[0] / 2., dt * wm[1] / 2., dt * wm[2] / 2.])?;
        let half_scalar = ((1. + increment[0]) / 2.).sqrt();
        let half = [
            half_scalar,
            increment[1] / (2. * half_scalar),
            increment[2] / (2. * half_scalar),
            increment[3] / (2. * half_scalar),
        ];
        let qmid = normalized(su2::hamilton(&s.quaternion, &half))?;
        let mut out = s;
        out.quaternion = normalized(su2::hamilton(&s.quaternion, &increment))?;
        out.omega_body = self.solve_inertia(h1)?;
        let force = add(
            rotation(qmid)?.apply(&wrench.force_body)?,
            wrench.force_inertial,
        );
        out.velocity = add(s.velocity, force.map(|v| (v / self.mass) * dt));
        out.position = add(s.position, scale(add(s.velocity, out.velocity), dt / 2.));
        Ok(StepResult {
            state: validated(out)?,
            iterations: iteration,
            residual,
        })
    }
}
