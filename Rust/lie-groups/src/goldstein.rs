//! Goldstein §§4.2–4.7, in column-vector, active right-handed notation.
//! R maps body components to inertial components; Goldstein's A_G = R^T.
//! Forward maps are scalar generic. Branch-sensitive inverse maps use f64.
use crate::{default_tolerance, su2, Complex, Error, Matrix, Real, Scalar, So3, So3Algebra, Su2};

/// Axis rotation from supplied cosine and sine; exact rings are supported.
/// The caller supplies c²+s²=1; this returns a matrix, not a checked group value.
pub fn axis_rotation<T: Scalar>(axis: usize, c: T, s: T) -> Result<Matrix<T, 3, 3>, Error> {
    let z = T::zero();
    let o = T::one();
    Ok(match axis {
        0 => Matrix::from_rows([
            [o, z.clone(), z.clone()],
            [z.clone(), c.clone(), -s.clone()],
            [z, s, c],
        ]),
        1 => Matrix::from_rows([
            [c.clone(), z.clone(), s.clone()],
            [z.clone(), o, z.clone()],
            [-s, z, c],
        ]),
        2 => Matrix::from_rows([
            [c.clone(), -s.clone(), z.clone()],
            [s, c, z.clone()],
            [z.clone(), z, o],
        ]),
        _ => {
            return Err(Error::IndexOutOfRange {
                index: axis,
                dimension: 3,
            })
        }
    })
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum EulerConvention {
    Zxz,
    Zyz,
    Zyx,
}
/// Explicit products: ZXZ/ZYZ use (phi,theta,psi); ZYX uses (yaw psi,pitch
/// theta,roll phi), R=Rz(psi)Ry(theta)Rx(phi). These are intrinsic moving-axis
/// rotations, equivalently the reverse-order extrinsic fixed-axis sequence.
pub fn euler_matrix<T: Scalar>(
    convention: EulerConvention,
    cosines: [T; 3],
    sines: [T; 3],
) -> Matrix<T, 3, 3> {
    let [c0, c1, c2] = cosines;
    let [s0, s1, s2] = sines;
    let (a, b, c) = match convention {
        EulerConvention::Zxz => (
            axis_rotation(2, c0, s0),
            axis_rotation(0, c1, s1),
            axis_rotation(2, c2, s2),
        ),
        EulerConvention::Zyz => (
            axis_rotation(2, c0, s0),
            axis_rotation(1, c1, s1),
            axis_rotation(2, c2, s2),
        ),
        EulerConvention::Zyx => (
            axis_rotation(2, c0, s0),
            axis_rotation(1, c1, s1),
            axis_rotation(0, c2, s2),
        ),
    };
    &(&a.expect("fixed axis") * &b.expect("fixed axis")) * &c.expect("fixed axis")
}
pub fn from_euler<T: Real>(convention: EulerConvention, angles: [T; 3]) -> Result<So3<T>, Error> {
    if angles.iter().any(|a| !a.is_finite()) {
        return Err(Error::NonFinite);
    }
    So3::try_from_matrix(
        euler_matrix(convention, angles.map(Real::cos), angles.map(Real::sin)),
        default_tolerance(),
    )
}
/// G'=B A B^T F' for orthogonal B; no orthogonality claim is made for generic
/// input matrices. For B in SO(3), Ad_B(hat(w))=hat(Bw).
pub fn orthogonal_conjugate<T: Scalar>(
    b: &Matrix<T, 3, 3>,
    a: &Matrix<T, 3, 3>,
) -> Matrix<T, 3, 3> {
    &(b * a) * &b.transpose()
}
/// Characteristic-polynomial residual for ANY 3x3 matrix, with no division:
/// A³-tr(A)A²+s₂(A)A-det(A)I, s₂=sum of principal 2x2 minors.
pub fn cayley_hamilton_residual<T: Scalar>(a: &Matrix<T, 3, 3>) -> Matrix<T, 3, 3> {
    let a2 = a * a;
    let mut s2 = T::zero();
    for i in 0..3 {
        for j in (i + 1)..3 {
            s2 = s2 + a[(i, i)].clone() * a[(j, j)].clone() - a[(i, j)].clone() * a[(j, i)].clone();
        }
    }
    &(&(&(a * &a2) - &a2.scaled(a.trace())) + &a.scaled(s2))
        - &Matrix::identity().scaled(a.determinant())
}
/// Goldstein K=[[alpha,beta],[gamma,delta]]=U(q)^dagger. Thus
/// alpha=w+iz, beta=y+ix, gamma=-beta*, delta=alpha*, and |alpha|²+|beta|²=1.
/// This is an SU(2) lift, retaining the distinction q versus -q.
/// Raw exact Cayley–Klein dictionary from scalar-first Hamilton coordinates.
/// Group membership requires sum(q_i²)=1; this polynomial map does not check it.
pub fn cayley_klein_parameters<T: Scalar>(q: &[T; 4]) -> [Complex<T>; 4] {
    let [w, x, y, z] = q.clone();
    [
        Complex::new(w.clone(), z.clone()),
        Complex::new(y.clone(), x.clone()),
        Complex::new(-y, x),
        Complex::new(w, -z),
    ]
}
/// K=U(q)^dagger for real coefficient symbols; no conjugation of formal
/// variables is needed. Membership is subject to the same unit-norm condition.
pub fn cayley_klein_matrix<T: Scalar>(q: &[T; 4]) -> Matrix<Complex<T>, 2, 2> {
    let [a, b, c, d] = cayley_klein_parameters(q);
    Matrix::from_rows([[a, b], [c, d]])
}
/// Checked unit SU(2) lift in Goldstein's Cayley–Klein notation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CayleyKlein<T> {
    lift: Su2<T>,
}
impl<T: Real> CayleyKlein<T> {
    pub fn from_lift(lift: Su2<T>) -> Self {
        Self { lift }
    }
    pub fn lift(&self) -> &Su2<T> {
        &self.lift
    }
    pub fn parameters(&self) -> [Complex<T>; 4] {
        cayley_klein_parameters(self.lift.quaternion())
    }
    pub fn matrix(&self) -> Matrix<Complex<T>, 2, 2> {
        self.lift.matrix().adjoint()
    }
    pub fn try_from_parameters(p: [Complex<T>; 4], tol: T) -> Result<Self, Error> {
        let [a, b, c, d] = p;
        let k = Matrix::from_rows([[a, b], [c, d]]);
        Ok(Self::from_lift(Su2::try_from_matrix(k.adjoint(), tol)?))
    }
    pub fn goldstein_matrix(&self) -> Result<Matrix<T, 3, 3>, Error> {
        Ok(self.lift.to_so3()?.matrix().transpose())
    }
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EulerAngles {
    pub angles: [f64; 3],
    pub singular: bool,
}
/// Canonical principal extraction. Proper Euler: theta in [0,pi]. ZYX:
/// theta in [-pi/2,pi/2]. At a pole the FINAL product factor's angle is fixed
/// to zero (psi for proper Euler, roll phi for ZYX). Singular is a numerical
/// chart classification, |sin(theta)| or |cos(theta)| <= 64 epsilon.
pub fn extract_euler(r: &So3<f64>, convention: EulerConvention) -> EulerAngles {
    let m = r.matrix();
    let tol = 64.0 * f64::EPSILON;
    let (angles, singular) = match convention {
        EulerConvention::Zxz | EulerConvention::Zyz => {
            let sin_theta = m[(0, 2)].hypot(m[(1, 2)]);
            let theta = sin_theta.atan2(m[(2, 2)]);
            if sin_theta > tol {
                let (phi, psi) = if convention == EulerConvention::Zxz {
                    (m[(0, 2)].atan2(-m[(1, 2)]), m[(2, 0)].atan2(m[(2, 1)]))
                } else {
                    (m[(1, 2)].atan2(m[(0, 2)]), m[(2, 1)].atan2(-m[(2, 0)]))
                };
                ([phi, theta, psi], false)
            } else {
                let sign = if convention == EulerConvention::Zyz && m[(2, 2)] < 0.0 {
                    -1.0
                } else {
                    1.0
                };
                (
                    [(sign * m[(1, 0)]).atan2(sign * m[(0, 0)]), theta, 0.0],
                    true,
                )
            }
        }
        EulerConvention::Zyx => {
            let cos_theta = m[(0, 0)].hypot(m[(1, 0)]);
            let theta = (-m[(2, 0)]).atan2(cos_theta);
            if cos_theta > tol {
                (
                    [
                        m[(1, 0)].atan2(m[(0, 0)]),
                        theta,
                        m[(2, 1)].atan2(m[(2, 2)]),
                    ],
                    false,
                )
            } else {
                ([(-m[(0, 1)]).atan2(m[(1, 1)]), theta, 0.0], true)
            }
        }
    };
    EulerAngles { angles, singular }
}
/// One principal lift of R (w >= 0). At pi, choose the largest diagonal pivot;
/// the sign of the axis is intrinsically ambiguous. Imported lifts themselves
/// are never canonicalized by this API. No normalization is performed: a matrix
/// imported at loose tolerance may yield a lift rejected at default tolerance.
pub fn lift_rotation(r: &So3<f64>) -> Result<Su2<f64>, Error> {
    let m = r.matrix();
    let t = m.trace();
    let mut q = [0.0; 4];
    if t > 0.0 {
        let s = 2.0 * (1.0 + t).max(0.0).sqrt();
        q = [
            s / 4.0,
            (m[(2, 1)] - m[(1, 2)]) / s,
            (m[(0, 2)] - m[(2, 0)]) / s,
            (m[(1, 0)] - m[(0, 1)]) / s,
        ];
    } else {
        let i = if m[(0, 0)] >= m[(1, 1)] && m[(0, 0)] >= m[(2, 2)] {
            0
        } else if m[(1, 1)] >= m[(2, 2)] {
            1
        } else {
            2
        };
        let j = (i + 1) % 3;
        let k = (i + 2) % 3;
        let s = 2.0 * (1.0 + m[(i, i)] - m[(j, j)] - m[(k, k)]).max(0.0).sqrt();
        q[0] = (m[(k, j)] - m[(j, k)]) / s;
        q[i + 1] = s / 4.0;
        q[j + 1] = (m[(j, i)] + m[(i, j)]) / s;
        q[k + 1] = (m[(k, i)] + m[(i, k)]) / s;
    }
    if q[0] < 0.0 {
        q = q.map(|v| -v)
    }
    Su2::try_from_unit_quaternion(q, default_tolerance())
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AxisAngle {
    pub axis: [f64; 3],
    pub angle: f64,
}
/// Principal angle in [0,pi], with axis [1,0,0] at identity. atan2 of the
/// quaternion vector norm avoids acos(trace) cancellation near zero and pi.
pub fn axis_angle(r: &So3<f64>) -> Result<AxisAngle, Error> {
    let q = *lift_rotation(r)?.quaternion();
    let n = q[1].hypot(q[2]).hypot(q[3]);
    Ok(AxisAngle {
        axis: if n == 0.0 {
            [1.0, 0.0, 0.0]
        } else {
            [q[1] / n, q[2] / n, q[3] / n]
        },
        angle: 2.0 * n.atan2(q[0]),
    })
}
/// Exact Rodrigues polynomial I+s hat(n)+(1-c)hat(n)^2. Caller supplies
/// unit n, c=cos(phi), s=sin(phi); Goldstein Eq.4.62 uses its transpose.
pub fn rodrigues<T: Scalar>(axis: [T; 3], c: T, s: T) -> Matrix<T, 3, 3> {
    let k = So3Algebra::new(axis).matrix();
    &(&Matrix::identity() + &k.scaled(s)) + &(&k * &k).scaled(T::one() - c)
}
/// Euler products lifted before projecting, so 2pi changes may change the lift.
pub fn euler_lift<T: Real>(convention: EulerConvention, angles: [T; 3]) -> Result<Su2<T>, Error> {
    if angles.iter().any(|x| !x.is_finite()) {
        return Err(Error::NonFinite);
    }
    let unit = |axis: usize, a: T| {
        let mut q = [
            (a * T::from_f64(0.5)).cos(),
            T::zero(),
            T::zero(),
            T::zero(),
        ];
        q[axis + 1] = (a * T::from_f64(0.5)).sin();
        q
    };
    let [p, t, s] = angles;
    let (a, b, c) = match convention {
        EulerConvention::Zxz => (unit(2, p), unit(0, t), unit(2, s)),
        EulerConvention::Zyz => (unit(2, p), unit(1, t), unit(2, s)),
        EulerConvention::Zyx => (unit(2, p), unit(1, t), unit(0, s)),
    };
    Su2::try_from_unit_quaternion(
        su2::hamilton(&su2::hamilton(&a, &b), &c),
        default_tolerance(),
    )
}
