//! SU(2) as unit quaternions: q = (a0, a1, a2, a3) represents a0 1 + i (a1 s1 + a2 s2 + a3 s3)
//! with Pauli matrices s_k, so (1/2) tr U = a0 and the matrix product is the quaternion product
//! below. Sampling uses `philox::Stream`.

use crate::philox::Stream;

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Quaternion {
    pub a: [f64; 4],
}

impl Quaternion {
    pub const IDENTITY: Quaternion = Quaternion { a: [1.0, 0.0, 0.0, 0.0] };
    pub const ZERO: Quaternion = Quaternion { a: [0.0; 4] };

    pub fn new(a0: f64, a1: f64, a2: f64, a3: f64) -> Self {
        Quaternion { a: [a0, a1, a2, a3] }
    }

    /// Matrix product U V. With U = a0 + i a.s and V = b0 + i b.s:
    /// U V = (a0 b0 - a.b) + i (a0 b + b0 a - a x b).s
    pub fn mul(self, o: Quaternion) -> Quaternion {
        let [a0, a1, a2, a3] = self.a;
        let [b0, b1, b2, b3] = o.a;
        Quaternion::new(
            a0 * b0 - a1 * b1 - a2 * b2 - a3 * b3,
            a0 * b1 + b0 * a1 - (a2 * b3 - a3 * b2),
            a0 * b2 + b0 * a2 - (a3 * b1 - a1 * b3),
            a0 * b3 + b0 * a3 - (a1 * b2 - a2 * b1),
        )
    }

    /// Hermitian conjugate (the inverse for a unit quaternion).
    pub fn conj(self) -> Quaternion {
        Quaternion::new(self.a[0], -self.a[1], -self.a[2], -self.a[3])
    }

    pub fn add(self, o: Quaternion) -> Quaternion {
        Quaternion::new(self.a[0] + o.a[0], self.a[1] + o.a[1], self.a[2] + o.a[2], self.a[3] + o.a[3])
    }

    pub fn scale(self, s: f64) -> Quaternion {
        Quaternion::new(self.a[0] * s, self.a[1] * s, self.a[2] * s, self.a[3] * s)
    }

    /// det of the 2x2 matrix; for sums of SU(2) elements this is |W|^2 and W = |W| V with V in SU(2).
    pub fn norm_sq(self) -> f64 {
        self.a.iter().map(|x| x * x).sum()
    }

    pub fn half_trace(self) -> f64 {
        self.a[0]
    }

    pub fn normalized(self) -> Quaternion {
        self.scale(1.0 / self.norm_sq().sqrt())
    }
}

/// Haar-random SU(2): a uniform point on S^3.
pub fn haar(stream: &mut Stream) -> Quaternion {
    loop {
        let q = Quaternion::new(stream.next_gaussian(), stream.next_gaussian(), stream.next_gaussian(), stream.next_gaussian());
        let n = q.norm_sq();
        if n > 1e-300 {
            return q.scale(1.0 / n.sqrt());
        }
    }
}

/// Uniform direction on S^2 scaled to radius r.
fn sphere2(stream: &mut Stream, r: f64) -> [f64; 3] {
    let cos_t = 2.0 * stream.next_f64_open01() - 1.0;
    let sin_t = (1.0 - cos_t * cos_t).max(0.0).sqrt();
    let phi = 2.0 * std::f64::consts::PI * stream.next_f64_open01();
    [r * sin_t * phi.cos(), r * sin_t * phi.sin(), r * cos_t]
}

fn heat_bath_impl(staple: Quaternion, beta: f64, stream: &mut Stream, accept_step: bool) -> Quaternion {
    let k = staple.norm_sq().sqrt();
    let a = beta * k;
    if a < 1e-12 {
        return haar(stream);
    }
    let v = staple.scale(1.0 / k);
    // Kennedy-Pendleton: x0 with density sqrt(1 - x0^2) exp(a x0) on [-1, 1].
    let x0 = loop {
        let r1 = stream.next_f64_open01();
        let r2 = stream.next_f64_open01();
        let r3 = stream.next_f64_open01();
        let c = (2.0 * std::f64::consts::PI * r2).cos();
        let lambda2 = -(r1.ln() + c * c * r3.ln()) / (2.0 * a);
        let r4 = stream.next_f64_open01();
        if !accept_step || r4 * r4 <= 1.0 - lambda2 {
            break (1.0 - 2.0 * lambda2).max(-1.0);
        }
    };
    let [x1, x2, x3] = sphere2(stream, (1.0 - x0 * x0).max(0.0).sqrt());
    // The local weight is exp(beta a0(U W)) = exp(a a0(U V)); with X = U V distributed as above, U = X V^dagger.
    Quaternion::new(x0, x1, x2, x3).mul(v.conj())
}

/// U distributed as exp(beta a0(U W)) dU (Haar), W the staple sum.
pub fn heat_bath(staple: Quaternion, beta: f64, stream: &mut Stream) -> Quaternion {
    heat_bath_impl(staple, beta, stream, true)
}

/// Negative control: the Kennedy-Pendleton proposal without its accept step. Wrong distribution.
pub fn heat_bath_broken(staple: Quaternion, beta: f64, stream: &mut Stream) -> Quaternion {
    heat_bath_impl(staple, beta, stream, false)
}

/// Over-relaxation U -> V^dagger U^dagger V^dagger with V = W/|W|; keeps a0(U W) and is an involution.
pub fn overrelax(u: Quaternion, staple: Quaternion) -> Quaternion {
    let n = staple.norm_sq();
    if n < 1e-300 {
        return u;
    }
    let vd = staple.scale(1.0 / n.sqrt()).conj();
    vd.mul(u.conj()).mul(vd)
}
