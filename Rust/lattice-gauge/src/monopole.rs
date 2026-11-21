//! Monopoles: maximal abelian gauge (MAG) and abelian projection for SU(2), and DeGrand-Toussaint monopole
//! currents for any U(1) link angles (abelian-projected SU(2), or the compact U(1) counter-model directly).
//!
//! MAG maximizes R = sum_{x,mu} (1/2) tr(s3 U s3 U^dagger) = sum (a0^2 + a3^2 - a1^2 - a2^2) over gauge
//! transformations. At site x the relevant part is (1/2) tr(s3 g X g^dagger) with the traceless Hermitian
//! X = sum_mu [U(x,mu) s3 U(x,mu)^dagger + U(x-mu,mu)^dagger s3 U(x-mu,mu)] = x.s, maximized by the g that rotates
//! x onto +z; over-relaxation applies g^omega. Abelian projection: with U = [[a0 + i a3, ...], ...], theta = arg(a0 + i a3).
//! DeGrand-Toussaint: theta_p = [theta_p] + 2 pi n_p with [theta_p] in (-pi, pi]; the monopole charge of a 3-cube is
//! m = -sum over its six faces (boundary orientation) of n_p, an integer, and sum of |m| over cubes / cube count is the
//! monopole density. Currents are conserved (zero divergence), which the tests check.

use crate::su2::Quaternion;
use crate::torus::Torus;
use crate::wilson::Lattice;
use std::f64::consts::PI;

/// Rotation action of a unit quaternion on 3-vectors: g (v.s) g^dagger = (R v).s, computed directly.
fn adjoint(g: Quaternion, v: [f64; 3]) -> [f64; 3] {
    let q = g.mul(Quaternion::new(0.0, v[0], v[1], v[2])).mul(g.conj());
    [q.a[1], q.a[2], q.a[3]]
}

/// Unit quaternion g with adjoint(g, v) along +z (v nonzero), raised to the power omega (same axis, angle x omega).
fn rotate_to_z(v: [f64; 3], omega: f64) -> Quaternion {
    let n = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
    if n < 1e-300 {
        return Quaternion::IDENTITY;
    }
    let u = [v[0] / n, v[1] / n, v[2] / n];
    let an = (u[0] * u[0] + u[1] * u[1]).sqrt();
    let axis = if an < 1e-14 { [1.0, 0.0, 0.0] } else { [u[1] / an, -u[0] / an, 0.0] };
    let half = u[2].clamp(-1.0, 1.0).acos() / 2.0;
    if half < 1e-15 {
        return Quaternion::IDENTITY;
    }
    // the sign of the half-angle depends on the product convention; take the one that maps u onto +z
    // (the wrong sign rotates by twice the angle, which can still leave a positive z component, so compare both)
    let sgn = {
        let gp = Quaternion::new(half.cos(), half.sin() * axis[0], half.sin() * axis[1], half.sin() * axis[2]);
        let gm = Quaternion::new(half.cos(), -half.sin() * axis[0], -half.sin() * axis[1], -half.sin() * axis[2]);
        if adjoint(gp, u)[2] >= adjoint(gm, u)[2] { 1.0 } else { -1.0 }
    };
    let h = sgn * half * omega;
    Quaternion::new(h.cos(), h.sin() * axis[0], h.sin() * axis[1], h.sin() * axis[2])
}

/// MAG functional per link: mean of a0^2 + a3^2 - a1^2 - a2^2 (1 when every link is diagonal).
pub fn mag_functional(l: &Lattice) -> f64 {
    l.links.iter().map(|q| q.a[0] * q.a[0] + q.a[3] * q.a[3] - q.a[1] * q.a[1] - q.a[2] * q.a[2]).sum::<f64>() / l.links.len() as f64
}

/// Iterate MAG over checkerboard site parities until the functional changes by less than `tol` per sweep.
/// Returns (functional, sweeps). The links are gauge transformed in place (all observables unchanged).
pub fn fix_mag(l: &mut Lattice, omega: f64, tol: f64, max_sweeps: usize) -> (f64, usize) {
    let d = l.torus.dimension();
    let mut last = mag_functional(l);
    let s3 = [0.0, 0.0, 1.0];
    for sweep in 1..=max_sweeps {
        for parity in 0..2 {
            let t = &l.torus;
            let mut g = vec![Quaternion::IDENTITY; t.site_count()];
            for x in 0..t.site_count() {
                if t.site_coordinates(x).iter().sum::<usize>() % 2 != parity {
                    continue;
                }
                let mut v = [0.0; 3];
                for mu in 0..d {
                    let a = adjoint(l.links[t.link_index(x, mu)], s3);
                    let xm = t.neighbor(x, mu, false);
                    let b = adjoint(l.links[t.link_index(xm, mu)].conj(), s3);
                    for k in 0..3 {
                        v[k] += a[k] + b[k];
                    }
                }
                g[x] = rotate_to_z(v, omega);
            }
            l.gauge_transform(&g);
        }
        let r = mag_functional(l);
        if (r - last).abs() < tol {
            return (r, sweep);
        }
        last = r;
    }
    (last, max_sweeps)
}

/// Abelian projection: theta(x,mu) = arg(a0 + i a3) for every link.
pub fn abelian_angles(l: &Lattice) -> Vec<f64> {
    l.links.iter().map(|q| q.a[3].atan2(q.a[0])).collect()
}

fn principal(x: f64) -> (f64, i64) {
    // x = p + 2 pi n with p in (-pi, pi]
    let n = ((x + PI) / (2.0 * PI)).ceil() as i64 - 1;
    let p = x - 2.0 * PI * n as f64;
    if p <= -PI { (p + 2.0 * PI, n - 1) } else { (p, n) }
}

/// Integer n_p of every plaquette (index Torus::plaquette_index) for U(1) link angles.
pub fn plaquette_integers(t: &Torus, theta: &[f64]) -> Vec<i64> {
    let d = t.dimension();
    let mut n = vec![0i64; t.plaquette_count()];
    for x in 0..t.site_count() {
        for mu in 0..d {
            for nu in mu + 1..d {
                let th = theta[t.link_index(x, mu)] + theta[t.link_index(t.neighbor(x, mu, true), nu)]
                    - theta[t.link_index(t.neighbor(x, nu, true), mu)] - theta[t.link_index(x, nu)];
                n[t.plaquette_index(x, mu, nu)] = principal(th).1;
            }
        }
    }
    n
}

/// Monopole charge of the 3-cube at x spanned by a < b < c: minus the boundary sum of plaquette integers,
/// boundary = +P_bc(x+a) - P_bc(x) - P_ac(x+b) + P_ac(x) + P_ab(x+c) - P_ab(x).
pub fn cube_charge(t: &Torus, n: &[i64], x: usize, a: usize, b: usize, c: usize) -> i64 {
    let p = |s: usize, i: usize, j: usize| n[t.plaquette_index(s, i, j)];
    -(p(t.neighbor(x, a, true), b, c) - p(x, b, c) - p(t.neighbor(x, b, true), a, c) + p(x, a, c) + p(t.neighbor(x, c, true), a, b) - p(x, a, b))
}

/// Monopole density: mean |m| over all 3-cubes (d = 4: four per site).
pub fn monopole_density(t: &Torus, theta: &[f64]) -> f64 {
    let n = plaquette_integers(t, theta);
    let d = t.dimension();
    let (mut s, mut count) = (0i64, 0usize);
    for x in 0..t.site_count() {
        for a in 0..d {
            for b in a + 1..d {
                for c in b + 1..d {
                    s += cube_charge(t, &n, x, a, b, c).abs();
                    count += 1;
                }
            }
        }
    }
    s as f64 / count as f64
}
