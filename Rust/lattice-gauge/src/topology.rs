//! Topological charge of SU(2) lattice fields in d = 4 (clover definition) and cooling.
//!
//! Clover C_{mu nu}(x): the four plaquette loops in the (mu, nu) plane that start and end at x, all with the same
//! orientation. For SU(2) its antihermitian part is i f.s with f = (1/4) vec(C) per leaf average, and
//! q(x) = (1/(2 pi^2)) (f_01.f_23 - f_02.f_13 + f_03.f_12), the lattice form of (1/(32 pi^2)) eps tr(F F) up to sign
//! (the sign does not affect chi_t = <Q^2>/V). After cooling, Q clusters near integers; that is the test of the
//! normalization. Reference values: Teper, hep-th/9812187, Table 15 (cooled, 15-20 sweeps).

use crate::su2::Quaternion;
use crate::wilson::{Lattice, Update};

fn clover(l: &Lattice, x: usize, mu: usize, nu: usize) -> Quaternion {
    let t = &l.torus;
    let u = |s: usize, m: usize| l.links[t.link_index(s, m)];
    let f = |s: usize, m: usize| t.neighbor(s, m, true);
    let b = |s: usize, m: usize| t.neighbor(s, m, false);
    let (xpm, xpn, xmm, xmn) = (f(x, mu), f(x, nu), b(x, mu), b(x, nu));
    let leaf1 = u(x, mu).mul(u(xpm, nu)).mul(u(xpn, mu).conj()).mul(u(x, nu).conj());
    let leaf2 = u(x, nu).mul(u(b(xpn, mu), mu).conj()).mul(u(xmm, nu).conj()).mul(u(xmm, mu));
    let leaf3 = u(xmm, mu).conj().mul(u(b(xmm, nu), nu).conj()).mul(u(b(xmm, nu), mu)).mul(u(xmn, nu));
    let leaf4 = u(xmn, nu).conj().mul(u(xmn, mu)).mul(u(f(xmn, mu), nu)).mul(u(x, mu).conj());
    leaf1.add(leaf2).add(leaf3).add(leaf4)
}

fn vec4(q: Quaternion) -> [f64; 3] {
    [q.a[1] / 4.0, q.a[2] / 4.0, q.a[3] / 4.0]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Topological charge density at x (d = 4).
pub fn charge_density(l: &Lattice, x: usize) -> f64 {
    let f01 = vec4(clover(l, x, 0, 1));
    let f02 = vec4(clover(l, x, 0, 2));
    let f03 = vec4(clover(l, x, 0, 3));
    let f12 = vec4(clover(l, x, 1, 2));
    let f13 = vec4(clover(l, x, 1, 3));
    let f23 = vec4(clover(l, x, 2, 3));
    (dot(f01, f23) - dot(f02, f13) + dot(f03, f12)) / (2.0 * std::f64::consts::PI * std::f64::consts::PI)
}

/// Total clover topological charge Q (d = 4).
pub fn charge(l: &Lattice) -> f64 {
    assert_eq!(l.torus.dimension(), 4, "topological charge needs d = 4");
    (0..l.torus.site_count()).map(|x| charge_density(l, x)).sum()
}

/// `n` cooling sweeps.
pub fn cool(l: &mut Lattice, n: usize) {
    for _ in 0..n {
        l.sweep_with(Update::Cool);
    }
}
