//! Centre vortices for SU(2): direct maximal centre gauge (DMCG), centre projection and P-vortices.
//!
//! DMCG maximizes R = sum_{x,mu} (tr U)^2 / 4 = sum a0^2 over gauge transformations; since (tr U)^2 = tr_adj U + 1 this
//! is maximizing sum tr O(U), O the 3x3 adjoint (rotation) matrix. At site x, with neighbours fixed, the relevant part is
//! tr(O(g) M), M = sum_mu O(U(x,mu)) + sum_mu O(U(x-mu,mu))^T, maximized over SO(3) by the polar factor of M with a
//! determinant correction (from the singular value decomposition, computed by Jacobi rotations). Centre projection:
//! Z(x,mu) = sign(a0); a plaquette with Z-product -1 is pierced by a P-vortex.

use crate::su2::Quaternion;
use crate::wilson::Lattice;

type M3 = [[f64; 3]; 3];

/// Adjoint matrix O(g): g (v.s) g^dagger = (O v).s
pub fn adjoint_matrix(g: Quaternion) -> M3 {
    let mut o = [[0.0; 3]; 3];
    for j in 0..3 {
        let mut e = [0.0; 3];
        e[j] = 1.0;
        let q = g.mul(Quaternion::new(0.0, e[0], e[1], e[2])).mul(g.conj());
        for i in 0..3 {
            o[i][j] = q.a[i + 1];
        }
    }
    o
}

fn matmul(a: &M3, b: &M3) -> M3 {
    let mut c = [[0.0; 3]; 3];
    for i in 0..3 { for j in 0..3 { for k in 0..3 { c[i][j] += a[i][k] * b[k][j]; } } }
    c
}

fn transpose(a: &M3) -> M3 {
    let mut t = [[0.0; 3]; 3];
    for i in 0..3 { for j in 0..3 { t[i][j] = a[j][i]; } }
    t
}

#[allow(dead_code)]
fn det(a: &M3) -> f64 {
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1]) - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
}

/// Symmetric eigen-decomposition by cyclic Jacobi rotations: returns (eigenvalues, V) with S = V diag V^T.
fn jacobi_eigen(mut s: M3) -> ([f64; 3], M3) {
    let mut v = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    for _ in 0..50 {
        let off = s[0][1].abs() + s[0][2].abs() + s[1][2].abs();
        if off < 1e-15 * (s[0][0].abs() + s[1][1].abs() + s[2][2].abs() + 1e-300) {
            break;
        }
        for (p, q) in [(0, 1), (0, 2), (1, 2)] {
            if s[p][q].abs() < 1e-300 { continue; }
            let theta = (s[q][q] - s[p][p]) / (2.0 * s[p][q]);
            let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
            let t = if theta == 0.0 { 1.0 } else { t };
            let c = 1.0 / (t * t + 1.0).sqrt();
            let sn = t * c;
            let mut r = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
            r[p][p] = c; r[q][q] = c; r[p][q] = sn; r[q][p] = -sn;
            s = matmul(&transpose(&r), &matmul(&s, &r));
            v = matmul(&v, &r);
        }
    }
    ([s[0][0], s[1][1], s[2][2]], v)
}

/// The rotation R in SO(3) maximizing tr(R M). With M = U diag(s1, s2, s3) V^T, U and V proper rotations,
/// s1 >= s2 >= |s3| (s3 signed), the maximum is s1 + s2 + s3 at R = V U^T.
pub fn best_rotation(m: &M3) -> M3 {
    let (ev, v) = jacobi_eigen(matmul(&transpose(m), m)); // M^T M = V S^2 V^T
    let mut order = [0usize, 1, 2];
    order.sort_by(|&a, &b| ev[b].partial_cmp(&ev[a]).unwrap());
    let vc = |i: usize| [v[0][i], v[1][i], v[2][i]];
    let (v0, v1) = (vc(order[0]), vc(order[1]));
    let v2 = cross(v0, v1); // proper V
    let apply = |x: [f64; 3]| [dot3(m[0], x), dot3(m[1], x), dot3(m[2], x)];
    let u0 = normalize(apply(v0));
    let w = apply(v1);
    let u1 = normalize(sub(w, scale(u0, dot3(u0, w))));
    let u2 = cross(u0, u1);
    // R = sum_k v_k u_k^T
    let mut r = [[0.0; 3]; 3];
    for (vk, uk) in [(v0, u0), (v1, u1), (v2, u2)] {
        for i in 0..3 { for j in 0..3 { r[i][j] += vk[i] * uk[j]; } }
    }
    r
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 { a[0] * b[0] + a[1] * b[1] + a[2] * b[2] }
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] { [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]] }
fn scale(a: [f64; 3], s: f64) -> [f64; 3] { [a[0] * s, a[1] * s, a[2] * s] }
fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] { [a[0] - b[0], a[1] - b[1], a[2] - b[2]] }
fn normalize(a: [f64; 3]) -> [f64; 3] {
    let n = dot3(a, a).sqrt();
    if n < 1e-300 { [1.0, 0.0, 0.0] } else { scale(a, 1.0 / n) }
}

/// Unit quaternion with adjoint matrix R (one of the two preimages).
pub fn quaternion_from_rotation(r: &M3) -> Quaternion {
    // Our adjoint convention may be R or R^T of the standard formula; fix by checking.
    let tr = r[0][0] + r[1][1] + r[2][2];
    let w = ((1.0 + tr).max(0.0)).sqrt() / 2.0;
    let candidates = if w > 1e-6 {
        let x = (r[2][1] - r[1][2]) / (4.0 * w);
        let y = (r[0][2] - r[2][0]) / (4.0 * w);
        let z = (r[1][0] - r[0][1]) / (4.0 * w);
        vec![Quaternion::new(w, x, y, z), Quaternion::new(w, -x, -y, -z)]
    } else {
        // angle near pi: pick the axis from the largest diagonal entry
        let i = (0..3).max_by(|&a, &b| r[a][a].partial_cmp(&r[b][b]).unwrap()).unwrap();
        let (j, k) = ((i + 1) % 3, (i + 2) % 3);
        let s = ((r[i][i] - r[j][j] - r[k][k] + 1.0).max(0.0)).sqrt() / 2.0;
        let mut a = [0.0; 3];
        a[i] = s;
        a[j] = (r[j][i] + r[i][j]) / (4.0 * s);
        a[k] = (r[k][i] + r[i][k]) / (4.0 * s);
        let w2 = (r[k][j] - r[j][k]) / (4.0 * s);
        vec![Quaternion::new(w2, a[0], a[1], a[2]), Quaternion::new(-w2, a[0], a[1], a[2]), Quaternion::new(w2, -a[0], -a[1], -a[2])]
    };
    candidates
        .into_iter()
        .min_by(|p, q| dist(&adjoint_matrix(p.normalized()), r).partial_cmp(&dist(&adjoint_matrix(q.normalized()), r)).unwrap())
        .unwrap()
        .normalized()
}

fn dist(a: &M3, b: &M3) -> f64 {
    let mut s = 0.0;
    for i in 0..3 { for j in 0..3 { s += (a[i][j] - b[i][j]).powi(2); } }
    s
}

/// DMCG functional: mean a0^2 per link (1/4 at random, 1 when every link is +-1).
pub fn mcg_functional(l: &Lattice) -> f64 {
    l.links.iter().map(|q| q.a[0] * q.a[0]).sum::<f64>() / l.links.len() as f64
}

/// Iterate DMCG over checkerboard site parities until the functional changes by less than `tol` per sweep.
pub fn fix_mcg(l: &mut Lattice, tol: f64, max_sweeps: usize) -> (f64, usize) {
    let d = l.torus.dimension();
    let mut last = mcg_functional(l);
    for sweep in 1..=max_sweeps {
        for parity in 0..2 {
            let t = &l.torus;
            let mut g = vec![Quaternion::IDENTITY; t.site_count()];
            for x in 0..t.site_count() {
                if t.site_coordinates(x).iter().sum::<usize>() % 2 != parity {
                    continue;
                }
                let mut m = [[0.0; 3]; 3];
                for mu in 0..d {
                    let a = adjoint_matrix(l.links[t.link_index(x, mu)]);
                    let b = transpose(&adjoint_matrix(l.links[t.link_index(t.neighbor(x, mu, false), mu)]));
                    for i in 0..3 { for j in 0..3 { m[i][j] += a[i][j] + b[i][j]; } }
                }
                g[x] = quaternion_from_rotation(&best_rotation(&m));
            }
            l.gauge_transform(&g);
        }
        let r = mcg_functional(l);
        if (r - last).abs() < tol {
            return (r, sweep);
        }
        last = r;
    }
    (last, max_sweeps)
}

/// Centre projection Z = sign(a0) per link.
pub fn centre_projection(l: &Lattice) -> Vec<i8> {
    l.links.iter().map(|q| if q.a[0] >= 0.0 { 1 } else { -1 }).collect()
}

/// P-vortex density: fraction of plaquettes whose centre-projected product is -1.
pub fn p_vortex_density(l: &Lattice, z: &[i8]) -> f64 {
    let t = &l.torus;
    let d = t.dimension();
    let mut pierced = 0usize;
    for x in 0..t.site_count() {
        for mu in 0..d {
            for nu in mu + 1..d {
                let p = z[t.link_index(x, mu)] * z[t.link_index(t.neighbor(x, mu, true), nu)] * z[t.link_index(t.neighbor(x, nu, true), mu)] * z[t.link_index(x, nu)];
                if p < 0 { pierced += 1; }
            }
        }
    }
    pierced as f64 / t.plaquette_count() as f64
}
