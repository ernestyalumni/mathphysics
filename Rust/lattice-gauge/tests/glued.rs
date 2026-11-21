//! Oracles: at lambda = 0 every block is an independent n x n torus (exact 2D plaquette); the gluing integral
//! reproduces the exact 2D square defect 4 log Z(n,n) - log Z(2n,2n).
use lattice_gauge::character::{square_defect_2d, torus_plaquette_2d};
use lattice_gauge::glued::GluedLattice;
use lattice_gauge::stats::jackknife;

#[test]
fn blocks_are_independent_tori_at_lambda_zero() {
    let (n, beta) = (2usize, 3.0);
    let mut l = GluedLattice::hot(2, n, beta, 0.0, 5).unwrap();
    for _ in 0..300 { l.update(2); }
    let p: Vec<f64> = (0..30000).map(|_| {
        l.update(2);
        // mean a0 over plaquettes with nonzero weight (bulk and blocked)
        let (mut s, mut c) = (0.0, 0.0);
        for q in 0..l.plaquettes.len() { if l.kind[q] != 1 { s += l.plaquette_value(q); c += 1.0; } }
        s / c
    }).collect();
    let (m, e) = jackknife(&[&p], 20, |v| v[0]);
    let exact = torus_plaquette_2d(beta, (n * n) as u64);
    assert!((m - exact).abs() < 5.0 * e, "{m} +- {e} vs {exact}");
}

#[test]
fn gluing_integral_matches_exact_square_defect() {
    let (n, beta, k) = (2usize, 2.4, 8usize);
    let h = 1.0 / k as f64;
    let (mut integral, mut var) = (0.0, 0.0);
    for i in 0..=k {
        let mut l = GluedLattice::hot(2, n, beta, i as f64 * h, 40 + i as u64).unwrap();
        for _ in 0..300 { l.update(2); }
        let s: Vec<f64> = (0..20000).map(|_| { l.update(2); l.boundary_difference() }).collect();
        let (m, e) = jackknife(&[&s], 20, |v| v[0]);
        let w = if i == 0 || i == k { 1.0 } else if i % 2 == 1 { 4.0 } else { 2.0 };
        integral += w * m;
        var += (w * e).powi(2);
    }
    let delta = -beta * integral * h / 3.0;
    let err = beta * var.sqrt() * h / 3.0;
    let exact = square_defect_2d(beta, n as u64);
    assert!((delta - exact).abs() < 5.0 * err + 0.01, "Delta {delta} +- {err} vs exact {exact}");
}
