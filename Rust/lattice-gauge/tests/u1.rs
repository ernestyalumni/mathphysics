//! Oracles: von Mises marginal E[cos x] = I_1(kappa)/I_0(kappa); exact 2D U(1) torus
//! <cos theta_p> = sum_n I_n^{V-1} I_n' / sum_n I_n^V; gauge invariance; thread-count independence.
use lattice_gauge::character::bessel_i_scaled;
use lattice_gauge::philox::Stream;
use lattice_gauge::stats::jackknife;
use lattice_gauge::u1::*;

#[test]
fn von_mises_marginal() {
    for &k in &[0.5, 2.0, 10.0, 60.0] {
        let mut s = Stream::new(8, 0, (k * 10.0) as u32, 0);
        let n = 300_000;
        let xs: Vec<f64> = (0..n).map(|_| von_mises(k, &mut s).cos()).collect();
        let m = xs.iter().sum::<f64>() / n as f64;
        let sd = (xs.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / n as f64).sqrt();
        let want = bessel_i_scaled(1, k) / bessel_i_scaled(0, k);
        assert!((m - want).abs() < 6.0 * sd / (n as f64).sqrt(), "kappa {k}: {m} vs {want}");
    }
}

#[test]
fn two_dimensional_exact() {
    for &(side, beta) in &[(2usize, 2.0f64), (4, 1.0), (4, 3.0)] {
        let exact = torus_plaquette_u1_2d(beta, (side * side) as u64);
        let mut l = U1Lattice::hot(&[side, side], beta, 60 + side as u64, 0).unwrap();
        for _ in 0..500 { l.update(2); }
        let p: Vec<f64> = (0..40000).map(|_| { l.update(2); l.plaquette() }).collect();
        let (m, e) = jackknife(&[&p], 20, |x| x[0]);
        assert!((m - exact).abs() < 5.0 * e, "L={side} beta={beta}: {m} +- {e} vs {exact}");
    }
}

#[test]
fn gauge_invariance_and_threads() {
    let mut a = U1Lattice::hot(&[4, 4, 4, 4], 1.1, 3, 0).unwrap();
    let mut b = U1Lattice::hot(&[4, 4, 4, 4], 1.1, 3, 0).unwrap();
    a.threads = 1;
    b.threads = 6;
    for _ in 0..5 { a.update(1); b.update(1); }
    assert!(a.theta == b.theta);
    let p = a.plaquette();
    let mut s = Stream::new(4, 0, 0, 0);
    let g: Vec<f64> = (0..a.torus.site_count()).map(|_| 6.0 * s.next_f64_open01()).collect();
    a.gauge_transform(&g);
    assert!((a.plaquette() - p).abs() < 1e-12);
}
