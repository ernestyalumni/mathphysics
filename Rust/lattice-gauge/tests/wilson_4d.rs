//! Oracles: gauge invariance of plaquette and Wilson loops; the strong-coupling limit
//! <a0(U_p)> = beta/4 + O(beta^3) in any dimension; cold start gives W = 1.
use lattice_gauge::philox::Stream;
use lattice_gauge::stats::jackknife;
use lattice_gauge::su2::haar;
use lattice_gauge::wilson::Lattice;

#[test]
fn gauge_invariance() {
    let mut l = Lattice::hot(&[4, 4, 4, 4], 2.3, 3, 0).unwrap();
    for _ in 0..5 { l.update(1); }
    let before = (l.plaquette(), l.wilson_loop(2, 2), l.wilson_loop(1, 3));
    let mut s = Stream::new(99, 0, 0, 0);
    let g: Vec<_> = (0..l.torus.site_count()).map(|_| haar(&mut s)).collect();
    l.gauge_transform(&g);
    let after = (l.plaquette(), l.wilson_loop(2, 2), l.wilson_loop(1, 3));
    assert!((before.0 - after.0).abs() < 1e-12 && (before.1 - after.1).abs() < 1e-12 && (before.2 - after.2).abs() < 1e-12);
}

#[test]
fn cold_start_and_unit_loop() {
    let l = Lattice::cold(&[4, 4, 4, 4], 2.3, 0, 0).unwrap();
    assert_eq!(l.plaquette(), 1.0);
    assert_eq!(l.wilson_loop(2, 3), 1.0);
    assert!((l.wilson_loop(1, 1) - l.plaquette()).abs() < 1e-15);
}

#[test]
fn strong_coupling_slope() {
    let beta = 0.1;
    let mut l = Lattice::hot(&[4, 4, 4, 4], beta, 5, 0).unwrap();
    for _ in 0..20 { l.update(1); }
    let p: Vec<f64> = (0..400).map(|_| { l.update(1); l.plaquette() }).collect();
    let (m, e) = jackknife(&[&p], 20, |v| v[0]);
    // beta/4 - beta^3/96 from the one-plaquette ratio I_2/I_1; higher orders are O(beta^5) in 4D.
    let want = beta / 4.0 - beta.powi(3) / 96.0;
    assert!((m - want).abs() < 5.0 * e + 1e-5, "{m} +- {e} vs {want}");
}

#[test]
fn loop_table_matches_single_loops_and_threads() {
    let mut l = Lattice::hot(&[4, 4, 6, 4], 2.4, 8, 0).unwrap();
    for _ in 0..5 { l.update(1); }
    l.threads = 1;
    let t1 = l.wilson_loop_table(3);
    l.threads = 7;
    let t7 = l.wilson_loop_table(3);
    assert_eq!(t1, t7, "loop table depends on thread count");
    for r in 1..=3 {
        for t in 1..=3 {
            assert!((t1[r - 1][t - 1] - l.wilson_loop(r, t)).abs() < 1e-12, "W({r},{t})");
        }
    }
}
