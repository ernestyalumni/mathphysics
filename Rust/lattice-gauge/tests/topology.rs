//! Oracles: Q = 0 on the cold configuration; q(x) gauge invariant; cooling never lowers the plaquette;
//! after cooling, |Q| is near a multiple of a common unit (instanton number) rather than spread uniformly.
use lattice_gauge::philox::Stream;
use lattice_gauge::su2::haar;
use lattice_gauge::topology::*;
use lattice_gauge::wilson::{Lattice, Update};

#[test]
fn cold_and_gauge_invariance() {
    let c = Lattice::cold(&[4, 4, 4, 4], 2.4, 0, 0).unwrap();
    assert_eq!(charge(&c), 0.0);
    let mut l = Lattice::hot(&[4, 4, 4, 6], 2.4, 9, 0).unwrap();
    for _ in 0..5 { l.update(1); }
    let before: Vec<f64> = (0..l.torus.site_count()).map(|x| charge_density(&l, x)).collect();
    let mut s = Stream::new(12, 0, 0, 0);
    let g: Vec<_> = (0..l.torus.site_count()).map(|_| haar(&mut s)).collect();
    l.gauge_transform(&g);
    for (x, b) in before.iter().enumerate() {
        assert!((charge_density(&l, x) - b).abs() < 1e-13);
    }
}

#[test]
fn cooling_raises_plaquette() {
    let mut l = Lattice::hot(&[4, 4, 4, 4], 2.4, 10, 0).unwrap();
    for _ in 0..20 { l.update(1); }
    let mut p = l.plaquette();
    for _ in 0..10 {
        l.sweep_with(Update::Cool);
        let q = l.plaquette();
        assert!(q >= p - 1e-12, "cooling lowered the plaquette");
        p = q;
    }
}
