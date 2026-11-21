//! Oracles: monopole charges are invariant under U(1) gauge transformations with re-wrapping (tests the face-sign
//! pattern); total magnetic charge in every closed 3-slice is zero; cold configuration has none; MAG never lowers its
//! functional, leaves the plaquette unchanged, and is stationary at convergence.
use lattice_gauge::monopole::*;
use lattice_gauge::philox::Stream;
use lattice_gauge::torus::Torus;
use lattice_gauge::u1::U1Lattice;
use lattice_gauge::wilson::Lattice;

fn charges(t: &Torus, theta: &[f64]) -> Vec<i64> {
    let n = plaquette_integers(t, theta);
    let mut m = Vec::new();
    for x in 0..t.site_count() {
        for a in 0..4 { for b in a + 1..4 { for c in b + 1..4 { m.push(cube_charge(t, &n, x, a, b, c)); } } }
    }
    m
}

#[test]
fn charges_gauge_invariant_and_neutral() {
    let mut l = U1Lattice::hot(&[4, 4, 4, 6], 0.5, 5, 0).unwrap();
    for _ in 0..10 { l.update(1); }
    let m0 = charges(&l.torus, &l.theta);
    assert!(m0.iter().any(|&q| q != 0), "expected monopoles at strong coupling");
    let mut s = Stream::new(2, 0, 0, 0);
    let g: Vec<f64> = (0..l.torus.site_count()).map(|_| 20.0 * s.next_f64_open01()).collect();
    l.gauge_transform(&g);
    assert_eq!(m0, charges(&l.torus, &l.theta), "monopole charges changed under a gauge transformation");
    // closed 3-slices: cubes (0,1,2) at fixed x_3 sum to zero
    let t = &l.torus;
    let n = plaquette_integers(t, &l.theta);
    for x3 in 0..6 {
        let total: i64 = (0..t.site_count()).filter(|&x| t.site_coordinate(x, 3) == x3).map(|x| cube_charge(t, &n, x, 0, 1, 2)).sum();
        assert_eq!(total, 0, "net magnetic charge in slice x3={x3}");
    }
    let c = U1Lattice::cold(&[4, 4, 4, 4], 1.0, 0, 0).unwrap();
    assert_eq!(monopole_density(&c.torus, &c.theta), 0.0);
}

#[test]
fn mag_gauge_fixing() {
    let mut l = Lattice::hot(&[4, 4, 4, 4], 2.4, 6, 0).unwrap();
    for _ in 0..10 { l.update(1); }
    let p = l.plaquette();
    let r0 = mag_functional(&l);
    let mut last = r0;
    for _ in 0..30 {
        let (r, _) = fix_mag(&mut l, 1.0, 0.0, 1);
        assert!(r >= last - 1e-12, "MAG functional decreased: {last} -> {r}");
        last = r;
    }
    assert!(last > r0 + 0.1, "MAG barely moved: {r0} -> {last}");
    assert!((l.plaquette() - p).abs() < 1e-12);
    let (r, sweeps) = fix_mag(&mut l, 1.7, 1e-12, 5000);
    assert!(sweeps < 5000 && r >= last - 1e-12, "MAG did not converge");
}
