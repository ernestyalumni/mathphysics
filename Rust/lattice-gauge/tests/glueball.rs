//! Oracles: APE smearing keeps links in SU(2), leaves time links alone, and is gauge covariant, so
//! the slice operators are gauge invariant; on a cold configuration every operator is 1.
use lattice_gauge::glueball::*;
use lattice_gauge::philox::Stream;
use lattice_gauge::su2::haar;
use lattice_gauge::wilson::Lattice;

#[test]
fn smearing_and_invariance() {
    let mut l = Lattice::hot(&[4, 4, 4, 6], 2.3, 5, 0).unwrap();
    for _ in 0..3 { l.update(1); }
    let t = 3;
    let s = ape_step(&l.torus, &l.links, t, 0.5);
    for site in 0..l.torus.site_count() {
        let i = l.torus.link_index(site, t);
        assert_eq!(s[i], l.links[i], "time link changed");
    }
    assert!(s.iter().all(|q| (q.norm_sq() - 1.0).abs() < 1e-13));
    let before = smeared_slice_operators(&l.torus, &l.links, t, 0.5, &[0, 2, 5]);
    let mut st = Stream::new(31, 0, 0, 0);
    let g: Vec<_> = (0..l.torus.site_count()).map(|_| haar(&mut st)).collect();
    l.gauge_transform(&g);
    let after = smeared_slice_operators(&l.torus, &l.links, t, 0.5, &[0, 2, 5]);
    for (a, b) in before.iter().flatten().zip(after.iter().flatten()) {
        assert!((a - b).abs() < 1e-12, "slice operator not gauge invariant");
    }
}

#[test]
fn cold_configuration() {
    let l = Lattice::cold(&[4, 4, 4, 4], 2.3, 0, 0).unwrap();
    let ops = smeared_slice_operators(&l.torus, &l.links, 3, 0.5, &[0, 3]);
    assert!(ops.iter().flatten().all(|&x| (x - 1.0).abs() < 1e-14));
}

#[test]
fn correlator_at_zero_is_variance() {
    let series = vec![vec![1.0, 2.0, 3.0, 4.0], vec![0.0, 1.0, 0.0, 1.0]];
    let all: Vec<f64> = series.iter().flatten().cloned().collect();
    let mean = all.iter().sum::<f64>() / all.len() as f64;
    let var = all.iter().map(|x| (x - mean) * (x - mean)).sum::<f64>() / all.len() as f64;
    assert!((correlator(&series, mean)[0] - var).abs() < 1e-15);
}
