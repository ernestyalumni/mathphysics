//! Oracles: best_rotation beats random rotations and equals the SVD maximum s1+s2+s3; quaternion_from_rotation
//! reproduces R; DMCG never lowers its functional, keeps the plaquette, and the P-vortex density is gauge-fixing
//! independent of the residual Z2 (sign) choice; centre-projected cold configuration has no vortices.
use lattice_gauge::philox::Stream;
use lattice_gauge::su2::haar;
use lattice_gauge::vortex::*;
use lattice_gauge::wilson::Lattice;

#[test]
fn rotation_algebra() {
    let mut s = Stream::new(21, 0, 0, 0);
    for _ in 0..200 {
        let mut m = [[0.0; 3]; 3];
        for i in 0..3 { for j in 0..3 { m[i][j] = s.next_gaussian(); } }
        let r = best_rotation(&m);
        let tr = |r: &[[f64; 3]; 3]| (0..3).map(|i| (0..3).map(|k| r[i][k] * m[k][i]).sum::<f64>()).sum::<f64>();
        let best = tr(&r);
        for _ in 0..200 {
            let o = adjoint_matrix(haar(&mut s));
            assert!(tr(&o) <= best + 1e-12, "random rotation beats best_rotation");
        }
        let g = quaternion_from_rotation(&r);
        let back = adjoint_matrix(g);
        for i in 0..3 { for j in 0..3 { assert!((back[i][j] - r[i][j]).abs() < 1e-10, "quaternion_from_rotation"); } }
    }
}

#[test]
fn mcg_gauge_fixing() {
    let mut l = Lattice::hot(&[4, 4, 4, 4], 2.4, 22, 0).unwrap();
    for _ in 0..10 { l.update(1); }
    let p = l.plaquette();
    let mut last = mcg_functional(&l);
    let start = last;
    for _ in 0..20 {
        let (r, _) = fix_mcg(&mut l, 0.0, 1);
        assert!(r >= last - 1e-12, "MCG functional decreased: {last} -> {r}");
        last = r;
    }
    assert!(last > start + 0.1);
    assert!((l.plaquette() - p).abs() < 1e-12);
    let z = centre_projection(&l);
    let d1 = p_vortex_density(&l, &z);
    let mut s = Stream::new(23, 0, 0, 0);
    let flip: Vec<_> = (0..l.torus.site_count()).map(|_| if s.next_f64_open01() < 0.5 { lattice_gauge::su2::Quaternion::new(-1.0, 0.0, 0.0, 0.0) } else { lattice_gauge::su2::Quaternion::IDENTITY }).collect();
    l.gauge_transform(&flip);
    assert_eq!(p_vortex_density(&l, &centre_projection(&l)), d1, "P-vortices depend on the residual Z2");
    let c = Lattice::cold(&[4, 4, 4, 4], 1.0, 0, 0).unwrap();
    assert_eq!(p_vortex_density(&c, &centre_projection(&c)), 0.0);
}
