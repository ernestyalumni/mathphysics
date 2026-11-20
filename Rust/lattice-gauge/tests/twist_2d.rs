//! Oracle: exact twisted 2D torus from the character expansion, Z_tw = sum_n (-1)^{n-1} b_n^V
//! (character::torus_twist_2d), against Monte Carlo of the twisted action <z_p a0(U_p)>.
use lattice_gauge::character::{torus_plaquette_2d, torus_twist_2d};
use lattice_gauge::stats::jackknife;
use lattice_gauge::wilson::Lattice;

#[test]
fn twisted_action_matches_exact() {
    for &(side, beta) in &[(2usize, 4.0f64), (4, 3.0), (4, 6.0)] {
        let v = (side * side) as u64;
        let (ratio, exact) = torus_twist_2d(beta, v);
        let mut l = Lattice::hot(&[side, side], beta, 50 + side as u64, 0).unwrap();
        l.set_twist(0, 1);
        for _ in 0..500 { l.update(2); }
        let p: Vec<f64> = (0..40000).map(|_| { l.update(2); l.action_density() }).collect();
        let (m, e) = jackknife(&[&p], 20, |x| x[0]);
        eprintln!("L={side} beta={beta}: twisted MC {m:.5}({e:.5}) exact {exact:.5}; untwisted {:.5}; Z_tw/Z = {ratio:.4}", torus_plaquette_2d(beta, v));
        assert!((m - exact).abs() < 5.0 * e, "L={side} beta={beta}: {m} +- {e} vs {exact}");
    }
}

#[test]
fn twist_flips_one_plaquette_per_plane() {
    let mut l = Lattice::cold(&[4, 4, 4, 4], 1.0, 0, 0).unwrap();
    l.set_twist(0, 1);
    let z = l.plaquette_sign.as_ref().unwrap();
    assert_eq!(z.iter().filter(|&&s| s < 0.0).count(), 16);
    // cold start: sum_p z_p / N_P = 1 - 2 * 16 / 1536
    assert!((l.action_density() - (1.0 - 32.0 / 1536.0)).abs() < 1e-15);
}
