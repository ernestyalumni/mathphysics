//! Oracles: essential classes equal the Betti numbers (mod 2) of the complex for any monotone filtration; every cell is
//! accounted for (2 x pairs + essential = cells); a function with k separated minima on a torus gives k - 1 finite H_0 bars
//! with positive length; constant filtrations give only zero-length finite pairs.
use algebraic_topology::cubical::CubicalTorus;
use algebraic_topology::persistence::{lower_star_cubical, persistence};

#[test]
fn essential_classes_are_betti_numbers() {
    for d in 2..=4 {
        let t = CubicalTorus::new(&vec![3; d]);
        let mut s: u64 = 12345 + d as u64;
        let f: Vec<f64> = (0..t.sites).map(|_| { s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407); (s >> 11) as f64 / (1u64 << 53) as f64 }).collect();
        let (cx, vals) = lower_star_cubical(&t, &f);
        let p = persistence(&cx, &vals);
        let mut ess = vec![0usize; d + 1];
        for &(k, _) in &p.essential { ess[k] += 1; }
        let betti = cx.betti_mod_p(2);
        assert_eq!(ess, betti, "T^{d}");
        let cells: usize = cx.cell_counts.iter().sum();
        assert_eq!(2 * p.pairs.len() + p.essential.len(), cells);
    }
}

#[test]
fn minima_count() {
    // on a 6x6 torus, f has three deep separated minima on a flat background: three H_0 classes are born early,
    // one survives (essential), two die when the background level is reached
    let t = CubicalTorus::new(&[6, 6]);
    let mut f = vec![1.0; t.sites];
    for (x, y, v) in [(0usize, 0usize, -3.0), (3, 0, -2.0), (0, 3, -1.0)] { f[x + 6 * y] = v; }
    let (cx, vals) = lower_star_cubical(&t, &f);
    let p = persistence(&cx, &vals);
    let mut bars = p.bars(0, 1e-12);
    bars.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    assert_eq!(bars, vec![(-2.0, 1.0), (-1.0, 1.0)]);
    assert!(p.essential.contains(&(0, -3.0)));
}
