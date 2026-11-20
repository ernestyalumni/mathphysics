//! 2x2 torus at beta = 4: separates the exact torus formula (with higher representations) from
//! the infinite-volume value I_2/I_1.
use lattice_gauge::character::{plaquette_2d_infinite, torus_plaquette_2d};
use lattice_gauge::stats::jackknife;
use lattice_gauge::wilson::Lattice;
fn main() {
    for &(side, beta) in &[(2usize, 4.0f64), (2, 8.0), (4, 8.0)] {
        let mut l = Lattice::hot(&[side, side], beta, 77, 0).unwrap();
        for _ in 0..500 { l.update(2); }
        let p: Vec<f64> = (0..200_000).map(|_| { l.update(2); l.plaquette() }).collect();
        let (m, e) = jackknife(&[&p], 50, |v| v[0]);
        let v = (side * side) as u64;
        println!("L={side} beta={beta}: MC {m:.6} +- {e:.6} | torus formula {:.6} | infinite volume {:.6}", torus_plaquette_2d(beta, v), plaquette_2d_infinite(beta));
    }
}
