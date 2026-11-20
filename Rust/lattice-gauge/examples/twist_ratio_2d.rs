//! Z_tw / Z on a 2D torus by thermodynamic integration over beta (the roadmap T08 method), against the
//! exact character-expansion value. Usage: twist_ratio_2d <side> <beta_max> <intervals (even)> <sweeps> <seed>
use lattice_gauge::character::torus_twist_2d;
use lattice_gauge::stats::jackknife;
use lattice_gauge::wilson::Lattice;
fn mean_err(side: usize, beta: f64, twisted: bool, sweeps: usize, seed: u64) -> (f64, f64) {
    let mut l = Lattice::hot(&[side, side], beta, seed, twisted as u32).unwrap();
    if twisted { l.set_twist(0, 1); }
    for _ in 0..500 { l.update(2); }
    let p: Vec<f64> = (0..sweeps).map(|_| { l.update(2); l.action_density() }).collect();
    jackknife(&[&p], 20, |x| x[0])
}
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let side: usize = a[1].parse().unwrap();
    let bmax: f64 = a[2].parse().unwrap();
    let n: usize = a[3].parse().unwrap();
    let sweeps: usize = a[4].parse().unwrap();
    let seed: u64 = a[5].parse().unwrap();
    let v = (side * side) as f64;
    let h = bmax / n as f64;
    let (mut integral, mut var) = (0.0, 0.0);
    for i in 0..=n {
        let beta = i as f64 * h;
        let w = if i == 0 || i == n { 1.0 } else if i % 2 == 1 { 4.0 } else { 2.0 };
        let (mt, et) = mean_err(side, beta, true, sweeps, seed + 2 * i as u64);
        let (m0, e0) = mean_err(side, beta, false, sweeps, seed + 2 * i as u64 + 1);
        integral += w * v * (mt - m0);
        var += (w * v).powi(2) * (et * et + e0 * e0);
    }
    integral *= h / 3.0;
    let err = var.sqrt() * h / 3.0;
    let exact = torus_twist_2d(bmax, v as u64).0;
    println!("side={side} beta={bmax} intervals={n} sweeps={sweeps}: ln(Z_tw/Z) = {integral:.5} +- {err:.5} (MC, Simpson) vs exact {:.5}; Z_tw/Z = {:.5} vs {:.5}",
        exact.ln(), integral.exp(), exact);
}
