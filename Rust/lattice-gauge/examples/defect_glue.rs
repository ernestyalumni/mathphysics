//! Delta_beta(n) by block gluing (see src/glued.rs): Simpson over lambda in [0, 1].
//! Usage: defect_glue <d> <block n> <beta> <intervals (even)> <therm> <sweeps> <seed>
//! In d = 2 the exact square defect is printed for comparison.
use lattice_gauge::character::square_defect_2d;
use lattice_gauge::glued::GluedLattice;
use lattice_gauge::stats::{jackknife, tau_int};
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let d: usize = a[1].parse().unwrap();
    let n: usize = a[2].parse().unwrap();
    let beta: f64 = a[3].parse().unwrap();
    let k: usize = a[4].parse().unwrap();
    let (therm, sweeps): (usize, usize) = (a[5].parse().unwrap(), a[6].parse().unwrap());
    let seed: u64 = a[7].parse().unwrap();
    let h = 1.0 / k as f64;
    let (mut integral, mut var) = (0.0, 0.0);
    for i in 0..=k {
        let lambda = i as f64 * h;
        let mut l = GluedLattice::hot(d, n, beta, lambda, seed + i as u64).unwrap();
        for _ in 0..therm { l.update(2); }
        let s: Vec<f64> = (0..sweeps).map(|_| { l.update(2); l.boundary_difference() }).collect();
        let (m, e) = jackknife(&[&s], 20, |v| v[0]);
        let w = if i == 0 || i == k { 1.0 } else if i % 2 == 1 { 4.0 } else { 2.0 };
        integral += w * m;
        var += (w * e).powi(2);
        eprintln!("lambda={lambda:.3} <boundary difference>={m:.5}({:.0}) tau={:.2}", e * 1e5, tau_int(&s, 6.0));
    }
    let delta = -beta * integral * h / 3.0;
    let err = beta * var.sqrt() * h / 3.0;
    let exact = if d == 2 { format!(" exact {:.5}", square_defect_2d(beta, n as u64)) } else { String::new() };
    println!("glue d={d} n={n} beta={beta} intervals={k} sweeps={sweeps}: Delta = {delta:.5} +- {err:.5}{exact}");
}
