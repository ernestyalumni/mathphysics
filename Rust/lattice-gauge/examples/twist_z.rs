//! Twist free energy by integration over the twist strength z of the stack plaquettes (YM-N02 estimator):
//! ln(Z_tw/Z) = -beta int_{-1}^{1} <sum_stack a0(U_p)>_z dz (Simpson). In 2D compared with the exact value.
//! Usage: twist_z <dims> <beta> <intervals (even)> <therm> <sweeps> <seed> [threads]
use lattice_gauge::character::torus_twist_2d;
use lattice_gauge::stats::jackknife;
use lattice_gauge::wilson::Lattice;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let dims: Vec<usize> = a[1].split(',').map(|x| x.parse().unwrap()).collect();
    let beta: f64 = a[2].parse().unwrap();
    let n: usize = a[3].parse().unwrap();
    let (therm, sweeps): (usize, usize) = (a[4].parse().unwrap(), a[5].parse().unwrap());
    let seed: u64 = a[6].parse().unwrap();
    let threads: usize = a.get(7).map(|x| x.parse().unwrap()).unwrap_or(1);
    let h = 2.0 / n as f64;
    let (mut integral, mut var) = (0.0, 0.0);
    for i in 0..=n {
        let z = 1.0 - h * i as f64;
        let mut l = Lattice::hot(&dims, beta, seed + i as u64, 0).unwrap();
        l.threads = threads;
        l.set_twist_strength(0, 1, z);
        let stack = l.twist_stack(0, 1);
        for _ in 0..therm { l.update(2); }
        let s: Vec<f64> = (0..sweeps).map(|k| { l.update(2); if k % 500 == 499 { l.reunitarize(); } l.plaquette_sum(&stack) }).collect();
        let (m, e) = jackknife(&[&s], 20, |v| v[0]);
        let w = if i == 0 || i == n { 1.0 } else if i % 2 == 1 { 4.0 } else { 2.0 };
        integral += w * m;
        var += (w * e).powi(2);
        eprintln!("z={z:+.3} <stack a0 sum>={m:.5}({:.0})", e * 1e5);
    }
    let ln_ratio = -beta * integral * h / 3.0;
    let err = beta * var.sqrt() * h / 3.0;
    let exact = if dims.len() == 2 { format!(" exact {:.5}", torus_twist_2d(beta, (dims[0] * dims[1]) as u64).0.ln()) } else { String::new() };
    println!("twist dims={dims:?} beta={beta} intervals={n} sweeps={sweeps}: ln(Z_tw/Z) = {ln_ratio:.5} +- {err:.5}{exact}; Z_tw/Z = {:.5}", ln_ratio.exp());
}
