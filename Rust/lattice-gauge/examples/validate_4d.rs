//! YM-N01 validation runs in d = 4. Usage:
//!   validate_4d <side> <beta> <therm> <meas> <overrelax> <seed>
//! Prints one line of key=value results: plaquette <a0(U_p)> with jackknife error, the integrated
//! autocorrelation time of the plaquette, the pseudo specific heat C_V = beta^2 N_P Var(P) per
//! configuration (Engels and Scheideler, hep-lat/9607041, eq. 6) with jackknife error, and timing.
use lattice_gauge::stats::{jackknife, tau_int};
use lattice_gauge::wilson::Lattice;
use std::time::Instant;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let side: usize = a[1].parse().unwrap();
    let beta: f64 = a[2].parse().unwrap();
    let therm: usize = a[3].parse().unwrap();
    let meas: usize = a[4].parse().unwrap();
    let or: usize = a[5].parse().unwrap();
    let seed: u64 = a[6].parse().unwrap();
    let mut l = Lattice::hot(&[side; 4], beta, seed, 0).unwrap();
    let t0 = Instant::now();
    for i in 0..therm {
        l.update(or);
        if i % 100 == 99 {
            l.reunitarize();
        }
    }
    let t1 = Instant::now();
    let mut p = Vec::with_capacity(meas);
    for i in 0..meas {
        l.update(or);
        if i % 100 == 99 {
            l.reunitarize();
        }
        p.push(l.plaquette());
    }
    let t2 = Instant::now();
    let np = l.torus.plaquette_count() as f64;
    let p2: Vec<f64> = p.iter().map(|x| x * x).collect();
    let bins = 20;
    let (mean, err) = jackknife(&[&p], bins, |m| m[0]);
    let (cv, cv_err) = jackknife(&[&p, &p2], bins, |m| beta * beta * np * (m[1] - m[0] * m[0]));
    let tau = tau_int(&p, 6.0);
    let sweep_ms = (t2 - t1).as_secs_f64() * 1e3 / (meas as f64 * (1 + or) as f64);
    println!(
        "side={side} beta={beta} therm={therm} meas={meas} or={or} seed={seed} threads={} plaquette={mean:.6} err={err:.6} tau_int={tau:.2} cv={cv:.4} cv_err={cv_err:.4} ms_per_sweep={sweep_ms:.3} therm_s={:.1}",
        l.threads,
        (t1 - t0).as_secs_f64()
    );
}
