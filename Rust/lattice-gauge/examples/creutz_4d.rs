//! Creutz ratios in d = 4 (YM-N01). Usage:
//!   creutz_4d <side> <beta> <therm> <meas> <every> <overrelax> <seed> <rmax> <out.tsv>
//! Measures all W(r,t), r,t <= rmax, every `every` updates; writes the time series to out.tsv
//! (columns: plaquette, then W(1,1), W(1,2), ..., W(rmax,rmax)) and prints chi(R,R) and chi(R,R-1)
//! with jackknife errors. Reference for a sqrt(sigma): Teper, hep-th/9812187, Table 10.
use lattice_gauge::stats::{jackknife, tau_int};
use lattice_gauge::wilson::Lattice;
use std::io::Write;
use std::time::Instant;

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let side: usize = a[1].parse().unwrap();
    let beta: f64 = a[2].parse().unwrap();
    let therm: usize = a[3].parse().unwrap();
    let meas: usize = a[4].parse().unwrap();
    let every: usize = a[5].parse().unwrap();
    let or: usize = a[6].parse().unwrap();
    let seed: u64 = a[7].parse().unwrap();
    let rmax: usize = a[8].parse().unwrap();
    let out = &a[9];
    let mut l = Lattice::hot(&[side; 4], beta, seed, 0).unwrap();
    let t0 = Instant::now();
    for i in 0..therm {
        l.update(or);
        if i % 100 == 99 {
            l.reunitarize();
        }
    }
    let mut f = std::io::BufWriter::new(std::fs::File::create(out).unwrap());
    let mut series: Vec<Vec<f64>> = vec![Vec::with_capacity(meas); 1 + rmax * rmax];
    for m in 0..meas {
        for _ in 0..every {
            l.update(or);
        }
        if m % 20 == 19 {
            l.reunitarize();
        }
        let w = l.wilson_loop_table(rmax);
        let mut row = vec![l.plaquette()];
        for r in 0..rmax {
            row.extend_from_slice(&w[r]);
        }
        writeln!(f, "{}", row.iter().map(|x| format!("{x:.12e}")).collect::<Vec<_>>().join("\t")).unwrap();
        for (s, x) in series.iter_mut().zip(row) {
            s.push(x);
        }
    }
    f.flush().unwrap();
    let idx = |r: usize, t: usize| 1 + (r - 1) * rmax + (t - 1);
    let bins = 20;
    let (p, pe) = jackknife(&[&series[0]], bins, |m| m[0]);
    println!(
        "side={side} beta={beta} therm={therm} meas={meas} every={every} or={or} seed={seed} plaquette={p:.6}({pe:.6}) tau_plaq={:.2} minutes={:.1}",
        tau_int(&series[0], 6.0),
        t0.elapsed().as_secs_f64() / 60.0
    );
    for r in 2..=rmax {
        for t in [r - 1, r] {
            if t < 2 {
                continue;
            }
            let cols = [&series[idx(r, t)][..], &series[idx(r - 1, t - 1)][..], &series[idx(r, t - 1)][..], &series[idx(r - 1, t)][..]];
            let (chi, err) = jackknife(&cols, bins, |m| -((m[0] * m[1]) / (m[2] * m[3])).ln());
            println!("  chi({r},{t}) = {chi:.5} +- {err:.5}");
        }
    }
}
