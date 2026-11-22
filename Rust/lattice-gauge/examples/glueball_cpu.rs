//! Lattice 0++ glueball operator streams with CPU updates (no CUDA), in the format of `glueball_gpu`.
//! Usage: glueball_cpu <L> <T> <beta> <therm> <meas> <every> <seed> <out.tsv>
//! Lattice L^3 x T, time = direction 3. Each measurement writes four rows `level<TAB>O(0) .. O(T-1)`: the
//! zero-momentum spatial-plaquette operator on each time slice after 0, 4, 8, 16 APE steps (alpha 0.5).
//! Each update is one heat-bath sweep and two over-relaxation sweeps. Analyse the streams with `glueball_fit`.
use lattice_gauge::glueball::smeared_slice_operators;
use lattice_gauge::wilson::Lattice;
use std::io::Write;
use std::time::Instant;

const STEPS: [usize; 4] = [0, 4, 8, 16];

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() != 9 {
        eprintln!("usage: glueball_cpu <L> <T> <beta> <therm> <meas> <every> <seed> <out.tsv>");
        std::process::exit(2);
    }
    let (l, nt): (usize, usize) = (a[1].parse().unwrap(), a[2].parse().unwrap());
    let beta: f64 = a[3].parse().unwrap();
    let (therm, meas, every): (usize, usize, usize) = (a[4].parse().unwrap(), a[5].parse().unwrap(), a[6].parse().unwrap());
    let seed: u64 = a[7].parse().unwrap();
    let mut lattice = Lattice::hot(&[l, l, l, nt], beta, seed, 0).unwrap();
    let t0 = Instant::now();
    for _ in 0..therm {
        lattice.update(2);
    }
    let mut f = std::io::BufWriter::new(std::fs::File::create(&a[8]).unwrap());
    for m in 0..meas {
        for _ in 0..every {
            lattice.update(2);
        }
        if m % 50 == 49 {
            lattice.reunitarize();
        }
        let o = smeared_slice_operators(&lattice.torus, &lattice.links, 3, 0.5, &STEPS);
        for (level, row) in o.into_iter().enumerate() {
            writeln!(f, "{level}\t{}", row.iter().map(|x| format!("{x:.12e}")).collect::<Vec<_>>().join("\t")).unwrap();
        }
        if m % 1000 == 999 {
            f.flush().unwrap();
        }
    }
    f.flush().unwrap();
    println!(
        "glueball_cpu L={l} T={nt} beta={beta} therm={therm} meas={meas} every={every} seed={seed} minutes={:.1}",
        t0.elapsed().as_secs_f64() / 60.0
    );
}
