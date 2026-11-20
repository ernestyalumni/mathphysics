//! Lattice 0++ glueball mass (the lattice mass gap) with GPU updates and CPU measurements.
//! Usage: glueball_gpu <L> <T> <beta> <therm> <meas> <every> <seed> <out.tsv>
//! Lattice L^3 x T, time = direction 3. Operators: spatial plaquettes after 0, 4, 8, 16 APE steps
//! (alpha 0.5). Prints effective masses m_eff(t) = -ln C(t+1)/C(t) with jackknife errors.
//! Reference: Teper, hep-th/9812187, Table 10 (am_0+ = 1.216(23) at beta 2.3, 0.959(20) at 2.4).
use lattice_gauge::cuda::GpuLattice;
use lattice_gauge::glueball::smeared_slice_operators;
use lattice_gauge::stats::jackknife;
use lattice_gauge::wilson::Lattice;
use std::io::Write;
use std::time::Instant;

const STEPS: [usize; 4] = [0, 4, 8, 16];

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let (l, nt): (usize, usize) = (a[1].parse().unwrap(), a[2].parse().unwrap());
    let beta: f64 = a[3].parse().unwrap();
    let (therm, meas, every): (usize, usize, usize) = (a[4].parse().unwrap(), a[5].parse().unwrap(), a[6].parse().unwrap());
    let seed: u64 = a[7].parse().unwrap();
    let mut cpu = Lattice::hot(&[l, l, l, nt], beta, seed, 0).unwrap();
    let mut gpu = GpuLattice::from_cpu(&cpu).unwrap();
    let t0 = Instant::now();
    for _ in 0..therm { gpu.update(2).unwrap(); }
    let mut f = std::io::BufWriter::new(std::fs::File::create(&a[8]).unwrap());
    // ops[level][config][t]
    let mut ops: Vec<Vec<Vec<f64>>> = vec![Vec::with_capacity(meas); STEPS.len()];
    for m in 0..meas {
        for _ in 0..every { gpu.update(2).unwrap(); }
        cpu.links = gpu.download().unwrap();
        if m % 50 == 49 { cpu.reunitarize(); gpu.upload(&cpu.links).unwrap(); }
        let o = smeared_slice_operators(&cpu.torus, &cpu.links, 3, 0.5, &STEPS);
        for (lev, row) in o.into_iter().enumerate() {
            writeln!(f, "{lev}\t{}", row.iter().map(|x| format!("{x:.12e}")).collect::<Vec<_>>().join("\t")).unwrap();
            ops[lev].push(row);
        }
    }
    f.flush().unwrap();
    println!("glueball L={l} T={nt} beta={beta} therm={therm} meas={meas} every={every} seed={seed} minutes={:.1}", t0.elapsed().as_secs_f64() / 60.0);
    let tmax = 4.min(nt / 2 - 1);
    for (lev, cfgs) in ops.iter().enumerate() {
        // per-configuration A_t = mean_t0 O(t0) O(t0+t), B = mean_t0 O(t0)
        let mut cols: Vec<Vec<f64>> = vec![Vec::with_capacity(cfgs.len()); tmax + 2];
        for o in cfgs {
            for t in 0..=tmax {
                cols[t].push((0..nt).map(|t0| o[t0] * o[(t0 + t) % nt]).sum::<f64>() / nt as f64);
            }
            cols[tmax + 1].push(o.iter().sum::<f64>() / nt as f64);
        }
        let refs: Vec<&[f64]> = cols.iter().map(|c| &c[..]).collect();
        let mut line = format!("  APE steps {:>2}:", STEPS[lev]);
        for t in 0..tmax {
            let (m, e) = jackknife(&refs, 20, |v| {
                let b2 = v[tmax + 1] * v[tmax + 1];
                -((v[t + 1] - b2) / (v[t] - b2)).ln()
            });
            line += &format!("  m_eff({t})={m:.3}({:.0})", e * 1000.0);
        }
        println!("{line}");
    }
}
