//! P-vortex density and centre-projected Creutz ratios after direct maximal centre gauge, against beta.
//! Usage: vortex_scan <side> <betas> <therm> <configs> <every> <seed> <threads>
use lattice_gauge::stats::jackknife;
use lattice_gauge::su2::Quaternion;
use lattice_gauge::vortex::{centre_projection, fix_mcg, p_vortex_density};
use lattice_gauge::wilson::Lattice;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let side: usize = a[1].parse().unwrap();
    let betas: Vec<f64> = a[2].split(',').map(|x| x.parse().unwrap()).collect();
    let (therm, configs, every): (usize, usize, usize) = (a[3].parse().unwrap(), a[4].parse().unwrap(), a[5].parse().unwrap());
    let seed: u64 = a[6].parse().unwrap();
    let threads: usize = a[7].parse().unwrap();
    let rmax = 3;
    for (i, &beta) in betas.iter().enumerate() {
        let mut l = Lattice::hot(&[side; 4], beta, seed + i as u64, 0).unwrap();
        l.threads = threads;
        for k in 0..therm { l.update(2); if k % 100 == 99 { l.reunitarize(); } }
        let mut dens = Vec::new();
        let mut w: Vec<Vec<f64>> = vec![Vec::new(); rmax * rmax];
        let mut mcg = 0.0;
        for _ in 0..configs {
            for _ in 0..every { l.update(2); }
            l.reunitarize();
            let mut g = Lattice::cold(&[side; 4], beta, seed, 1).unwrap();
            g.threads = threads;
            g.links = l.links.clone();
            let (r, _) = fix_mcg(&mut g, 1e-9, 2000);
            mcg += r / configs as f64;
            let z = centre_projection(&g);
            dens.push(p_vortex_density(&g, &z));
            let mut c = Lattice::cold(&[side; 4], beta, seed, 2).unwrap();
            c.threads = threads;
            c.links = z.iter().map(|&s| Quaternion::new(s as f64, 0.0, 0.0, 0.0)).collect();
            let t = c.wilson_loop_table(rmax);
            for r in 0..rmax { for s in 0..rmax { w[r * rmax + s].push(t[r][s]); } }
        }
        let (p, pe) = jackknife(&[&dens], 20, |v| v[0]);
        let mut line = format!("vortex side={side} beta={beta:.3} configs={configs} mcg_functional={mcg:.4} p_vortex_density={p:.5}({:.0})", pe * 1e5);
        for r in 2..=rmax {
            let idx = |a: usize, b: usize| (a - 1) * rmax + (b - 1);
            let cols = [&w[idx(r, r)][..], &w[idx(r - 1, r - 1)][..], &w[idx(r, r - 1)][..], &w[idx(r - 1, r)][..]];
            let (chi, e) = jackknife(&cols, 20, |m| -((m[0] * m[1]) / (m[2] * m[3])).ln());
            line += &format!(" chi_cp({r},{r})={chi:.4}({:.0})", e * 1e4);
        }
        line += &format!(" [random-vortex estimate -ln(1-2p)={:.4}]", -(1.0 - 2.0 * p).ln());
        println!("{line}");
    }
}
