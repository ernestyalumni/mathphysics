//! Monopole density against beta. SU(2): maximal abelian gauge + abelian projection + DeGrand-Toussaint.
//! U(1): DeGrand-Toussaint on the links directly (the counter-model: monopoles should disappear in the
//! Coulomb phase, beta >~ 1.0). Usage: monopole_scan <su2|u1> <side> <betas, comma-separated> <therm> <configs> <every> <seed>
use lattice_gauge::monopole::{abelian_angles, fix_mag, monopole_density, worldline_topology};
use lattice_gauge::stats::jackknife;
use lattice_gauge::u1::U1Lattice;
use lattice_gauge::wilson::Lattice;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let group = a[1].as_str();
    let side: usize = a[2].parse().unwrap();
    let betas: Vec<f64> = a[3].split(',').map(|x| x.parse().unwrap()).collect();
    let (therm, configs, every): (usize, usize, usize) = (a[4].parse().unwrap(), a[5].parse().unwrap(), a[6].parse().unwrap());
    let seed: u64 = a[7].parse().unwrap();
    for (i, &beta) in betas.iter().enumerate() {
        let (mut rho, mut plaq, mut magf) = (Vec::new(), Vec::new(), Vec::new());
        let (mut share, mut clusters, mut wound) = (Vec::new(), Vec::new(), 0usize);
        if group == "u1" {
            let mut l = U1Lattice::hot(&[side; 4], beta, seed + i as u64, 0).unwrap();
            for _ in 0..therm { l.update(2); }
            for _ in 0..configs {
                for _ in 0..every { l.update(2); }
                rho.push(monopole_density(&l.torus, &l.theta));
                plaq.push(l.plaquette());
                let (_, c, sh, w) = worldline_topology(&l.torus, &l.theta);
                share.push(sh); clusters.push(c as f64); if w.iter().any(|&x| x != 0) { wound += 1; }
            }
        } else {
            let mut l = Lattice::hot(&[side; 4], beta, seed + i as u64, 0).unwrap();
            for k in 0..therm { l.update(2); if k % 100 == 99 { l.reunitarize(); } }
            for _ in 0..configs {
                for _ in 0..every { l.update(2); }
                l.reunitarize();
                let mut g = Lattice::cold(&[side; 4], beta, seed, 1).unwrap();
                g.links = l.links.clone();
                let (r, _) = fix_mag(&mut g, 1.7, 1e-10, 3000);
                let ang = abelian_angles(&g);
                rho.push(monopole_density(&g.torus, &ang));
                let (_, c, sh, w) = worldline_topology(&g.torus, &ang);
                share.push(sh); clusters.push(c as f64); if w.iter().any(|&x| x != 0) { wound += 1; }
                plaq.push(l.plaquette());
                magf.push(r);
            }
        }
        let (m, e) = jackknife(&[&rho], 20, |v| v[0]);
        let (p, pe) = jackknife(&[&plaq], 20, |v| v[0]);
        let mag = if magf.is_empty() { String::new() } else { format!(" mag_functional={:.4}", magf.iter().sum::<f64>() / magf.len() as f64) };
        let mean = |v: &Vec<f64>| v.iter().sum::<f64>() / v.len() as f64;
        println!("{group} side={side} beta={beta:.3} configs={configs}: monopole density per cube = {m:.5}({:.0}) plaquette={p:.5}({:.0}){mag} worldlines: clusters={:.1} largest_share={:.3} configs_with_winding={wound}/{configs}", e * 1e5, pe * 1e5, mean(&clusters), mean(&share));
    }
}
