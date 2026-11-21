//! Cooled topological charge densities for persistence analysis. Usage: dump_charge <side> <beta> <therm> <configs> <every> <cool> <seed> <out.tsv>
//! Writes `# sides L,L,L,L`, then per configuration a line `Q <value>` followed by q(x) for every site in index order (one line, tab-separated).
use lattice_gauge::topology::{charge, charge_density, cool};
use lattice_gauge::wilson::Lattice;
use std::io::Write;
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let side: usize = a[1].parse().unwrap();
    let beta: f64 = a[2].parse().unwrap();
    let (therm, configs, every, ncool): (usize, usize, usize, usize) = (a[3].parse().unwrap(), a[4].parse().unwrap(), a[5].parse().unwrap(), a[6].parse().unwrap());
    let seed: u64 = a[7].parse().unwrap();
    let mut f = std::io::BufWriter::new(std::fs::File::create(&a[8]).unwrap());
    writeln!(f, "# sides {}", vec![side.to_string(); 4].join(",")).unwrap();
    let mut l = Lattice::hot(&[side; 4], beta, seed, 0).unwrap();
    for i in 0..therm { l.update(2); if i % 100 == 99 { l.reunitarize(); } }
    for _ in 0..configs {
        for _ in 0..every { l.update(2); }
        l.reunitarize();
        let mut k = Lattice::cold(&[side; 4], beta, seed, 1).unwrap();
        k.links = l.links.clone();
        cool(&mut k, ncool);
        writeln!(f, "Q\t{:.6}", charge(&k)).unwrap();
        let q: Vec<String> = (0..k.torus.site_count()).map(|x| format!("{:.6e}", charge_density(&k, x))).collect();
        writeln!(f, "{}", q.join("\t")).unwrap();
    }
}
