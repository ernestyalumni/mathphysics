//! Topology of P-vortex surfaces from a lattice-gauge vortex_scan dump (`# sides L,L,L,L`, then `config site mu nu`).
//! Each pierced plaquette (x; mu, nu) is dual to the square spanned by the complementary directions (rho < sigma) with
//! corners at the dual sites x - a e_rho - b e_sigma, (a, b) = (0,0), (1,0), (1,1), (0,1) (dual site labelled by the
//! original site index). Prints per configuration: components, the largest component's share of the vortex area, its
//! Euler characteristic, closedness, branching edges and orientability; and totals over all components.
use algebraic_topology::surface::{analyze, chain_complex};
use std::collections::BTreeMap;
fn main() {
    let path = std::env::args().nth(1).expect("dump file");
    let text = std::fs::read_to_string(&path).unwrap();
    let mut lines = text.lines();
    let sides: Vec<usize> = lines.next().unwrap().trim_start_matches("# sides ").split(',').map(|x| x.parse().unwrap()).collect();
    assert!(sides.len() == 4 && sides.iter().all(|&n| n >= 3));
    let mut strides = [1usize; 4];
    for i in 1..4 { strides[i] = strides[i - 1] * sides[i - 1]; }
    let back = |s: usize, dir: usize| { let c = (s / strides[dir]) % sides[dir]; if c == 0 { s + (sides[dir] - 1) * strides[dir] } else { s - strides[dir] } };
    let mut configs: BTreeMap<usize, Vec<[usize; 4]>> = BTreeMap::new();
    for l in lines {
        let v: Vec<usize> = l.split('\t').map(|x| x.parse().unwrap()).collect();
        let (cfg, x, mu, nu) = (v[0], v[1], v[2], v[3]);
        let comp: Vec<usize> = (0..4).filter(|&d| d != mu && d != nu).collect();
        let (r, s) = (comp[0], comp[1]);
        let a = x;
        let b = back(x, r);
        let c = back(b, s);
        let d = back(x, s);
        configs.entry(cfg).or_default().push([a, b, c, d]);
    }
    let mut all_closed = true;
    for (cfg, sq) in &configs {
        let comps = analyze(sq);
        let total: usize = comps.iter().map(|c| c.squares).sum();
        let big = &comps[0];
        all_closed &= comps.iter().all(|c| c.closed);
        let chi_sum: i64 = comps.iter().map(|c| c.euler).sum();
        let small_orientable = comps.iter().skip(1).filter(|c| c.orientable == Some(true)).count();
        let small_nonorientable = comps.iter().skip(1).filter(|c| c.orientable == Some(false)).count();
        println!(
            "config {cfg}: plaquettes {total}, components {}, largest share {:.3}, largest chi {} branching {} orientable {:?}; sum chi {chi_sum}; small components orientable/non-orientable/branching {}/{}/{}",
            comps.len(), big.squares as f64 / total as f64, big.euler, big.branching_edges, big.orientable,
            small_orientable, small_nonorientable, comps.len() - 1 - small_orientable - small_nonorientable
        );
        if std::env::var("HOMOLOGY").is_ok() {
            let t0 = std::time::Instant::now();
            let cc = chain_complex(sq);
            let b2 = cc.betti_mod_p(2)[2];
            let bq = cc.betti_mod_p(1_000_003)[2];
            // each small closed orientable component contributes exactly 1 to b2 over Q and over Z2
            println!("   homology: b2(Z2) = {b2}, b2(Q) = {bq}; giant component b2(Q) = {} (after removing {small_orientable} small orientable components) [{:.1}s]",
                bq as i64 - small_orientable as i64, t0.elapsed().as_secs_f64());
        }
    }
    println!("all components closed (Bianchi identity check): {all_closed}");
}
