//! Torelon and 0++ glueball operator streams from the same configurations.
//! Usage: torelon_scan <L> <T> <beta> <therm> <meas> <every> <seed> <out-prefix> [gpu]
//! Lattice L^3 x T, time = direction 3. After 0, 4, 8, 16 APE steps (alpha 0.5, spatial links only) each measurement writes
//!   <out-prefix>-plaq.tsv: the slice average of a0 over spatial plaquettes (as in `glueball_gpu`);
//!   <out-prefix>-tor.tsv:  the slice average, over the three spatial directions j and the L^2 transverse positions, of
//!                          a0 of the Polyakov loop winding once around direction j (zero-momentum torelon operator).
//! Rows are `level<TAB>O(0) .. O(T-1)`, four rows per measurement; analyse each file with `glueball_fit`.
//! A flux tube winding a spatial direction of length L has energy E_T(L) ~ sigma L - pi/(3L) (Nambu-Goto, leading
//! terms); the plaquette operator couples to torelon pairs with energy near 2 E_T(L).
//! With `gpu` (and the `cuda` feature) updates run on the GPU; measurements always run on the CPU.
use lattice_gauge::glueball::{ape_step, slice_plaquettes};
use lattice_gauge::su2::Quaternion;
use lattice_gauge::torus::Torus;
use lattice_gauge::wilson::Lattice;
use std::io::Write;
use std::time::Instant;

const STEPS: [usize; 4] = [0, 4, 8, 16];

fn slice_polyakov(torus: &Torus, links: &[Quaternion], time_dir: usize) -> Vec<f64> {
    let d = torus.dimension();
    let nt = torus.side_lengths()[time_dir];
    let mut sum = vec![0.0; nt];
    let mut count = vec![0usize; nt];
    for site in 0..torus.site_count() {
        let t = torus.site_coordinate(site, time_dir);
        for j in 0..d {
            if j == time_dir || torus.site_coordinate(site, j) != 0 {
                continue;
            }
            let mut p = Quaternion::IDENTITY;
            let mut s = site;
            for _ in 0..torus.side_lengths()[j] {
                p = p.mul(links[torus.link_index(s, j)]);
                s = torus.neighbor(s, j, true);
            }
            sum[t] += p.half_trace();
            count[t] += 1;
        }
    }
    sum.iter().zip(&count).map(|(s, c)| s / *c as f64).collect()
}

fn write_row(f: &mut impl Write, level: usize, row: &[f64]) {
    writeln!(f, "{level}\t{}", row.iter().map(|x| format!("{x:.12e}")).collect::<Vec<_>>().join("\t")).unwrap();
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    if a.len() < 9 {
        eprintln!("usage: torelon_scan <L> <T> <beta> <therm> <meas> <every> <seed> <out-prefix> [gpu]");
        std::process::exit(2);
    }
    let (l, nt): (usize, usize) = (a[1].parse().unwrap(), a[2].parse().unwrap());
    let beta: f64 = a[3].parse().unwrap();
    let (therm, meas, every): (usize, usize, usize) = (a[4].parse().unwrap(), a[5].parse().unwrap(), a[6].parse().unwrap());
    let seed: u64 = a[7].parse().unwrap();
    let use_gpu = a.get(9).map(|s| s == "gpu").unwrap_or(false);
    let mut lattice = Lattice::hot(&[l, l, l, nt], beta, seed, 0).unwrap();
    #[cfg(feature = "cuda")]
    let mut gpu = if use_gpu { Some(lattice_gauge::cuda::GpuLattice::from_cpu(&lattice).unwrap()) } else { None };
    #[cfg(not(feature = "cuda"))]
    if use_gpu {
        eprintln!("built without the cuda feature");
        std::process::exit(2);
    }
    let t0 = Instant::now();
    let mut fp = std::io::BufWriter::new(std::fs::File::create(format!("{}-plaq.tsv", a[8])).unwrap());
    let mut ft = std::io::BufWriter::new(std::fs::File::create(format!("{}-tor.tsv", a[8])).unwrap());
    for m in 0..therm + meas * every {
        #[cfg(feature = "cuda")]
        if let Some(g) = gpu.as_mut() {
            g.update(2).unwrap();
        } else {
            lattice.update(2);
        }
        #[cfg(not(feature = "cuda"))]
        lattice.update(2);
        if m < therm || (m - therm) % every != every - 1 {
            continue;
        }
        let k = (m - therm) / every;
        #[cfg(feature = "cuda")]
        if let Some(g) = gpu.as_mut() {
            lattice.links = g.download().unwrap();
            if k % 50 == 49 {
                lattice.reunitarize();
                g.upload(&lattice.links).unwrap();
            }
        } else if k % 50 == 49 {
            lattice.reunitarize();
        }
        #[cfg(not(feature = "cuda"))]
        if k % 50 == 49 {
            lattice.reunitarize();
        }
        let mut links = lattice.links.clone();
        let mut done = 0;
        for (level, &steps) in STEPS.iter().enumerate() {
            while done < steps {
                links = ape_step(&lattice.torus, &links, 3, 0.5);
                done += 1;
            }
            write_row(&mut fp, level, &slice_plaquettes(&lattice.torus, &links, 3));
            write_row(&mut ft, level, &slice_polyakov(&lattice.torus, &links, 3));
        }
        if k % 1000 == 999 {
            fp.flush().unwrap();
            ft.flush().unwrap();
        }
    }
    fp.flush().unwrap();
    ft.flush().unwrap();
    println!(
        "torelon_scan L={l} T={nt} beta={beta} therm={therm} meas={meas} every={every} seed={seed} gpu={use_gpu} minutes={:.1}",
        t0.elapsed().as_secs_f64() / 60.0
    );
}
