//! Delta_beta(n) = N_P(2n) int [<a0>_{n^4} - <a0>_{(2n)^4}] dbeta from replica scans, integrated from beta0 = the first
//! common grid point (Simpson), assuming Delta(beta0) ~ 0; the printed integrand lets one check that assumption.
//! Usage: defect_window <small scans, comma-separated> <large scans, comma-separated>
use std::collections::BTreeMap;
fn avg(paths: &str) -> (BTreeMap<i64, (f64, f64)>, f64) {
    let mut acc: BTreeMap<i64, Vec<(f64, f64)>> = BTreeMap::new();
    let mut np = 0.0;
    for p in paths.split(',') {
        for line in std::fs::read_to_string(p).unwrap().lines() {
            let c: Vec<f64> = line.split('\t').map(|x| x.parse().unwrap()).collect();
            acc.entry((c[0] * 1000.0).round() as i64).or_default().push((c[1], c[2]));
            np = c[4];
        }
    }
    let m = acc.into_iter().map(|(b, v)| {
        let k = v.len() as f64;
        let mean = v.iter().map(|x| x.0).sum::<f64>() / k;
        let err = v.iter().map(|x| x.1 * x.1).sum::<f64>().sqrt() / k;
        (b, (mean, err))
    }).collect();
    (m, np)
}
fn main() {
    let a: Vec<String> = std::env::args().collect();
    let (small, _) = avg(&a[1]);
    let (large, np_large) = avg(&a[2]);
    let grid: Vec<i64> = large.keys().filter(|b| small.contains_key(b)).cloned().collect();
    let h = (grid[1] - grid[0]) as f64 / 1000.0;
    println!("beta\tintegrand N_P(2n)(<a0>_n - <a0>_2n)\terr\tDelta-Delta(beta0)\terr");
    let f: Vec<(f64, f64)> = grid.iter().map(|b| {
        let (ms, es) = small[b]; let (ml, el) = large[b];
        (np_large * (ms - ml), np_large * (es * es + el * el).sqrt())
    }).collect();
    for i in 0..grid.len() {
        let (d, e) = if i >= 2 && i % 2 == 0 {
            let (mut v, mut var) = (0.0, 0.0);
            for k in 0..=i { let w = if k == 0 || k == i { 1.0 } else if k % 2 == 1 { 4.0 } else { 2.0 }; v += w * f[k].0; var += (w * f[k].1).powi(2); }
            (format!("{:.4}", v * h / 3.0), format!("{:.4}", var.sqrt() * h / 3.0))
        } else { ("".into(), "".into()) };
        println!("{:.3}\t{:.4}\t{:.4}\t{d}\t{e}", grid[i] as f64 / 1000.0, f[i].0, f[i].1);
    }
}
