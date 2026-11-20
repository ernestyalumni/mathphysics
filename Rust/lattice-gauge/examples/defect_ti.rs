//! Defect and purities of yang-mills-lattice.tex (ym:defect, ym:four-purities) from plaquette scans
//! by thermodynamic integration (Simpson in beta; independent errors per volume).
//! Usage:
//!   defect_ti 2d <n> <scan n x n> <scan 2n x 2n>
//!   defect_ti 4d <n> <scan n,n,n,n> <scan 2n,n,n,n> <scan 2n,2n,n,n> <scan 2n,2n,2n,n> <scan 2n,2n,2n,2n>
//! All scans must share the beta grid (even number of intervals).
use lattice_gauge::character::square_defect_2d;

struct Scan { beta: Vec<f64>, mean: Vec<f64>, err: Vec<f64>, np: f64 }

fn read(path: &str) -> Scan {
    let mut s = Scan { beta: vec![], mean: vec![], err: vec![], np: 0.0 };
    for line in std::fs::read_to_string(path).unwrap().lines() {
        let c: Vec<f64> = line.split('\t').map(|x| x.parse().unwrap()).collect();
        s.beta.push(c[0]); s.mean.push(c[1]); s.err.push(c[2]); s.np = c[4];
    }
    s
}

/// (log Z, error) at grid index i (even) by Simpson's rule from beta = 0.
fn log_z(s: &Scan, i: usize) -> (f64, f64) {
    if i == 0 { return (0.0, 0.0); }
    let h = s.beta[1] - s.beta[0];
    let (mut v, mut var) = (0.0, 0.0);
    for k in 0..=i {
        let w = if k == 0 || k == i { 1.0 } else if k % 2 == 1 { 4.0 } else { 2.0 };
        v += w * s.mean[k];
        var += (w * s.err[k]).powi(2);
    }
    (s.np * v * h / 3.0, s.np * var.sqrt() * h / 3.0)
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let n: f64 = a[2].parse().unwrap();
    let scans: Vec<Scan> = a[3..].iter().map(|p| read(p)).collect();
    let m = scans[0].beta.len();
    assert!(scans.iter().all(|s| s.beta == scans[0].beta), "scans must share the beta grid");
    let ln2 = 2f64.ln();
    if a[1] == "2d" {
        println!("beta\tdefect\terr\texact");
        for i in (2..m).step_by(2) {
            let (l0, e0) = log_z(&scans[0], i);
            let (l1, e1) = log_z(&scans[1], i);
            let d = 4.0 * l0 - l1;
            let e = (16.0 * e0 * e0 + e1 * e1).sqrt();
            println!("{:.3}\t{:.5}\t{:.5}\t{:.5}", scans[0].beta[i], d, e, square_defect_2d(scans[0].beta[i], n as u64));
        }
    } else {
        println!("beta\tDelta\terr\tlogP1\tlogP2\tlogP3\tlogP4\tP4>=e^-Delta\tgap_bound(2n)");
        for i in (2..m).step_by(2) {
            let l: Vec<(f64, f64)> = scans.iter().map(|s| log_z(s, i)).collect();
            let delta = 16.0 * l[0].0 - l[4].0;
            let derr = (256.0 * l[0].1 * l[0].1 + l[4].1 * l[4].1).sqrt();
            let lp: Vec<(f64, f64)> = (1..5).map(|k| (l[k].0 - 2.0 * l[k - 1].0, (l[k].1.powi(2) + 4.0 * l[k - 1].1.powi(2)).sqrt())).collect();
            let check = lp[3].0 >= -delta - 3.0 * (lp[3].1 + derr);
            let bound = if delta < ln2 { format!("{:.4}", (1.0 / (delta.exp() - 1.0)).ln() / n) } else { "-".into() };
            println!("{:.3}\t{:.4}\t{:.4}\t{}\t{}\t{}",
                scans[0].beta[i], delta, derr,
                lp.iter().map(|(v, e)| format!("{v:.4}({:.0})", e * 1e4)).collect::<Vec<_>>().join("\t"),
                if check { "ok" } else { "VIOLATED" }, bound);
        }
    }
}
