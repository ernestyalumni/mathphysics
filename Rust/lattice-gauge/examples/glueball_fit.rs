//! Variational analysis and correlated fits of 0++ glueball operator streams (from `glueball_gpu` or `glueball_cpu`).
//! Usage: glueball_fit <T> <bins per replica> <replica.tsv>[:<replica.tsv>...]
//!
//! For each measurement and smearing levels i, j the correlator matrix is
//!   C_ij(t) = <A_ij(t)> - <B_i><B_j>,  A_ij(t) = (1/T) sum_t0 [O_i(t0) O_j(t0+t) + O_j(t0) O_i(t0+t)] / 2,  B_i = (1/T) sum_t0 O_i(t0).
//! Errors: jackknife over contiguous bins within each replica (all replicas pooled), applied to the whole analysis.
//!
//! Variational step: solve C(t1) v = lambda C(t0) v with (t0, t1) = (0, 1) on the full sample (Cholesky of C(0), then
//! Jacobi on the symmetric reduced matrix) and keep the vector of the largest lambda. The projected correlator
//! c(t) = v^T C(t) v is analysed with v held fixed.
//!
//! Output, all for the projected correlator: effective masses -ln c(t+1)/c(t) and the periodic (cosh) effective mass; then
//! one-state fits c(t) = A [exp(-m t) + exp(-m (T - t))] over every window [t_min, t_max] with 1 <= t_min < t_max <= T/2,
//! correlated chi^2 with the jackknife covariance, m and its jackknife error (refitted in every jackknife sample).
//! Windows are all printed: choosing one afterwards is a decision to report, not a property of the data.
//! Pre-registered for YM-N05 (2026-11-21, before the new runs): primary window [2,5], check window [3,6].
use std::io::BufRead;

const LEVELS: usize = 4;

fn read_replica(path: &str, nt: usize) -> Vec<Vec<Vec<f64>>> {
    let f = std::io::BufReader::new(std::fs::File::open(path).unwrap_or_else(|e| panic!("{path}: {e}")));
    let mut out: Vec<Vec<Vec<f64>>> = Vec::new();
    let mut cur: Vec<Vec<f64>> = Vec::new();
    for line in f.lines() {
        let line = line.unwrap();
        let mut it = line.split('\t');
        let level: usize = it.next().unwrap().parse().unwrap();
        let row: Vec<f64> = it.map(|x| x.parse().unwrap()).collect();
        assert_eq!(row.len(), nt, "{path}: row of length {} for T = {nt}", row.len());
        assert_eq!(level, cur.len(), "{path}: levels out of order");
        cur.push(row);
        if cur.len() == LEVELS {
            out.push(std::mem::take(&mut cur));
        }
    }
    out
}

/// Primary per-measurement quantities: A_ij(t) for i <= j, t = 0..=tmax, then B_i.
fn primaries(o: &[Vec<f64>], nt: usize, tmax: usize) -> Vec<f64> {
    let mut v = Vec::with_capacity(LEVELS * (LEVELS + 1) / 2 * (tmax + 1) + LEVELS);
    for i in 0..LEVELS {
        for j in i..LEVELS {
            for t in 0..=tmax {
                let s: f64 = (0..nt).map(|t0| o[i][t0] * o[j][(t0 + t) % nt] + o[j][t0] * o[i][(t0 + t) % nt]).sum();
                v.push(s / (2.0 * nt as f64));
            }
        }
    }
    for i in 0..LEVELS {
        v.push(o[i].iter().sum::<f64>() / nt as f64);
    }
    v
}

/// Correlator matrices C[t][i][j] from primary means.
fn matrices(p: &[f64], tmax: usize) -> Vec<[[f64; LEVELS]; LEVELS]> {
    let b = &p[p.len() - LEVELS..];
    let mut c = vec![[[0.0; LEVELS]; LEVELS]; tmax + 1];
    let mut k = 0;
    for i in 0..LEVELS {
        for j in i..LEVELS {
            for t in 0..=tmax {
                let x = p[k] - b[i] * b[j];
                c[t][i][j] = x;
                c[t][j][i] = x;
                k += 1;
            }
        }
    }
    c
}

fn cholesky(a: &[[f64; LEVELS]; LEVELS]) -> Option<[[f64; LEVELS]; LEVELS]> {
    let mut l = [[0.0; LEVELS]; LEVELS];
    for i in 0..LEVELS {
        for j in 0..=i {
            let s: f64 = a[i][j] - (0..j).map(|k| l[i][k] * l[j][k]).sum::<f64>();
            if i == j {
                if s <= 0.0 {
                    return None;
                }
                l[i][i] = s.sqrt();
            } else {
                l[i][j] = s / l[j][j];
            }
        }
    }
    Some(l)
}

/// Eigen-decomposition of a symmetric matrix by cyclic Jacobi rotations: (eigenvalues, eigenvectors as columns).
fn jacobi(mut a: [[f64; LEVELS]; LEVELS]) -> ([f64; LEVELS], [[f64; LEVELS]; LEVELS]) {
    let mut v = [[0.0; LEVELS]; LEVELS];
    for i in 0..LEVELS {
        v[i][i] = 1.0;
    }
    for _ in 0..100 {
        let off: f64 = (0..LEVELS).flat_map(|i| (0..LEVELS).map(move |j| (i, j))).filter(|(i, j)| i != j).map(|(i, j)| a[i][j] * a[i][j]).sum();
        if off < 1e-30 {
            break;
        }
        for p in 0..LEVELS {
            for q in p + 1..LEVELS {
                if a[p][q].abs() < 1e-300 {
                    continue;
                }
                let theta = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                for k in 0..LEVELS {
                    let (akp, akq) = (a[k][p], a[k][q]);
                    a[k][p] = c * akp - s * akq;
                    a[k][q] = s * akp + c * akq;
                }
                for k in 0..LEVELS {
                    let (apk, aqk) = (a[p][k], a[q][k]);
                    a[p][k] = c * apk - s * aqk;
                    a[q][k] = s * apk + c * aqk;
                }
                for k in 0..LEVELS {
                    let (vkp, vkq) = (v[k][p], v[k][q]);
                    v[k][p] = c * vkp - s * vkq;
                    v[k][q] = s * vkp + c * vkq;
                }
            }
        }
    }
    let mut e = [0.0; LEVELS];
    for i in 0..LEVELS {
        e[i] = a[i][i];
    }
    (e, v)
}

/// Generalized eigenvector of C(1) v = lambda C(0) v for the largest lambda (normalized v^T C(0) v = 1).
fn variational_vector(c0: &[[f64; LEVELS]; LEVELS], c1: &[[f64; LEVELS]; LEVELS]) -> ([f64; LEVELS], f64) {
    let l = cholesky(c0).expect("C(0) is not positive definite");
    // M = L^{-1} C1 L^{-T}
    let solve_lower = |b: [f64; LEVELS]| -> [f64; LEVELS] {
        let mut x = [0.0; LEVELS];
        for i in 0..LEVELS {
            x[i] = (b[i] - (0..i).map(|k| l[i][k] * x[k]).sum::<f64>()) / l[i][i];
        }
        x
    };
    let mut y = [[0.0; LEVELS]; LEVELS]; // y = L^{-1} C1, column by column
    for j in 0..LEVELS {
        let col = solve_lower([c1[0][j], c1[1][j], c1[2][j], c1[3][j]]);
        for i in 0..LEVELS {
            y[i][j] = col[i];
        }
    }
    let mut m = [[0.0; LEVELS]; LEVELS]; // m = y L^{-T} = (L^{-1} y^T)^T
    for i in 0..LEVELS {
        let row = solve_lower(y[i]);
        for j in 0..LEVELS {
            m[i][j] = row[j];
        }
    }
    for i in 0..LEVELS {
        for j in i + 1..LEVELS {
            let s = 0.5 * (m[i][j] + m[j][i]);
            m[i][j] = s;
            m[j][i] = s;
        }
    }
    let (e, w) = jacobi(m);
    let k = (0..LEVELS).max_by(|&a, &b| e[a].partial_cmp(&e[b]).unwrap()).unwrap();
    // v = L^{-T} w_k
    let mut v = [0.0; LEVELS];
    for i in (0..LEVELS).rev() {
        v[i] = (w[i][k] - (i + 1..LEVELS).map(|r| l[r][i] * v[r]).sum::<f64>()) / l[i][i];
    }
    (v, e[k])
}

fn project(c: &[[f64; LEVELS]; LEVELS], v: &[f64; LEVELS]) -> f64 {
    (0..LEVELS).map(|i| (0..LEVELS).map(|j| v[i] * c[i][j] * v[j]).sum::<f64>()).sum()
}

/// Mass with cosh(m (t + 1 - T/2)) / cosh(m (t - T/2)) = r, by bisection; NaN if no solution.
fn cosh_mass(r: f64, t: usize, nt: usize) -> f64 {
    let h = nt as f64 / 2.0;
    let g = |m: f64| (m * (t as f64 + 1.0 - h)).cosh() / (m * (t as f64 - h)).cosh() - r;
    let (mut lo, mut hi) = (1e-6, 10.0);
    if !(r > 0.0 && r < 1.0) || g(lo) * g(hi) > 0.0 {
        return f64::NAN;
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if g(lo) * g(mid) <= 0.0 { hi = mid } else { lo = mid }
    }
    0.5 * (lo + hi)
}

fn invert(a: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = a.len();
    let mut m: Vec<Vec<f64>> = a.iter().enumerate().map(|(i, r)| { let mut r = r.clone(); r.extend((0..n).map(|j| if i == j { 1.0 } else { 0.0 })); r }).collect();
    for col in 0..n {
        let piv = (col..n).max_by(|&x, &y| m[x][col].abs().partial_cmp(&m[y][col].abs()).unwrap())?;
        if m[piv][col].abs() < 1e-300 {
            return None;
        }
        m.swap(col, piv);
        let d = m[col][col];
        for k in 0..2 * n { m[col][k] /= d; }
        for r in 0..n {
            if r != col {
                let f = m[r][col];
                for k in 0..2 * n { m[r][k] -= f * m[col][k]; }
            }
        }
    }
    Some(m.into_iter().map(|r| r[n..].to_vec()).collect())
}

/// Correlated one-state cosh fit over times ts with inverse covariance w: (m, A, chi^2).
fn cosh_fit(c: &[f64], ts: &[usize], nt: usize, w: &[Vec<f64>]) -> (f64, f64, f64) {
    let basis = |m: f64| -> Vec<f64> { ts.iter().map(|&t| (-m * t as f64).exp() + (-m * (nt - t) as f64).exp()).collect() };
    let chi = |m: f64| -> (f64, f64) {
        let f = basis(m);
        let n = ts.len();
        let fwf: f64 = (0..n).map(|i| (0..n).map(|j| f[i] * w[i][j] * f[j]).sum::<f64>()).sum();
        let fwc: f64 = (0..n).map(|i| (0..n).map(|j| f[i] * w[i][j] * c[j]).sum::<f64>()).sum();
        let a = fwc / fwf;
        let r: Vec<f64> = (0..n).map(|i| c[i] - a * f[i]).collect();
        ((0..n).map(|i| (0..n).map(|j| r[i] * w[i][j] * r[j]).sum::<f64>()).sum(), a)
    };
    let mut best = (f64::INFINITY, 0.05);
    let mut m = 0.02;
    while m < 4.0 {
        let (x, _) = chi(m);
        if x < best.0 { best = (x, m) }
        m += 0.01;
    }
    let (mut lo, mut hi) = (best.1 - 0.01, best.1 + 0.01);
    for _ in 0..80 {
        let (a, b) = (lo + 0.382 * (hi - lo), lo + 0.618 * (hi - lo));
        if chi(a).0 < chi(b).0 { hi = b } else { lo = a }
    }
    let m = 0.5 * (lo + hi);
    let (x, a) = chi(m);
    (m, a, x)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 {
        eprintln!("usage: glueball_fit <T> <bins per replica> <replica.tsv>[:<replica.tsv>...]");
        std::process::exit(2);
    }
    let nt: usize = args[1].parse().unwrap();
    let bins_per: usize = args[2].parse().unwrap();
    let tmax = nt / 2;
    // bin sums of the primaries over every replica
    let mut bins: Vec<(Vec<f64>, usize)> = Vec::new();
    let mut total = 0usize;
    for path in args[3].split(':') {
        let cfgs = read_replica(path, nt);
        let size = cfgs.len() / bins_per;
        assert!(size >= 1, "{path}: fewer measurements than bins");
        for b in 0..bins_per {
            let mut s: Vec<f64> = Vec::new();
            for o in &cfgs[b * size..(b + 1) * size] {
                let p = primaries(o, nt, tmax);
                if s.is_empty() { s = p } else { for (x, y) in s.iter_mut().zip(p) { *x += y } }
            }
            bins.push((s, size));
            total += size;
        }
        eprintln!("{path}: {} measurements, {bins_per} bins of {size}", cfgs.len());
    }
    let dim = bins[0].0.len();
    let sum_all: Vec<f64> = (0..dim).map(|k| bins.iter().map(|b| b.0[k]).sum()).collect();
    let means = |skip: Option<usize>| -> Vec<f64> {
        let (s, n): (Vec<f64>, usize) = match skip {
            None => (sum_all.clone(), total),
            Some(b) => ((0..dim).map(|k| sum_all[k] - bins[b].0[k]).collect(), total - bins[b].1),
        };
        s.into_iter().map(|x| x / n as f64).collect()
    };
    let nb = bins.len();
    let full = matrices(&means(None), tmax);
    let (v, lambda) = variational_vector(&full[0], &full[1]);
    let jk: Vec<Vec<[[f64; LEVELS]; LEVELS]>> = (0..nb).map(|b| matrices(&means(Some(b)), tmax)).collect();
    let corr = |c: &Vec<[[f64; LEVELS]; LEVELS]>| -> Vec<f64> { c.iter().map(|m| project(m, &v)).collect() };
    let c_full = corr(&full);
    let c_jk: Vec<Vec<f64>> = jk.iter().map(corr).collect();
    let jk_err = |f: &dyn Fn(&[f64]) -> f64| -> (f64, f64) {
        let est = f(&c_full);
        let vals: Vec<f64> = c_jk.iter().map(|c| f(c)).collect();
        let m = vals.iter().sum::<f64>() / nb as f64;
        let var = vals.iter().map(|x| (x - m) * (x - m)).sum::<f64>() * (nb as f64 - 1.0) / nb as f64;
        (est, var.sqrt())
    };
    println!("T={nt} replicas={} bins={nb} measurements={total}", args[3].split(':').count());
    println!("variational vector (levels 0,4,8,16 APE; C(0)-normalized): [{}], lambda(0,1) = {lambda:.5} -> m = {:.4}",
        v.iter().map(|x| format!("{x:.4}")).collect::<Vec<_>>().join(", "), -lambda.ln());
    let mut line = String::from("effective mass   ");
    let mut line2 = String::from("cosh eff. mass   ");
    for t in 0..tmax.min(8) {
        let (m, e) = jk_err(&|c: &[f64]| -(c[t + 1] / c[t]).ln());
        line += &format!(" t={t}: {m:.3}({:.0})", e * 1000.0);
        let (m, e) = jk_err(&|c: &[f64]| cosh_mass(c[t + 1] / c[t], t, nt));
        line2 += &format!(" t={t}: {m:.3}({:.0})", e * 1000.0);
    }
    println!("{line}\n{line2}");
    println!("one-state cosh fits, correlated (jackknife covariance):");
    for tmin in 1..tmax {
        for tend in tmin + 1..=tmax.min(tmin + 4) {
            let ts: Vec<usize> = (tmin..=tend).collect();
            let n = ts.len();
            let mut cov = vec![vec![0.0; n]; n];
            let cbar: Vec<f64> = ts.iter().map(|&t| c_jk.iter().map(|c| c[t]).sum::<f64>() / nb as f64).collect();
            for c in &c_jk {
                for i in 0..n {
                    for j in 0..n {
                        cov[i][j] += (c[ts[i]] - cbar[i]) * (c[ts[j]] - cbar[j]);
                    }
                }
            }
            for r in cov.iter_mut() { for x in r.iter_mut() { *x *= (nb as f64 - 1.0) / nb as f64; } }
            let Some(w) = invert(&cov) else { println!("  [{tmin},{tend}]: singular covariance"); continue };
            let pick = |c: &[f64]| -> Vec<f64> { ts.iter().map(|&t| c[t]).collect() };
            let (m, _, chi2) = cosh_fit(&pick(&c_full), &ts, nt, &w);
            let ms: Vec<f64> = c_jk.iter().map(|c| cosh_fit(&pick(c), &ts, nt, &w).0).collect();
            let mbar = ms.iter().sum::<f64>() / nb as f64;
            let err = (ms.iter().map(|x| (x - mbar) * (x - mbar)).sum::<f64>() * (nb as f64 - 1.0) / nb as f64).sqrt();
            let dof = n - 2;
            let chi_s = if dof > 0 { format!("{:.2}", chi2 / dof as f64) } else { "-".into() };
            let flag = if m < 0.03 || m > 3.98 { "  (fit at the edge of the search range: noise)" } else { "" };
            println!("  [{tmin},{tend}]: am = {m:.4}({:.0})  chi2/dof = {chi_s} (dof {dof}){flag}", err * 10000.0);
        }
    }
}
