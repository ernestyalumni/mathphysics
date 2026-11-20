//! Zero-momentum 0++ glueball operators and correlators: the lattice mass gap.
//!
//! One direction is Euclidean time (`time_dir`). Spatial links are APE-smeared within each time
//! slice, U' = P[(1 - alpha) U + (alpha / (2(d-2))) sum of spatial staples], where for SU(2) the
//! projection P of a sum of SU(2) elements is division by its norm. The operator on slice t is the
//! average of a0 over spatial plaquettes of the smeared links. From the time series of these slice
//! operators, C(t) = <O(t0) O(t0 + t)> - <O>^2 averaged over t0, and the effective mass
//! m_eff(t) = -ln(C(t+1)/C(t)), which bounds the lowest 0++ mass from above as t grows.

use crate::su2::Quaternion;
use crate::torus::Torus;

/// One APE smearing step of all spatial links (links in `time_dir` are left unchanged).
pub fn ape_step(torus: &Torus, links: &[Quaternion], time_dir: usize, alpha: f64) -> Vec<Quaternion> {
    let d = torus.dimension();
    let n_spatial = (d - 1) as f64;
    let mut out = links.to_vec();
    let link = |s: usize, m: usize| links[torus.link_index(s, m)];
    for site in 0..torus.site_count() {
        for mu in 0..d {
            if mu == time_dir {
                continue;
            }
            let x_mu = torus.neighbor(site, mu, true);
            let mut w = Quaternion::ZERO;
            for nu in 0..d {
                if nu == mu || nu == time_dir {
                    continue;
                }
                let x_nu = torus.neighbor(site, nu, true);
                let x_mnu = torus.neighbor(site, nu, false);
                let x_mu_mnu = torus.neighbor(x_mu, nu, false);
                // paths from x to x+mu: up nu, along mu, down nu; and down, along, up
                let up = link(site, nu).mul(link(x_nu, mu)).mul(link(x_mu, nu).conj());
                let down = link(x_mnu, nu).conj().mul(link(x_mnu, mu)).mul(link(x_mu_mnu, nu));
                w = w.add(up).add(down);
            }
            let i = torus.link_index(site, mu);
            let s = links[i].scale(1.0 - alpha).add(w.scale(alpha / (2.0 * (n_spatial - 1.0))));
            out[i] = s.normalized();
        }
    }
    out
}

/// Slice operators O(t): average of a0 over spatial plaquettes in each time slice.
pub fn slice_plaquettes(torus: &Torus, links: &[Quaternion], time_dir: usize) -> Vec<f64> {
    let d = torus.dimension();
    let nt = torus.side_lengths()[time_dir];
    let mut sum = vec![0.0; nt];
    let mut count = vec![0usize; nt];
    let link = |s: usize, m: usize| links[torus.link_index(s, m)];
    for site in 0..torus.site_count() {
        let t = torus.site_coordinate(site, time_dir);
        for mu in 0..d {
            if mu == time_dir {
                continue;
            }
            let x_mu = torus.neighbor(site, mu, true);
            for nu in mu + 1..d {
                if nu == time_dir {
                    continue;
                }
                let x_nu = torus.neighbor(site, nu, true);
                let p = link(site, mu).mul(link(x_mu, nu)).mul(link(x_nu, mu).conj()).mul(link(site, nu).conj());
                sum[t] += p.half_trace();
                count[t] += 1;
            }
        }
    }
    sum.iter().zip(count).map(|(s, c)| s / c as f64).collect()
}

/// Slice operators at several smearing levels: entry [level][t], with `steps[level]` cumulative APE steps.
pub fn smeared_slice_operators(torus: &Torus, links: &[Quaternion], time_dir: usize, alpha: f64, steps: &[usize]) -> Vec<Vec<f64>> {
    let mut current = links.to_vec();
    let mut done = 0;
    let mut out = Vec::new();
    for &s in steps {
        while done < s {
            current = ape_step(torus, &current, time_dir, alpha);
            done += 1;
        }
        out.push(slice_plaquettes(torus, &current, time_dir));
    }
    out
}

/// Connected correlator C(t), t = 0..nt/2, from per-configuration slice operators (configs x nt),
/// given the mean of O over all configurations and slices.
pub fn correlator(series: &[Vec<f64>], mean: f64) -> Vec<f64> {
    let nt = series[0].len();
    let mut c = vec![0.0; nt / 2 + 1];
    for o in series {
        for (t, ct) in c.iter_mut().enumerate() {
            let mut s = 0.0;
            for t0 in 0..nt {
                s += (o[t0] - mean) * (o[(t0 + t) % nt] - mean);
            }
            *ct += s / nt as f64;
        }
    }
    c.iter().map(|x| x / series.len() as f64).collect()
}
