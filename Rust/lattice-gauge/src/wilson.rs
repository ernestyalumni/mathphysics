//! SU(2) Wilson lattice gauge theory on a periodic torus: S = beta sum_p (1 - a0(U_p)) with
//! U_p = U(x,mu) U(x+mu,nu) U(x+nu,mu)^dagger U(x,nu)^dagger and a0 = (1/2) tr.
//!
//! A sweep updates, for each direction mu and each parity of the coordinate sum of x, all links
//! (x, mu) of that class at once. Two such links share no plaquette when every side length is even,
//! so the update is exact in parallel. Each link update draws from Stream::new(seed, sweep, link,
//! replica), so the result is independent of the number of threads.

use crate::philox::Stream;
use crate::su2::{self, Quaternion};
use crate::torus::Torus;

pub struct Lattice {
    pub torus: Torus,
    pub links: Vec<Quaternion>,
    pub beta: f64,
    pub seed: u64,
    pub replica: u32,
    pub sweep: u32,
    pub threads: usize,
    classes: Vec<Vec<usize>>, // link indices per (direction, parity)
    /// Optional coupling sign z_p = +-1 per plaquette (Torus::plaquette_index): S = beta sum_p (1 - z_p a0(U_p)).
    /// A centre twist sets z_p = -1 on one stack of plaquettes.
    pub plaquette_sign: Option<Vec<f64>>,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Update {
    HeatBath,
    /// Negative control: Kennedy-Pendleton without its accept step.
    HeatBathBroken,
    OverRelax,
    /// Cooling: U -> W^dagger/|W|, the link that maximizes the local a0(U W). Deterministic; used to
    /// expose instantons before measuring the topological charge.
    Cool,
}

impl Lattice {
    /// Cold start (all links 1). Every side must be even and at least 2.
    pub fn cold(side_lengths: &[usize], beta: f64, seed: u64, replica: u32) -> Result<Self, String> {
        if side_lengths.iter().any(|&n| n % 2 != 0) {
            return Err("every side length must be even for the checkerboard update".into());
        }
        let torus = Torus::new(side_lengths).ok_or("invalid torus")?;
        let d = torus.dimension();
        let mut classes = vec![Vec::new(); 2 * d];
        for site in 0..torus.site_count() {
            let parity = torus.site_coordinates(site).iter().sum::<usize>() % 2;
            for mu in 0..d {
                classes[2 * mu + parity].push(torus.link_index(site, mu));
            }
        }
        let threads = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
        let links = vec![Quaternion::IDENTITY; torus.link_count()];
        Ok(Lattice { torus, links, beta, seed, replica, sweep: 0, threads, classes, plaquette_sign: None })
    }

    /// Hot start: Haar-random links from the stream with sweep = u32::MAX.
    pub fn hot(side_lengths: &[usize], beta: f64, seed: u64, replica: u32) -> Result<Self, String> {
        let mut l = Lattice::cold(side_lengths, beta, seed, replica)?;
        for (i, u) in l.links.iter_mut().enumerate() {
            *u = su2::haar(&mut Stream::new(seed, u32::MAX, i as u32, replica));
        }
        Ok(l)
    }

    /// Link indices of each update class, in sweep order: class 2 mu + parity.
    pub fn update_classes(&self) -> &[Vec<usize>] {
        &self.classes
    }

    fn link(&self, site: usize, mu: usize) -> Quaternion {
        self.links[self.torus.link_index(site, mu)]
    }

    /// Staple sum W for link (site, mu): the local action is beta a0(U(site,mu) W).
    pub fn staple(&self, site: usize, mu: usize) -> Quaternion {
        let t = &self.torus;
        let x_mu = t.neighbor(site, mu, true);
        let mut w = Quaternion::ZERO;
        for nu in 0..t.dimension() {
            if nu == mu {
                continue;
            }
            let x_nu = t.neighbor(site, nu, true);
            let x_mnu = t.neighbor(site, nu, false);
            let x_mu_mnu = t.neighbor(x_mu, nu, false);
            // upper: U(x+mu,nu) U(x+nu,mu)^dagger U(x,nu)^dagger
            let mut up = self.link(x_mu, nu).mul(self.link(x_nu, mu).conj()).mul(self.link(site, nu).conj());
            // lower: U(x+mu-nu,nu)^dagger U(x-nu,mu)^dagger U(x-nu,nu)
            let mut down = self.link(x_mu_mnu, nu).conj().mul(self.link(x_mnu, mu).conj()).mul(self.link(x_mnu, nu));
            if let Some(z) = &self.plaquette_sign {
                let (a, b) = if mu < nu { (mu, nu) } else { (nu, mu) };
                up = up.scale(z[t.plaquette_index(site, a, b)]);
                down = down.scale(z[t.plaquette_index(x_mnu, a, b)]);
            }
            w = w.add(up).add(down);
        }
        w
    }

    /// One sweep of the given update over all links; increments the sweep counter.
    pub fn sweep_with(&mut self, update: Update) {
        let d = self.torus.dimension();
        for class in 0..2 * d {
            let indices = &self.classes[class];
            // Spawning threads costs more than updating a few thousand links; results do not depend on it.
            let threads = self.threads.max(1).min((indices.len() / 2048).max(1));
            let chunk = (indices.len() + threads - 1) / threads;
            let this = &*self;
            let new_values: Vec<(usize, Quaternion)> = std::thread::scope(|s| {
                let handles: Vec<_> = indices
                    .chunks(chunk.max(1))
                    .map(|part| {
                        s.spawn(move || {
                            part.iter()
                                .map(|&link| {
                                    let site = link / d;
                                    let mu = link % d;
                                    let w = this.staple(site, mu);
                                    let mut stream = Stream::new(this.seed, this.sweep, link as u32, this.replica);
                                    let u = match update {
                                        Update::HeatBath => su2::heat_bath(w, this.beta, &mut stream),
                                        Update::HeatBathBroken => su2::heat_bath_broken(w, this.beta, &mut stream),
                                        Update::OverRelax => su2::overrelax(this.links[link], w),
                                        Update::Cool => {
                                            let n = w.norm_sq();
                                            if n < 1e-300 { this.links[link] } else { w.conj().scale(1.0 / n.sqrt()) }
                                        }
                                    };
                                    (link, u)
                                })
                                .collect::<Vec<_>>()
                        })
                    })
                    .collect();
                handles.into_iter().flat_map(|h| h.join().expect("update thread panicked")).collect()
            });
            for (link, u) in new_values {
                self.links[link] = u;
            }
        }
        self.sweep = self.sweep.wrapping_add(1);
    }

    /// One heat-bath sweep followed by `overrelax` over-relaxation sweeps.
    pub fn update(&mut self, overrelax: usize) {
        self.sweep_with(Update::HeatBath);
        for _ in 0..overrelax {
            self.sweep_with(Update::OverRelax);
        }
    }

    /// Re-unitarize all links (guards against rounding drift in long runs).
    pub fn reunitarize(&mut self) {
        for u in self.links.iter_mut() {
            *u = u.normalized();
        }
    }

    /// Average of a0(U_p) over all plaquettes.
    pub fn plaquette(&self) -> f64 {
        let t = &self.torus;
        let d = t.dimension();
        let mut sum = 0.0;
        for site in 0..t.site_count() {
            for mu in 0..d {
                let x_mu = t.neighbor(site, mu, true);
                for nu in mu + 1..d {
                    let x_nu = t.neighbor(site, nu, true);
                    let p = self.link(site, mu).mul(self.link(x_mu, nu)).mul(self.link(x_nu, mu).conj()).mul(self.link(site, nu).conj());
                    sum += p.half_trace();
                }
            }
        }
        sum / t.plaquette_count() as f64
    }

    /// Average of z_p a0(U_p): (1/N_P) d log Z / d beta. Equals `plaquette` without signs.
    pub fn action_density(&self) -> f64 {
        let Some(z) = &self.plaquette_sign else { return self.plaquette() };
        let t = &self.torus;
        let d = t.dimension();
        let mut sum = 0.0;
        for site in 0..t.site_count() {
            for mu in 0..d {
                let x_mu = t.neighbor(site, mu, true);
                for nu in mu + 1..d {
                    let x_nu = t.neighbor(site, nu, true);
                    let p = self.link(site, mu).mul(self.link(x_mu, nu)).mul(self.link(x_nu, mu).conj()).mul(self.link(site, nu).conj());
                    sum += z[t.plaquette_index(site, mu, nu)] * p.half_trace();
                }
            }
        }
        sum / t.plaquette_count() as f64
    }

    /// 't Hooft centre twist: z_p = -1 for the plaquettes in plane (first, second) at the sites with
    /// coordinates x_first = 0 and x_second = 0 (one plaquette in every (first, second) plane).
    pub fn set_twist(&mut self, first: usize, second: usize) {
        self.set_twist_strength(first, second, -1.0);
    }

    /// Plaquettes of the twist stack (same set as `set_twist`).
    pub fn twist_stack(&self, first: usize, second: usize) -> Vec<usize> {
        let t = &self.torus;
        (0..t.site_count())
            .filter(|&s| t.site_coordinate(s, first) == 0 && t.site_coordinate(s, second) == 0)
            .map(|s| t.plaquette_index(s, first.min(second), first.max(second)))
            .collect()
    }

    /// Coupling z on the twist stack (z = 1: untwisted, z = -1: centre twist), 1 elsewhere. Interpolating z gives
    /// ln(Z_tw/Z) = -beta int_{-1}^{1} <sum_stack a0(U_p)>_z dz, an integrand over L^{d-2} plaquettes only.
    pub fn set_twist_strength(&mut self, first: usize, second: usize, z: f64) {
        let mut sign = vec![1.0; self.torus.plaquette_count()];
        for p in self.twist_stack(first, second) {
            sign[p] = z;
        }
        self.plaquette_sign = Some(sign);
    }

    /// Sum of a0(U_p) over the given plaquettes (by Torus::plaquette_index).
    pub fn plaquette_sum(&self, plaquettes: &[usize]) -> f64 {
        let t = &self.torus;
        let planes = t.plane_count();
        plaquettes
            .iter()
            .map(|&p| {
                let site = p / planes;
                let (mu, nu) = t.plane_directions(p % planes);
                let x_mu = t.neighbor(site, mu, true);
                let x_nu = t.neighbor(site, nu, true);
                self.link(site, mu).mul(self.link(x_mu, nu)).mul(self.link(x_nu, mu).conj()).mul(self.link(site, nu).conj()).half_trace()
            })
            .sum()
    }

    /// Product of links along `steps` steps from `site` in direction `mu` (forward) and the end site.
    fn line(&self, mut site: usize, mu: usize, steps: usize) -> (Quaternion, usize) {
        let mut q = Quaternion::IDENTITY;
        for _ in 0..steps {
            q = q.mul(self.link(site, mu));
            site = self.torus.neighbor(site, mu, true);
        }
        (q, site)
    }

    /// Rectangular Wilson loop <a0> of extent r in one direction and t in the other, averaged over all
    /// sites, all planes mu < nu, and both orientations (r along mu or along nu).
    pub fn wilson_loop(&self, r: usize, t: usize) -> f64 {
        let tor = &self.torus;
        let d = tor.dimension();
        let mut sum = 0.0;
        let mut count = 0usize;
        for site in 0..tor.site_count() {
            for mu in 0..d {
                for nu in mu + 1..d {
                    for &(a, b, la, lb) in &[(mu, nu, r, t), (nu, mu, r, t)] {
                        let (q1, x1) = self.line(site, a, la);
                        let (q2, _) = self.line(x1, b, lb);
                        let (q3, x3) = self.line(site, b, lb);
                        let (q4, _) = self.line(x3, a, la);
                        sum += q1.mul(q2).mul(q4.conj()).mul(q3.conj()).half_trace();
                        count += 1;
                    }
                }
            }
        }
        sum / count as f64
    }

    /// All rectangular Wilson loops W(r, t) for 1 <= r, t <= rmax, averaged as in `wilson_loop`;
    /// entry [r-1][t-1]. Straight-line products L_a(x, l) for l <= rmax are built once per site and
    /// direction, so each loop costs three products. Threads split the sites; the sum is reduced in a
    /// fixed order, so the result does not depend on the thread count.
    pub fn wilson_loop_table(&self, rmax: usize) -> Vec<Vec<f64>> {
        let tor = &self.torus;
        let d = tor.dimension();
        let n = tor.site_count();
        // lines[(site * d + a) * rmax + (l - 1)] = U(x,a) U(x+a,a) ... (l links)
        let mut lines = vec![Quaternion::IDENTITY; n * d * rmax];
        for site in 0..n {
            for a in 0..d {
                let mut q = Quaternion::IDENTITY;
                let mut y = site;
                for l in 0..rmax {
                    q = q.mul(self.link(y, a));
                    y = tor.neighbor(y, a, true);
                    lines[(site * d + a) * rmax + l] = q;
                }
            }
        }
        let shift = |mut y: usize, a: usize, l: usize| {
            for _ in 0..l {
                y = tor.neighbor(y, a, true);
            }
            y
        };
        // Fixed-size chunks summed in chunk order: the result is independent of the thread count.
        const CHUNK: usize = 512;
        let chunks = (n + CHUNK - 1) / CHUNK;
        let threads = self.threads.max(1).min(chunks);
        let lines = &lines;
        let mut partials: Vec<(usize, Vec<f64>)> = std::thread::scope(|s| {
            let handles: Vec<_> = (0..threads)
                .map(|tid| {
                    s.spawn(move || {
                        let mut out = Vec::new();
                        for c in (tid..chunks).step_by(threads) {
                        let mut acc = vec![0.0; rmax * rmax];
                        for site in c * CHUNK..((c + 1) * CHUNK).min(n) {
                            for mu in 0..d {
                                for nu in mu + 1..d {
                                    for &(a, b) in &[(mu, nu), (nu, mu)] {
                                        for r in 1..=rmax {
                                            let xr = shift(site, a, r);
                                            let la = lines[(site * d + a) * rmax + r - 1];
                                            for t in 1..=rmax {
                                                let xt = shift(site, b, t);
                                                let lb_end = lines[(xr * d + b) * rmax + t - 1];
                                                let la_top = lines[(xt * d + a) * rmax + r - 1];
                                                let lb = lines[(site * d + b) * rmax + t - 1];
                                                acc[(r - 1) * rmax + t - 1] += la.mul(lb_end).mul(la_top.conj()).mul(lb.conj()).half_trace();
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        out.push((c, acc));
                        }
                        out
                    })
                })
                .collect();
            handles.into_iter().flat_map(|h| h.join().expect("loop thread panicked")).collect()
        });
        partials.sort_by_key(|p| p.0);
        let count = (n * d * (d - 1)) as f64; // sites x planes x 2 orientations
        (0..rmax)
            .map(|r| (0..rmax).map(|t| partials.iter().map(|p| p.1[r * rmax + t]).sum::<f64>() / count).collect())
            .collect()
    }

    /// Gauge transformation U(x,mu) -> g(x) U(x,mu) g(x+mu)^dagger (for invariance tests).
    pub fn gauge_transform(&mut self, g: &[Quaternion]) {
        let d = self.torus.dimension();
        for site in 0..self.torus.site_count() {
            for mu in 0..d {
                let i = self.torus.link_index(site, mu);
                let y = self.torus.neighbor(site, mu, true);
                self.links[i] = g[site].mul(self.links[i]).mul(g[y].conj());
            }
        }
    }
}

/// Creutz ratio chi(R,T) = -ln( W(R,T) W(R-1,T-1) / (W(R,T-1) W(R-1,T)) ) from a function W.
pub fn creutz_ratio(w: impl Fn(usize, usize) -> f64, r: usize, t: usize) -> f64 {
    -((w(r, t) * w(r - 1, t - 1)) / (w(r, t - 1) * w(r - 1, t))).ln()
}
