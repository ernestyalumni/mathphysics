//! Compact U(1) Wilson lattice gauge theory, the counter-model of the program: S = beta sum_p (1 - cos theta_p),
//! theta_p = theta(x,mu) + theta(x+mu,nu) - theta(x+nu,mu) - theta(x,nu). In 4D it has a massless Coulomb phase
//! above beta ~ 1.0, so any gap argument that also works here is wrong.
//!
//! Same conventions as `wilson`: links by Torus::link_index, checkerboard classes (direction, parity), one Philox
//! stream per link update keyed by (seed, sweep, link, replica). With the staple sum A = sum_j e^{i alpha_j} = k e^{i phi},
//! the local weight of theta is exp(beta k cos(theta + phi)).

use crate::philox::Stream;
use crate::torus::Torus;
use std::f64::consts::PI;

/// Von Mises sample on (-pi, pi] with density proportional to exp(kappa cos x) (Best and Fisher, Appl. Stat. 28 (1979) 152).
pub fn von_mises(kappa: f64, s: &mut Stream) -> f64 {
    if kappa < 1e-8 {
        return PI * (2.0 * s.next_f64_open01() - 1.0);
    }
    let tau = 1.0 + (1.0 + 4.0 * kappa * kappa).sqrt();
    let rho = (tau - (2.0 * tau).sqrt()) / (2.0 * kappa);
    let r = (1.0 + rho * rho) / (2.0 * rho);
    loop {
        let u1 = s.next_f64_open01();
        let u2 = s.next_f64_open01();
        let u3 = s.next_f64_open01();
        let z = (PI * u1).cos();
        let f = (1.0 + r * z) / (r + z);
        let c = kappa * (r - f);
        if c * (2.0 - c) - u2 > 0.0 || (c / u2).ln() + 1.0 - c >= 0.0 {
            let x = f.clamp(-1.0, 1.0).acos();
            return if u3 > 0.5 { x } else { -x };
        }
    }
}

fn wrap(x: f64) -> f64 {
    let y = x.rem_euclid(2.0 * PI);
    if y > PI { y - 2.0 * PI } else { y }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum U1Update {
    HeatBath,
    OverRelax,
}

pub struct U1Lattice {
    pub torus: Torus,
    pub theta: Vec<f64>,
    pub beta: f64,
    pub seed: u64,
    pub replica: u32,
    pub sweep: u32,
    pub threads: usize,
    classes: Vec<Vec<usize>>,
}

impl U1Lattice {
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
        let theta = vec![0.0; torus.link_count()];
        Ok(U1Lattice { torus, theta, beta, seed, replica, sweep: 0, threads, classes })
    }

    pub fn hot(side_lengths: &[usize], beta: f64, seed: u64, replica: u32) -> Result<Self, String> {
        let mut l = U1Lattice::cold(side_lengths, beta, seed, replica)?;
        for (i, t) in l.theta.iter_mut().enumerate() {
            *t = PI * (2.0 * Stream::new(seed, u32::MAX, i as u32, replica).next_f64_open01() - 1.0);
        }
        Ok(l)
    }

    fn th(&self, site: usize, mu: usize) -> f64 {
        self.theta[self.torus.link_index(site, mu)]
    }

    /// Staple sum A = sum_j e^{i alpha_j} with local action beta Re(e^{i theta} A) = beta k cos(theta + phi).
    pub fn staple(&self, site: usize, mu: usize) -> (f64, f64) {
        let t = &self.torus;
        let x_mu = t.neighbor(site, mu, true);
        let (mut re, mut im) = (0.0, 0.0);
        for nu in 0..t.dimension() {
            if nu == mu {
                continue;
            }
            let x_nu = t.neighbor(site, nu, true);
            let x_mnu = t.neighbor(site, nu, false);
            let x_mu_mnu = t.neighbor(x_mu, nu, false);
            // plaquette (x, mu, nu): theta(x,mu) + theta(x+mu,nu) - theta(x+nu,mu) - theta(x,nu)
            let a_up = self.th(x_mu, nu) - self.th(x_nu, mu) - self.th(site, nu);
            // plaquette (x-nu, mu, nu): theta(x-nu,mu) + theta(x+mu-nu,nu) - theta(x,mu) - theta(x-nu,nu); as a function
            // of theta(x,mu) this is cos(theta(x,mu) - (theta(x-nu,mu) + theta(x+mu-nu,nu) - theta(x-nu,nu)))
            let a_down = -(self.th(x_mnu, mu) + self.th(x_mu_mnu, nu) - self.th(x_mnu, nu));
            re += a_up.cos() + a_down.cos();
            im += a_up.sin() + a_down.sin();
        }
        (re, im)
    }

    pub fn sweep_with(&mut self, update: U1Update) {
        let d = self.torus.dimension();
        for class in 0..2 * d {
            let indices = &self.classes[class];
            let threads = self.threads.max(1).min((indices.len() / 2048).max(1));
            let chunk = (indices.len() + threads - 1) / threads;
            let this = &*self;
            let new_values: Vec<(usize, f64)> = std::thread::scope(|s| {
                let handles: Vec<_> = indices
                    .chunks(chunk.max(1))
                    .map(|part| {
                        s.spawn(move || {
                            part.iter()
                                .map(|&link| {
                                    let (re, im) = this.staple(link / d, link % d);
                                    let k = (re * re + im * im).sqrt();
                                    let phi = im.atan2(re);
                                    let t = match update {
                                        U1Update::HeatBath => {
                                            let mut st = Stream::new(this.seed, this.sweep, link as u32, this.replica);
                                            wrap(von_mises(this.beta * k, &mut st) - phi)
                                        }
                                        U1Update::OverRelax => wrap(-2.0 * phi - this.theta[link]),
                                    };
                                    (link, t)
                                })
                                .collect::<Vec<_>>()
                        })
                    })
                    .collect();
                handles.into_iter().flat_map(|h| h.join().expect("update thread panicked")).collect()
            });
            for (link, t) in new_values {
                self.theta[link] = t;
            }
        }
        self.sweep = self.sweep.wrapping_add(1);
    }

    pub fn update(&mut self, overrelax: usize) {
        self.sweep_with(U1Update::HeatBath);
        for _ in 0..overrelax {
            self.sweep_with(U1Update::OverRelax);
        }
    }

    /// Plaquette angle theta_p for plaquette (site, mu < nu).
    pub fn plaquette_angle(&self, site: usize, mu: usize, nu: usize) -> f64 {
        let t = &self.torus;
        self.th(site, mu) + self.th(t.neighbor(site, mu, true), nu) - self.th(t.neighbor(site, nu, true), mu) - self.th(site, nu)
    }

    /// Average of cos theta_p over all plaquettes.
    pub fn plaquette(&self) -> f64 {
        let t = &self.torus;
        let d = t.dimension();
        let mut sum = 0.0;
        for site in 0..t.site_count() {
            for mu in 0..d {
                for nu in mu + 1..d {
                    sum += self.plaquette_angle(site, mu, nu).cos();
                }
            }
        }
        sum / t.plaquette_count() as f64
    }

    /// Gauge transformation theta(x,mu) -> theta(x,mu) + g(x) - g(x+mu).
    pub fn gauge_transform(&mut self, g: &[f64]) {
        let d = self.torus.dimension();
        for site in 0..self.torus.site_count() {
            for mu in 0..d {
                let i = self.torus.link_index(site, mu);
                let y = self.torus.neighbor(site, mu, true);
                self.theta[i] = wrap(self.theta[i] + g[site] - g[y]);
            }
        }
    }
}

/// Exact 2D U(1) torus: e^{beta cos x} = sum_n I_n(beta) e^{i n x}, so Z = sum_{n in Z} I_n(beta)^V and
/// <cos theta_p> = sum_n I_n^{V-1} I_n' / Z with I_n' = (I_{n-1} + I_{n+1})/2 (I_{-n} = I_n).
pub fn torus_plaquette_u1_2d(beta: f64, plaquettes: u64) -> f64 {
    use crate::character::bessel_i_scaled;
    let i0 = bessel_i_scaled(0, beta);
    let v = plaquettes as f64;
    let (mut num, mut den) = (0.0, 0.0);
    let mut n = 0u32;
    loop {
        let i_n = bessel_i_scaled(n, beta);
        let w = (i_n / i0).powf(v);
        if n > 0 && (i_n == 0.0 || w < 1e-18 * den) {
            break;
        }
        let i_prev = bessel_i_scaled(if n == 0 { 1 } else { n - 1 }, beta);
        let deriv_ratio = (i_prev + bessel_i_scaled(n + 1, beta)) / (2.0 * i_n);
        let mult = if n == 0 { 1.0 } else { 2.0 }; // n and -n
        num += mult * w * deriv_ratio;
        den += mult * w;
        n += 1;
        if n > 2000 {
            break;
        }
    }
    num / den
}
