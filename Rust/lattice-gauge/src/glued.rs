//! Block gluing: an estimator of the hypercubic defect Delta_beta(n) = 2^d log Z(n^d) - log Z((2n)^d) by integration
//! over a gluing parameter lambda, in a single family of ensembles.
//!
//! The (2n)^d torus is cut into 2^d blocks of side n. Every plaquette that crosses a block boundary appears twice:
//! glued (the ordinary plaquette, weight lambda) and blocked (the plaquette closed by wrapping within its block,
//! weight 1 - lambda). Both use the same link variables: a link leaving a block either enters the next block (glued) or
//! wraps to the start of its own block (blocked). At lambda = 1 the system is the (2n)^d torus, at lambda = 0 it is 2^d
//! independent n^d tori, so with S_lambda = -beta sum_p w_p a0(U_p),
//!   Delta = log Z(0) - log Z(1) = -beta int_0^1 < sum_boundary [a0(glued) - a0(blocked)] >_lambda dlambda.
//! The integrand is a difference measured in one ensemble, over boundary plaquettes only.
//!
//! Plaquettes are general: four links with orientation signs and a weight; staples come from incidence lists. With n even,
//! the checkerboard classes (direction, site parity) of `wilson` stay independent (a wrapping plaquette joins sites n - 1
//! apart, of opposite parity), and updates use the same Philox streams.

use crate::philox::Stream;
use crate::su2::{self, Quaternion};
use crate::torus::Torus;

#[derive(Clone, Copy, Debug)]
pub struct Plaquette {
    pub links: [usize; 4],
    pub forward: [bool; 4], // false: the link enters as its conjugate
}

pub struct GluedLattice {
    pub torus: Torus,
    pub block: usize,
    pub links: Vec<Quaternion>,
    pub beta: f64,
    pub lambda: f64,
    pub seed: u64,
    pub sweep: u32,
    pub plaquettes: Vec<Plaquette>,
    /// 0 = bulk (weight 1), 1 = glued boundary (weight lambda), 2 = blocked boundary (weight 1 - lambda)
    pub kind: Vec<u8>,
    incidence: Vec<Vec<(usize, usize)>>, // link -> (plaquette, position)
    classes: Vec<Vec<usize>>,
}

impl GluedLattice {
    /// Torus of side 2*block in each of d directions; block must be even and >= 2.
    pub fn hot(d: usize, block: usize, beta: f64, lambda: f64, seed: u64) -> Result<Self, String> {
        if block % 2 != 0 || block < 2 {
            return Err("block side must be even and at least 2".into());
        }
        let torus = Torus::new(&vec![2 * block; d]).ok_or("invalid torus")?;
        let n = block;
        let inblock = |s: usize, dir: usize| -> usize {
            // forward neighbour within the block (wraps at the block edge)
            let c = torus.site_coordinate(s, dir);
            if c % n == n - 1 {
                let mut t = s;
                for _ in 0..n - 1 { t = torus.neighbor(t, dir, false); }
                t
            } else {
                torus.neighbor(s, dir, true)
            }
        };
        let mut plaquettes = Vec::new();
        let mut kind = Vec::new();
        for x in 0..torus.site_count() {
            for mu in 0..d {
                for nu in mu + 1..d {
                    let crosses = torus.site_coordinate(x, mu) % n == n - 1 || torus.site_coordinate(x, nu) % n == n - 1;
                    let mk = |fwd: &dyn Fn(usize, usize) -> usize| Plaquette {
                        links: [torus.link_index(x, mu), torus.link_index(fwd(x, mu), nu), torus.link_index(fwd(x, nu), mu), torus.link_index(x, nu)],
                        forward: [true, true, false, false],
                    };
                    let glued = mk(&|s, dir| torus.neighbor(s, dir, true));
                    if crosses {
                        plaquettes.push(glued);
                        kind.push(1);
                        plaquettes.push(mk(&inblock));
                        kind.push(2);
                    } else {
                        plaquettes.push(glued);
                        kind.push(0);
                    }
                }
            }
        }
        let mut incidence = vec![Vec::new(); torus.link_count()];
        for (p, q) in plaquettes.iter().enumerate() {
            for (i, &l) in q.links.iter().enumerate() {
                incidence[l].push((p, i));
            }
        }
        let mut classes = vec![Vec::new(); 2 * d];
        for s in 0..torus.site_count() {
            let parity = torus.site_coordinates(s).iter().sum::<usize>() % 2;
            for mu in 0..d {
                classes[2 * mu + parity].push(torus.link_index(s, mu));
            }
        }
        let links = (0..torus.link_count()).map(|i| su2::haar(&mut Stream::new(seed, u32::MAX, i as u32, 0))).collect();
        Ok(GluedLattice { torus, block, links, beta, lambda, seed, sweep: 0, plaquettes, kind, incidence, classes })
    }

    fn weight(&self, p: usize) -> f64 {
        match self.kind[p] { 0 => 1.0, 1 => self.lambda, _ => 1.0 - self.lambda }
    }

    /// a0 of plaquette p.
    pub fn plaquette_value(&self, p: usize) -> f64 {
        let q = &self.plaquettes[p];
        let mut u = Quaternion::IDENTITY;
        for i in 0..4 {
            let l = self.links[q.links[i]];
            u = u.mul(if q.forward[i] { l } else { l.conj() });
        }
        u.half_trace()
    }

    /// Staple W of a link: the local action is beta a0(U W).
    pub fn staple(&self, link: usize) -> Quaternion {
        let mut w = Quaternion::ZERO;
        for &(p, i) in &self.incidence[link] {
            let wt = self.weight(p);
            if wt == 0.0 { continue; }
            let q = &self.plaquettes[p];
            // C = product of the other three links in cyclic order after position i
            let mut c = Quaternion::IDENTITY;
            for k in 1..4 {
                let j = (i + k) % 4;
                let l = self.links[q.links[j]];
                c = c.mul(if q.forward[j] { l } else { l.conj() });
            }
            // U_p = U_l^{s} C up to cyclic order; a0(U_l C) for s = +1, a0(U_l^dagger C) = a0(U_l C^dagger) for s = -1
            let contrib = if q.forward[i] { c } else { c.conj() };
            w = w.add(contrib.scale(wt));
        }
        w
    }

    /// One heat-bath sweep followed by `overrelax` over-relaxation sweeps (sequential within each class).
    pub fn update(&mut self, overrelax: usize) {
        for pass in 0..=overrelax {
            for class in 0..self.classes.len() {
                for idx in 0..self.classes[class].len() {
                    let l = self.classes[class][idx];
                    let w = self.staple(l);
                    self.links[l] = if pass == 0 {
                        su2::heat_bath(w, self.beta, &mut Stream::new(self.seed, self.sweep, l as u32, 0))
                    } else {
                        su2::overrelax(self.links[l], w)
                    };
                }
            }
            self.sweep = self.sweep.wrapping_add(1);
        }
    }

    /// The integrand sum over boundary plaquettes of a0(glued) - a0(blocked).
    pub fn boundary_difference(&self) -> f64 {
        (0..self.plaquettes.len())
            .map(|p| match self.kind[p] { 1 => self.plaquette_value(p), 2 => -self.plaquette_value(p), _ => 0.0 })
            .sum()
    }
}
