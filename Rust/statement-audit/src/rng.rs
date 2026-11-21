//! Deterministic SplitMix64 and generators of small exact numbers.
use crate::num::{self, q, qi};
use crate::value::{real, Mat, V};
use qft_symbolics::exact::BigQ;
use std::rc::Rc;

pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Rng {
        Rng(seed ^ 0x9E37_79B9_7F4A_7C15)
    }
    pub fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
    pub fn below(&mut self, n: u64) -> u64 {
        self.next() % n.max(1)
    }
    pub fn coin(&mut self) -> bool {
        self.below(2) == 0
    }
    /// A rational in [0, 1], with 0, 1 and simple fractions over-represented.
    pub fn unit_rational(&mut self) -> BigQ {
        match self.below(8) {
            0 => qi(0),
            1 => qi(1),
            2 => q(1, 2),
            _ => {
                let d = 1 + self.below(12) as i64;
                q(self.below(d as u64 + 1) as i64, d)
            }
        }
    }
    /// A rational in (0, 1].
    pub fn unit_rational_pos(&mut self) -> BigQ {
        let d = 1 + self.below(12) as i64;
        q(1 + self.below(d as u64) as i64, d)
    }
    /// A positive rational: simple fractions, occasionally tiny or large.
    pub fn pos_rational(&mut self) -> BigQ {
        match self.below(10) {
            0 => q(1, 1 << (1 + self.below(30))),
            1 => qi(10 + self.below(90) as i64),
            _ => q(1 + self.below(9) as i64, 1 + self.below(6) as i64),
        }
    }
    /// A rational of either sign, with special values over-represented.
    pub fn rational(&mut self) -> BigQ {
        match self.below(12) {
            0 => qi(0),
            1 => qi(1),
            2 => q(1, 4),
            3 => q(693_147, 1_000_000), // just below log 2
            4 => qi(2),
            5 => num::neg(&self.pos_rational()),
            _ => self.pos_rational(),
        }
    }
    pub fn real(&mut self) -> V {
        real(self.rational())
    }
    pub fn small_int(&mut self) -> i64 {
        self.below(7) as i64 - 3
    }
    pub fn matrix(&mut self, n: usize, complex: bool) -> V {
        let e = (0..n * n)
            .map(|_| {
                if complex {
                    V::Cx(qi(self.small_int()), qi(self.small_int()))
                } else if self.below(5) == 0 {
                    real(q(self.small_int(), 2))
                } else {
                    real(qi(self.small_int()))
                }
            })
            .collect();
        V::Mat(Rc::new(Mat { rows: n, cols: n, e }))
    }
    /// A finitely supported sequence of reals.
    pub fn sequence(&mut self, len: usize, nonnegative: bool) -> V {
        let e = (0..len)
            .map(|_| {
                let x = self.rational();
                if nonnegative { real(num::abs(&x)) } else { real(x) }
            })
            .collect();
        V::Seq(Rc::new(e))
    }
}
