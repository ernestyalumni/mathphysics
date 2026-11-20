//! Philox4x32-10 counter-based random numbers (Salmon, Moraes, Dror, Shaw, SC'11), bit-identical
//! to Random123's `philox4x32` with 10 rounds.
//!
//! Stream convention for this crate, which GPU code must follow to reproduce a run exactly:
//! key = (seed low 32 bits, seed high 32 bits); counter = (block, sweep, link index, replica).
//! Each block yields four u32 in order; `block` then increments. A link update draws only from
//! its own stream, so the result does not depend on thread count or update order within a parity.

const M0: u32 = 0xD251_1F53;
const M1: u32 = 0xCD9E_8D57;
const W0: u32 = 0x9E37_79B9;
const W1: u32 = 0xBB67_AE85;

fn mulhilo(a: u32, b: u32) -> (u32, u32) {
    let p = (a as u64) * (b as u64);
    ((p >> 32) as u32, p as u32)
}

/// One Philox4x32-10 block.
pub fn philox4x32(counter: [u32; 4], key: [u32; 2]) -> [u32; 4] {
    let mut c = counter;
    let mut k = key;
    for round in 0..10 {
        if round > 0 {
            k[0] = k[0].wrapping_add(W0);
            k[1] = k[1].wrapping_add(W1);
        }
        let (hi0, lo0) = mulhilo(M0, c[0]);
        let (hi1, lo1) = mulhilo(M1, c[2]);
        c = [hi1 ^ c[1] ^ k[0], lo1, hi0 ^ c[3] ^ k[1], lo0];
    }
    c
}

/// A reproducible stream identified by (seed, sweep, link, replica).
#[derive(Clone, Debug)]
pub struct Stream {
    key: [u32; 2],
    counter: [u32; 4],
    buffer: [u32; 4],
    used: usize,
}

impl Stream {
    pub fn new(seed: u64, sweep: u32, link: u32, replica: u32) -> Self {
        Stream { key: [seed as u32, (seed >> 32) as u32], counter: [0, sweep, link, replica], buffer: [0; 4], used: 4 }
    }

    pub fn next_u32(&mut self) -> u32 {
        if self.used == 4 {
            self.buffer = philox4x32(self.counter, self.key);
            self.counter[0] = self.counter[0].wrapping_add(1);
            self.used = 0;
        }
        let x = self.buffer[self.used];
        self.used += 1;
        x
    }

    /// Uniform in the open interval (0, 1), with 53 random bits; never 0 or 1.
    pub fn next_f64_open01(&mut self) -> f64 {
        let hi = self.next_u32() as u64;
        let lo = self.next_u32() as u64;
        let bits = ((hi << 32) | lo) >> 11;
        (bits as f64 + 0.5) * (1.0 / (1u64 << 53) as f64)
    }

    /// Standard normal by Box-Muller (one value per call; the second is discarded for simplicity).
    pub fn next_gaussian(&mut self) -> f64 {
        let r1 = self.next_f64_open01();
        let r2 = self.next_f64_open01();
        (-2.0 * r1.ln()).sqrt() * (2.0 * std::f64::consts::PI * r2).cos()
    }
}
