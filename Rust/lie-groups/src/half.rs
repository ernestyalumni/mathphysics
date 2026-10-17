//! Two-byte IEEE binary16 storage with round-to-nearest, ties-to-even conversion.
//! Arithmetic/transcendentals are SOFTWARE promoted to f32, then rounded back.
//! No hardware-half throughput or correctly-rounded transcendentals are claimed.
use std::{
    cmp::Ordering,
    fmt,
    ops::{Add, Div, Mul, Neg, Sub},
};
#[repr(transparent)]
#[derive(Clone, Copy, Default)]
pub struct Half(u16);
fn rounded(value: u64, shift: u32) -> u64 {
    let quotient = value >> shift;
    let rest = value & ((1_u64 << shift) - 1);
    let midpoint = 1_u64 << (shift - 1);
    quotient + u64::from(rest > midpoint || (rest == midpoint && quotient & 1 == 1))
}
impl Half {
    pub const fn from_bits(bits: u16) -> Self {
        Self(bits)
    }
    pub const fn to_bits(self) -> u16 {
        self.0
    }
    /// Direct f64 conversion avoids double rounding through f32.
    pub fn from_f64(value: f64) -> Self {
        let bits = value.to_bits();
        let sign = ((bits >> 48) & 0x8000) as u16;
        let exponent = ((bits >> 52) & 0x7ff) as i32;
        let fraction = bits & ((1_u64 << 52) - 1);
        if exponent == 0x7ff {
            return Self(
                sign | if fraction == 0 {
                    0x7c00
                } else {
                    0x7e00 | ((fraction >> 42) as u16 & 0x01ff)
                },
            );
        }
        if exponent == 0 {
            return Self(sign);
        }
        let mut e = exponent - 1023;
        if e > 15 {
            return Self(sign | 0x7c00);
        }
        if e < -25 {
            return Self(sign);
        }
        let mantissa = (1_u64 << 52) | fraction;
        if e < -14 {
            return Self(sign | rounded(mantissa, (28 - e) as u32) as u16);
        }
        let mut m = rounded(mantissa, 42);
        if m == 2048 {
            m = 1024;
            e += 1;
        }
        if e > 15 {
            Self(sign | 0x7c00)
        } else {
            Self(sign | (((e + 15) as u16) << 10) | (m as u16 - 1024))
        }
    }
    pub fn to_f32(self) -> f32 {
        let sign = ((self.0 as u32) & 0x8000) << 16;
        let e = (self.0 >> 10) & 31;
        let m = self.0 & 1023;
        if e == 0 {
            let value = (m as f32) * (1.0 / 16_777_216.0);
            if sign != 0 {
                -value
            } else {
                value
            }
        } else {
            let exponent = if e == 31 { 255 } else { u32::from(e) + 112 };
            f32::from_bits(sign | (exponent << 23) | (u32::from(m) << 13))
        }
    }
    pub fn to_f64(self) -> f64 {
        f64::from(self.to_f32())
    }
}
impl fmt::Debug for Half {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Half({})", self.to_f32())
    }
}
impl fmt::Display for Half {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.to_f32())
    }
}
impl PartialEq for Half {
    fn eq(&self, b: &Self) -> bool {
        self.to_f32() == b.to_f32()
    }
}
impl PartialOrd for Half {
    fn partial_cmp(&self, b: &Self) -> Option<Ordering> {
        self.to_f32().partial_cmp(&b.to_f32())
    }
}
macro_rules! operation {($trait:ident,$method:ident,$op:tt)=>{
    impl $trait for Half {type Output=Self;fn $method(self,b:Self)->Self{Self::from_f64(f64::from(self.to_f32() $op b.to_f32()))}}
};}
operation!(Add,add,+);
operation!(Sub,sub,-);
operation!(Mul,mul,*);
operation!(Div,div,/);
impl Neg for Half {
    type Output = Self;
    fn neg(self) -> Self {
        Self(self.0 ^ 0x8000)
    }
}
