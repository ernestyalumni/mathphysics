//! Exact averages of monomials over unit spheres, and rational multiples of powers of pi.
//!
//! Conventions: the sphere S^(d-1) is the unit sphere of R^d with the rotation-invariant
//! probability measure. Angular averages are exact rationals. Areas and momentum-shell
//! coefficients are exact rationals times an integer power of pi ([`PiMonomial`]).
//! Nothing here is floating point.
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero};

/// The number `coefficient * pi^pi_exponent` with an exact rational coefficient.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PiMonomial {
    pub coefficient: BigRational,
    pub pi_exponent: i32,
}

/// The exact rational `numerator / denominator`, or `None` when the denominator is zero.
pub fn exact_rational(numerator: i64, denominator: i64) -> Option<BigRational> {
    if denominator == 0 {
        return None;
    }
    Some(BigRational::new(
        BigInt::from(numerator),
        BigInt::from(denominator),
    ))
}

impl PiMonomial {
    pub fn new(coefficient: BigRational, pi_exponent: i32) -> Self {
        Self {
            coefficient,
            pi_exponent,
        }
    }
    pub fn multiply(&self, other: &Self) -> Self {
        Self::new(
            &self.coefficient * &other.coefficient,
            self.pi_exponent + other.pi_exponent,
        )
    }
    pub fn scale(&self, factor: &BigRational) -> Self {
        Self::new(&self.coefficient * factor, self.pi_exponent)
    }
    /// `None` when the divisor is zero.
    pub fn divide(&self, other: &Self) -> Option<Self> {
        if other.coefficient.is_zero() {
            return None;
        }
        Some(Self::new(
            &self.coefficient / &other.coefficient,
            self.pi_exponent - other.pi_exponent,
        ))
    }
}

/// Average of `n_1^e_1 * ... * n_d^e_d` over the unit sphere of R^d, where
/// `d = exponents.len()` and `e_i = exponents[i]`.
///
/// Formula: zero if some exponent is odd; otherwise, with `e_i = 2 a_i` and `k = sum a_i`,
/// `prod_i (2 a_i - 1)!! / (d (d+2) ... (d+2k-2))`. Returns `None` for `d = 0`.
pub fn sphere_monomial_average(exponents: &[u32]) -> Option<BigRational> {
    let dimension = exponents.len();
    if dimension == 0 {
        return None;
    }
    if exponents.iter().any(|exponent| exponent % 2 == 1) {
        return Some(BigRational::zero());
    }
    let mut numerator = BigInt::one();
    let mut half_degree: u64 = 0;
    for exponent in exponents {
        let half = u64::from(exponent / 2);
        half_degree += half;
        let mut odd_factor: u64 = 1;
        while odd_factor < 2 * half {
            numerator *= BigInt::from(odd_factor);
            odd_factor += 2;
        }
    }
    let mut denominator = BigInt::one();
    for step in 0..half_degree {
        denominator *= BigInt::from(dimension as u64 + 2 * step);
    }
    Some(BigRational::new(numerator, denominator))
}

fn factorial(value: u64) -> BigInt {
    (1..=value).fold(BigInt::one(), |product, factor| product * BigInt::from(factor))
}

/// Area of the unit sphere S^(d-1) in R^d, `2 pi^(d/2) / Gamma(d/2)`, for `d >= 1`.
///
/// Even `d = 2m`: `2 pi^m / (m-1)!`. Odd `d = 2m+1`: `2 * 4^m * m! * pi^m / (2m)!`.
/// Examples: `d = 1` gives 2 (two points), `d = 2` gives `2 pi`, `d = 3` gives `4 pi`,
/// `d = 4` gives `2 pi^2`.
pub fn unit_sphere_area(dimension: u32) -> Option<PiMonomial> {
    if dimension == 0 {
        return None;
    }
    let half = u64::from(dimension / 2);
    let coefficient = if dimension % 2 == 0 {
        BigRational::new(BigInt::from(2), factorial(half - 1))
    } else {
        BigRational::new(
            BigInt::from(2) * BigInt::from(4).pow(half as u32) * factorial(half),
            factorial(2 * half),
        )
    };
    Some(PiMonomial::new(coefficient, half as i32))
}

/// Coefficient `c` in `int_{1/L <= |k| <= 1} d^dk / ((2 pi)^d |k|^d) = c * log L`:
/// the area of S^(d-1) divided by `(2 pi)^d`. For `d = 4` it is `1 / (8 pi^2)`.
pub fn logarithmic_shell_coefficient(dimension: u32) -> Option<PiMonomial> {
    let area = unit_sphere_area(dimension)?;
    let two_to_dimension = BigRational::from_integer(BigInt::from(2).pow(dimension));
    Some(PiMonomial::new(
        area.coefficient / two_to_dimension,
        area.pi_exponent - dimension as i32,
    ))
}

/// Slope `b_*` of the inverse squared bare coupling of four-dimensional pure Yang-Mills
/// theory per unit of `log(1/a)`: `gauge_loop_coefficient * casimir / (8 pi^2)`.
///
/// Normalization: the action is `(1 / (2 g^2)) * sum_{mu<nu} |F_{mu nu}|^2` in an invariant
/// inner product, and `casimir` is defined by `Killing = -casimir * inner product`. The
/// product `g^2 * casimir` does not depend on the scale of the inner product. The textbook
/// one-loop input is `gauge_loop_coefficient = 11/3`; this function checks arithmetic, it
/// does not derive that number. With 11/3 the result is `11 * casimir / (24 pi^2)`.
pub fn one_loop_inverse_coupling_slope(
    casimir: &BigRational,
    gauge_loop_coefficient: &BigRational,
) -> PiMonomial {
    let shell = logarithmic_shell_coefficient(4).expect("dimension four is positive");
    shell.scale(&(casimir * gauge_loop_coefficient))
}
