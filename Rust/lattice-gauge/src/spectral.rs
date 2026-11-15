//! Purity of the trace of a nonnegative transfer operator, in exact integer arithmetic.
//!
//! For eigenvalues `lambda_1 >= lambda_2 >= ... >= 0` and `Z(n) = sum_j lambda_j^n`, put
//! `s_n = sum_{j>=2} (lambda_j / lambda_1)^n` and `P(n) = Z(2n) / Z(n)^2`. The mono-notes
//! (Proposition `ym:prop:purity`) prove `s_n <= 1/P(n) - 1 <= 2 s_n + s_n^2` and
//! `s_{2n} <= s_n^2`. This module evaluates both sides exactly for integer spectra, after
//! clearing denominators.

/// `Z(n) = sum_j lambda_j^n`; `None` on overflow.
pub fn trace_of_power(spectrum: &[u64], exponent: u32) -> Option<u128> {
    spectrum.iter().try_fold(0u128, |total, &eigenvalue| {
        total.checked_add(u128::from(eigenvalue).checked_pow(exponent)?)
    })
}

/// The three inequalities of the purity proposition, each as an exact integer comparison.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PurityBounds {
    /// `s_n <= 1/P(n) - 1`.
    pub excited_sum_below_defect: bool,
    /// `1/P(n) - 1 <= 2 s_n + s_n^2`.
    pub defect_below_quadratic: bool,
    /// `s_{2n} <= s_n^2`.
    pub doubling_improves: bool,
}

/// Evaluates the bounds for a spectrum of nonnegative integers with a positive maximum.
/// `None` if the spectrum is empty, its maximum is zero, or an intermediate overflows.
///
/// With `t = lambda_1^n`, `e = sum_{j>=2} lambda_j^n` (so `s_n = e / t`),
/// `a = Z(n) = t + e`, `b = Z(2n)`: the inequalities are
/// `e b <= t (a^2 - b)`, `(a^2 - b) t^2 <= b (2 e t + e^2)`,
/// `(b - t^2) <= e^2` (both sides of the last divided by `t^2`).
pub fn purity_bounds(spectrum: &[u64], exponent: u32) -> Option<PurityBounds> {
    let largest = *spectrum.iter().max()?;
    if largest == 0 || exponent == 0 {
        return None;
    }
    let top = u128::from(largest).checked_pow(exponent)?;
    let trace = trace_of_power(spectrum, exponent)?;
    let trace_doubled = trace_of_power(spectrum, 2 * exponent)?;
    let excited = trace - top;
    let defect_numerator = trace.checked_mul(trace)?.checked_sub(trace_doubled)?;
    let top_squared = top.checked_mul(top)?;
    Some(PurityBounds {
        excited_sum_below_defect: excited.checked_mul(trace_doubled)? <= top.checked_mul(defect_numerator)?,
        defect_below_quadratic: defect_numerator.checked_mul(top_squared)?
            <= trace_doubled.checked_mul(
                excited.checked_mul(top)?.checked_mul(2)?.checked_add(excited.checked_mul(excited)?)?,
            )?,
        doubling_improves: trace_doubled - top_squared <= excited.checked_mul(excited)?,
    })
}
