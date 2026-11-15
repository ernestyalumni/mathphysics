//! Exact enumeration of Z2 lattice gauge theory on small tori.
//!
//! Links carry `sigma = +1 or -1`. A plaquette has value `sigma_p`, the product of its four
//! links, and a sign `s_p = -1` if it is twisted, `+1` otherwise. A plaquette is *satisfied*
//! if `s_p sigma_p = +1`. With `P` plaquettes and `j` satisfied ones, the Wilson energy is
//! `sum_p (1 - s_p sigma_p) = 2 (P - j)`, so the twisted partition function is
//! `Z[twist](beta) = sum_j N_j exp(-2 beta (P - j))`, where `N_j` is the histogram returned
//! here. The histogram is the complete exact answer for every beta.
use crate::torus::Torus;

/// Largest number of links enumerated (2^28 configurations).
pub const MAXIMUM_ENUMERATED_LINKS: usize = 28;

/// `N_j` for `j = 0..=P`: the number of link configurations with exactly `j` satisfied
/// plaquettes. `None` if `twisted` has the wrong length or there are more than
/// [`MAXIMUM_ENUMERATED_LINKS`] links. Uses a Gray code: one link flip per step.
pub fn satisfied_plaquette_histogram(torus: &Torus, twisted: &[bool]) -> Option<Vec<u128>> {
    let plaquette_count = torus.plaquette_count();
    let link_count = torus.link_count();
    if twisted.len() != plaquette_count || link_count > MAXIMUM_ENUMERATED_LINKS {
        return None;
    }
    let coboundaries: Vec<Vec<usize>> = (0..link_count)
        .map(|link| torus.link_coboundary(link))
        .collect();
    // unsatisfied[p] is true when s_p sigma_p = -1. All links start at +1.
    let mut unsatisfied = twisted.to_vec();
    let mut satisfied_count = unsatisfied.iter().filter(|&&flag| !flag).count();
    let mut histogram = vec![0u128; plaquette_count + 1];
    histogram[satisfied_count] += 1;
    for step in 1u64..(1u64 << link_count) {
        let link = step.trailing_zeros() as usize;
        for &plaquette in &coboundaries[link] {
            if unsatisfied[plaquette] {
                satisfied_count += 1;
            } else {
                satisfied_count -= 1;
            }
            unsatisfied[plaquette] = !unsatisfied[plaquette];
        }
        histogram[satisfied_count] += 1;
    }
    Some(histogram)
}

/// `binomial(n, k)` in `u128`; `None` on overflow.
pub fn binomial(n: u32, k: u32) -> Option<u128> {
    if k > n {
        return Some(0);
    }
    let mut value: u128 = 1;
    for step in 0..k.min(n - k) {
        value = value.checked_mul(u128::from(n - step))? / u128::from(step + 1);
    }
    Some(value)
}

/// Closed form in two dimensions on an `L x L` torus (`P = L^2` plaquettes, `2 L^2` links),
/// from the character expansion: only the empty set and the set of all plaquettes survive
/// the sum over links, so with `x = e^beta`,
/// `sum_U prod_p x^{s_p sigma_p} = 2^{2P - P} [ (x + 1/x)^P + sign (x - 1/x)^P ]`,
/// `sign = prod_p s_p`. Hence `N_j = 2^P binomial(P, j) (1 + sign (-1)^{P - j})`.
pub fn two_dimensional_histogram_closed_form(side: u32, odd_number_of_twists: bool) -> Option<Vec<u128>> {
    let plaquettes = side.checked_mul(side)?;
    let prefactor = 1u128.checked_shl(plaquettes)?;
    let mut histogram = Vec::with_capacity(plaquettes as usize + 1);
    for satisfied in 0..=plaquettes {
        let parity_is_even = (plaquettes - satisfied) % 2 == 0;
        let bracket: u128 = if parity_is_even != odd_number_of_twists { 2 } else { 0 };
        histogram.push(prefactor.checked_mul(binomial(plaquettes, satisfied)?)?.checked_mul(bracket)?);
    }
    Some(histogram)
}
