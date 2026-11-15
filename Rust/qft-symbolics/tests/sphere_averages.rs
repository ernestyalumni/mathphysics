//! Oracles: the SymPy integrals in docs/plans/yang-mills-mass-gap/ledger/evidence/check_oneloop.py
//! (<n1^2> = 1/4, <n1^2 n0^2> = 1/24 on S^3, area 2 pi^2, slope 11 C / (24 pi^2)) and the
//! tensor identity <n_a n_b n_c n_d> = (dd + dd + dd) / (d (d + 2)).
use num_bigint::BigInt;
use num_rational::BigRational;
use num_traits::{One, Zero};
use qft_symbolics::sphere::{
    exact_rational, logarithmic_shell_coefficient, one_loop_inverse_coupling_slope,
    sphere_monomial_average, unit_sphere_area, PiMonomial,
};

fn q(numerator: i64, denominator: i64) -> BigRational {
    exact_rational(numerator, denominator).unwrap()
}

#[test]
fn averages_on_the_three_sphere_match_the_symbolic_integrals() {
    assert_eq!(sphere_monomial_average(&[2, 0, 0, 0]).unwrap(), q(1, 4));
    assert_eq!(sphere_monomial_average(&[2, 2, 0, 0]).unwrap(), q(1, 24));
    assert_eq!(sphere_monomial_average(&[4, 0, 0, 0]).unwrap(), q(1, 8));
    assert_eq!(sphere_monomial_average(&[1, 1, 0, 0]).unwrap(), BigRational::zero());
    assert_eq!(sphere_monomial_average(&[0, 0, 0, 0]).unwrap(), BigRational::one());
}

#[test]
fn negative_control_a_wrong_quartic_average_is_rejected() {
    assert_ne!(sphere_monomial_average(&[2, 2, 0, 0]).unwrap(), q(1, 12));
    assert_ne!(sphere_monomial_average(&[2, 2, 0, 0]).unwrap(), q(1, 16));
    assert!(sphere_monomial_average(&[]).is_none());
    assert!(exact_rational(1, 0).is_none());
}

#[test]
fn inserting_the_unit_constraint_does_not_change_an_average() {
    // sum_i <n_i^2 f> = <f> because sum_i n_i^2 = 1 on the sphere.
    for dimension in 1..=7usize {
        for base in [vec![0u32; dimension], {
            let mut exponents = vec![0u32; dimension];
            exponents[0] = 2;
            exponents[dimension - 1] += 4;
            exponents
        }] {
            let mut total = BigRational::zero();
            for index in 0..dimension {
                let mut exponents = base.clone();
                exponents[index] += 2;
                total += sphere_monomial_average(&exponents).unwrap();
            }
            assert_eq!(total, sphere_monomial_average(&base).unwrap());
        }
    }
}

#[test]
fn sphere_areas_and_the_four_dimensional_shell() {
    assert_eq!(unit_sphere_area(1).unwrap(), PiMonomial::new(q(2, 1), 0));
    assert_eq!(unit_sphere_area(2).unwrap(), PiMonomial::new(q(2, 1), 1));
    assert_eq!(unit_sphere_area(3).unwrap(), PiMonomial::new(q(4, 1), 1));
    assert_eq!(unit_sphere_area(4).unwrap(), PiMonomial::new(q(2, 1), 2));
    assert_eq!(unit_sphere_area(5).unwrap(), PiMonomial::new(q(8, 3), 2));
    assert_eq!(unit_sphere_area(6).unwrap(), PiMonomial::new(q(1, 1), 3));
    assert!(unit_sphere_area(0).is_none());
    assert_eq!(
        logarithmic_shell_coefficient(4).unwrap(),
        PiMonomial::new(q(1, 8), -2)
    );
}

#[test]
fn one_loop_slope_is_eleven_casimir_over_twenty_four_pi_squared() {
    let eleven_thirds = q(11, 3);
    // Casimir constants 2 N^2 of su(N) in the normalization <X,Y> = -(1/N) Re tr(XY).
    for (size, expected) in [(2i64, q(11, 3)), (3, q(33, 4))] {
        let casimir = BigRational::from_integer(BigInt::from(2 * size * size));
        let slope = one_loop_inverse_coupling_slope(&casimir, &eleven_thirds);
        assert_eq!(slope, PiMonomial::new(expected, -2));
        assert_eq!(slope.coefficient, &casimir * q(11, 24));
    }
    let slope = one_loop_inverse_coupling_slope(&q(8, 1), &eleven_thirds);
    assert_ne!(slope, PiMonomial::new(q(11, 6), -2), "the b_0 = b_*/2 value must not pass");
    let unit = PiMonomial::new(q(1, 1), 0);
    assert_eq!(slope.multiply(&unit), slope);
    assert_eq!(slope.divide(&slope).unwrap(), unit);
    assert!(slope.divide(&PiMonomial::new(BigRational::zero(), 0)).is_none());
}
