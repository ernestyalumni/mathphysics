//! Reusable Lorentz and operator-calculus checks, independent of textbook files.
use qft_symbolics::{
    checked,
    computation::{Budget, Limits},
    differential::DifferentialOperator,
    exact::{BigC, BigQ},
    ordered::OrderedPolynomial,
    polynomial::BigPoly,
    relativity::{self as r, Axis, Rapidity, Tolerance, UnitAxis},
};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut b = Budget::new(Limits {
        max_operations: 10_000_000_000,
        ..Limits::default()
    });
    println!("1. RAPIDITY AND THE BOOST SIGN");
    let eta = Rapidity::from_beta(0.6)?;
    let matrix = r::boost_f64(eta, UnitAxis::axial(Axis::Z))?;
    let residual = r::metric_residual_f64(
        &matrix,
        Tolerance {
            absolute: 1e-12,
            relative: 1e-12,
        },
    )?;
    assert!(residual.within_tolerance);
    println!(
        "   beta=0.6 -> eta={:.12}; positive tz/zt entries={:.12}",
        eta.eta(),
        matrix[0][3]
    );
    println!(
        "   max |Lambda^T g Lambda-g|={:.3e}; numerical tolerance check, not exact proof",
        residual.absolute
    );
    assert!(Rapidity::from_beta(1.0).is_err());
    let composed = eta.collinear_add(eta)?;
    assert!((composed.beta()? - (0.6 + 0.6) / (1.0 + 0.6 * 0.6)).abs() < 1e-14);
    println!("   Collinear rapidities add; the corresponding velocities obey Einstein addition.\n");

    println!("2. EXACT COMPONENTS AND RAISING/LOWERING");
    let u = BigQ::parse("1", "3", &mut b)?;
    let boost = r::exact_half_velocity_boost(&u, Axis::Z, &mut b)?;
    let v = [BigC::from(2), BigC::from(0), BigC::from(0), BigC::from(3)];
    let transformed = r::apply(&boost, &v, &mut b)?;
    assert_eq!(boost[0][0].re, BigQ::parse("5", "4", &mut b)?);
    assert_eq!(boost[0][3].re, BigQ::parse("3", "4", &mut b)?);
    assert_eq!(transformed[0].re, BigQ::parse("19", "4", &mut b)?);
    assert_eq!(transformed[3].re, BigQ::parse("21", "4", &mut b)?);
    assert_eq!(r::raise(&r::lower(&v, &mut b)?, &mut b)?, v);
    assert_eq!(r::dot(&transformed, &transformed, &mut b)?, BigC::from(5));
    let g = r::metric::<BigC>(&mut b)?;
    let gram = checked::mul(
        &checked::mul(&checked::transpose(&boost, &mut b)?, &g, &mut b)?,
        &boost,
        &mut b,
    )?;
    assert_eq!(gram, g);
    println!("   u=tanh(eta/2)=1/3 gives cosh=5/4, sinh=3/4 exactly.");
    println!("   (2,0,0,3) -> (19/4,0,0,21/4); Minkowski squared norm stays 5.\n");

    println!("3. SYMBOLIC BOOST CONDITION");
    let symbolic = r::symbolic_axial_boost("c_eta", "s_eta", Axis::Z, &mut b)?;
    assert!(!symbolic.condition.is_zero());
    println!(
        "   c and s stay symbolic; required zero polynomial: {}",
        symbolic.condition
    );
    println!("   Their names alone do not establish c^2-s^2=1.\n");

    println!("4. ORDERED ALGEBRA AND FIRST-ORDER CONJUGATION");
    let a = OrderedPolynomial::generator("A", &mut b)?;
    let c = OrderedPolynomial::generator("B", &mut b)?;
    let d = OrderedPolynomial::generator("C", &mut b)?;
    let jacobi = a
        .commutator(&c, &mut b)?
        .commutator(&d, &mut b)?
        .try_add(&c.commutator(&d, &mut b)?.commutator(&a, &mut b)?, &mut b)?
        .try_add(&d.commutator(&a, &mut b)?.commutator(&c, &mut b)?, &mut b)?;
    assert!(jacobi.terms().is_empty());
    assert_ne!(a.try_mul(&c, &mut b)?, c.try_mul(&a, &mut b)?);
    let i = BigPoly::try_constant(BigC::new(0.into(), 1.into()), &mut b)?;
    let eps = BigPoly::try_var("eps", &mut b)?.try_mul(&i, &mut b)?;
    let one = OrderedPolynomial::scalar(BigPoly::try_constant(1, &mut b)?, &mut b)?;
    let small = c.try_scale(&eps, &mut b)?;
    let conjugated = one
        .try_sub(&small, &mut b)?
        .try_mul(&a, &mut b)?
        .try_mul(&one.try_add(&small, &mut b)?, &mut b)?
        .truncate(&["eps"], 1, &mut b)?;
    assert_eq!(
        conjugated.coefficient("eps", 1, &mut b)?,
        a.commutator(&c, &mut b)?.try_scale(&i, &mut b)?
    );
    println!("   Jacobi expands to zero with noncommuting A,B,C.");
    println!("   (1-i eps B) A (1+i eps B) = A+i eps[A,B] + O(eps^2).\n");

    println!("5. ORBITAL DIFFERENTIAL OPERATORS");
    let coordinates = || ["t", "x", "y", "z"].map(str::to_owned);
    let l01 = DifferentialOperator::orbital_lorentz(coordinates(), 0, 1, 1.into(), &mut b)?;
    let l12 = DifferentialOperator::orbital_lorentz(coordinates(), 1, 2, 1.into(), &mut b)?;
    let l02 = DifferentialOperator::orbital_lorentz(coordinates(), 0, 2, 1.into(), &mut b)?;
    let minus_i = BigPoly::try_constant(BigC::new(0.into(), (-1).into()), &mut b)?;
    assert!(l01
        .commutator(&l12, &mut b)?
        .try_sub(&l02.try_scale(&minus_i, &mut b)?, &mut b)?
        .terms()
        .is_empty());
    let p0 = DifferentialOperator::momentum(coordinates(), 0, 1.into(), &mut b)?;
    let p3 = DifferentialOperator::momentum(coordinates(), 3, 1.into(), &mut b)?;
    assert!(p0.commutator(&p3, &mut b)?.terms().is_empty());
    println!("   Leibniz expansion verifies [L01,L12]=-i L02 and [P0,P3]=0 (hbar=1).");
    println!("   These are exact differential-operator normal forms, not a general quantum-domain theorem.");
    println!("   Field commutators leave a possible commuting term undetected; its absence is not inferred.");
    println!("\nAll positive checks and rejection controls met their declared expectations.");
    Ok(())
}
