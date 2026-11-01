//! Public exact-arithmetic client; no book, server, frontend or model required.
use qft_symbolics::{
    clifford,
    computation::{Budget, Limits},
    exact::{BigC, BigQ},
    polynomial::BigPoly,
    scalar::Q,
};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut budget = Budget::new(Limits {
        max_operations: 100_000_000,
        ..Limits::default()
    });
    let huge = BigQ::from(2).checked_pow(256, &mut budget)?;
    println!("2^256 = {huge}");
    let third = BigQ::parse(&huge.to_string(), "3", &mut budget)?;
    println!(
        "2^256/3 + 2/3 = {}",
        third.checked_add(&BigQ::parse("2", "3", &mut budget)?, &mut budget)?
    );
    let z = BigC::new(huge.clone(), huge.clone());
    println!(
        "(2^256 + i 2^256) / itself = {}",
        z.checked_div(&z, &mut budget)?
    );
    let factor = BigPoly::try_constant(BigC::new(huge, 0.into()), &mut budget)?;
    let trace = clifford::try_trace_exact(&["fresh_mu", "fresh_nu"], &mut budget)?;
    println!(
        "Big coefficient times the same shared gamma-trace engine = {}",
        factor.try_mul(&trace, &mut budget)?
    );
    let mut small = Budget::default();
    println!(
        "Checked small-backend overflow = {:?}",
        Q::new(i128::MAX, 1)
            .checked_add(&Q::ONE, &mut small)
            .unwrap_err()
    );
    let mut limited = Budget::new(Limits {
        max_operations: 0,
        ..Limits::default()
    });
    println!(
        "Explicit exhausted budget = {:?}",
        clifford::try_trace_exact(&["a", "b"], &mut limited).unwrap_err()
    );
    println!("Exact arithmetic within declared resource limits; not a formal proof.");
    Ok(())
}
