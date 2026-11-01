//! A direct public-library consumer, with no gateway or serialization dependency.
use lie_groups::su2_observation::{observe, Chart};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let norm = 38f64.sqrt();
    let axis = [2. / norm, 3. / norm, -5. / norm];
    let observation = observe(axis, 1.37, 0.42)?;
    assert_eq!(observation.axis, axis);
    assert!(matches!(observation.stereographic, Chart::Finite(_)));
    assert!(observation.checks.unitarity_error < 1e-12);
    println!(
        "coordinates={:?} bloch={:?}",
        observation.coordinates, observation.bloch
    );
    Ok(())
}
