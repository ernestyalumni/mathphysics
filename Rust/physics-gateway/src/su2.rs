//! Strict JSON adapter for the reusable typed SU(2) library.
use lie_groups::su2_observation::{observe, validate_axis, validate_t, Chart, Slice};
use serde_json::{json, Value};
fn number(value: &Value) -> Result<f64, &'static str> {
    value
        .as_f64()
        .filter(|x| x.is_finite())
        .ok_or("finite_number_required")
}

pub fn handle_bytes(body: &[u8]) -> Result<Value, &'static str> {
    let input: Value = serde_json::from_slice(body).map_err(|_| "invalid_json")?;
    let object = input.as_object().ok_or("object_required")?;
    if object.len() != 3
        || !["axis", "t", "slice"]
            .iter()
            .all(|k| object.contains_key(*k))
    {
        return Err("expected_axis_t_slice");
    }
    let values = input["axis"]
        .as_array()
        .filter(|a| a.len() == 3)
        .ok_or("axis_requires_three_numbers")?;
    let axis = [
        number(&values[0])?,
        number(&values[1])?,
        number(&values[2])?,
    ];

    validate_axis(axis).map_err(|e| e.code())?;
    let t = number(&input["t"])?;
    validate_t(t).map_err(|e| e.code())?;
    let c = number(&input["slice"])?;
    let o = observe(axis, t, c).map_err(|e| e.code())?;
    let stereographic = match o.stereographic {
        Chart::Finite(point) => json!({"status":"finite","point":point}),
        Chart::Pole => json!({"status":"pole","point":null}),
    };
    let slice = match o.slice {
        Slice::Sphere { coordinate, radius } => {
            json!({"coordinate":coordinate,"status":"sphere","radius":radius})
        }
        Slice::Point { coordinate } => {
            json!({"coordinate":coordinate,"status":"point","radius":0.0})
        }
        Slice::Empty { coordinate } => {
            json!({"coordinate":coordinate,"status":"empty","radius":null})
        }
    };
    Ok(json!({
        "convention":"pauli-plus-i-v1","t":o.t,"axis":o.axis,
        "coordinates":o.coordinates,"alpha":o.alpha,"beta":o.beta,"matrix":o.matrix,
        "stereographic":stereographic,"slice":slice,"bloch":o.bloch,
        "checks":{"unitarity_error":o.checks.unitarity_error,"determinant_error":o.checks.determinant_error},
        "evidence_kind":"numerical","formal_proof":false
    }))
}
