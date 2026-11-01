use physics_contracts::{
    valid_item_id, validate_json, Arithmetic, Bundle, Checker, EvidenceKind, Numerical, Outcome,
    MAX_INPUT_BYTES, MAX_RECORDS,
};
use serde_json::{json, Value};
fn value() -> Value {
    serde_json::from_str(include_str!("../examples/valid-bundle.json")).unwrap()
}
fn check(value: &Value) -> physics_contracts::ValidationReport {
    validate_json(&serde_json::to_vec(value).unwrap())
}
fn rejects(v: &Value, code: &str) {
    let r = check(v);
    assert_eq!(r.status, "invalid_input", "{r:?}");
    assert!(r.errors.iter().any(|e| e.code == code), "{r:?}");
}
#[test]
fn namespaced_ids_do_not_confuse_equations_exercises_or_other_material() {
    for kind in [
        "equation",
        "exercise",
        "prose",
        "diagram",
        "example",
        "definition",
        "theorem",
    ] {
        assert!(valid_item_id(&format!("srednicki:ch035:{kind}:35.1")));
    }
    for id in [
        "srednicki:ch098:exercise:1",
        "srednicki:ch000:exercise:1",
        "srednicki:ch35:exercise:1",
        "srednicki:ch035:unknown:1",
        "Srednicki:ch098:exercise:1",
        "srednicki:ch035:exercise:",
        "srednicki:ch035:exercise:1:extra",
        "srednicki:ché5:exercise:1",
    ] {
        assert!(!valid_item_id(id), "{id}");
    }
    let mut v = value();
    let mut item = v["items"][0].clone();
    item["id"] = json!("srednicki:ch035:exercise:35.1");
    item["evidence_ids"] = json!([]);
    v["items"].as_array_mut().unwrap().push(item);
    assert_eq!(check(&v).status, "record_valid");
}
#[test]
fn valid_metadata_never_claims_proof_or_file_verification() {
    let r = check(&value());
    assert!(r.errors.is_empty());
    assert_eq!(r.status, "record_valid");
    assert!(r.assurance.contains("no mathematical proof"));
    assert!(r.assurance.contains("digest verification"));
    assert!(r.assurance.contains("artifact replay"));
}
#[test]
fn unknown_and_duplicate_json_fields_and_versions_rejected() {
    let mut v = value();
    v["items"][0]["verified"] = json!(true);
    rejects(&v, "malformed_json");
    let mut v = value();
    v["items"][0]["sources"][0]["locator"]["extra"] = json!(1);
    rejects(&v, "malformed_json");
    let mut v = value();
    v["schema_version"] = json!("v999");
    rejects(&v, "unsupported_version");
    let input = serde_json::to_string(&value()).unwrap().replacen(
        "\"schema_version\":",
        "\"schema_version\":\"wrong\",\"schema_version\":",
        1,
    );
    assert_eq!(validate_json(input.as_bytes()).status, "invalid_input");
}
#[test]
fn missing_source_invalid_digest_and_page_ranges_rejected() {
    let mut v = value();
    v["items"][0]["sources"] = json!([]);
    rejects(&v, "source_count");
    for bad in ["abc", &"A".repeat(64), &"g".repeat(64)] {
        let mut v = value();
        v["items"][0]["sources"][0]["sha256"] = json!(bad);
        rejects(&v, "invalid_source");
    }
    let mut v = value();
    v["items"][0]["sources"][0]["locator"]["pdf_page"] = json!(0);
    rejects(&v, "invalid_locator");
    let mut v = value();
    v["items"][0]["sources"][0]["locator"] = json!({"kind":"text","line_start":10,"line_end":9});
    rejects(&v, "invalid_locator");
}
#[test]
fn unknown_and_duplicate_assumptions_are_not_ignored() {
    let mut v = value();
    v["items"][0]["assumption_refs"] = json!(["missing"]);
    rejects(&v, "unknown_reference");
    let mut v = value();
    v["items"][0]["assumption_refs"] = json!(["commutative", "commutative"]);
    rejects(&v, "duplicate_reference");
    let mut v = value();
    v["items"][0]["convention_ref"] = json!("missing");
    rejects(&v, "unknown_convention");
    let mut v = value();
    let a = v["assumptions"][0].clone();
    v["assumptions"].as_array_mut().unwrap().push(a);
    rejects(&v, "duplicate_id");
}
#[test]
fn evidence_must_bind_statement_domain_assumptions_and_conventions() {
    for key in ["statement", "domain"] {
        let mut v = value();
        v["evidence"][0][key] = json!("different");
        rejects(&v, "statement_mismatch");
    }
    let mut v = value();
    v["evidence"][0]["assumption_refs"] = json!([]);
    rejects(&v, "assumption_mismatch");
    let mut v = value();
    v["conventions"]
        .as_array_mut()
        .unwrap()
        .push(json!({"id":"other","statement":"other convention"}));
    v["evidence"][0]["convention_ref"] = json!("other");
    rejects(&v, "assumption_mismatch");
}
#[test]
fn orphan_missing_and_cross_item_evidence_links_rejected() {
    let mut v = value();
    v["evidence"][0]["item_id"] = json!("srednicki:ch001:equation:1.1");
    rejects(&v, "unknown_item");
    let mut v = value();
    v["items"][0]["evidence_ids"] = json!([]);
    rejects(&v, "missing_backlink");
    let mut v = value();
    v["items"][0]["evidence_ids"] = json!(["missing"]);
    rejects(&v, "unknown_reference");
    let mut v = value();
    let mut i = v["items"][0].clone();
    i["id"] = json!("srednicki:ch035:exercise:35.1");
    v["items"].as_array_mut().unwrap().push(i);
    rejects(&v, "evidence_link_mismatch");
}
#[test]
fn numeric_agreement_cannot_be_promoted_to_established_exact_claim() {
    let mut v = value();
    v["evidence"][0]["kind"] = json!("numerical");
    v["evidence"][0]["arithmetic"] = json!("approximate");
    v["evidence"][0]["outcome"] = json!("established");
    v["evidence"][0]["numerical"] = json!({"absolute_tolerance":1e-10,"relative_tolerance":1e-10,"precision_bits":53,"error_note":"finite points only"});
    rejects(&v, "numerical_overclaim");
    v["evidence"][0]["outcome"] = json!("computed");
    assert_eq!(check(&v).status, "record_valid");
    v["evidence"][0]["numerical"]["absolute_tolerance"] = json!(-1.);
    rejects(&v, "invalid_numerics");
}
#[test]
fn nan_and_infinity_are_rejected_through_typed_library_api() {
    for invalid in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        let mut b: Bundle = serde_json::from_value(value()).unwrap();
        let e = &mut b.evidence[0];
        e.kind = EvidenceKind::Numerical;
        e.arithmetic = Arithmetic::Approximate;
        e.numerical = Some(Numerical {
            absolute_tolerance: invalid,
            relative_tolerance: 1e-6,
            precision_bits: 53,
            error_note: "finite witness".into(),
        });
        assert!(physics_contracts::validate(&b)
            .errors
            .iter()
            .any(|x| x.code == "invalid_numerics"));
    }
}
#[test]
fn external_agreement_and_source_derivations_cannot_claim_established() {
    for kind in ["external_comparison", "source_derivation"] {
        let mut v = value();
        v["evidence"][0]["kind"] = json!(kind);
        v["evidence"][0]["outcome"] = json!("established");
        rejects(&v, "evidence_overclaim");
    }
}
#[test]
fn formal_metadata_requires_exact_arithmetic_and_linked_checker() {
    let mut b: Bundle = serde_json::from_value(value()).unwrap();
    let e = &mut b.evidence[0];
    e.kind = EvidenceKind::FormalKernel;
    e.arithmetic = Arithmetic::Exact;
    e.outcome = Outcome::Established;
    assert!(physics_contracts::validate(&b)
        .errors
        .iter()
        .any(|x| x.code == "missing_checker"));
    b.evidence[0].checker = Some(Checker {
        name: "example-checker".into(),
        version: "0.1".into(),
        implementation_sha256: "0".repeat(64),
        artifact_sha256: "1".repeat(64),
    });
    assert!(physics_contracts::validate(&b)
        .errors
        .iter()
        .any(|x| x.code == "missing_checker"));
    b.evidence[0].checker.as_mut().unwrap().artifact_sha256 = "0".repeat(64);
    let r = physics_contracts::validate(&b);
    assert!(r.errors.is_empty());
    assert_eq!(r.status, "record_valid");
    assert!(r.assurance.contains("no mathematical proof"));
    b.evidence[0].arithmetic = Arithmetic::Approximate;
    assert!(physics_contracts::validate(&b)
        .errors
        .iter()
        .any(|x| x.code == "exactness_mismatch"));
}
#[test]
fn failures_require_reasons_and_positive_declared_budgets() {
    let mut v = value();
    v["evidence"][0]["outcome"] = json!("inconclusive");
    v["evidence"][0]["reason"] = Value::Null;
    rejects(&v, "missing_reason");
    for key in ["wall_clock_seconds", "memory_bytes", "operation_limit"] {
        let mut v = value();
        v["evidence"][0]["resources"][key] = json!(0);
        rejects(&v, "invalid_budget");
    }
}
#[test]
fn byte_and_record_limits_fail_closed() {
    assert_eq!(
        validate_json(&vec![b' '; MAX_INPUT_BYTES + 1]).errors[0].code,
        "resource_limit"
    );
    let mut b: Bundle = serde_json::from_value(value()).unwrap();
    b.items[0].evidence_ids = vec!["example-trace".into(); MAX_RECORDS + 1];
    assert_eq!(
        physics_contracts::validate(&b).errors[0].code,
        "resource_limit"
    );
    for raw in [b"null".as_slice(), b"{}".as_slice(), b"not json".as_slice()] {
        assert_eq!(validate_json(raw).status, "invalid_input");
    }
}
#[test]
fn typed_roundtrip_preserves_all_evidence_and_source_states() {
    let b: Bundle = serde_json::from_value(value()).unwrap();
    let encoded = serde_json::to_vec(&b).unwrap();
    assert_eq!(serde_json::from_slice::<Value>(&encoded).unwrap(), value());
    assert!(validate_json(&encoded).errors.is_empty());
}

#[test]
fn front_and_back_matter_have_truthful_item_locations() {
    for id in [
        "srednicki:frontmatter:prose:preface",
        "srednicki:backmatter:prose:bibliography",
    ] {
        assert!(valid_item_id(id));
    }
    assert!(!valid_item_id("srednicki:ch000:prose:preface"));
}
#[test]
fn approximate_external_comparisons_require_a_numerical_error_contract() {
    let mut v = value();
    v["evidence"][0]["kind"] = json!("external_comparison");
    v["evidence"][0]["arithmetic"] = json!("approximate");
    rejects(&v, "invalid_numerics");
    v["evidence"][0]["numerical"] = json!({"absolute_tolerance":1e-10,"relative_tolerance":0.,"precision_bits":53,"error_note":"finite comparison, not proof"});
    assert_eq!(check(&v).status, "record_valid");
}
#[test]
fn one_artifact_path_cannot_claim_two_contents() {
    let mut v = value();
    let mut a = v["evidence"][0]["artifacts"][0].clone();
    a["sha256"] = json!("1".repeat(64));
    v["evidence"][0]["artifacts"]
        .as_array_mut()
        .unwrap()
        .push(a);
    rejects(&v, "artifact_identity_conflict");
    let mut v = value();
    let a = v["evidence"][0]["artifacts"][0].clone();
    v["evidence"][0]["artifacts"]
        .as_array_mut()
        .unwrap()
        .push(a);
    rejects(&v, "duplicate_artifact");
}
