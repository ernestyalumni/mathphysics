//! Metadata contracts, not a theorem prover or an evidence-artifact checker.
//! A valid record can still describe an incorrect mathematical result.
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};

pub const SCHEMA_VERSION: &str = "mathphysics.contracts.v1";
pub const MAX_INPUT_BYTES: usize = 1_048_576;
pub const MAX_RECORDS: usize = 10_000;
pub const ASSURANCE: &str = "metadata_validation_only; no mathematical proof, source certification, digest verification, computation or artifact replay performed";

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Bundle {
    pub schema_version: String,
    pub assumptions: Vec<Definition>,
    pub conventions: Vec<Definition>,
    pub items: Vec<MathItem>,
    pub evidence: Vec<Evidence>,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Definition {
    pub id: String,
    pub statement: String,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct MathItem {
    pub id: String,
    pub statement: String,
    pub domain: String,
    pub sources: Vec<Source>,
    pub assumption_refs: Vec<String>,
    pub convention_ref: String,
    pub implementation: Implementation,
    pub evidence_ids: Vec<String>,
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Implementation {
    Unassessed,
    Specified,
    Implemented,
    ResearchGap,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Source {
    pub path: String,
    pub sha256: String,
    pub locator: Locator,
    pub review: SourceReview,
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum SourceReview {
    Candidate,
    TextChecked,
    VisuallyChecked,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
pub enum Locator {
    Pdf {
        pdf_page: u32,
        printed_page: Option<String>,
    },
    Text {
        line_start: u32,
        line_end: u32,
    },
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Evidence {
    pub id: String,
    pub item_id: String,
    pub statement: String,
    pub domain: String,
    pub assumption_refs: Vec<String>,
    pub convention_ref: String,
    pub kind: EvidenceKind,
    pub outcome: Outcome,
    pub arithmetic: Arithmetic,
    pub operation: Operation,
    pub artifacts: Vec<Artifact>,
    pub checker: Option<Checker>,
    pub numerical: Option<Numerical>,
    pub reason: Option<String>,
    pub resources: Resources,
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum EvidenceKind {
    ExactAlgebra,
    Numerical,
    FormalKernel,
    CertificateReplay,
    ExternalComparison,
    SourceDerivation,
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Outcome {
    Established,
    Computed,
    Counterexample,
    Inconclusive,
    Unsupported,
    InvalidInput,
}

#[derive(Clone, Copy, Debug, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Arithmetic {
    Exact,
    Approximate,
    NotApplicable,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Operation {
    pub name: String,
    pub input_sha256: String,
    pub implementation_sha256: String,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Artifact {
    pub path: String,
    pub sha256: String,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Checker {
    pub name: String,
    pub version: String,
    pub implementation_sha256: String,
    pub artifact_sha256: String,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Numerical {
    pub absolute_tolerance: f64,
    pub relative_tolerance: f64,
    pub precision_bits: u32,
    pub error_note: String,
}

#[derive(Clone, Debug, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct Resources {
    pub wall_clock_seconds: u64,
    pub memory_bytes: u64,
    pub operation_limit: u64,
}

#[derive(Clone, Debug, Deserialize, Serialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct ValidationError {
    pub path: String,
    pub code: String,
    pub message: String,
}

#[derive(Clone, Debug, Serialize)]
pub struct ValidationReport {
    pub schema_version: &'static str,
    pub status: &'static str,
    pub assurance: &'static str,
    pub item_count: usize,
    pub evidence_count: usize,
    pub errors: Vec<ValidationError>,
}

fn error(errors: &mut Vec<ValidationError>, path: &str, code: &str, message: &str) {
    if errors.len() < 128 {
        errors.push(ValidationError {
            path: path.into(),
            code: code.into(),
            message: message.into(),
        });
    }
}
fn text(value: &str) -> bool {
    !value.trim().is_empty() && value.len() <= 16_384 && !value.contains('\0')
}
fn name(value: &str) -> bool {
    !value.is_empty()
        && value.len() <= 128
        && value
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b"-_.".contains(&b))
}
fn digest(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
/// Namespaces distinguish equations, exercises and unnumbered material.
/// Chapter identifiers are three digits; frontmatter/backmatter are explicit scopes.
/// Srednicki chapters are restricted to 1..=97.
pub fn valid_item_id(value: &str) -> bool {
    if value.len() > 256 {
        return false;
    }
    let parts: Vec<_> = value.split(':').collect();
    if parts.len() != 4
        || !name(parts[0])
        || !name(parts[3])
        || !parts[0]
            .bytes()
            .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit() || b == b'-')
        || !parts[0].as_bytes()[0].is_ascii_lowercase()
    {
        return false;
    }
    let location = parts[1];
    if !matches!(location, "frontmatter" | "backmatter") {
        if location.len() != 5
            || !location.starts_with("ch")
            || !location.as_bytes()[2..].iter().all(u8::is_ascii_digit)
        {
            return false;
        }
        let Ok(number) = location[2..].parse::<u16>() else {
            return false;
        };
        if number == 0 || (parts[0] == "srednicki" && number > 97) {
            return false;
        }
    }
    matches!(
        parts[2],
        "equation" | "exercise" | "prose" | "diagram" | "example" | "definition" | "theorem"
    )
}
fn check_refs(
    values: &[String],
    known: &BTreeSet<&str>,
    path: &str,
    errors: &mut Vec<ValidationError>,
) {
    let mut seen = BTreeSet::new();
    if values.len() > 256 {
        error(errors, path, "resource_limit", "too many references");
    }
    for value in values {
        if !known.contains(value.as_str()) {
            error(errors, path, "unknown_reference", "undeclared reference");
        }
        if !seen.insert(value) {
            error(errors, path, "duplicate_reference", "duplicate reference");
        }
    }
}
fn definitions<'a>(
    values: &'a [Definition],
    path: &str,
    errors: &mut Vec<ValidationError>,
) -> BTreeSet<&'a str> {
    let mut known = BTreeSet::new();
    for (index, value) in values.iter().enumerate() {
        let p = format!("{path}[{index}]");
        if !name(&value.id) || !text(&value.statement) {
            error(
                errors,
                &p,
                "invalid_definition",
                "definition requires an identifier and statement",
            );
        }
        if !known.insert(value.id.as_str()) {
            error(
                errors,
                &p,
                "duplicate_id",
                "duplicate definition identifier",
            );
        }
    }
    known
}
fn check_source(source: &Source, path: &str, errors: &mut Vec<ValidationError>) {
    if !text(&source.path) || !digest(&source.sha256) {
        error(
            errors,
            path,
            "invalid_source",
            "source requires a path and lowercase SHA-256 digest",
        );
    }
    match &source.locator {
        Locator::Pdf {
            pdf_page,
            printed_page,
        } => {
            if *pdf_page == 0 || printed_page.as_ref().is_some_and(|p| !text(p)) {
                error(
                    errors,
                    path,
                    "invalid_locator",
                    "PDF pages are one-based; printed labels must be nonempty when supplied",
                );
            }
        }
        Locator::Text {
            line_start,
            line_end,
        } => {
            if *line_start == 0 || line_end < line_start {
                error(
                    errors,
                    path,
                    "invalid_locator",
                    "text lines require a nonempty one-based range",
                );
            }
        }
    }
}

/// Check metadata consistency only. Declared artifact hashes and evidence outcomes
/// are not authenticated or replayed. Resource fields are declarations, not OS limits.
pub fn validate(bundle: &Bundle) -> ValidationReport {
    let mut errors = Vec::new();
    let mut total = bundle
        .assumptions
        .len()
        .saturating_add(bundle.conventions.len())
        .saturating_add(bundle.items.len())
        .saturating_add(bundle.evidence.len());
    if total > MAX_RECORDS {
        return record_limit(bundle);
    }
    for item in &bundle.items {
        total = total
            .saturating_add(item.sources.len())
            .saturating_add(item.assumption_refs.len())
            .saturating_add(item.evidence_ids.len());
        if total > MAX_RECORDS {
            return record_limit(bundle);
        }
    }
    for evidence in &bundle.evidence {
        total = total
            .saturating_add(evidence.artifacts.len())
            .saturating_add(evidence.assumption_refs.len());
        if total > MAX_RECORDS {
            return record_limit(bundle);
        }
    }
    if bundle.schema_version != SCHEMA_VERSION {
        error(
            &mut errors,
            "schema_version",
            "unsupported_version",
            "unsupported contract version",
        );
    }
    let assumptions = definitions(&bundle.assumptions, "assumptions", &mut errors);
    let conventions = definitions(&bundle.conventions, "conventions", &mut errors);
    let mut items = BTreeMap::new();
    let evidence_ids: BTreeSet<_> = bundle.evidence.iter().map(|e| e.id.as_str()).collect();
    let mut backlinks: BTreeMap<&str, BTreeSet<&str>> = BTreeMap::new();
    for item in &bundle.items {
        for evidence_id in &item.evidence_ids {
            backlinks
                .entry(evidence_id.as_str())
                .or_default()
                .insert(item.id.as_str());
        }
    }
    for (index, item) in bundle.items.iter().enumerate() {
        let p = format!("items[{index}]");
        if !valid_item_id(&item.id) {
            error(
                &mut errors,
                &p,
                "invalid_item_id",
                "invalid namespaced chapter item identifier",
            );
        }
        if items.insert(item.id.as_str(), item).is_some() {
            error(&mut errors, &p, "duplicate_id", "duplicate item identifier");
        }
        if !text(&item.statement) || !text(&item.domain) {
            error(
                &mut errors,
                &p,
                "missing_domain_statement",
                "explicit statement and domain are required",
            );
        }
        if item.sources.is_empty() || item.sources.len() > 256 {
            error(
                &mut errors,
                &p,
                "source_count",
                "one to 256 sources are required",
            );
        }
        for (i, source) in item.sources.iter().enumerate() {
            check_source(source, &format!("{p}.sources[{i}]"), &mut errors);
        }
        check_refs(
            &item.assumption_refs,
            &assumptions,
            &format!("{p}.assumption_refs"),
            &mut errors,
        );
        if !conventions.contains(item.convention_ref.as_str()) {
            error(
                &mut errors,
                &p,
                "unknown_convention",
                "undeclared convention",
            );
        }
        check_refs(
            &item.evidence_ids,
            &evidence_ids,
            &format!("{p}.evidence_ids"),
            &mut errors,
        );
    }
    let mut seen = BTreeSet::new();
    for (index, evidence) in bundle.evidence.iter().enumerate() {
        let p = format!("evidence[{index}]");
        if !name(&evidence.id) || !seen.insert(evidence.id.as_str()) {
            error(
                &mut errors,
                &p,
                "invalid_evidence_id",
                "invalid or duplicate evidence identifier",
            );
        }
        check_refs(
            &evidence.assumption_refs,
            &assumptions,
            &format!("{p}.assumption_refs"),
            &mut errors,
        );
        if !conventions.contains(evidence.convention_ref.as_str()) {
            error(
                &mut errors,
                &p,
                "unknown_convention",
                "undeclared convention",
            );
        }
        if let Some(item) = items.get(evidence.item_id.as_str()) {
            if evidence.statement != item.statement || evidence.domain != item.domain {
                error(
                    &mut errors,
                    &p,
                    "statement_mismatch",
                    "evidence statement and domain must exactly match the referenced item",
                );
            }
            if evidence.convention_ref != item.convention_ref
                || evidence.assumption_refs.iter().collect::<BTreeSet<_>>()
                    != item.assumption_refs.iter().collect::<BTreeSet<_>>()
            {
                error(
                    &mut errors,
                    &p,
                    "assumption_mismatch",
                    "evidence must carry the item's complete assumptions and convention",
                );
            }
            if !item.evidence_ids.contains(&evidence.id) {
                error(
                    &mut errors,
                    &p,
                    "missing_backlink",
                    "referenced item must link this evidence",
                );
            }
        } else {
            error(
                &mut errors,
                &p,
                "unknown_item",
                "evidence references an unknown item",
            );
        }
        if backlinks
            .get(evidence.id.as_str())
            .is_some_and(|linked| linked.iter().any(|id| *id != evidence.item_id))
        {
            error(
                &mut errors,
                &p,
                "evidence_link_mismatch",
                "evidence may only be linked to its declared item",
            );
        }
        if !name(&evidence.operation.name)
            || !digest(&evidence.operation.input_sha256)
            || !digest(&evidence.operation.implementation_sha256)
        {
            error(
                &mut errors,
                &p,
                "invalid_operation",
                "operation requires a name and input/implementation digests",
            );
        }
        if evidence.artifacts.is_empty() || evidence.artifacts.len() > 256 {
            error(
                &mut errors,
                &p,
                "artifact_count",
                "one to 256 artifacts are required",
            );
        }
        let mut artifact_paths = BTreeMap::new();
        for artifact in &evidence.artifacts {
            if !text(&artifact.path) || !digest(&artifact.sha256) {
                error(
                    &mut errors,
                    &p,
                    "invalid_artifact",
                    "artifact requires a path and digest",
                );
            }
            if let Some(previous) =
                artifact_paths.insert(artifact.path.as_str(), artifact.sha256.as_str())
            {
                let code = if previous == artifact.sha256 {
                    "duplicate_artifact"
                } else {
                    "artifact_identity_conflict"
                };
                error(
                    &mut errors,
                    &p,
                    code,
                    "an artifact path must identify a unique record in an evidence snapshot",
                );
            }
        }
        if evidence.resources.wall_clock_seconds == 0
            || evidence.resources.memory_bytes == 0
            || evidence.resources.operation_limit == 0
        {
            error(
                &mut errors,
                &p,
                "invalid_budget",
                "all resource declarations must be positive",
            );
        }
        let failed = matches!(
            evidence.outcome,
            Outcome::Inconclusive | Outcome::Unsupported | Outcome::InvalidInput
        );
        if failed && !evidence.reason.as_ref().is_some_and(|x| text(x)) {
            error(
                &mut errors,
                &p,
                "missing_reason",
                "unsuccessful attempts require a reason",
            );
        }
        if evidence.reason.as_ref().is_some_and(|x| !text(x)) {
            error(
                &mut errors,
                &p,
                "invalid_reason",
                "reason must be nonempty when supplied",
            );
        }
        match evidence.kind {
            EvidenceKind::Numerical => {
                if evidence.arithmetic != Arithmetic::Approximate
                    || evidence.outcome == Outcome::Established
                {
                    error(&mut errors, &p, "numerical_overclaim", "numerical evidence must be approximate and cannot establish an exact claim");
                }
            }
            EvidenceKind::ExactAlgebra
            | EvidenceKind::FormalKernel
            | EvidenceKind::CertificateReplay => {
                if evidence.arithmetic != Arithmetic::Exact || evidence.numerical.is_some() {
                    error(
                        &mut errors,
                        &p,
                        "exactness_mismatch",
                        "exact evidence cannot use approximate arithmetic or numerical tolerance",
                    );
                }
            }
            EvidenceKind::ExternalComparison | EvidenceKind::SourceDerivation => {
                if evidence.outcome == Outcome::Established {
                    error(
                        &mut errors,
                        &p,
                        "evidence_overclaim",
                        "comparison or source derivation alone does not establish a proof",
                    );
                }
            }
        }
        if evidence.arithmetic == Arithmetic::Approximate {
            match &evidence.numerical {
                Some(n) if n.absolute_tolerance.is_finite() && n.relative_tolerance.is_finite()
                    && n.absolute_tolerance >= 0. && n.relative_tolerance >= 0.
                    && (n.absolute_tolerance > 0. || n.relative_tolerance > 0.)
                    && n.precision_bits > 0 && text(&n.error_note) => {},
                _ => error(&mut errors, &p, "invalid_numerics", "all approximate evidence requires finite positive error budget, precision and error explanation"),
            }
            if evidence.outcome == Outcome::Established {
                error(
                    &mut errors,
                    &p,
                    "approximate_overclaim",
                    "approximate arithmetic cannot establish an exact claim",
                );
            }
        } else if evidence.numerical.is_some() {
            error(
                &mut errors,
                &p,
                "kind_mismatch",
                "numerical error budgets require approximate arithmetic",
            );
        }
        if matches!(
            evidence.kind,
            EvidenceKind::FormalKernel | EvidenceKind::CertificateReplay
        ) {
            match &evidence.checker {
                Some(c) if text(&c.name) && text(&c.version) && digest(&c.implementation_sha256)
                    && digest(&c.artifact_sha256) && evidence.artifacts.iter().any(|a|a.sha256 == c.artifact_sha256) => {},
                _ => error(&mut errors, &p, "missing_checker", "formal/certificate evidence needs a checker version/digest linked to a declared artifact"),
            }
        } else if evidence.checker.is_some() {
            error(
                &mut errors,
                &p,
                "checker_kind_mismatch",
                "checker metadata belongs to formal or replay evidence only",
            );
        }
    }
    report(bundle, errors)
}
fn record_limit(bundle: &Bundle) -> ValidationReport {
    report(
        bundle,
        vec![ValidationError {
            path: "$".into(),
            code: "resource_limit".into(),
            message: "bundle exceeds record budget".into(),
        }],
    )
}
fn report(bundle: &Bundle, errors: Vec<ValidationError>) -> ValidationReport {
    ValidationReport {
        schema_version: SCHEMA_VERSION,
        status: if errors.is_empty() {
            "record_valid"
        } else {
            "invalid_input"
        },
        assurance: ASSURANCE,
        item_count: bundle.items.len(),
        evidence_count: bundle.evidence.len(),
        errors,
    }
}
/// Decode a bounded JSON bundle; unknown and duplicate fields are rejected.
pub fn validate_json(input: &[u8]) -> ValidationReport {
    let invalid = |code: &str, message: &str| ValidationReport {
        schema_version: SCHEMA_VERSION,
        status: "invalid_input",
        assurance: ASSURANCE,
        item_count: 0,
        evidence_count: 0,
        errors: vec![ValidationError {
            path: "$".into(),
            code: code.into(),
            message: message.into(),
        }],
    };
    if input.len() > MAX_INPUT_BYTES {
        return invalid("resource_limit", "input exceeds one MiB");
    }
    match serde_json::from_slice::<Bundle>(input) {
        Ok(bundle) => validate(&bundle),
        Err(e) => invalid(
            "malformed_json",
            &format!(
                "invalid contract JSON at line {}, column {}",
                e.line(),
                e.column()
            ),
        ),
    }
}
