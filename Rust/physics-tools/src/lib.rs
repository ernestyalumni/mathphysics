//! Bounded native bridge, not a general expression engine or proof checker.
use serde::Deserialize;
use serde_json::{json, Value};
use std::io::Read;
use std::panic::{catch_unwind, UnwindSafe};

pub const INPUT_LIMIT: usize = 1_048_576;
pub const MAX_LABELS: usize = 12;
pub const MAX_LABEL_BYTES: usize = 64;
pub const MAX_PAIRINGS: u64 = 10_395;
pub const DOMAIN: &str = "formal_symmetric_h";

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Request {
    command: String,
    #[serde(default)]
    labels: Option<Vec<String>>,
    #[serde(default)]
    domain: Option<String>,
    #[serde(default)]
    budget: Option<u64>,
}

/// JSON response and process exit code. Only computed/discovery responses exit zero.
#[derive(Debug)]
pub struct Outcome {
    pub response: Value,
    pub exit_code: u8,
}

pub fn error(status: &str, code: &str, message: &str, exit_code: u8) -> Outcome {
    Outcome {
        response: json!({"schema_version":"physics-tools.v1", "status":status,
            "evidence_kind":null,"formal_proof":false,
            "error":{"code":code,"message":message}}),
        exit_code,
    }
}

fn success(status: &str, result: Value) -> Outcome {
    Outcome {
        response: json!({"schema_version":"physics-tools.v1", "status":status,
            "formal_proof":false,"evidence_kind": if status == "computed" {Some("exact_algebra")} else {None},
            "result":result}),
        exit_code: 0,
    }
}

pub fn discovery() -> Value {
    json!({"version":env!("CARGO_PKG_VERSION"),"runtime":"native_rust",
        "backend":"qft-symbolics 0.1.0", "transport":"single JSON stdin/stdout CLI",
        "mcp":false,"operations":["discover","gamma_trace"],
        "input_limit_bytes":INPUT_LIMIT,"max_labels":MAX_LABELS,"label_bytes":[1,MAX_LABEL_BYTES],
        "label_grammar":"[A-Za-z_][A-Za-z0-9_]*", "domain":DOMAIN,
        "budget":{"unit":"pairing_leaves_upper_bound","minimum":1,"maximum":MAX_PAIRINGS,"required":true},
        "kernel_bounds":{"arithmetic":"checked i128 rational intermediates; no float fallback",
            "max_polynomial_terms":8192,"max_polynomial_degree":64,"panic":"inconclusive"},
        "convention":{"anticommutator":"2h","srednicki":"h=-eta; eta=(-,+,+,+)",
            "trace_identity":4,"gamma5":false,"labels":"free, never implicitly summed; repeated labels identify the same formal slot"},
        "request_schema":{"discover":{"command":"discover"},"gamma_trace":{"command":"gamma_trace",
            "labels":"array of identifiers (empty word allowed)","domain":DOMAIN,"budget":"positive integer <=10395"},
            "unknown_fields":"rejected","duplicate_fields":"rejected","trailing_json":"rejected"},
        "statuses":{"success":["discovered","computed"],"errors":["invalid_input","unsupported","over_budget","inconclusive"]},
        "exit_codes":{"success":0,"invalid_input":2,"unsupported":3,"over_budget":4,"inconclusive":5},
        "limitations":["no general F09/F10 expression/proof engine","seed F12/F13 bridge; not complete F13 acceptance",
            "no formal proof","no contractions or implicit sums","no arbitrary execution, shell, providers, Python or C++"]})
}

fn valid_label(s: &str) -> bool {
    let mut bytes = s.bytes();
    (1..=MAX_LABEL_BYTES).contains(&s.len())
        && bytes
            .next()
            .is_some_and(|b| b.is_ascii_alphabetic() || b == b'_')
        && bytes.all(|b| b.is_ascii_alphanumeric() || b == b'_')
}

/// Maximum number of pairing leaves, checked even though input length is capped.
pub fn pairing_cost(n: usize) -> Option<u64> {
    if n > MAX_LABELS {
        return None;
    }
    if n % 2 == 1 {
        return Some(1);
    }
    (1..n)
        .step_by(2)
        .try_fold(1u64, |a, b| a.checked_mul(b as u64))
}

/// Execute a kernel closure behind the same panic boundary used in production.
/// The default Rust panic hook can still write diagnostic text to stderr.
fn guarded<F: FnOnce() -> Value + UnwindSafe>(f: F) -> Outcome {
    match catch_unwind(f) {
        Ok(v) => success("computed", v),
        Err(_) => error(
            "inconclusive",
            "kernel_panic",
            "Kernel arithmetic or expansion bound failed",
            5,
        ),
    }
}

// Trace produces integral real coefficients only. Convert its bounded canonical
// display into structured monomials without inspecting/changing private kernel data.
// Strict identifier grammar excludes all display separators in caller labels.
fn polynomial_terms(display: &str) -> Option<Vec<Value>> {
    if display == "0" {
        return Some(vec![]);
    }
    display
        .split(" + ")
        .map(|term| {
            let mut parts = term.split('*');
            let coefficient = parts.next()?.parse::<i128>().ok()?.to_string();
            let factors: Vec<_> = parts.map(str::to_owned).collect();
            if factors
                .iter()
                .any(|v| !v.starts_with("h(") || !v.ends_with(')'))
            {
                return None;
            }
            Some(json!({"coefficient":coefficient,"factors":factors}))
        })
        .collect()
}

pub fn handle_bytes(input: &[u8]) -> Outcome {
    if input.len() > INPUT_LIMIT {
        return error("over_budget", "input_limit", "JSON input exceeds 1 MiB", 4);
    }
    let r: Request = match serde_json::from_slice(input) {
        Ok(r) => r,
        Err(_) => {
            return error(
                "invalid_input",
                "malformed_json",
                "Expected one strict request object; unknown/duplicate fields forbidden",
                2,
            )
        }
    };
    match r.command.as_str() {
        "discover" => {
            // Reject even null optional fields on discovery (presence differs from absence).
            let keys: Value = serde_json::from_slice(input).expect("already parsed");
            if keys.as_object().is_none_or(|m| m.len() != 1) {
                return error(
                    "invalid_input",
                    "discovery_fields",
                    "Discovery accepts command only",
                    2,
                );
            }
            return success("discovered", discovery());
        }
        "gamma_trace" => (),
        _ => {
            return error(
                "unsupported",
                "operation",
                "Supported commands are discover and gamma_trace",
                3,
            )
        }
    }
    let (Some(labels), Some(domain), Some(budget)) = (r.labels, r.domain, r.budget) else {
        return error(
            "invalid_input",
            "missing_fields",
            "gamma_trace requires labels, domain and budget",
            2,
        );
    };
    if domain.trim().is_empty() {
        return error(
            "invalid_input",
            "empty_domain",
            "A nonempty domain is required",
            2,
        );
    }
    if domain != DOMAIN {
        return error(
            "unsupported",
            "domain",
            "Only formal_symmetric_h is supported",
            3,
        );
    }
    if labels.len() > MAX_LABELS {
        return error(
            "over_budget",
            "word_length",
            "At most 12 labels are supported",
            4,
        );
    }
    if !labels.iter().all(|s| valid_label(s)) {
        return error(
            "invalid_input",
            "label",
            "Labels must be 1..64 byte ASCII identifiers",
            2,
        );
    }
    if budget == 0 {
        return error("invalid_input", "zero_budget", "Budget must be positive", 2);
    }
    let cost = pairing_cost(labels.len()).expect("validated length");
    if budget > MAX_PAIRINGS || budget < cost {
        return error(
            "over_budget",
            "pairing_budget",
            "Budget must cover pairing bound and be <=10395",
            4,
        );
    }
    guarded(|| {
        let refs: Vec<_> = labels.iter().map(String::as_str).collect();
        let p = qft_symbolics::clifford::trace(&refs);
        let display = p.to_string();
        let terms = polynomial_terms(&display).expect("trace serialization invariant");
        json!({"operation":"gamma_trace","domain":domain,"labels":labels,
            "backend":"qft-symbolics 0.1.0","arithmetic":"checked_i128_exact",
            "convention":{"anticommutator":"2h","srednicki":"h=-eta; eta=(-,+,+,+)","trace_identity":4},
            "assumptions":["four-component Dirac trace without gamma5", "h is formal symmetric",
                "labels are free and never imply sums", "repeated labels do not impose metric normalization or dimension contractions"],
            "pairing_bound":cost,"budget":budget,
            "polynomial":{"coefficient_encoding":"decimal integer strings","factor_encoding":"h(byte_length:label;byte_length:label)",
                "terms":terms,"display":display}})
    })
}

/// Reads at most limit+1 bytes; never buffers an unbounded stdin stream.
pub fn handle_reader(reader: impl Read) -> Outcome {
    let mut bytes = Vec::new();
    match reader
        .take((INPUT_LIMIT + 1) as u64)
        .read_to_end(&mut bytes)
    {
        Ok(_) => handle_bytes(&bytes),
        Err(_) => error("invalid_input", "input_io", "Cannot read request", 2),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn panic_is_inconclusive() {
        let r = guarded(|| panic!("simulated checked overflow"));
        assert_eq!(r.response["status"], "inconclusive");
        assert_eq!(r.exit_code, 5);
        assert!(r.response["evidence_kind"].is_null());
    }
    #[test]
    fn checked_cost() {
        assert_eq!(pairing_cost(12), Some(10395));
        assert_eq!(pairing_cost(usize::MAX), None);
    }
}
