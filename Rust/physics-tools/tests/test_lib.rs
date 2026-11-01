use physics_tools::{handle_bytes, handle_reader, DOMAIN, INPUT_LIMIT};
use qft_symbolics::{matrices as m, scalar::C};
use serde_json::{json, Value};
use std::collections::BTreeMap;

fn request(labels: &[&str]) -> Value {
    json!({"command":"gamma_trace","labels":labels,"domain":DOMAIN,"budget":10395})
}
fn run(v: Value) -> physics_tools::Outcome {
    handle_bytes(&serde_json::to_vec(&v).unwrap())
}
fn metric(a: &str, b: &str) -> String {
    let (a, b) = if a <= b { (a, b) } else { (b, a) };
    format!("h({}:{};{}:{})", a.len(), a, b.len(), b)
}
fn evaluate(v: &Value, values: &BTreeMap<String, i128>) -> i128 {
    v["result"]["polynomial"]["terms"]
        .as_array()
        .unwrap()
        .iter()
        .map(|t| {
            t["factors"].as_array().unwrap().iter().fold(
                t["coefficient"].as_str().unwrap().parse::<i128>().unwrap(),
                |a, f| a * values[f.as_str().unwrap()],
            )
        })
        .sum()
}
#[test]
fn independent_zero_two_four_formula_with_new_labels() {
    let labels = ["new_time", "space_zz", "NeverSeen", "index_42"];
    let mut values = BTreeMap::new();
    for (i, a) in labels.iter().enumerate() {
        for (j, b) in labels.iter().enumerate() {
            values.insert(metric(a, b), ((i + j + 1) * (i + j + 1)) as i128 - 3);
        }
    }
    let h = |i: usize, j: usize| values[&metric(labels[i], labels[j])];
    for (word, expected) in [
        (&labels[..0], 4),
        (&labels[..2], 4 * h(0, 1)),
        (
            &labels[..4],
            4 * (h(0, 1) * h(2, 3) - h(0, 2) * h(1, 3) + h(0, 3) * h(1, 2)),
        ),
    ] {
        let out = run(request(word));
        assert_eq!(out.exit_code, 0);
        assert_eq!(out.response["status"], "computed");
        assert_eq!(out.response["evidence_kind"], "exact_algebra");
        assert_eq!(evaluate(&out.response, &values), expected);
    }
}
#[test]
fn explicit_matrices_exhaustive_two_four_and_repeated_slots() {
    let names = ["t", "x", "y", "z"];
    let mut values = BTreeMap::new();
    for i in 0..4 {
        for j in 0..4 {
            values.insert(
                metric(names[i], names[j]),
                if i == j { -m::ETA[i] as i128 } else { 0 },
            );
        }
    }
    for n in [0, 2, 4] {
        for code in 0..4usize.pow(n) {
            let mut k = code;
            let mut word = Vec::new();
            let mut product = m::identity::<4>();
            for _ in 0..n {
                let i = k % 4;
                k /= 4;
                word.push(names[i]);
                product = m::mul(product, m::gamma(i));
            }
            let out = run(request(&word));
            let actual = evaluate(&out.response, &values);
            assert_eq!(
                C::from(qft_symbolics::scalar::Q::new(actual, 1)),
                m::trace(product),
                "{word:?}"
            );
        }
    }
    // Repeated free slot is h(a,a), not an implicit four-dimensional sum or unit norm.
    values.insert(metric("a", "a"), 7);
    let out = run(request(&["a", "a"]));
    assert_eq!(evaluate(&out.response, &values), 28);
    let out = run(request(&["a", "a", "a", "a"]));
    assert_eq!(evaluate(&out.response, &values), 196);
}
#[test]
fn strict_parser_and_labels() {
    for bytes in [br#"{"command":"discover","command":"discover"}"#.as_slice(),
        br#"{"command":"discover","x":1}"#,br#"{"command":"discover"} {}"#,
        br#"{"command":"gamma_trace","labels":null,"labels":[],"domain":"formal_symmetric_h","budget":1}"#,
        br#"{"command":"gamma_trace","labels":[],"domain":"formal_symmetric_h","budget":1,"budget":1}"#,
        br#"[]"#,br#"{"command":"discover","labels":null}"#] {
        assert_eq!(handle_bytes(bytes).response["status"],"invalid_input","{bytes:?}");
    }
    for s in [
        "",
        "bad*label",
        "h(a,b)",
        "a b",
        "1bad",
        "μ",
        &"a".repeat(65),
    ] {
        assert_eq!(run(request(&[s])).response["error"]["code"], "label");
    }
    assert_eq!(run(request(&[&"a".repeat(64)])).exit_code, 0);
}
#[test]
fn budgets_domains_odd_and_discovery() {
    assert_eq!(
        run(json!({"command":"discover"})).response["result"]["mcp"],
        false
    );
    assert_eq!(
        run(json!({"command":"exec"})).response["status"],
        "unsupported"
    );
    for (key, value, status) in [
        ("domain", json!(" "), "invalid_input"),
        ("domain", json!("euclidean"), "unsupported"),
        ("budget", json!(0), "invalid_input"),
        ("budget", json!(10396), "over_budget"),
        ("budget", json!(2), "over_budget"),
        ("budget", json!(-1), "invalid_input"),
        ("budget", json!(1.5), "invalid_input"),
    ] {
        let mut v = request(&["a", "b", "c", "d"]);
        v[key] = value;
        let out = run(v);
        assert_ne!(out.exit_code, 0);
        assert_eq!(out.response["status"], status);
    }
    assert_eq!(
        run(request(&["a"; 13])).response["error"]["code"],
        "word_length"
    );
    assert_eq!(run(request(&["a"; 12])).exit_code, 0);
    assert_eq!(
        run(request(&["a"; 3])).response["result"]["polynomial"]["terms"],
        json!([])
    );
    let mut v = request(&[]);
    v.as_object_mut().unwrap().remove("domain");
    assert_eq!(run(v).response["error"]["code"], "missing_fields");
}
#[test]
fn bounded_reader_and_exact_limit() {
    let mut bytes = serde_json::to_vec(&json!({"command":"discover"})).unwrap();
    bytes.resize(INPUT_LIMIT, b' ');
    assert_eq!(handle_reader(bytes.as_slice()).exit_code, 0);
    bytes.push(b' ');
    assert_eq!(
        handle_reader(bytes.as_slice()).response["status"],
        "over_budget"
    );
    struct Infinite;
    impl std::io::Read for Infinite {
        fn read(&mut self, b: &mut [u8]) -> std::io::Result<usize> {
            b.fill(b' ');
            Ok(b.len())
        }
    }
    assert_eq!(
        handle_reader(Infinite).response["error"]["code"],
        "input_limit"
    );
}
