use serde_json::{json, Value};
use std::{
    io::Write,
    process::{Command, Stdio},
};
fn cli(input: &[u8], args: &[&str]) -> (i32, Value) {
    let mut child = Command::new(env!("CARGO_BIN_EXE_physics-tools"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    // Separate writer avoids pipe deadlock for large requests.
    let mut stdin = child.stdin.take().unwrap();
    let bytes = input.to_vec();
    let writer = std::thread::spawn(move || {
        let _ = stdin.write_all(&bytes);
    });
    let out = child.wait_with_output().unwrap();
    writer.join().unwrap();
    (
        out.status.code().unwrap(),
        serde_json::from_slice(&out.stdout).unwrap(),
    )
}
#[test]
fn discovery_and_actual_supplied_words() {
    let (code, v) = cli(br#"{"command":"discover"}"#, &[]);
    assert_eq!(code, 0);
    assert_eq!(v["result"]["runtime"], "native_rust");
    for n in [0, 2, 4] {
        let input = json!({"command":"gamma_trace","labels":vec!["brandNew";n],"domain":"formal_symmetric_h","budget":3});
        let (code, v) = cli(&serde_json::to_vec(&input).unwrap(), &[]);
        assert_eq!(code, 0);
        assert_eq!(v["status"], "computed");
        let terms = v["result"]["polynomial"]["terms"].as_array().unwrap();
        assert_eq!(terms.len(), 1);
        assert_eq!(terms[0]["coefficient"], "4");
        assert_eq!(terms[0]["factors"].as_array().unwrap().len(), n / 2);
    }
}
#[test]
fn distinct_errors_and_nonzero_exits() {
    for (input,code,status) in [(br#"{"command":"exec"}"#.as_slice(),3,"unsupported"),(br#"{"command":"discover","bad":1}"#,2,"invalid_input"),(br#"{"command":"gamma_trace","labels":["a","b","c","d"],"domain":"formal_symmetric_h","budget":1}"#,4,"over_budget"),(br#"{"command":"discover","command":"discover"}"#,2,"invalid_input")] {
        let (c,v)=cli(input,&[]);assert_eq!(c,code);assert_eq!(v["status"],status);assert_eq!(v["formal_proof"],false);
    }
    let (c, v) = cli(br#"{"command":"discover"}"#, &["shell"]);
    assert_eq!(c, 2);
    assert_eq!(v["error"]["code"], "arguments");
    let (c, v) = cli(&vec![b' '; physics_tools::INPUT_LIMIT + 1], &[]);
    assert_eq!(c, 4);
    assert_eq!(v["error"]["code"], "input_limit");
}

#[test]
fn real_kernel_expansion_failure_is_inconclusive() {
    // 11!! = 10395 distinct monomials exceeds the unchanged kernel's 8192 term bound.
    let labels: Vec<_> = (0..12).map(|i| format!("fresh_{i}")).collect();
    let input = json!({"command":"gamma_trace","labels":labels,"domain":"formal_symmetric_h","budget":10395});
    let (code, v) = cli(&serde_json::to_vec(&input).unwrap(), &[]);
    assert_eq!(code, 5);
    assert_eq!(v["status"], "inconclusive");
    assert_eq!(v["error"]["code"], "kernel_panic");
    assert!(v["evidence_kind"].is_null());
}

#[test]
fn cli_label_length_and_domain_boundaries() {
    let base = json!({"command":"gamma_trace","labels":["newLabel"],"domain":"formal_symmetric_h","budget":10395});
    for (key, value, code, error) in [
        ("labels", json!([""]), 2, "label"),
        ("labels", json!(["a*b"]), 2, "label"),
        ("labels", json!(["a".repeat(65)]), 2, "label"),
        ("labels", json!(vec!["free"; 13]), 4, "word_length"),
        ("domain", json!(""), 2, "empty_domain"),
        ("domain", json!("summed_indices"), 3, "domain"),
        ("budget", json!(0), 2, "zero_budget"),
    ] {
        let mut input = base.clone();
        input[key] = value;
        let (c, v) = cli(&serde_json::to_vec(&input).unwrap(), &[]);
        assert_eq!(c, code);
        assert_eq!(v["error"]["code"], error);
    }
}
