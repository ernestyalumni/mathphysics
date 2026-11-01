use serde_json::Value;
use std::process::Command;
fn run(args: &[&str]) -> (i32, Value) {
    let out = Command::new(env!("CARGO_BIN_EXE_qft-symbolics"))
        .args(args)
        .output()
        .unwrap();
    (
        out.status.code().unwrap(),
        serde_json::from_slice(&out.stdout).unwrap(),
    )
}
#[test]
fn large_rationals_and_zero_errors_are_real_cli_results() {
    let (code, v) = run(&[
        "rational",
        "add",
        "340282366920938463463374607431768211456",
        "3",
        "2",
        "3",
    ]);
    assert_eq!(code, 0);
    assert_eq!(v["result"], "113427455640312821154458202477256070486");
    assert_eq!(v["formal_proof"], false);
    for args in [
        ["rational", "div", "1", "2", "0", "1"],
        ["rational", "add", "1", "0", "1", "2"],
    ] {
        let (code, v) = run(&args);
        assert_eq!(code, 2);
        assert_eq!(v["status"], "INVALID_INPUT");
        assert_eq!(v["formal_proof"], false);
    }
}
#[test]
fn checked_trace_reports_budgets_without_panicking() {
    let (code, v) = run(&["gamma-trace-checked", "u", "v"]);
    assert_eq!(code, 0);
    assert_eq!(v["result"], "4*h(1:u;1:v)");
    let args = std::iter::once("gamma-trace-checked")
        .chain(std::iter::repeat_n("x", 13))
        .collect::<Vec<_>>();
    let (code, v) = run(&args);
    assert_eq!(code, 2);
    assert_eq!(v["status"], "INCONCLUSIVE");
    assert_eq!(v["error"]["code"], "resource_limit");
}
