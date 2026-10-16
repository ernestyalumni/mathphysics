use std::process::Command;
fn cli(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_qft-symbolics"))
        .args(args)
        .output()
        .unwrap()
}
#[test]
fn actual_cli_positive_and_negative_controls() {
    let good = cli(&["schouten"]);
    assert_eq!(good.status.code(), Some(0));
    assert!(String::from_utf8(good.stdout)
        .unwrap()
        .contains("\"status\":\"PASS\""));
    let bad = cli(&["schouten", "--wrong-sign"]);
    assert_eq!(bad.status.code(), Some(1));
    assert!(String::from_utf8(bad.stdout)
        .unwrap()
        .contains("\"status\":\"FAIL\""));
}
#[test]
fn actual_cli_budget_and_invalid_input_never_pass() {
    for args in [
        vec!["gamma-trace"; 14],
        vec!["smga-fixture", "2", "0"],
        vec!["unknown"],
    ] {
        let row = cli(&args);
        assert_eq!(row.status.code(), Some(2));
        assert!(String::from_utf8(row.stdout)
            .unwrap()
            .contains("\"status\":\"INCONCLUSIVE\""));
    }
}
#[test]
fn actual_cli_json_escapes_untrusted_labels() {
    let row = cli(&["gamma-trace", "a\"\\", "b"]);
    assert_eq!(row.status.code(), Some(0));
    let output = String::from_utf8(row.stdout).unwrap();
    assert!(output.contains("\\\""));
    assert!(output.contains("\\\\"));
}
