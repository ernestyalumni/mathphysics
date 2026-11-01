use std::{
    io::Write,
    process::{Command, Stdio},
};
fn invoke(args: &[&str], input: &[u8]) -> std::process::Output {
    let mut c = Command::new(env!("CARGO_BIN_EXE_physics-contracts"))
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    c.stdin.take().unwrap().write_all(input).unwrap();
    c.wait_with_output().unwrap()
}
#[test]
fn real_cli_validates_stdin_and_rejects_overclaims_with_structured_errors() {
    let good = invoke(
        &["validate"],
        include_bytes!("../examples/valid-bundle.json"),
    );
    assert!(good.status.success());
    let v: serde_json::Value = serde_json::from_slice(&good.stdout).unwrap();
    assert_eq!(v["status"], "record_valid");
    assert!(v["assurance"]
        .as_str()
        .unwrap()
        .contains("metadata_validation_only"));
    let bad = invoke(
        &["validate", "-"],
        include_bytes!("../examples/invalid-numerical-proof.json"),
    );
    assert_eq!(bad.status.code(), Some(2));
    let v: serde_json::Value = serde_json::from_slice(&bad.stdout).unwrap();
    assert_eq!(v["status"], "invalid_input");
}
#[test]
fn missing_files_and_bad_commands_have_no_false_success() {
    for args in [
        &["validate", "/no/such/physics-record"][..],
        &["prove"][..],
        &["validate", "-", "extra"][..],
    ] {
        let out = invoke(args, b"");
        assert_eq!(out.status.code(), Some(2));
        let v: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
        assert_eq!(v["status"], "invalid_input");
    }
}
