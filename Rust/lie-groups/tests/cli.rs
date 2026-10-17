use std::process::Command;
#[test]
fn cli_distinguishes_verified_wrong_and_unsupported() {
    for (args, code, status) in [
        (vec!["verify"], 0, "PASS"),
        (vec!["verify", "--wrong-sign"], 1, "FAIL"),
        (vec!["missing-command"], 2, "INCONCLUSIVE"),
        (vec!["symbolic", "extra"], 2, "INCONCLUSIVE"),
    ] {
        let result = Command::new(env!("CARGO_BIN_EXE_lie-groups"))
            .args(args)
            .output()
            .unwrap();
        assert_eq!(result.status.code(), Some(code));
        assert!(String::from_utf8(result.stdout)
            .unwrap()
            .contains(&format!("\"status\":\"{status}\"")));
    }
}
