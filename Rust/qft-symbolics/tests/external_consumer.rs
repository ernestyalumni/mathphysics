//! A separate Cargo package qualifies public APIs without chapter/HTTP/UI dependencies.
use std::{path::PathBuf, process::Command};
#[test]
fn independent_rust_consumer_runs_public_libraries() {
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/consumers/research-client/Cargo.toml");
    let target = std::env::temp_dir().join(format!(
        "mathphysics-external-client-{}",
        std::process::id()
    ));
    let output = Command::new(env!("CARGO"))
        .args(["run", "--offline", "--locked", "--quiet", "--manifest-path"])
        .arg(manifest)
        .env("CARGO_TARGET_DIR", target)
        .output()
        .expect("run independent Cargo client");
    assert!(
        output.status.success(),
        "consumer stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(String::from_utf8_lossy(&output.stdout).contains("research_consumer_ok"));
}
