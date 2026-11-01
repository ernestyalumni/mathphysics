use physics_contracts::{validate_json, ASSURANCE, MAX_INPUT_BYTES, SCHEMA_VERSION};
use std::{
    env,
    fs::File,
    io::{self, Read},
    process,
};
fn main() {
    let args: Vec<_> = env::args().skip(1).collect();
    let result: Result<Vec<u8>, String> = (|| {
        if args.first().map(String::as_str) != Some("validate") || !(1..=2).contains(&args.len()) {
            return Err("usage: physics-contracts validate [FILE|-]; stdin is the default".into());
        }
        let mut bytes = Vec::new();
        match args.get(1).map(String::as_str) {
            None | Some("-") => io::stdin()
                .lock()
                .take((MAX_INPUT_BYTES + 1) as u64)
                .read_to_end(&mut bytes),
            Some(path) => File::open(path)
                .map_err(|_| "input file cannot be opened")?
                .take((MAX_INPUT_BYTES + 1) as u64)
                .read_to_end(&mut bytes),
        }
        .map_err(|_| "input cannot be read")?;
        Ok(bytes)
    })();
    match result {
        Ok(bytes) => {
            let report = validate_json(&bytes);
            println!(
                "{}",
                serde_json::to_string(&report).expect("report uses serializable fields")
            );
            process::exit(if report.errors.is_empty() { 0 } else { 2 });
        }
        Err(message) => {
            println!(
                "{}",
                serde_json::json!({"schema_version":SCHEMA_VERSION,"status":"invalid_input","assurance":ASSURANCE,"errors":[{"path":"$","code":"input_error","message":message}]})
            );
            process::exit(2);
        }
    }
}
