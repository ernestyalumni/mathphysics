//! Single-request JSON stdin/stdout CLI; diagnostics belong on stderr.
use std::io::{self, Write};
use std::process::ExitCode;
fn main() -> ExitCode {
    let result = if std::env::args_os().len() != 1 {
        physics_tools::error(
            "invalid_input",
            "arguments",
            "Supply one JSON request on stdin; no CLI arguments",
            2,
        )
    } else {
        physics_tools::handle_reader(io::stdin().lock())
    };
    let mut stdout = io::stdout().lock();
    if serde_json::to_writer(&mut stdout, &result.response).is_err() || writeln!(stdout).is_err() {
        return ExitCode::from(5);
    }
    ExitCode::from(result.exit_code)
}
