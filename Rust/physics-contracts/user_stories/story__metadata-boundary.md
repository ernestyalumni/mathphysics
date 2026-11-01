<!-- pdd-story-prompts: prompts/lib_rust.prompt, prompts/main_rust.prompt -->
# Mathematical evidence metadata boundary

## Story
As a mathematical physics researcher using AI tools, I want malformed, inconsistent or overstated evidence records rejected while valid metadata remains explicitly distinguished from mathematical proof, so that an agent cannot treat bookkeeping success as verified physics.

## Accepted source
F01 in the approved platform plan; the user authorized implementation and subagents on 2026-10-31. The bounded request is preserved in `../docs/accepted-request.txt`.

## Acceptance evidence
The independent tests in `tests/test_lib.rs` and `tests/test_main.rs` exercise unknown references, broken links, numerical overclaims, source locators, formal-record declarations and input/resource failures. The JSON schemas are a structural companion; the Rust library checks cross-record consistency. A fully well-formed but false proof declaration must never cause this tool to report a checked mathematical proof.

No automatically generated story contract or semantic-model validation is claimed. This story preserves the approved independent acceptance meaning; a future generated contract must remain derived output.
