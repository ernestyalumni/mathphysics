# Mathematical item and evidence metadata schemas

Four sibling Draft 2020-12 schemas define bundle, item, evidence and shared types for Rust/physics-contracts. Register the four $id resources locally; the identifying URLs are not hosted endpoints. JSON shapes alone cannot check registry membership, complete assumption equality, bidirectional evidence links or aggregate limits; use the Rust CLI for those checks.

Record validation is not proof, digest authentication, source certification or artifact replay. Source-review and outcome fields are declarations. Evidence records cannot promote approximate or external-comparison results into established exact claims.

Independent schema tests passed in the first implementation wave, including 432 evidence-kind/arithmetic/outcome/budget combinations. The existing Python JSON Schema checker was used as a development oracle only; the runtime validator is entirely Rust. Reproduction scripts, fixtures and logs are retained outside the source repository in the wave evidence directory recorded by docs/verification/PLATFORM-WAVE1.md.
