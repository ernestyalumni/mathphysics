# lattice-gauge: prompt-driven development status

Read `../../PDD.md` and the Monoclaw playbooks it names first.

Four units are in this crate's prompt graph (`architecture.json`, `.pddrc`): `torus`, `cochain`, `z2_gauge`, `spectral`. They implement the first slice of the Yang–Mills tooling roadmap (`../../docs/plans/yang-mills-mass-gap/05-TOOLING-ROADMAP.md`, units T05 and the start of T10).

**How they were produced.** Intent `yang-mills-construction-audit-rust-lattice-gauge-d3c6da77` was approved by Ernest on 2026-11-15. `pdd intent apply` recorded the intent in `docs/PRODUCT_INTENT.md` and `docs/intents/`, then entered its provider-backed architecture workflow, which had completed 4 of 19 steps when a 25-minute limit stopped it. It generated no prompts and no code. The units here were therefore written by hand, prompt first, as in `Rust/qft-symbolics`. No regeneration from the prompts has been attempted; do not claim one.

**Rules for changing this crate.** Edit the prompt, then the code, then the independent test. Tests in `tests/` are oracles and are not generated output. No crates.io dependency without Ernest's approval; the crate currently has none.
