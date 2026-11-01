# physics-gateway

Local Rust HTTP companion for the shared TypeScript workbench. Source is integrated at
`mathphysics/Rust/physics-gateway`.
No frontend or desktop packaging is included. No API writes files or launches
models/harnesses. See `prompts/gateway_rust.prompt` for the endpoint contract.

Build/run (use an existing built frontend directory and explicit approved roots):

```sh
export CARGO_TARGET_DIR="${XDG_CACHE_HOME:-$HOME/.cache}/mathphysics-workbench/rust-target"
cargo build --offline --locked
"$CARGO_TARGET_DIR/debug/physics-gateway" \
  --static-dir /absolute/path/to/workbench/dist \
  --root papers=/absolute/path/to/papers --port 8876
```

Open http://127.0.0.1:8876. Serve the frontend through this gateway: other origins
are denied, including a separate frontend dev-server origin. Requests from
localhost must use a matching http://localhost:8876 Origin. Roots labels derive
from IDs. Files and root labels are read-only; configure only intended roots.
The static directory is trusted built application code. API source text is JSON,
and HTML files cannot be fetched through the source-text endpoint.

Endpoints: `/api/roots`, `/api/files?root=ID&path=REL&q=TEXT`,
`/api/text?root=ID&path=REL`, `/api/pdf?root=ID&path=REL`, POST `/api/tool`,
`/api/harnesses`, `/api/health`. Paths and search parameters are URL encoded.
All non-stream responses are JSON. Errors use `{error:{code,message}}`.
Harness IDs: claude-code, codex, grok-build, hermes-agent, pi; installed executables
only mean `not_connected`, never an active integration. Detection aliases are
claude, codex, grok, hermes, pi. Text kinds are tex/markdown/text; files also expose
directory/pdf/file. PDF GET/HEAD supports one byte range, including suffix ranges.
Text .tex pairing looks for same-stem lowercase `.pdf`; no metadata route added.
Modified times use integer seconds since Unix epoch (zero if unavailable).

Limits: four workers, 1MiB request bodies, 2MiB source text, 200 results,
10000 visited directory entries and approximately500ms scan. Hidden entries are
visible. Search avoids common build caches and canonical directory cycles.
Configured roots should not be concurrently modified by untrusted processes:
canonicalization and subsequent filesystem access are separate operations.
This is a local development companion, with Host/Origin defenses; it is not an
internet-facing authenticated multiuser server. Slow client socket hardening
and OS-handle-based confinement against hostile concurrent filesystem renames
are not implemented.

Verification:

```sh
cargo fmt --check
cargo clippy --offline --locked --all-targets -- -D warnings
cargo test --offline --locked
```

HTTP integration tests spawn/terminate the binary on ephemeral loopback ports,
use temporary fixtures and stdTcpStream, and require permission to bind local
sockets. Tests passed after socket sandbox escalation. No model/provider calls.
The actual Clifford kernel is exercised; exact algebra is not formal proof.

PDD provenance: request.txt, intent-plan.json (read-only planner against actual
integration target), raw planner stdout in intent-plan.raw.txt, versioned prompt,
.pddrc and architecture.json. Implementation and tests were manually authored.
Full PDD intent apply, regeneration, semantic story contract and sync remain
UNVERIFIED. The integrated dependency uses `../physics-tools`; Cargo.lock uses
a path dependency. Build output is exclusively outside the source tree.

Approved parsed-source extension: `.mmd` is markdown; `.jsonl`, `.yaml`, `.yml`,
and `.csv` are text through the same bounded UTF8 JSON endpoint. Static and PDF
responses use `X-Frame-Options: SAMEORIGIN`, permitting workbench iframe previews.
Extension verification used a separate external target directory.


From repository root, `./run_workbench.sh [--port 8876]` builds the TypeScript
workbench and gateway with offline cached dependencies and serves read-only
`mathphysics` and optional `books` roots. `--build-only` compiles without starting a server.
Frontend source must first be integrated at `TypeScript/physics-workbench`.
All frontend snapshots, npm dependencies and Rust output live in
`${XDG_CACHE_HOME:-$HOME/.cache}/mathphysics-workbench`. Snapshots remain so a
running server can continue to serve its original build. Remove old snapshots
only after their servers stop. No source-tree build files are created.

Overrides: `WORKBENCH_CACHE_DIR` (must be outside repository),
`WORKBENCH_FRONTEND_DIR`, `WORKBENCH_BOOK_ROOT`, `WORKBENCH_PORT`,
`WORKBENCH_NPM_OFFLINE` (default `1`), and optional
`WORKBENCH_NODE_MODULES` (absolute installed dependency directory matching the
frontend lockfile). Without that dependency override, `npm ci --offline` uses
the npm cache; missing cached packages fail rather than fetching the network.
The default books root is only the approved workspace corpus
`Data/Public/books/Physics/Srednicki-QuantumFieldTheory`, resolved from the
physical repository workspace. If it is missing, startup prints a notice and
serves `mathphysics` alone. An explicit `WORKBENCH_BOOK_ROOT` selects another
approved corpus; an absent override is also optional.

For a first-use cache miss, explicitly opt into online npm dependency bootstrap:
`WORKBENCH_NPM_OFFLINE=0 ./run_workbench.sh --build-only`. This runs the local
user's `npm ci` inside the external cache snapshot, with network access permitted
only by this opt-in. Default launches never silently retry online. Alternatively
provide preinstalled dependencies using `WORKBENCH_NODE_MODULES`. Rust builds
remain offline in both modes.

For frontend test fixtures, the parent owns the fixture generator and test
source. Run the gateway directly with a third explicitly approved read-only
root; no fixture root is enabled by the ordinary launcher:

```sh
"$CARGO_TARGET_DIR/debug/physics-gateway" \
  --static-dir /absolute/path/to/cached/build/dist \
  --root mathphysics=/absolute/path/to/mathphysics \
  --root books=/absolute/path/to/Data/Public/books/Physics/Srednicki-QuantumFieldTheory \
  --root fixtures=/absolute/path/to/generated/fixtures --port 8877
```

Tests can use `/api/files?root=fixtures` and
`/api/text?root=fixtures&path=example.mmd`. Generate fixtures before starting the
gateway. The gateway never creates or modifies them. Select an unused port;
this command does not replace, stop or restart any existing service.


Integration verification (2026-10-31): actual repository gateway passed all five
HTTP integration tests, offline locked Clippy with warnings denied, and format
check with target `/tmp/mathphysics-gateway-integrated-target`. Root launcher
passed `bash -n` and `--build-only` using the actual integrated frontend, default
offline npm installation, and external cache
`/tmp/mathphysics-workbench-launcher-check`. TypeScript checking, Vite build and
actual integrated Rust gateway build passed. Existing servers were not stopped
or restarted. Vite reported its existing large-bundle advisory. No frontend
source or fixture generator was changed by gateway integration.
