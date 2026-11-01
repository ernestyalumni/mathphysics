# Reusable Rust mathematical physics research library

## Story

As a mathematical physics researcher using mathphysics, I want Srednicki’s Quantum Field Theory to drive a reusable, extensible Rust library of mathematical and physical capabilities, so that the same tested operations can reproduce textbook calculations and support active QFT research, other books and papers—especially versioned arXiv papers—without rewriting the mathematics for every subsection or source.

The book is a requirements and validation corpus, not the public API hierarchy. Progress means reusable capabilities with explicit limits and independent evidence, not a growing count of chapter-specific functions. The goal is production research use; that status must be earned rather than inferred from Rust implementation or passing textbook fixtures.

## Requirements

- Production mathematics MUST be organized by coherent capabilities, such as coefficient arithmetic, expressions, tensor/index algebra, Clifford algebra, spinors, group representations, kinematics, diagrams, amplitudes and verification. MUST NOT generate a separate production solver/module/function for every chapter, subsection, equation, exercise or paper when it duplicates a shared capability.
- A new source item MUST first map to existing public APIs. If something is missing, specify the smallest mathematically justified reusable extension, inspect existing ownership and consumers, and implement it once at the correct library boundary. Avoid both duplicated solvers and speculative all-purpose abstractions unsupported by actual calculations.
- Textbook/paper item IDs, expected results, source locations and provenance MUST live in corpus data, examples and independent acceptance fixtures. Thin source-specific consumers MAY assemble inputs, translate declared conventions and call shared APIs; they MUST NOT hide a second implementation of the algebra. Chapter-numbered test names are useful and are not themselves a duplication defect.
- Source citations and named algorithms or convention profiles MAY appear in documented APIs where they identify real mathematics. For example, BCFW, Parke–Taylor or a Srednicki convention adapter is not automatically forbidden. The prohibited dependency is a solver that knows a subsection/answer instead of operating on mathematical inputs.
- Public Rust APIs MUST accept caller-supplied mathematical objects and expose domains, assumptions, conventions, exact versus approximate behavior, resource limits and structured failures. No hardcoded expected answers, silent normalization, fixed author symbols or branch-by-problem-number. Unsupported input and exhaustion MUST NOT become verified results.
- CLI, HTTP/MCP, AI harnesses and both frontends MUST be clients of the libraries. They MUST NOT independently reimplement production physics. Transport, rendering, textbook loading and model access MUST NOT be required for a Rust client to perform a supported calculation.
- New books/papers MUST be onboardable through versioned source records and thin input/expectation fixtures, without a new source-specific engine. A paper identifier/version or content digest, exact claim, conventions, assumptions and evidence strength MUST remain available. Bibliographic metadata and matching output do not prove a paper’s claims.
- Production readiness MUST require documented/stable-enough interfaces, dependency/convention compatibility, reproducible builds, error/resource contracts, meaningful independent tests and workload-appropriate performance evidence. Use idiomatic Rust ownership, composition, traits/generics where justified, module visibility and Result-based boundaries. Do not literally port C++ inheritance, header organization or unsafe ownership patterns.
- Existing HrdwCCppCUDA and Propulsion code SHOULD inform reusable interfaces, implementation/test separation, numerical policies and independent validation; prefer recently maintained examples judged by Git history. Inspect actual code and tradeoffs. Do not copy code blindly, infer correctness from recency, or mandate a framework solely because a reference uses it.
- Independent test/oracle mathematics MAY intentionally repeat a derivation using a different representation. MUST NOT eliminate that independence in the name of deduplication. Production code reuse and independent verification serve different purposes.
- Preserve existing correct libraries, callers, notebooks, exploration state and source manuscripts. Any extraction, convention unification or API replacement MUST have characterization tests and a compatibility/migration plan. A textbook convention adapter must not silently change signs when moved into a library.

## Independent acceptance scenarios

1. **External Rust consumer:** A small crate outside the source/corpus tree depends on a public mathphysics library and computes a supported Clifford/tensor/spinor/group operation using fresh caller-supplied inputs. It needs no Srednicki files, chapter number, HTTP service, frontend, Python/C++ runtime or AI account.
2. **One capability, several consumers:** The same implementation is invoked by a Srednicki fixture, an independent non-Srednicki fixture and a CLI/agent or interactive client. Distinct symbol names, boundary cases and unsupported domains expose memorized answers or hidden assumptions. A genuinely source-reviewed paper case is identified as such; a synthetic research-style fixture is not mislabeled as an arXiv reproduction.
3. **No subsection duplication:** For a new item requiring an already supported operation, the change adds only source data and a thin fixture/example. For a new operation, reviewers can identify its single owning module and show how multiple inputs or consumers use it. Chapter coverage cannot be increased by cloning an algorithm.
4. **Independent mathematical validation:** Expected identities are established separately from the implementation path, with negative controls, assumption/convention changes, singular inputs and numerical/error budgets where applicable. Round trips or agreement between wrappers sharing one kernel are not independent proof.
5. **Compatibility and failures:** Existing public behavior stays characterized across extraction/refactoring. Checked-range overflow, unsupported regularization, ambiguous conventions and exhausted budgets produce explicit failures/inconclusive outcomes. No float fallback or default convention silently changes the mathematical question.
6. **Research claim/evidence distinction:** Reproducing an equation, finding a numerical witness, checking an identity in a declared algebra, replaying a certificate and obtaining a formal proof remain distinct. Production-capable tooling must not label unproved paper claims formally proved.

## Explainable verification acceptance

A researcher MUST be able to run a repository-owned walkthrough that names each
mathematical claim, invokes the public Rust operation, explains how an independent
expected result is obtained, and demonstrates a meaningful negative control.
A passing test, an exact symbolic identity, an evaluated witness, a conditionally
replayed derivation and a formal proof MUST remain visibly distinct. Resource
limits, pole restrictions, ordered products and convention translations MUST NOT
be hidden by a success message. Test names/counts MUST NOT imply chapter completion.
The mono notes MUST remain independent of implementation details; mappings and
explanations belong in implementation documentation near the owning library.

## Lorentz and operator-calculus acceptance

Reusable rapidity MUST support declared numerical tolerances, beta/eta conversion,
collinear composition and inverse boosts in Srednicki's (-+++) convention with
positive off-diagonal entries. Exact rational parameterizations MUST be distinguished
from an assertion that rapidity itself is rational. Symbolic conditions MUST remain
explicit. Metric/index identities and independent expected components must catch
wrong signatures or boost signs.

Chapter-specific fixtures MAY test first-order coefficient extraction, ordered
commutators, associative Jacobi identities and differential product rules. They MUST
state group-law/covariance assumptions and distinguish representation identities
from general unbounded-operator proofs. Unitarity alone MUST NOT imply commuting
translations, and a term invisible to field commutators MUST NOT be erased or
inferred absent. Those mathematical limitations remain part of the result.

## Current boundary

The current qft-symbolics and lie-groups crates already expose reusable mathematical modules. Existing named Srednicki checks remain valid source-specific acceptance cases. The first checked arithmetic, typed structural-expression and independently replayed scalar-algebra slices are implemented. They do not establish a full production CAS, arbitrary operator/tensor rewriting, a general parser, a completed arXiv prover or formal QFT proof. See the owning library documentation and execution evidence.

See [reuse architecture and reference audit](../docs/plans/rust-physics-platform/REUSABLE-LIBRARY.md). The source manuscripts and reference repositories are unchanged.
