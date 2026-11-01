> Superseded integration: Explore now launches independently. See [standalone separation](SEPARATION.md). Earlier shared-shell behavior/evidence below is historical.

# Game exploration delivery plan

Status: first playable SU(2) browser slice implemented. See [implementation and evidence](IMPLEMENTATION.md). Remaining destinations and RTS are proposed. This plan adds work alongside existing Srednicki and notebook tasks. The user requested a new story, not abandonment of either workstream.

## Experience and style

Treat `experience = notebook | explore | strategy (future)` independently of `style = notebook | pixel`. The first two styles already exist. Experiences can differ in layout, controls and useful capabilities; they share durable mathematical objects and source/evidence links. Avoid an immediate 2×3 mode explosion in navigation: one clearly named Explore entry is sufficient for the first prototype. The storage migration and transitions need tests before shipping.

First concept: a small top-down observatory with an SU(2) workbench, a contextual tool belt (inspect, manipulate, slice, animate, compare, record) and a journal. Stardew Valley suggests welcoming places and persistent experiments; Zelda suggests contextual tools and discovery; Super Mario World suggests a legible destination map, short interaction loops and clear feedback. These are design interpretations of the user's references, not evidence of measured usability. An avatar and traversal are optional interaction experiments. Offer direct access to each laboratory.

Future strategy mode may explore a command map of experiments, agent assignments, task queues and dependencies, inspired by StarCraft/Age of Empires. It is a separate proposal, not a copied toolbar or part of the first game's acceptance. GTA's broad open-world reference remains in the story for later exploration.

## Proposed SU(2) conventions and first lab

Freeze these in a future Rust prompt only after checking existing Lie-group interfaces; a story does not silently alter those libraries.

Let σ₁,σ₂,σ₃ be the usual Pauli matrices and define

```
U = a0 I + i(a1 σ1 + a2 σ2 + a3 σ3)
a0²+a1²+a2²+a3² = 1
U = [[alpha, beta], [-conj(beta), conj(alpha)]]
alpha = a0+i a3, beta = a2+i a1
```

For a unit axis n and dimensionless rotation angle t in radians:

```
U(t) = exp(-i t n·σ/2)
a0 = cos(t/2), a = -n sin(t/2)
U(0)=I, U(pi)=-i n·σ, U(2pi)=-I, U(4pi)=I
```

This plus-i coordinate representation has a minus-cross-product composition rule. Do not directly reuse scalar-first Hamilton quaternion multiplication without an explicit tested conversion. Independently check `(i σ1)(i σ2)=-i σ3`, the reverse order, and matrix multiplication order. Keep U and -U distinct; an SO(3) sign canonicalizer would destroy this lab's central distinction.

Use one versioned experiment state for complex-coordinate selection, matrix, path, t and views. Display alpha and beta in linked complex planes, including the coupled norm constraint; these are coordinate views, not a claim that SU(2) is itself a single complex plane or that an arbitrary chosen field is holomorphic. Later complex-analysis tools (phase/magnitude coloring, analytic functions, branch cuts) need their own explicitly named functions/domains and legends.

Proposed views:

- A 3D stereographic chart x=a/(1-a0) of S³ excluding I; show the missing pole and off-chart status. Near-pole visual limits must be labeled, not silently treated as finite mathematical coordinates.
- A coordinate slice a0=c with radius sqrt(1-c²), -1≤c≤1. At ±1 it degenerates to a point; outside that range there is no slice. Distinguish a selected state lying off the slice from a missing state.
- Optional Bloch sphere S² with ρ↦UρU†, clearly distinguished from SU(2)≅S³ and the double cover of SO(3). U and -U produce the same Bloch action while retaining distinct group coordinates. For n=z, a +x Bloch vector initially moves toward +y.

Reject non-finite values and invalid axes/matrices. An explicit user action may normalize an axis with its changed values visible; never hide invalid inputs through automatic projection. Declare tolerances separately for mathematical invariants and GPU visualization error. A 4π path endpoint must remain inspectable before looping. Compute state from absolute t, not accumulated frame increments. Physical time needs a separately declared generator, units and angular rate.

## Work packages

| ID | Exclusive proposed scope | Depends on | Deliverable / acceptance |
|---|---|---|---|
| GX01 | Experience-state contract and acceptance fixtures | current workbench | Versioned experiment data, style/experience distinction, migration, transitions preserve notebook drafts. |
| GX02 | Rust SU(2) operations and bounded gateway endpoint | GX01 | Audit existing LieGroups, add only missing operations; independent Pauli fixtures, invariants, composition, path checkpoints, invalid-input tests. |
| GX03 | TypeScript linked mathematical views | GX01, GX02 | Complex planes, labeled 3D chart, slices and t controls; poles/branches, stale-result handling, keyboard readouts. |
| GX04 | Exploration world and contextual interactions | GX01 | Small original map, direct lab entry, tool belt, optional discovery goals; no forced platforming or access unlocks. Can prototype in parallel with GX02/03 using clearly labeled fixtures. |
| GX05 | Experiment journal and notebook/source integration | GX01, GX03 | Pin/reopen/compare exact state and evidence; source/PDF links; migration and interrupted-save checks. |
| GX06 | Independent math, browser and usability review | GX02–GX05 | Verify numerical truth separately from visual/game progress, replay/frame-rate independence, keyboard/reduced-motion/sound controls; observe users and revise. |
| GX07 | Later RTS experience design spike | GX06 | Separate user review of experiment/agent command UX; no implication of implemented live agent connections. |

Proposed paths must be frozen by the integrator after inspecting active work. Do not delegate overlapping writes to main.ts/model.ts, gateway routes or existing LieGroups. A suitable first parallel round is Rust mathematics (GX02), graphics (GX03) and world interaction (GX04), after GX01 freezes their interface. Independent validation must use fixtures beyond implementation round trips. Budget geometry resolution, frame time, memory and cancellation before increasing complexity. Record actual target-device performance before choosing GPU compute.

## Reviewed interactive-demo references

Inspected public source/documentation 2026-10-31; references, not installed dependencies or imported assets:

- [MathBox](https://github.com/unconed/mathbox), [primitives](https://github.com/unconed/mathbox/blob/master/docs/primitives.md): linked mathematical graphing and stereographic/4D view primitives. Candidate inspiration for projections and animated coordinate systems; MIT shown at repository root. Recheck a pinned revision and complete asset notices before adoption.
- [CindyJS](https://github.com/CindyJS/CindyJS): interactive mathematical content; study direct manipulation of mathematical objects rather than building a conventional form for every action. No reuse decision yet.
- [PhET Joist](https://github.com/phetsims/joist): simulation navigation and framework integration. Study focused simulation screens, discoverability and accessible interaction. License/build compatibility must be checked for specific reused files.
- [Three.js examples](https://github.com/mrdoob/three.js/tree/dev/examples): candidate rendering interaction examples; inspect a pinned clipping/slicing example before adoption. Not a commitment to a new engine dependency.
- [Mathigon Polypad API documentation](https://github.com/mathigon/documentation/blob/main/polypad/index.md): manipulatives, sliders and serialized interactive canvases. [Mathigon's licensing FAQ](https://id.mathigon.org/faqs) distinguishes non-commercial use and licensing for other uses; do not assume source visibility grants unrestricted product reuse. Design reference only.

No third-party renderer selected, code downloaded into the product or game assets copied in this story change. No CUDA dependency added. These web references support choosing candidates to inspect, not claims that any demo already implements this SU(2) experience.

## PDD and verification status

The exact request and successful read-only plan are preserved here. Intent ID: `game-like-mathematical-physics-exploration-d51a8097`. The planner identified the existing PDD subproject but extracted no concrete constraints/examples and asked a generic scope question. The story preserves the complete request; agent inspection resolves scope to the existing workbench plus proposed Rust math/gateway extensions. No need to ask the user to choose paths or story slugs.

The initial story-only turn was followed by explicit implementation approval. The first runtime slice now has manual prompt-first implementation and independent browser/mathematical checks in IMPLEMENTATION.md. No intent apply, generated story contract, provider generation or whole-project regeneration is claimed. Future behavior changes begin in the corresponding prompt.

Independent mathematical design review identified complex-pair signs, reversed quaternion composition, path sign, projection poles, SU(2)/Bloch distinction and clock-driven playback as concrete risks. The proposed formulas and acceptance scenarios account for them; this is design review, not executable proof or runtime verification.

The subsequent scene-first redesign supersedes the default worksheet. See [current scene-first evidence](SCENE-FIRST.md).

Latest interaction/art iteration: [playful contemporary pixel polish](PLAYFUL-PIXEL.md). Reference libraries are inspiration, not constraints.
