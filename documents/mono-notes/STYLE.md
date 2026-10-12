# Mathematical exposition in these notes

Ernest's `LaTeXandpdfs/Propulsion.tex`, especially the Runge–Kutta discussion,
is a guide to the exposition. These conventions were made explicit by Ernest
on 2026-09-12 and will be refined with him.

- Preserve his sequence, meaningful equation labels, wording and mathematical
  work when incorporating an existing passage almost verbatim. Add connected
  explanations and the missing calculations. Record necessary mathematical
  corrections with their reasons; do not silently change the meaning.
- Keep textbook exercises/problems and their solutions inline in the relevant
  shared topic, rather than creating a separate solutions manuscript.
- Reproduce exercise/problem statements and supplied hints verbatim, including
  the source notation and punctuation. Repair OCR artifacts and encode the
  mathematics correctly in LaTeX; reserve rewording, notation changes, and
  clarification for the solution. Cite the book, exercise identifier, and
  printed page, and distinguish any editorial correction from the quotation.
- Use the shared unnumbered `exercise` environment (amsthm, definition style),
  with the optional heading in this approved order: book/volume, original
  exercise/problem identifier, printed page, and bibliography citation.
  Use the Hairer example as the default template for future insertions:
  `\begin{exercise}[Hairer I, II.1.1, p.~141 \cite{HNW1993}]`.
  Follow it with `\paragraph{Solution.}` for an unfinished working solution.
  A completed proof may use `\begin{proof}[Solution]` and `\end{proof}`.
- Treat notation translation as part of the mathematics. State each author's
  symbols, their roles, and the explicit map to our symbols. Keep the original
  notation long enough to make the comparison checkable. Name the point at
  which our notation takes over.
- Write the mass ratio explicitly as `\frac{M_f}{M(0)}` in our
  notation, or `\frac{m_f}{m_0}` when using Sutton's mass symbols.
  Use `\frac{M(0)}{M_f}` for the reciprocal. Ernest prefers the
  explicit fractions and rejects both `MR` and `r_m` as shorthand
  (2026-10-07). Show the actual masses for stages and subsystem boxes;
  keep `M(t)`, `M_f`, `m_p(t)` and
  `\zeta = 1-\frac{M_f}{M(0)}`. Define each source author's convention
  when translating it. Simplify reciprocal fractions and split long
  displays at meaningful equalities so both PDF editions stay readable.
- Prefer explicit partial-derivative fractions, including in prose:
  `\frac{\partial f^r}{\partial t}` and
  `\frac{\partial^2 f^r}{\partial t^2}`. Do not routinely replace them
  with derivative subscripts such as `f_t^r`, `f_{tt}^r`, or operators
  such as `\partial_t`. Reserve derivative shorthand for sustained algebra
  extending over several pages, when it materially improves readability.
  Even then, define each form explicitly **before** first use with `\frac`,
  and repeat the definitions at new derivations or after a long gap.
  A distant definition or equation reference alone is not a sufficient reminder.
- Keep component reminders to one short sentence, for example:
  `$f^r$ is component $r$ of $\mathbf f$.` Repeat only what the nearby
  calculation needs; do not reintroduce powers, covariance, metric raising,
  and stage indices in a paragraph each time.
- Assume the reader knows Einstein summation convention. Use it without
  explanations or demonstrations expanding contractions into explicit sums.
  State an index range when the coordinate space changes.
- In vector-valued derivative expressions, show contractions with scalar
  coordinate/component factors:
  `\frac{\partial\mathbf f}{\partial x^j}f^j`, rather than
  `\frac{\partial\mathbf f}{\partial\mathbf x}\mathbf f`.
  This keeps the vector-valued output in bold and makes the contraction
  explicit. Apply the same convention to other vectors and higher derivatives.
- State the variables held fixed and the evaluation point of partial derivatives.
- Choose readable shorthand after defining it: each abbreviation must expand
  to a specified expression. Distinguish collisions such as stage abscissa
  versus output weight, derivative versus increment, stage index versus state
  index, and exact solution versus numerical state.
- Organize shared physics topics around the physical question, definitions
  and derivation. Combine complementary textbook presentations; place their
  equation numbers and notation translations as parallel source references.
  Introduce the principal variables on separate lines with definitions and
  units when setting up a physical model.
- State the equation, spaces, domains and assumptions before deriving a
  result. Identify additional assumptions at the precise specialization:
  autonomous/nonautonomous, constant matrix/nonlinear field, exact starting
  data/computed state, fixed interval/long-time claim.
- Show the algebra explicitly: apply the product and chain rules term by term,
  show intermediate substitutions, rename dummy indices visibly, and identify
  the symmetry or cancellation used. Explain numerical factors such as 2, 3,
  and factorial denominators rather than merely asserting them.
- For a physical balance, first display the equation in named rates or
  inventories. Then display each component formula with its origin, and
  substitute one component at a time in an explicit chain of equations.
  Put an author's equation beside our translated equation when comparing
  notation. Prose should explain the displayed mathematics, not replace it.
- Keep time-dependent inputs explicit in transient models. State exactly
  where an input is held fixed for an equilibrium or stability calculation,
  and distinguish a moving quasi-steady target from the actual solution.
  Justify neglecting a term by its magnitude relative to retained terms;
  a small state variable alone does not establish a small time derivative.
- Keep commentary such as “add and subtract,” “iterate,” and “sum the
  geometric series,” and display the corresponding work. Show the inserted
  zero before regrouping and applying an inequality; write the first few
  recurrence substitutions, the induction step where needed, and the terms
  that cancel in a sum. Handle exceptional cases before dividing by a
  potentially zero quantity.
- Define a proof's maps by explicit formulas and identify their arguments.
  State each assumed estimate separately, explain the role of its constants,
  and distinguish hypotheses from derived conclusions. Keep dependencies
  visible when uniformity matters: for example, what is fixed as h varies.
  Brief notation reminders do not justify omitting the mathematical work.
- Give reusable equations stable labels and cite those labels where they are
  actually substituted. Identify external book volume, section, printed page,
  and equation number where available. Every citation needs a real bibliography
  entry. A second printing is not a new edition.
- Complete a TODO with the missing argument. A checked coefficient table is
  useful evidence after a derivation; it does not replace the derivation.
- Keep the prose to the point. Explain why the displayed step is being taken
  and how it follows from the previous one. Avoid rhetorical transitions,
  vague assurances, or extra generalizations without hypotheses.
- Separate general proofs from symbolic examples, finite numerical checks,
  and physical validation, at the point where each appears. Do not add
  closing "Scope of what is claimed" sections, transcription-change logs or
  other meta-commentary about the notes themselves; state assumptions where
  they are used and keep the text on the physics and mathematics
  (Ernest, 2026-10-06).
- Preserve readable mathematics in both portrait and screen PDFs. Split long
  displays at meaningful equalities or terms, and add equation references
  where page breaks separate an argument. Do not hide overflow by shrinking
  an entire derivation into unreadable text.


These conventions are inherited from Propulsion/documents/notes/STYLE.md.
For new work here, the shared topic is the source for both PDF editions.
Historical imported notes retain their existing exposition; these conventions
guide additions and deliberate revisions, rather than silent rewriting.
