# TOOLS.md - MathPhysics Agent

## Project-Specific Notes

**Repository**
- Root: `workspace2/repos/mathphysics`
- Main document: `amplitudes/99-master.tex`
- Compile command: `cd amplitudes && pdflatex 99-master.tex && pdflatex 99-master.tex`

**Cadabra2**
- Docker image: `cadabra2-ubuntu:24.04`
- Example scripts live in Monoclaw at `Python/Cadabra2/Srednicki/` and `Python/Cadabra2/spinors/`
- Use for symbolic spinor and amplitude calculations

**LaTeX Workflow**
- Always add bibliography block when using `\cite{}`
- Priority order: `.tex` > `.mmd` > plain text

**Memory**
- Long-term: `../../MEMORY.md` (main session only)
- Daily: `memory/YYYY-MM-DD.md`

**Git Rules**
- Never push directly to master
- Use `feat/`, `fix/`, `chore/` branches
- Ernest performs final merge to master

Add any new physics-specific tools or conventions here.
