# TASK-008 — Public research note and product-facing README

Status: TODO
Priority: P1
Depends on: TASK-002, TASK-003
Estimated size: 2–3 hours
Branch suggestion: `feat/task-008-research-note`

## Goal

Make the verified reproduction understandable to a physicist or scientific-AI researcher in a 90-second skim. This is a research note, not a novelty paper.

## Deliverables

1. Refresh the repo-root `README.md` with:
   - one-sentence product/research thesis;
   - exact quick-start command;
   - architecture diagram in text or Mermaid;
   - current verified scope and limitations;
   - links to generated result reports.
2. Add a short section to the existing amplitudes LaTeX notes covering:
   - the two independent computation paths;
   - the tested multiplicity range;
   - truth label: **reproduction**;
   - how this becomes a discovery loop later.
3. Add a “What is not yet built” section: no unattended discovery, no originality claim, no external users unless evidence is added.

## Definition of done

A technical reader can identify what was computed, how it was verified, how to rerun it, and what remains hypothetical without opening more than three files.

## Out of scope

Full paper restructuring, new derivations, generated PDFs, marketing superlatives, or claims about replacing researchers.
