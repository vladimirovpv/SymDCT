# `framed_workflow` — manual overlay history

File-specific edit log for `framed_workflow.svg` / `framed_workflow_with_title_arrows.svg`.
General technique lives in `public/docs/diagram/mermaid-to-pdf-skill.md`; this
file is instance history only.

- Added top title text: `SymDCT Pipeline Workflow`.
- Added custom elbow connectors:
  - `YAML -> S1A`
  - `POSCAR -> S1A`
  - `ENERG -> S4A`
- Merged `YAML` and `POSCAR` through an invisible junction node before `S1A`.
- Adjusted arrow thickness/endpoints to match native Mermaid arrows and avoid
  frame overlap.
- Increased top canvas margin via `viewBox` adjustment (`0 -45 W H`) to
  improve title spacing.
- Replaced the yellow cluster background (`#ffffde`, stroke `#aaaa33`) by a muted
  grey-blue (fill `#DCE5EE`, stroke `#7F97B0`) in the SVG. The PDFs were
  recolored in place by replacing the two color operators in the page content
  stream, so layout and fonts of the existing figure are unchanged.
