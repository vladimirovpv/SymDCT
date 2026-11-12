# Jump-network core

`src.jump_network` contains the lower-level graph and symbolic-matrix building
blocks used mainly by Stage 2.

This package is not the user-facing Stage-1 CLI. Instead, it provides the core
objects that interpret jump graphs and turn them into symbolic transition
matrices.

## Main components

- `core.py` — graph-oriented jump-network container logic
- `neighbors.py` — neighbour discovery helper used by the geometry workflow
- `symbolic_context.py` — shared symbolic variables and conventions
- `transition_matrix.py` — symbolic transition-matrix assembly

## When users should care

Most users do not run this package directly. It becomes relevant when:

- debugging Stage-2 symbolic matrix construction
- tracing how graph edges become matrix entries
- understanding neighbour-search assumptions used during network assembly

## Related documentation

- Stage-1 CLI: `src/jump_network_builder/README.md`
- Stage-2 CLI: `src/matrix_builder/README.md`
- Pipeline contracts: `docs/contracts.md`
