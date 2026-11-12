# Symmetry context

`src.symmetry_context` contains the symmetry-handling layer shared by Stages 1
and 2.

Its job is to:

- refine structures into a symmetry-aware representation
- label equivalent jump sites
- convert transport axes between conventional and structure-cell bases
- compute axis-dependent site-equivalence classes for matrix reduction

## Main components

- `jump_symmetry_context.py` — façade used by the pipeline stages
- `symmetry_refiner.py` — structure refinement via pymatgen/spglib
- `node_labeler.py` — symmetry-based node labeling
- `axis_equivalence.py` — axis-dependent equivalence analysis
- `smatrix_handler.py` — conventional/primitive basis conversions

## When users should care

Most users do not invoke this package directly. It is most useful when:

- checking why two sites were grouped together
- understanding how a user-specified axis was interpreted
- debugging conventional-cell vs primitive-cell behaviour

## Related documentation

- Stage-1 geometry stage: `src/jump_network_builder/README.md`
- Stage-2 matrix stage: `src/matrix_builder/README.md`
- Example workflows: `examples/README.md`
