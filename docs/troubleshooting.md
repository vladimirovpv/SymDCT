# Troubleshooting index

This page collects the most common user-facing warnings and errors across all
stages of the public pipeline.

## Stage 1 - jump-network builder

### `WARNING: <N> neighbors found for bond A-B without a specified cutoff.`

- Meaning: the pair was found within the global neighbour-search radius but
  `cutoffs.yaml` has no explicit cutoff for that pair.
- Fix: add the missing pair to `cutoffs.yaml` if you want strict control.

### `ERROR: No neighbors found within the cutoff distance <max_jump_length>. Try to increase the cutoff.`

- Meaning: the global neighbour search found no candidate jumps.
- Fix: increase `max_jump_length`, check `fname`, and confirm the `matrix`
  species selection did not remove the intended jump sites.

### `FileNotFoundError: POSCAR file '<path>' not found. Check the 'fname' parameter in <yaml>.`

- Meaning: the structure path in the main input YAML could not be resolved.
- Fix: correct `fname` or move the POSCAR.

## Stage 2 - matrix builder

### `Stage 2 requires jump_structure.node_link.json...`

- Meaning: Stage 2 was given `jump_structure.adjlist` instead of the canonical
  node-link JSON.
- Fix: use `jump_structure.node_link.json`.

### `FileNotFoundError: Graph file not found: <path>`

- Meaning: the Stage-1 graph artifact is missing.
- Fix: run Stage 1 first for the same example.

### `Axis must be a 3-component [h, k, l] vector...`

- Meaning: an `axes:` entry in the main YAML is malformed.
- Fix: use plain integer Miller vectors such as `[0, 0, 1]`.

## Stage 3 - symbolic solver

### `Run matrix_builder to generate tc_axis_payload_<h>_<k>_<l>.pickle files.`

- Meaning: Stage 3 cannot find Stage-2 payload pickles.
- Fix: run Stage 2 first.

### `interplanar_distance_squared is required for solver input`

- Meaning: the Stage-2 payload is incomplete or incompatible.
- Fix: regenerate the payloads with the current Stage-2 code.

### `Denominator is zero; transition matrix is singular for t=0`

- Meaning: the symbolic limit needed by the solver failed.
- Fix: inspect the Stage-2 reduced matrix and verify the physical model and
  axis reduction.

## Stage 4 - numeric evaluator

### `dc expression file not found: <path>`

- Meaning: the specified Stage-3 symbolic output file is missing.
- Fix: run Stage 3 first or point to the correct `dc_symbolic_*.txt`.

### `WARNING: Energetics template/production mismatch ...`

- Meaning: `energetics.yaml` no longer matches the current
  `energetics_template.yaml`.
- Fix: synchronize the production energetics file before trusting Stage-4
  results.

### `Automatic synchronization failed: ...`

- Meaning: `--auto-sync` or `DC_AUTO_SYNC=1` was used but the repair step
  failed.
- Fix: run `src.common_utils.synchronizer` manually and resolve placeholders.

## Full pipeline via `main.py`

### Stage 4 skipped due to missing energetics

Typical message:

```text
Production energetics file not found: <path>; skipping Stage 4 (numeric evaluation).
```

- Meaning: the end-to-end pipeline completed Stages 1-3 but found no
  `input/energetics.yaml`.
- Fix: add `energetics.yaml` if numeric evaluation is required.
