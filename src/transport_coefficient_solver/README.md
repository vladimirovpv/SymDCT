# Transport coefficient solver (Stage 3)

`src.transport_coefficient_solver` solves the symbolic diffusion-coefficient expressions from the payloads exported by Stage 2.

For the full workflow, see `README.md`. For the Stage-2 payload contract, see `docs/contracts.md`.

## Inputs

Stage 3 consumes per-axis payloads in the example `output/` directory:

- `tc_axis_payload_<h>_<k>_<l>.pickle`

These are produced by `src.matrix_builder`.

## Typical CLI

```bash
uv run python -m src.transport_coefficient_solver examples/hcp_vac/output
```

You can also pass a YAML path; the solver resolves `<yaml>/../output` automatically.

## Outputs

The stage writes:

- `dc_symbolic_<h>_<k>_<l>.txt`
- `dc_symbolic_<h>_<k>_<l>.m`

## Python API

- `solve_from_output_dir(output_dir, ...)`
- `solve_from_matrix_result(matrix_result, ...)`
- `solve_from_payloads(payloads, ...)`

This path skips disk serialization between Stage 2 and Stage 3 when you are orchestrating directly in Python.

## Warnings and errors

### Missing Stage-2 payloads

```text
Run matrix_builder to generate tc_axis_payload_<h>_<k>_<l>.pickle files.
```

- Meaning: the solver did not find the per-axis payload pickles exported by Stage 2.
- Recommended action: run the matrix builder first for the same example.

### Invalid payload contract

Typical messages include:

```text
Unsupported payload contract type
axis_hkl_structure must be a 3-component sequence
transition_matrix must be a sympy Matrix
labels must be a sequence
```

- Meaning: a payload file exists but does not match the Stage-3 input schema.
- Recommended action: regenerate Stage-2 outputs with the current code rather than mixing artifacts from older or external versions.

### Missing transport-distance metadata

```text
interplanar_distance_squared is required for solver input
```

- Meaning: the payload is incomplete and lacks the axis-specific transport distance required for the diffusion-coefficient prefactor.
- Recommended action: regenerate the payloads via Stage 2.

### Singular symbolic system

```text
Denominator is zero; transition matrix is singular for t=0
```

- Meaning: the reduced symbolic system cannot be evaluated in the required limit.
- Recommended action: inspect the Stage-2 matrix for the affected axis and confirm that the jump network and site-equivalence reduction are physically meaningful.
