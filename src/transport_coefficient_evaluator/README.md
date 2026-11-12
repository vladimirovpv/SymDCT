# Transport coefficient evaluator (Stage 4)

`src.transport_coefficient_evaluator` performs the numeric evaluation step of the public diffusion pipeline: it reads one symbolic diffusion-coefficient expression from Stage 3, substitutes energies and attempt frequencies from `energetics.yaml`, and writes an Arrhenius plot.

For the full workflow and stage map, see `README.md`. For energetics synchronization guidance, see `src/common_utils/README.md`.

## Inputs

Stage 4 needs:

- one symbolic expression file from Stage 3: `dc_symbolic_<h>_<k>_<l>.txt`
- one production energetics file: `input/energetics.yaml`
- optionally the generated template: `output/energetics_template.yaml`

If the template path is omitted, the evaluator looks for it automatically in `../output/energetics_template.yaml` next to the production file.

## Typical CLI

```bash
uv run python -m src.transport_coefficient_evaluator.evaluator \
    examples/hcp_vac/output/dc_symbolic_0_0_1.txt \
    examples/hcp_vac/input/energetics.yaml
```

You may also pass the template path explicitly:

```bash
uv run python -m src.transport_coefficient_evaluator.evaluator \
    examples/hcp_vac/output/dc_symbolic_0_0_1.txt \
    examples/hcp_vac/input/energetics.yaml \
    examples/hcp_vac/output/energetics_template.yaml
```

Useful options:

- `--length-tolerance <float>` — override the jump-length matching tolerance
- `--auto-sync` — attempt automatic template/production synchronization before evaluation

## Outputs

For each symbolic expression file evaluated, Stage 4 writes:

- `arrhenius_plot_<example>_<h>_<k>_<l>_freq<model>.png`

into the same `output/` directory as the input `dc_symbolic_*.txt` file.

## Warnings and error messages

### Missing symbolic diffusion-coefficient file

```text
dc expression file not found: <path>
```

- Meaning: the requested Stage-3 output file does not exist.
- Recommended action: run Stage 3 first or point the evaluator at the correct `dc_symbolic_<h>_<k>_<l>.txt` file.

### Template/production mismatch

Typical message:

```text
WARNING: Energetics template/production mismatch ...
```

- Meaning: the Stage-1-generated `energetics_template.yaml` no longer matches the hand-edited `energetics.yaml`.
- Typical causes: changed cutoffs, changed jump geometry, regenerated Stage-1 outputs, or manual edits that removed or added jump entries.
- Recommended action: inspect the common utilities documentation and re-synchronize the production energetics file before trusting numeric results.

### Automatic synchronization failed

```text
Automatic synchronization failed: ...
Please run the synchronizer manually and fix placeholders. See src/common_utils/README.md for details.
```

- Meaning: `--auto-sync` or `DC_AUTO_SYNC=1` was used, but the repair step failed.
- Recommended action: run the synchronizer manually, inspect the generated file, and fill any placeholder `energy_sp` values before retrying.
