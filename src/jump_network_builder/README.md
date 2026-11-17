# Jump-network builder (Stage 1)

`src.jump_network_builder` constructs the geometric jump network from a POSCAR and a YAML config. It is the first stage in the public pipeline and generates the graph and site artifacts used downstream.

For the full workflow and stage map, see `README.md`. For the YAML and artifact contracts, see `docs/contracts.md`.

## Inputs

- main input YAML
- POSCAR referenced by `fname`
- optional `cutoffs.yaml`

If no explicit cutoffs file is provided, the stage auto-discovers a local `input/cutoffs.yaml` next to the main YAML. If that file is absent, it falls back to `max_jump_length`.

## Jump classes

Every directed jump belongs to a jump class, and each class carries one rate symbol `W<n>`. The key `jump_classification` selects the rule.

- `orbit` (default): jumps related by a space-group operation of the jump structure form one class. The operations come from spglib, with the node labels acting as atom types. The reverse of a jump belongs to the class that the symmetry assigns to it, so a polar structure can give a jump and its reverse different classes.
- `triplet`: jumps with equal start label, end label and length (within a tolerance) form one class. This is the rule of version 0.1.0.

Both partitions are computed. If they differ, the log names the affected jumps. If the neighbour selection does not respect the symmetry (for example after a cutoff that keeps only some of the symmetry-equivalent neighbours), a warning reports the number of missing images.

## Typical CLI

```bash
uv run python -m src.jump_network_builder examples/hcp_vac/input/hcp_vac.yaml
```

Optional explicit cutoffs override:

```bash
uv run python -m src.jump_network_builder \
    examples/hcp_vac/input/hcp_vac.yaml \
    examples/hcp_vac/input/cutoffs.yaml
```

## Outputs

Stage 1 writes into `examples/<example>/output/`:

- `jump_structure.node_link.json` — canonical Stage-2 graph input
- `jump_structure.adjlist` — inspection artifact
- `jump_sites.poscar` — refined jump-site structure
- `energetics_template.yaml` — template for Stage 4

## Stage-local notes

1. resolves `fname` and loads the structure
2. optionally removes matrix species from the jump-site set
3. performs a global neighbour search
4. applies pair-specific cutoffs from `cutoffs.yaml`
5. exports graph and template artifacts for the next stages

## Warnings and error messages

### Missing pair-specific cutoff

```text
WARNING: <N> neighbors found for bond A-B without a specified cutoff.
```

- Meaning: the pair was found within the global search radius, but `cutoffs.yaml` has no explicit entry for it.
- Behaviour: those neighbours are still kept.
- Recommended action: add the missing pair key to `cutoffs.yaml` if you want stricter geometric control.

### No neighbours found within the global cutoff

```text
ERROR: No neighbors found within the cutoff distance <max_jump_length>. Try to increase the cutoff.
```

- Meaning: the global neighbour search found no candidate jumps.
- Recommended action: increase `max_jump_length`, check the POSCAR, and verify that `matrix` does not remove the intended jump species.

### YAML parsing failure

- Meaning: the input YAML or `cutoffs.yaml` is invalid.
- Recommended action: fix the YAML syntax and rerun.

### Missing POSCAR referenced by `fname`

```text
FileNotFoundError: Structure file '<path>' not found. Check the 'fname' parameter in <yaml>.
```

- Meaning: the structure path in the main YAML could not be resolved.
- Recommended action: fix `fname` or move the POSCAR.
