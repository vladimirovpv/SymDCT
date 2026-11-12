# Matrix builder (Stage 2)

`src.matrix_builder` turns the Stage-1 jump graph into the full symbolic transition matrix and the axis-reduced matrices used by the solver.

For the full workflow, see `README.md`. For the input/output contracts, see `docs/contracts.md`.

## Inputs

- main YAML (especially `axes`)
- `jump_structure.node_link.json` from Stage 1
- `jump_sites.poscar` from Stage 1

The stage expects these artifacts in the matching example `output/` directory.

## Typical CLI

Stage 1 must have been run for the same input file, because Stage 2 reads its jump graph.

```bash
uv run python -m src.jump_network_builder examples/hcp_vac/input/hcp_vac.yaml
uv run python -m src.matrix_builder examples/hcp_vac/input/hcp_vac.yaml
```

## Outputs

Stage 2 writes into `examples/<example>/output/`:

- `tp_matrix_axis_<h>_<k>_<l>.m`
- `tc_axis_payload_<h>_<k>_<l>.pickle`

It also prints the full transition-matrix shape and a per-axis summary.

## Notes

- This stage must run after Stage 1 for the same example.
- The main orchestration helpers are `build_matrices(...)` and `export_axis_matrices_to_mathematica(...)`.

## Warnings and error messages

### Wrong Stage-1 graph artifact

```text
Stage 2 requires jump_structure.node_link.json. jump_structure.adjlist does not preserve the node/edge metadata required for matrix construction.
```

- Meaning: Stage 2 was pointed at the adjacency-list export instead of the canonical node-link JSON.
- Recommended action: use `jump_structure.node_link.json`.

### Missing Stage-1 outputs

```text
FileNotFoundError: Graph file not found: <path>
FileNotFoundError: jump_sites.poscar ...
```

- Meaning: Stage 1 has not been run yet, or its outputs are missing.
- Recommended action: rerun Stage 1 for the same example.

### Invalid axis definition

Typical messages:

```text
Axis must be a 3-component [h, k, l] vector, got: ...
Axis component at index <i> must be integer-valued, got: ...
Axis definitions must be plain [h, k, l] vectors; the screw parameter format [[h, k, l], n] is not supported.
```

- Meaning: the `axes:` entry in the main YAML is malformed.
- Recommended action: use plain Miller vectors such as `[1, 0, 0]`.

### Incomplete serialized edge metadata

```text
Serialized graph edge (...) misses attributes: [...]
```

- Meaning: the Stage-1 graph artifact is incomplete or incompatible.
- Recommended action: regenerate the Stage-1 outputs with the current code.

## Merging equivalent sites: `site_equivalence` and `merge_classes`

For every transport axis, sites are merged into equivalence classes (Ishioka
and Koiwa, Philos. Mag. A 52, 267 (1985), Sec. 3). The classes are the orbits
of the subgroup of space-group operations `x' = W x + w` (fractional
coordinates) with `hkl W = hkl`. This condition keeps the transport direction.
Inversion and other operations that reverse it are excluded. The xi component
of every jump is then unchanged, and the Fourier phases of symmetry-equivalent
jumps map onto each other. Two cases occur:

- `hkl . w` is an integer: the operation keeps the site in its plane (rotations
  about xi, mirror planes containing xi, glide planes with the glide vector in
  the plane). These are the "equivalent sublattices" of the paper.
- `hkl . w` is not an integer: a screw axis or a glide plane with a component
  along xi. It carries a site into another plane and shortens the period of the
  reduced one-dimensional lattice to `a/n` (case (b) of the paper, e.g. the
  6_3 screw axis along c in hcp).

Classes are numbered by increasing plane coordinate `s = (hkl . frac) mod 1`
(smallest value among the members). The labels of the reduced matrix are the
node label followed by this number, for example `O0`, `O1`.

Merging is exact if the partition into classes is *lumpable*
(equitable): for every site `j` of a class `d`, the sum of `A_ij(k)` over the
sites `i` of any class `c` is the same for all `j` in `d` (`R A(k) = A'(k) R`
along `k = hkl * t`, see `lumpability.py`). Classes defined by the operations
above satisfy this condition if symmetry-equivalent jumps have equal rates.
The program tests the condition symbolically for every axis along which sites
are merged. A passed test is reported in one line. A failed test prints a
warning and the relations between the rates that would make the merge exact.

Example (hcp vacancy, two basis sites A and B, one species for both; see
`src/matrix_builder/tests/test_hcp_vacancy.py`): they are merged along
`[2, -1, 0]` and `[0, 0, 1]`, and for every direction in a plane spanned by c
and an a axis, because a glide plane or a screw axis that keeps the direction
relates them. They stay separate along `[1, 0, 0]` (the `<10-10>` direction),
where no such operation exists. The merged matrix element would be complex
(`Im A12 = O(k^3)`) and the merge would not be lumpable. D is the same along
`[2, -1, 0]` and `[1, 0, 0]`.

Two optional YAML keys control the merge:

| Key | Values | Meaning |
|-----|--------|---------|
| `site_equivalence` | `axis_orbit` (default), `plane_orbit`, `orbit`, `wyckoff` | `axis_orbit`: sites related by operations with `hkl W = hkl`, including screw axes and glides along xi (a class may span several planes). `plane_orbit`: the same, but only operations that keep each site in its plane, a conservative choice with larger matrices (hcp O/T along c: 6x6 instead of 3x3, same D). `orbit`: classical spglib orbit within one plane; it is weaker, because two sites of one orbit can be related only by an operation that reverses the direction (e.g. inversion about a point in the plane). `wyckoff`: the bare Wyckoff letter within one plane, which also merges sites with different labels at the same position (e.g. orientation pseudo-species). |
| `merge_classes` | `check` (default), `force` | `check`: if a proposed merge is not lumpable, classes are split to the coarsest lumpable partition and the rate relations that would allow the merge are printed. `force`: classes are merged anyway; the result is valid only if the printed relations hold (impose them in the final formula or via `tweak_transition_rates`). |

Example (host-host diplon, three orientation pseudo-species at one site):
`site_equivalence: wyckoff` together with `merge_classes: force` reproduces the
single merged class; the program prints `W3 = W0, W6 = W0` as the condition
for the merge to be exact.
