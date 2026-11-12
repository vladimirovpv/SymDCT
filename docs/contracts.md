# YAML and artifact contracts

This document describes the main user-facing file contracts for the public
pipeline.

## 1. Main input YAML

Purpose: drive Stage 1 directly and Stages 2-3 indirectly through the same
configuration source.

Typical location:

```text
examples/<example>/input/<example>.yaml
```

Common keys:

| Key | Required | Type | Meaning |
|---|---|---|---|
| `fname` | yes | string | Path to the POSCAR, usually relative to the YAML file |
| `axes` | usually | list of 3-component lists | Miller directions for transport analysis |
| `matrix` | optional | list of strings | Species treated as host/matrix atoms and removed from the jump-site set |
| `max_jump_length` | recommended | float | Global neighbour-search radius used when no better value is derived from `cutoffs.yaml` |
| `tweak_transition_rates` | optional | mapping | Manual remapping/scaling hints for symbolic jump rates |

Minimal example:

```yaml
fname: hcp_vac.poscar
matrix: [Ti]
axes:
  - [0, 0, 1]
  - [1, 0, 0]
max_jump_length: 3.0
```

Common mistakes:

- `fname` points to a missing file
- `axes` contains non-integer or malformed entries
- `matrix` removes the wrong species and leaves no intended jump substructure

## 2. `cutoffs.yaml`

Purpose: define pair-specific geometric cutoffs for Stage 1.

Typical location:

```text
examples/<example>/input/cutoffs.yaml
```

Common keys:

| Key | Required | Type | Meaning |
|---|---|---|---|
| `units.length` | recommended | string | Unit label for all cutoff values |
| `bonds.<TYPE1-TYPE2>` | yes | float | Maximum bond length for that pair |

Minimal example:

```yaml
units:
  length: angstrom

bonds:
  A-A: 1.4
  A-B: 2.2
  B-B: 1.8
```

Semantics:

- bond keys are order-independent: `A-B` and `B-A` mean the same thing
- the largest listed cutoff is used for the initial global neighbour search
- omitted pairs are kept if found, but generate a warning

## 3. `energetics_template.yaml`

Purpose: Stage-1-generated template listing every discovered jump and the jump
index mapping needed later by Stages 3 and 4.

Typical location:

```text
examples/<example>/output/energetics_template.yaml
```

Important sections:

| Key | Required | Meaning |
|---|---|---|
| `sites` | yes | Site labels discovered from geometry |
| `jumps` | yes | Jump families grouped by pair key and representative length |
| `jumps.<PAIR>[].jump_id` | yes | Forward/backward global jump indices |
| `jumps.<PAIR>[].length` | yes | Representative jump length |

Semantics:

- this file is generated, not hand-maintained
- it is the geometry/reference contract for the numeric energetics file
- if Stage 1 is rerun after geometry or cutoff changes, this file may change

## 4. `energetics.yaml`

Purpose: user-authored numeric energetics for Stage 4.

Typical location:

```text
examples/<example>/input/energetics.yaml
```

Important sections:

| Key | Required | Meaning |
|---|---|---|
| `settings.length_tolerance` | recommended | Tolerance for matching jump lengths |
| `settings.units.*` | recommended | Unit metadata |
| `frequencies.default_frequency` | recommended | Fallback attempt frequency |
| `frequencies.<TYPE>` | optional | Per-site-type frequency override |
| `sites.<TYPE>` | yes | Site energy for each node label |
| `jumps.<PAIR>[].length` | yes | Representative jump length |
| `jumps.<PAIR>[].energy_sp` | yes | Saddle-point energy |
| `jumps.<PAIR>[].jump_id` | yes | Forward/backward global jump ids |

Minimal skeleton:

```yaml
settings:
  length_tolerance: 0.05
  units:
    energy: eV
    length: angstrom
    frequency: THz

frequencies:
  default_frequency: 1.0

sites:
  A: -3.2
  B: -3.4

jumps:
  A-B:
    - jump_id:
        forward: 0
        backward: 1
      length: 2.14
      energy_sp: 0.7
```

Common mistakes:

- keeping an outdated copy after `energetics_template.yaml` changed
- missing `energy_sp` values for new jumps
- mismatched lengths or pair keys relative to the template

## 5. Stage artifacts

The core pipeline artifacts are:

| Artifact | Producer | Consumer |
|---|---|---|
| `jump_structure.node_link.json` | Stage 1 | Stage 2 |
| `jump_sites.poscar` | Stage 1 | Stage 2 |
| `energetics_template.yaml` | Stage 1 | Stage 4 sync/mismatch checks |
| `tc_axis_payload_<h>_<k>_<l>.pickle` | Stage 2 | Stage 3 |
| `dc_symbolic_<h>_<k>_<l>.txt` | Stage 3 | Stage 4 |

The `.adjlist` graph export is for inspection only and is not the canonical
Stage-2 input artifact.
