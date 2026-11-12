# Example walkthroughs

This page summarizes the bundled user-facing examples and the main artifacts
you should expect after each stage.

Each example folder has an `input/` folder with the files needed for a run.
A run writes its artifacts into `output/`, which is not tracked by git.
Reference results of the examples are stored in `assets/`. Compare the files in
`output/` with them to check a run.

## Common artifact progression

For any example, the normal artifact flow is:

1. **Stage 1** → `jump_structure.node_link.json`, `jump_sites.poscar`, `energetics_template.yaml`
2. **Stage 2** → `tp_matrix_axis_<h>_<k>_<l>.m`, `tc_axis_payload_<h>_<k>_<l>.pickle`
3. **Stage 3** → `dc_symbolic_<h>_<k>_<l>.txt`, `dc_symbolic_<h>_<k>_<l>.m`
4. **Stage 4** → `arrhenius_plot_<example>_<h>_<k>_<l>_freq<model>.png`

## `hcp_vac`

Scenario: vacancy-like diffusion example on an hcp structure with multiple
transport axes and `tweak_transition_rates` enabled.

Main input:

- `examples/hcp_vac/input/hcp_vac.yaml`

Notable settings:

- three axes: `[2, -1, 0]`, `[-1, 2, 0]`, `[0, 0, 1]`
- `if_calc_diffcoeff: True`
- label-based `tweak_transition_rates`

Expected artifacts:

- **Stage 1:** jump graph, jump-site POSCAR, energetics template
- **Stage 2:** one matrix export and one payload pickle for each of the three axes
- **Stage 3:** one symbolic diffusion-coefficient text/mathematica pair per axis
- **Stage 4:** one Arrhenius plot per symbolic diffusion-coefficient file, once
  `input/energetics.yaml` is present and synchronized

## `bcc_int`

Scenario: interstitial diffusion example on a bcc-like structure with two
transport axes.

Main input:

- `examples/bcc_int/input/bcc_int.yaml`

Notable settings:

- matrix species filtered as `["Al"]`
- two axes: `[1, 1, 1]` and `[0, 1, 1]`
- no `tweak_transition_rates`

Expected artifacts:

- **Stage 1:** jump graph, jump-site POSCAR, energetics template
- **Stage 2:** axis-specific matrix exports and payload pickles for the two axes
- **Stage 3:** symbolic diffusion-coefficient files for both axes
- **Stage 4:** Arrhenius plots after providing the matching `energetics.yaml`

## `host-host_diplon`

Scenario: a multi-orientation-state test problem where coincident pseudo-species
represent a diplon defect with rotation and translation channels.

Main inputs:

- `examples/host-host_diplon/input/host-host_diplon.yaml`
- `examples/host-host_diplon/input/cutoffs.yaml`
- `examples/host-host_diplon/input/energetics.yaml`

Notable settings:

- orientation pseudo-species `Qx`, `Qy`, `Qz`
- small rotation cutoff plus ordinary translation cutoff
- one transport axis `[1, 1, 1]`

Expected artifacts:

- **Stage 1:** jump graph including zero-distance rotation neighbours and
  finite-distance translation neighbours, plus energetics template
- **Stage 2:** one reduced matrix/payload pair for the `[1, 1, 1]` axis
- **Stage 3:** one symbolic diffusion-coefficient expression showing the
  translation-controlled result
- **Stage 4:** one Arrhenius plot with the example's hand-filled energetics

For the physical interpretation of this example, see
`examples/host-host_diplon/README.md`.

## `hcp_int_OT`

Scenario: interstitial impurity on octahedral (O) and tetrahedral (T) sites of
an hcp lattice, Ishioka and Koiwa (1985), Sec. 4.1, with four jump types
`W_OO`, `W_TT`, `W_TO`, `W_OT` set by `cutoffs.yaml`.

Main input:

- `examples/hcp_int_OT/input/hcp_int_OT.yaml`

Notable settings:

- two axes: `[2, -1, 0]` (2x2 reduced matrix) and `[0, 0, 1]` (3x3 reduced
  matrix; the 6_3 screw axis along c halves the period)
- the result reproduces eqs. (29) and (31) of the paper
  (`src/matrix_builder/tests/test_hcp_interstitial_ot.py`)
