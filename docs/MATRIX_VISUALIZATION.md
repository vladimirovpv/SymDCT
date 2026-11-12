# Viewing Transition Matrices

The `matrix_builder` stage writes the transition-rate matrices in several formats.
This page lists the ways to inspect them.

A run writes its files into `examples/<name>/output/`, which is not tracked by git.
Reference results of the examples are stored in `examples/<name>/assets/`.
All commands are run from the root of the repository.
Stage 2 reads the jump graph written by Stage 1, so Stage 1 is run first.

## 1. **Quick Terminal Display (Recommended)**

Run Stage 2 and then the provided script, which prints the matrices as formatted tables:

```bash
uv run python -m src.jump_network_builder examples/hcp_vac/input/hcp_vac.yaml
uv run python -m src.matrix_builder examples/hcp_vac/input/hcp_vac.yaml
uv run python scripts/print_transition_matrices.py examples/hcp_vac/input/hcp_vac.yaml
```

This displays:
- Aligned table format with row/column labels
- Symbolic expressions (with W_AA, W_AB, cos(t), etc.)
- One matrix per diffusion axis
- Easy to read in terminal

**Example output:**
```
Diffusion Axis: (0, 0, 1)
Jump length factor: 1/2
Sites: ['A0', 'B0']

Transition Matrix for axis (0, 0, 1)
───────────────────────────────────
    |          A0 |          B0 |
───────────────────────────────────
 A0 |       -6*W3 | 6*W3*cos(t) |
 B0 | 6*W3*cos(t) |       -6*W3 |
───────────────────────────────────
```

## 2. **Mathematica Format (.m files)**

Stage 2 writes one file per axis, `tp_matrix_axis_<h>_<k>_<l>.m`.
The files appear in `examples/hcp_vac/output/` after running the stage on its own:

```bash
uv run python -m src.jump_network_builder examples/hcp_vac/input/hcp_vac.yaml
uv run python -m src.matrix_builder examples/hcp_vac/input/hcp_vac.yaml
```

Reference copies are in `examples/hcp_vac/assets/`.

```mathematica
(* View in Mathematica: *)
<< tp_matrix_axis_0_0_1.m
(* Then use: *)
MatrixForm[%]
```

**Example content:**
```mathematica
{{-6*W3, 6*W3*Cos[t]},
{6*W3*Cos[t], -6*W3}}
```

Stage 2 also writes `tc_axis_meta_<h>_<k>_<l>.m` next to each matrix file.
It holds the squared interplanar distance and the minimal jump length of the axis.

### Solving large matrices with the Wolfram Language

SymPy becomes slow for large reduced matrices. The script `scripts/solve_dc.wls`
evaluates the same expression as Stage 3 with the Wolfram Language
(Mathematica or the Wolfram Engine):

```bash
uv run python -m src.jump_network_builder examples/hcp_vac/input/hcp_vac.yaml
uv run python -m src.matrix_builder examples/hcp_vac/input/hcp_vac.yaml
wolframscript -file scripts/solve_dc.wls examples/hcp_vac/output/tp_matrix_axis_0_0_1.m
```

The result is written to `dc_symbolic_<h>_<k>_<l>_wl.m` in the same folder.
It is expressed in the jump rates `W0`, `W1`, ... . The meta file lists the
jump type of each rate under `"jumpLabels"`.

## 3. **Pickle Payload (.pickle files)**

The same Stage-2 command writes `tc_axis_payload_<h>_<k>_<l>.pickle`, one file
per axis. These files pass the reduced matrix from Stage 2 to Stage 3 when the
stages are run separately. `main.py` keeps the data in memory and writes no
pickle files. Pickle files are not part of the reference results.

For programmatic access in Python:

```python
import pickle
from pathlib import Path

pkl_file = Path("examples/hcp_vac/output/tc_axis_payload_0_0_1.pickle")
with open(pkl_file, "rb") as f:
    data = pickle.load(f)

print(data["transition_matrix"])              # sympy.Matrix, reduced matrix
print(data["labels"])                         # site-class labels, e.g. ['A0', 'B0']
print(data["axis_hkl_conventional"])          # axis, e.g. (0, 0, 1)
print(data["interplanar_distance_squared"])   # e.g. c**2/4
print(data["min_jump_length"])                # jump length factor
print(data["t_symbol"])                       # symbol of the phase variable
print(data["w_substitutions"])                # merged jump-rate symbols
print(data["jump_label_map"])                 # rate index -> jump type, e.g. {0: 'A->A'}
```

Only load pickle files that you created yourself.

## 4. **Symbolic diffusion coefficients (.txt and .m files)**

Stage 3 writes `dc_symbolic_<h>_<k>_<l>.txt` and `dc_symbolic_<h>_<k>_<l>.m`:

```bash
cat examples/hcp_vac/assets/dc_symbolic_0_0_1.txt
```

## 5. **Jump-network graph visualisation**

Render the full jump graph and the axis-specific quotient graphs from the
Stage-1 and Stage-2 files of a run:

```bash
uv run python scripts/visualize_jump_graphs.py examples/hcp_vac/input/hcp_vac.yaml
```

The images are written to:

```text
examples/hcp_vac/output/graph_visualization/
  full_jump_network.png
  quotient_axis_<h>_<k>_<l>.png
```

Reference images are in `examples/hcp_vac/assets/graph_visualization/`.

## Output File Naming

Matrices are saved with axis notation:
- `tp_matrix_axis_0_0_1.m` → axis (0, 0, 1)
- `tp_matrix_axis_2_-1_0.m` → axis (2, -1, 0)
- `tp_matrix_axis_-1_2_0.m` → axis (-1, 2, 0)

## Understanding the Symbols

- **W_AA** = Jump probability for A→A transitions (W0 in output)
- **W_AB** = Jump probability for A→B transitions (W1, W2, ... in output)
- **t** = Wave vector component (varies with axis)
- **cos(t), sin(t), exp(I*t), etc.** = Wave vector-dependent terms
- **I** = Imaginary unit

## Matrix Structure

In the hcp vacancy example with two site classes each transition matrix is 2×2 (one row/column per site class):

```
       [To A0]    [To B1]
From A0 [diag]    [off-diag]
From B1 [off-diag] [diag]
```

The diagonal elements represent net outflow rates (negative), while off-diagonal elements represent jump rates between sites.

## Interpreting Results with tweak_transition_rates

When using `tweak_transition_rates` (e.g., mapping "A->A: all" to "B->B"), the jump probabilities are redefined before matrix construction. The output matrices will reflect these customizations.

From the hcp_vac example:
```yaml
tweak_transition_rates:
  "A->A: all": "B->B"      # All A→A jumps behave like B→B jumps
  "A->B: 2.213": "B→A"     # A→B jumps behave like B→A jumps
```

This causes the matrix to have different W_AA and W_AB symbols in the output.

## Next Steps

Once matrices are generated, they can be:
1. **Analyzed symbolically** in Mathematica
2. **Evaluated numerically** with specific temperature/activation energy values
3. **Used for diffusion coefficient calculations** (Stage 3)
4. **Plotted** as functions of wave vector in Python using matplotlib
