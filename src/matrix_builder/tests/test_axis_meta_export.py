"""Stage 2 writes the scale factors needed by scripts/solve_dc.wls.

The file tc_axis_meta_<h>_<k>_<l>.m must reproduce, together with the matrix
file tp_matrix_axis_<h>_<k>_<l>.m, the diffusion coefficient of Stage 3.
"""

import re
from pathlib import Path

import sympy as sp
from ruamel.yaml import YAML
from sympy.parsing.mathematica import parse_mathematica

from src.common_utils.utils import load_cutoffs_config
from src.jump_network_builder.geometry import build_jump_network_from_poscar
from src.matrix_builder.matrix_builder import build_matrices, export_axis_matrices_to_mathematica
from src.symmetry_context.jump_symmetry_context import JumpNetworkSymmetryContext

EXAMPLE = Path(__file__).resolve().parents[3] / "examples" / "hcp_int_OT" / "input"


def _diffusion_coefficient_from_files(folder: Path, axis_label: str) -> sp.Expr:
    """Evaluate the expression of scripts/solve_dc.wls with SymPy."""
    matrix = sp.Matrix(parse_mathematica((folder / f"tp_matrix_axis_{axis_label}.m").read_text()))
    meta = (folder / f"tc_axis_meta_{axis_label}.m").read_text()
    d2 = parse_mathematica(re.search(r'"interplanarDistanceSquared" -> (.*),', meta).group(1))
    jump_length = parse_mathematica(re.search(r'"minJumpLength" -> (.*),', meta).group(1))
    t = next(s for s in matrix.free_symbols if s.name == "t")
    numerator = sp.diff(matrix.det(), t, 2).subs(t, 0)
    b = matrix.copy()
    b[0, :] = sp.ones(1, matrix.shape[0])
    expr = -sp.Rational(1, 2) * d2 * jump_length**2 * numerator / b.subs(t, 0).det()
    # compare by name: the parser returns symbols without assumptions
    return expr.subs({s: sp.Symbol(s.name) for s in expr.free_symbols})


def test_meta_file_reproduces_hcp_interstitial_coefficients(tmp_path):
    yaml_file = EXAMPLE / "hcp_int_OT.yaml"
    params = YAML(typ="safe").load(yaml_file.read_text(encoding="utf-8"))
    params["fname"] = str(EXAMPLE / "hcp_OT.poscar")
    cutoffs, max_len = load_cutoffs_config(yaml_file, params)
    graph, _, _ = build_jump_network_from_poscar(
        params["fname"], cutoffs, str(tmp_path), params["matrix"], max_len
    )
    ctx = JumpNetworkSymmetryContext(params["fname"], params["matrix"])
    result = build_matrices(graph=graph, symmetry_ctx=ctx, params=params)
    export_axis_matrices_to_mathematica(result, params, tmp_path)

    meta = (tmp_path / "tc_axis_meta_0_0_1.m").read_text()
    assert meta.startswith("<|") and meta.rstrip().endswith("|>")
    assert '"axis" -> {0, 0, 1}' in meta

    # W0: T->O, W1: T->T, W2: O->T, W3: O->O
    w0, w1, w2, w3, a, c = sp.symbols("W0 W1 W2 W3 a c")
    d_par = w0 * c**2 * (3 * w3 * w0 + 2 * w3 * w1 + 3 * w2 * w1) / (4 * (2 * w2 + w0) * (3 * w0 + 2 * w1))
    d_perp = w2 * w0 * a**2 / (2 * w2 + w0)
    assert sp.simplify(_diffusion_coefficient_from_files(tmp_path, "0_0_1") - d_par) == 0
    assert sp.simplify(_diffusion_coefficient_from_files(tmp_path, "2_-1_0") - d_perp) == 0
