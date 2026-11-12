"""bcc interstitial (T and O sites): reductions along <111>, <100> and <110>.

The primitive cell has nine sites.  The reduced matrices have sizes 2, 4 and 5
along the three directions, and all give the same D (isotropy of a rank-two
tensor in a cubic crystal).  Reference: Ishioka and Koiwa, Phil. Mag. A 52 (1985).
"""

from pathlib import Path

import pytest
import sympy as sp

from src.common_utils.utils import load_cutoffs_config
from src.jump_network_builder.geometry import build_jump_network_from_poscar
from src.matrix_builder.matrix_builder import build_matrices
from src.symmetry_context.jump_symmetry_context import JumpNetworkSymmetryContext
from src.transport_coefficient_solver.transport_coefficient_solver import solve_from_matrix_result

EXAMPLE = Path(__file__).resolve().parents[3] / "examples" / "bcc_int" / "input"
AXES = {(1, 1, 1): 2, (1, 0, 0): 4, (1, 1, 0): 5}


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    from ruamel.yaml import YAML

    tmp_path = tmp_path_factory.mktemp("bcc_int")
    yaml_file = EXAMPLE / "bcc_int.yaml"
    params = YAML(typ="safe").load(yaml_file.read_text(encoding="utf-8"))
    params["fname"] = str(EXAMPLE / "prim_int_bcc.poscar")
    params["axes"] = [list(a) for a in AXES]
    params.pop("tweak_transition_rates", None)
    cutoffs, max_len = load_cutoffs_config(yaml_file, params)
    graph, _, _ = build_jump_network_from_poscar(
        params["fname"], cutoffs, str(tmp_path), params["matrix"], max_len
    )
    ctx = JumpNetworkSymmetryContext(params["fname"], params["matrix"])
    result = build_matrices(graph=graph, symmetry_ctx=ctx, params=params)
    return result, solve_from_matrix_result(result)


@pytest.mark.parametrize("axis,size", list(AXES.items()))
def test_reduced_matrix_size(run, axis, size):
    result, _ = run
    assert result.axis_results[axis].reordered_tp_matrix.shape[0] == size


def test_same_D_along_all_directions(run):
    _, coeffs = run
    ref = coeffs[(1, 1, 1)]
    for axis in AXES:
        assert sp.simplify(coeffs[axis] - ref) == 0
