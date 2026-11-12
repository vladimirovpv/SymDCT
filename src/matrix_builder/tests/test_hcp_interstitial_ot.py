"""hcp interstitial on O and T sites (Ishioka and Koiwa 1985, Sec. 4.1).

The example reproduces their eqs. (29) and (31).  Along c the 2_1 screw axis
reduces the one-dimensional lattice to three sublattices (their Fig. 3).
"""

from pathlib import Path

import pytest
import sympy as sp

from src.common_utils.utils import load_cutoffs_config
from src.jump_network_builder.geometry import build_jump_network_from_poscar
from src.matrix_builder.matrix_builder import build_matrices
from src.symmetry_context.jump_symmetry_context import JumpNetworkSymmetryContext
from src.transport_coefficient_solver.transport_coefficient_solver import solve_from_matrix_result

EXAMPLE = Path(__file__).resolve().parents[3] / "examples" / "hcp_int_OT" / "input"


def _run(mode, tmp_path):
    from ruamel.yaml import YAML

    yaml_file = EXAMPLE / "hcp_int_OT.yaml"
    params = YAML(typ="safe").load(yaml_file.read_text(encoding="utf-8"))
    params["fname"] = str(EXAMPLE / "hcp_OT.poscar")
    params["site_equivalence"] = mode
    params["merge_classes"] = "check"
    cutoffs, max_len = load_cutoffs_config(yaml_file, params)
    graph, _, _ = build_jump_network_from_poscar(
        params["fname"], cutoffs, str(tmp_path), params["matrix"], max_len
    )
    ctx = JumpNetworkSymmetryContext(params["fname"], params["matrix"])
    result = build_matrices(graph=graph, symmetry_ctx=ctx, params=params)
    return result, solve_from_matrix_result(result)


def _symbols(expr):
    return {s.name: s for s in expr.free_symbols}


@pytest.fixture(scope="module")
def axis_orbit_run(tmp_path_factory):
    return _run("axis_orbit", tmp_path_factory.mktemp("hcp_ot_axis"))


def test_axis_orbit_matrix_sizes_match_paper(axis_orbit_run):
    result, _ = axis_orbit_run
    sizes = {ax: r.reordered_tp_matrix.shape[0] for ax, r in result.axis_results.items()}
    assert sizes[(2, -1, 0)] == 2
    assert sizes[(0, 0, 1)] == 3


def test_diffusion_coefficients_match_ishioka_koiwa(axis_orbit_run):
    _, coeffs = axis_orbit_run
    a, c = sp.symbols("a c", positive=True)
    d_perp = coeffs[(2, -1, 0)]
    d_par = coeffs[(0, 0, 1)]
    names = {**_symbols(d_perp), **_symbols(d_par)}
    w_oo, w_tt, w_to, w_ot = (names[f"W_{{{k}}}"] for k in ("OO", "TT", "TO", "OT"))
    a_sym = _symbols(d_perp)["a"]
    c_sym = _symbols(d_par)["c"]
    ref_perp = w_to * w_ot / (w_to + 2 * w_ot) * a_sym**2
    ref_par = (
        w_to * (3 * w_oo * w_to + 2 * w_oo * w_tt + 3 * w_tt * w_ot)
        / (4 * (2 * w_tt + 3 * w_to) * (w_to + 2 * w_ot))
        * c_sym**2
    )
    assert sp.simplify(d_perp - ref_perp) == 0
    assert sp.simplify(d_par - ref_par) == 0


def test_plane_orbit_gives_same_D_with_larger_matrix(axis_orbit_run, tmp_path):
    result, coeffs = _run("plane_orbit", tmp_path)
    assert result.axis_results[(0, 0, 1)].reordered_tp_matrix.shape[0] == 6
    _, ref = axis_orbit_run
    for axis in coeffs:
        assert sp.simplify(coeffs[axis] - ref[axis]) == 0
