"""hcp vacancy: merging of the two basis sites along different directions.

The two sites are related by symmetry.  A symmetry operation that relates them
and keeps the transport direction fixed exists along <2-1-10> and <0001>, and
for every direction in the planes spanned by c and an a axis.  It does not exist
along <10-10>, where a merge of the two sites would not be lumpable
(Im A12 = O(k^3)), so the sites stay separate.
In the hkl notation of the input file these directions are (2, -1, 0),
(0, 0, 1) and (1, 0, 0).
"""

from pathlib import Path

import pytest
import sympy as sp

from src.common_utils.utils import load_cutoffs_config
from src.jump_network_builder.geometry import build_jump_network_from_poscar
from src.matrix_builder.matrix_builder import build_matrices
from src.symmetry_context.jump_symmetry_context import JumpNetworkSymmetryContext
from src.transport_coefficient_solver.transport_coefficient_solver import solve_from_matrix_result

EXAMPLE = Path(__file__).resolve().parents[3] / "examples" / "hcp_vac" / "input"
A_AXIS, C_AXIS, PRISM_AXIS = (2, -1, 0), (0, 0, 1), (1, 0, 0)


def _run(tmp_path, **options):
    from ruamel.yaml import YAML

    yaml_file = EXAMPLE / "hcp_vac.yaml"
    params = YAML(typ="safe").load(yaml_file.read_text(encoding="utf-8"))
    # One species for both sites, so that only the symmetry decides the classes.
    poscar = tmp_path / "hcp_one_species.poscar"
    poscar.write_text(
        (EXAMPLE / "hcp_Be.poscar").read_text(encoding="utf-8").replace("A B\n1 1", "Be\n2"),
        encoding="utf-8",
    )
    params["fname"] = str(poscar)
    params["axes"] = [list(A_AXIS), list(C_AXIS), list(PRISM_AXIS)]
    params.pop("tweak_transition_rates", None)
    params.update(options)
    cutoffs, max_len = load_cutoffs_config(yaml_file, params)
    graph, _, _ = build_jump_network_from_poscar(
        params["fname"], cutoffs, str(tmp_path), params["matrix"], max_len
    )
    ctx = JumpNetworkSymmetryContext(params["fname"], params["matrix"])
    result = build_matrices(graph=graph, symmetry_ctx=ctx, params=params)
    return result, solve_from_matrix_result(result)


def _size(result, axis):
    return result.axis_results[axis].reordered_tp_matrix.shape[0]


@pytest.fixture(scope="module")
def default_run(tmp_path_factory):
    return _run(tmp_path_factory.mktemp("hcp_vac_default"))


def test_default_merges_along_a_and_c(default_run):
    result, _ = default_run
    assert _size(result, A_AXIS) == 1
    assert _size(result, C_AXIS) == 1


def test_default_keeps_sites_separate_along_prismatic_direction(default_run):
    result, _ = default_run
    assert _size(result, PRISM_AXIS) == 2


def test_diffusion_coefficients(default_run):
    _, coeffs = default_run
    w = next(s for s in coeffs[A_AXIS].free_symbols if s.name.startswith("W_"))
    a = next(s for s in coeffs[A_AXIS].free_symbols if s.name == "a")
    c = next(s for s in coeffs[C_AXIS].free_symbols if s.name == "c")
    # W_a = W_b = W: D_11 = 2 (3 W_b + W_a) (a/2)^2 = 2 W a^2
    assert sp.simplify(coeffs[A_AXIS] - 2 * w * a**2) == 0
    assert sp.simplify(coeffs[C_AXIS] - sp.Rational(3, 4) * w * c**2) == 0
    # in-plane isotropy of the rank-two tensor
    assert sp.simplify(coeffs[PRISM_AXIS] - coeffs[A_AXIS]) == 0


def test_plane_orbit_merges_only_within_a_plane_with_same_D(default_run, tmp_path):
    # Along a the two sites lie in one plane.  Along c and <10-10> they lie in
    # different planes, which the conservative relation does not merge.
    result, coeffs = _run(tmp_path, site_equivalence="plane_orbit")
    assert _size(result, A_AXIS) == 1
    assert _size(result, C_AXIS) == 2
    assert _size(result, PRISM_AXIS) == 2
    _, ref = default_run
    for axis in coeffs:
        assert sp.simplify(coeffs[axis] - ref[axis]) == 0
