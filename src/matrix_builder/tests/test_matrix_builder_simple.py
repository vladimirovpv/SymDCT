"""Smoke tests for the simplified matrix_builder helpers.

These tests do *not* attempt to validate full physical correctness of
the transition matrices.  Instead they ensure that the new, factored
helpers

* :func:`build_full_transition_matrix`, and
* :func:`find_equiv_and_reduce_matrix_simple`

are wired correctly and produce matrices of sensible shape for a small
example system.
"""

from __future__ import annotations

import pymatgen as pmg
import sympy as sp
from pymatgen.core import Lattice, Structure

from src.jump_network_builder.geometry import build_jump_network_from_poscar
from src.matrix_builder.matrix_builder import (
    build_full_transition_matrix,
    build_matrices_from_yaml,
    export_axis_matrices_to_mathematica,
    MatrixStageJumpNetwork,
    print_axis_build_summary,
)
from src.transport_coefficient_solver.io import read_axis_payloads


def test_build_full_transition_matrix_uses_jumpnetwork(monkeypatch, tmp_path):
    """The wrapper should delegate directly to JumpNetwork.get_full_tp_matrix."""

    class DummyJumpNet:
        def __init__(self):
            self.called_with = None

        def get_full_tp_matrix(self, subst_w):  # noqa: D401 - simple dummy
            self.called_with = subst_w
            return sp.zeros(2, 2)

    dummy = DummyJumpNet()

    M = build_full_transition_matrix(dummy, {"foo": 1})

    assert M.shape == (2, 2)
    assert dummy.called_with == {"foo": 1}


def test_build_matrices_from_yaml_smoke(tmp_path):
    """Smoke test: hermetic Stage-1 -> Stage-2 matrix build on synthetic input."""

    lattice = Lattice.cubic(3.0)
    species = ["X", "X"]
    coords = [(0.0, 0.0, 0.0), (0.5, 0.5, 0.5)]
    struct = Structure(lattice, species, coords)
    poscar_file = tmp_path / "POSCAR"
    pmg.io.vasp.Poscar(struct).write_file(poscar_file)

    yaml_file = tmp_path / "input.yaml"
    yaml_file.write_text("fname: POSCAR\naxes:\n  - [1, 0, 0]\n")
    params = {
        "fname": str(poscar_file),
        "axes": [[1, 0, 0]],
        "tweak_transition_rates": {},
        "max_jump_length": 2.6,
    }
    out_dir = tmp_path / "output"
    out_dir.mkdir()

    build_jump_network_from_poscar(
        str(poscar_file),
        {},
        str(out_dir),
        crystal_matrix=None,
        max_jump_length=float(params["max_jump_length"]),
    )

    result = build_matrices_from_yaml(str(yaml_file), params, out_dir)
    export_axis_matrices_to_mathematica(result, params, out_dir)

    assert result.full_tp_matrix.shape[0] == result.full_tp_matrix.shape[1]
    assert result.full_tp_matrix.shape[0] > 0
    payloads = read_axis_payloads(out_dir)
    assert set(payloads.keys()) == set(result.axis_results.keys())
    first_axis = next(iter(payloads))
    assert payloads[first_axis].labels


def test_print_axis_build_summary_includes_interplanar_distance(capsys):
    class _FakeAxisResult:
        def __init__(self):
            self.reordered_tp_matrix = sp.zeros(2, 2)
            self.min_jump_length = None
            self.interplanar_distance_squared = sp.Symbol("d_sq")

    class _FakeResult:
        def __init__(self):
            self.axis_results = {(1, 0, 0): _FakeAxisResult()}

    print_axis_build_summary(_FakeResult())
    out = capsys.readouterr().out
    assert "interplanar_distance_squared=d_sq" in out
    assert "effective_transport_distance_squared=d_sq" in out


def test_compute_squared_interplanar_distance_uses_hex_metric():
    a, c = sp.symbols("a c", positive=True)

    fake = type(
        "_Fake",
        (),
        {
            "lat_const": sp.Matrix(
                [
                    [a, 0, 0],
                    [-a / 2, sp.sqrt(3) * a / 2, 0],
                    [0, 0, c],
                ]
            )
        },
    )()

    dist_sq = MatrixStageJumpNetwork.compute_squared_interplanar_distance(
        fake,
        sp.Matrix([2, -1, 0]),
    )
    # The convention is the true interplanar spacing d^2 = 1 / (h^T G^-1 h),
    # where G is the metric tensor built from the conventional hexagonal
    # basis. For the conventional hex basis used here and Miller indices
    # (2, -1, 0), this evaluates to a**2 / 4.
    assert sp.simplify(dist_sq - a**2 / 4) == 0


def test_compute_squared_interplanar_distance_uses_conventional_metric_even_when_structure_matrix_present():
    """Squared metric distance must be computed from the conventional metric.

    Historically, JumpNetwork.compute_squared_interplanar_distance used a
    metric built directly from conventional lattice constants. The Stage-2
    adapter briefly applied the conventional→primitive matrix here, which
    caused incorrect values and anisotropic DCs in cubic/bcc examples.

    After the fix, even if a conv2prim_matrix is available on the object, the
    helper must still use the conventional metric derived from ``lat_const``.

    For a cubic lattice with conventional lattice constant ``a`` and Miller
    indices (0, 1, 1), the squared interplanar distance is
    ``d^2 = 1 / (h^T G^-1 h) = a**2 / 2`` (with G = diag(a**2, a**2, a**2)).
    This test asserts that value and thus confirms that conv2prim_matrix does
    not alter the metric used for this computation.
    """

    a = sp.symbols("a", positive=True)

    fake = type(
        "_Fake",
        (),
        {
            "lat_const": sp.diag(a, a, a),
            # conv2prim_matrix is intentionally present but must be ignored by
            # compute_squared_interplanar_distance.
            "conv2prim_matrix": sp.Matrix(
                [
                    [sp.Rational(1, 2), sp.Rational(1, 2), sp.Rational(1, 2)],
                    [sp.Rational(1, 2), sp.Rational(-1, 2), sp.Rational(-1, 2)],
                    [sp.Rational(-1, 2), sp.Rational(1, 2), sp.Rational(-1, 2)],
                ]
            ),
        },
    )()

    dist_sq = MatrixStageJumpNetwork.compute_squared_interplanar_distance(
        fake,
        sp.Matrix([0, 1, 1]),
    )
    assert sp.simplify(dist_sq - a**2 / 2) == 0


def test_compute_squared_interplanar_distance_cubic_111_and_110():
    """For cubic lattices the current code gives d^2(111)=a^2/3 and d^2(110)=a^2/2.

    This explicitly exercises the convention used in the diffusion prefactor
    for simple cubic/bcc examples: we take the true interplanar spacing
    ``d^2 = 1 / (h^T G^-1 h)`` with G built from the conventional lattice
    basis.
    """

    a = sp.symbols("a", positive=True)

    fake = type(
        "_Fake",
        (),
        {
            "lat_const": sp.diag(a, a, a),
        },
    )()

    d2_111 = MatrixStageJumpNetwork.compute_squared_interplanar_distance(
        fake,
        sp.Matrix([1, 1, 1]),
    )
    d2_110 = MatrixStageJumpNetwork.compute_squared_interplanar_distance(
        fake,
        sp.Matrix([1, 1, 0]),
    )

    assert sp.simplify(d2_111 - a**2 / 3) == 0
    assert sp.simplify(d2_110 - a**2 / 2) == 0
