from __future__ import annotations

import sympy as sp

from src.transport_coefficient_solver.contracts import AxisMatrixPayload
from src.transport_coefficient_solver.solver import (
    compute_axis_diffusion_coefficient,
    solve_all_axes,
)


def test_compute_axis_diffusion_coefficient_basic_case() -> None:
    t = sp.Symbol("t")
    w = sp.Symbol("W1")
    matrix = sp.Matrix([[-w, w * sp.cos(t)], [w * sp.cos(t), -w]])
    base_payload = AxisMatrixPayload(
        axis_hkl_structure=(0, 0, 1),
        axis_hkl_conventional=(0, 0, 1),
        transition_matrix=matrix,
        labels=("A0", "B1"),
        interplanar_distance_squared=sp.Integer(1),
        t_symbol=t,
    )
    scaled_payload = AxisMatrixPayload(
        axis_hkl_structure=(0, 0, 1),
        axis_hkl_conventional=(0, 0, 1),
        transition_matrix=matrix,
        labels=("A0", "B1"),
        interplanar_distance_squared=sp.Integer(1),
        min_jump_length=sp.Rational(1, 2),
        t_symbol=t,
    )

    result = compute_axis_diffusion_coefficient(base_payload)
    scaled_result = compute_axis_diffusion_coefficient(scaled_payload)

    assert result.axis_hkl_structure == (0, 0, 1)
    assert sp.simplify(result.diffusion_coefficient - (w / 2)) == 0
    assert sp.simplify(scaled_result.diffusion_coefficient - result.diffusion_coefficient / 4) == 0


def test_solve_all_axes_with_substitutions() -> None:
    t = sp.Symbol("t")
    w = sp.Symbol("W1")
    matrix = sp.Matrix([[-w, w * sp.cos(t)], [w * sp.cos(t), -w]])
    # Deliberately distinct from axis_hkl_structure to prove results are
    # keyed by axis_hkl_conventional, not axis_hkl_structure.
    payload = AxisMatrixPayload(
        axis_hkl_structure=(1, 0, 0),
        axis_hkl_conventional=(2, 0, 0),
        transition_matrix=matrix,
        labels=("A0", "B1"),
        interplanar_distance_squared=sp.Integer(4),
        t_symbol=t,
        w_substitutions={w: sp.Symbol("WAB")},
    )

    results = solve_all_axes([payload], apply_w_substitutions=True)

    assert (1, 0, 0) not in results
    assert (2, 0, 0) in results
    assert sp.simplify(results[(2, 0, 0)] - (2 * sp.Symbol("WAB"))) == 0
