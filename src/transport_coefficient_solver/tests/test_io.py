from __future__ import annotations

import sympy as sp

from src.transport_coefficient_solver.contracts import AxisMatrixPayload
from src.transport_coefficient_solver.io import (
    read_axis_payload,
    write_axis_payload,
    write_diffusion_coefficients,
)


def test_axis_payload_write_and_read(tmp_path) -> None:
    payload = AxisMatrixPayload(
        axis_hkl_structure=(1, 0, 0),
        axis_hkl_conventional=(2, -1, 0),
        transition_matrix=sp.eye(2),
        labels=("A0", "B1"),
        interplanar_distance_squared=sp.Integer(3),
        t_symbol=sp.Symbol("t"),
    )
    payload_file = tmp_path / "tc_axis_payload_2_-1_0.pickle"
    write_axis_payload(payload, payload_file)

    restored = read_axis_payload(payload_file)
    assert restored.axis_hkl_structure == (1, 0, 0)
    assert restored.axis_hkl_conventional == (2, -1, 0)
    assert restored.transition_matrix == sp.eye(2)
    assert restored.labels == ("A0", "B1")


def test_write_diffusion_coefficients_outputs_text_and_mathematica(tmp_path) -> None:
    coeffs = {(0, 0, 1): sp.Symbol("W1") / 2}
    files = write_diffusion_coefficients(coeffs, tmp_path, prefix="dc_symbolic")
    assert len(files) == 2
    assert (tmp_path / "dc_symbolic_0_0_1.txt").exists()
    assert (tmp_path / "dc_symbolic_0_0_1.m").exists()

