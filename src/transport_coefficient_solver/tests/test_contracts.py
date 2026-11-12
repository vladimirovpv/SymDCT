from __future__ import annotations

import pytest
import sympy as sp

from src.transport_coefficient_solver.contracts import AxisMatrixPayload


def test_axis_matrix_payload_pickle_roundtrip() -> None:
    payload = AxisMatrixPayload(
        axis_hkl_structure=(1, 0, 0),
        axis_hkl_conventional=(0, 0, 1),
        transition_matrix=sp.Matrix([[-sp.Symbol("W1"), sp.Symbol("W1")], [sp.Symbol("W1"), -sp.Symbol("W1")]]),
        labels=("A0", "B1"),
        interplanar_distance_squared=sp.Symbol("c") ** 2,
        t_symbol=sp.Symbol("t"),
        w_substitutions={sp.Symbol("W1"): sp.Symbol("WAB")},
    )

    serialized = payload.to_pickle_payload()
    restored = AxisMatrixPayload.from_pickle_payload(serialized)

    assert restored.axis_hkl_structure == (1, 0, 0)
    assert restored.axis_hkl_conventional == (0, 0, 1)
    assert restored.labels == ("A0", "B1")
    assert restored.transition_matrix == payload.transition_matrix
    assert restored.min_jump_length is None
    assert restored.interplanar_distance_squared == sp.Symbol("c") ** 2
    assert restored.t_symbol == sp.Symbol("t")
    assert restored.w_substitutions == {sp.Symbol("W1"): sp.Symbol("WAB")}
    assert restored.jump_label_display_map is None


def test_axis_matrix_payload_rejects_label_size_mismatch() -> None:
    bad_payload = {
        "contract": "transport-coefficient-axis-payload",
        "contract_version": 1,
        "axis_hkl_structure": [1, 0, 0],
        "transition_matrix": sp.eye(2),
        "labels": ["A0"],
    }
    with pytest.raises(ValueError, match="labels length must match"):
        AxisMatrixPayload.from_pickle_payload(bad_payload)
