from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from src.transport_coefficient_solver.adapters import payloads_from_matrix_build_result
from src.transport_coefficient_solver.contracts import AxisMatrixPayload
from src.transport_coefficient_solver.io import write_axis_payloads
from src.transport_coefficient_solver.transport_coefficient_solver import (
    solve_from_matrix_result,
    solve_from_output_dir,
    solve_from_payloads,
)


def _example_payload(axis: tuple[int, int, int]) -> AxisMatrixPayload:
    t = sp.Symbol("t")
    w = sp.Symbol("W1")
    matrix = sp.Matrix([[-w, w * sp.cos(t)], [w * sp.cos(t), -w]])
    return AxisMatrixPayload(
        axis_hkl_structure=axis,
        axis_hkl_conventional=axis,
        transition_matrix=matrix,
        labels=("A0", "B1"),
        interplanar_distance_squared=sp.Integer(1),
        t_symbol=t,
    )


def test_disk_and_in_memory_paths_produce_same_result(tmp_path) -> None:
    payload = _example_payload((0, 0, 1))
    expected = solve_from_payloads({(0, 0, 1): payload})

    write_axis_payloads({(0, 0, 1): payload}, tmp_path)
    disk_result = solve_from_output_dir(tmp_path, export=False)

    assert disk_result == expected


@dataclass
class _FakeAxisResult:
    hkl_structure: tuple[int, int, int]
    reordered_tp_matrix: sp.Matrix
    labels: tuple[str, ...]
    min_jump_length: sp.Expr | None
    interplanar_distance_squared: sp.Expr


class _FakeJumpProbabilities:
    def get_jump_probabilities(self) -> dict[int, sp.Expr]:
        return {1: sp.Symbol("W1")}


class _FakeJumpNetwork:
    def __init__(self) -> None:
        self.jump_probabilities = _FakeJumpProbabilities()

    def get_t(self) -> sp.Symbol:
        return sp.Symbol("t")


@dataclass
class _FakeMatrixBuildResult:
    jump_network: _FakeJumpNetwork
    axis_results: dict[tuple[int, int, int], _FakeAxisResult]


def test_solve_from_matrix_result_uses_adapter_contract() -> None:
    payload = _example_payload((1, 0, 0))
    fake_result = _FakeMatrixBuildResult(
        jump_network=_FakeJumpNetwork(),
        axis_results={
            (1, 0, 0): _FakeAxisResult(
                hkl_structure=(1, 0, 0),
                reordered_tp_matrix=payload.transition_matrix,
                labels=payload.labels,
                min_jump_length=payload.min_jump_length,
                interplanar_distance_squared=payload.interplanar_distance_squared,
            )
        },
    )

    adapted = payloads_from_matrix_build_result(fake_result)
    result = solve_from_matrix_result(fake_result)

    assert (1, 0, 0) in result
    assert result == solve_from_payloads(adapted)
