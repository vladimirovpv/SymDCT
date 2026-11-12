from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from src.transport_coefficient_solver.adapters import payloads_from_matrix_build_result


@dataclass
class _FakeAxisResult:
    hkl_structure: tuple[int, int, int]
    reordered_tp_matrix: sp.Matrix
    labels: tuple[str, ...]
    min_jump_length: sp.Expr | None
    interplanar_distance_squared: sp.Expr


class _FakeJumpProbabilities:
    def get_jump_probabilities(self) -> dict[int, sp.Expr]:
        return {0: sp.Symbol("W0"), 1: 0}


class _FakeJumpNetwork:
    def __init__(self) -> None:
        self.jump_probabilities = _FakeJumpProbabilities()

    def get_t(self) -> sp.Symbol:
        return sp.Symbol("t")


@dataclass
class _FakeMatrixBuildResult:
    jump_network: _FakeJumpNetwork
    axis_results: dict[tuple[int, int, int], _FakeAxisResult]


def test_payloads_from_matrix_build_result_adapts_axis_data() -> None:
    # axis_results is keyed by hkl_conventional; hkl_structure is
    # deliberately different here to prove the adapter threads both
    # through correctly rather than conflating them.
    fake = _FakeMatrixBuildResult(
        jump_network=_FakeJumpNetwork(),
        axis_results={
            (0, 0, 1): _FakeAxisResult(
                hkl_structure=(1, 0, 0),
                reordered_tp_matrix=sp.eye(2),
                labels=("A0", "B1"),
                min_jump_length=sp.Rational(1, 2),
                interplanar_distance_squared=sp.Symbol("c") ** 2,
            )
        },
    )

    payloads = payloads_from_matrix_build_result(fake)

    payload = payloads[(0, 0, 1)]
    assert payload.axis_hkl_structure == (1, 0, 0)
    assert payload.axis_hkl_conventional == (0, 0, 1)
    assert payload.transition_matrix == sp.eye(2)
    assert payload.labels == ("A0", "B1")
    assert payload.min_jump_length == sp.Rational(1, 2)
    assert payload.interplanar_distance_squared == sp.Symbol("c") ** 2
    assert payload.t_symbol == sp.Symbol("t")
    assert payload.w_substitutions == {sp.Symbol("W0", real=True): sp.Symbol("W0"), sp.Symbol("W1", real=True): 0}
