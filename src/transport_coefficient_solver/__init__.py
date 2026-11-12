"""Transport‑coefficient solver stage (Stage 3 of the public pipeline).

This package hosts the Stage‑3 logic of the public diffusion pipeline:
taking a full transition‑probability matrix and symmetry / graph context
from :mod:`jump_network_builder`, and computing direction‑dependent
transport coefficients symbolically.
"""

from __future__ import annotations

from .adapters import payloads_from_matrix_build_result, payloads_from_output_dir
from .contracts import AxisDiffusionResult, AxisMatrixPayload, PICKLE_CONTRACT_VERSION
from .io import read_axis_payloads, write_axis_payloads, write_diffusion_coefficients
from .solver import compute_axis_diffusion_coefficient, solve_all_axes
from .transport_coefficient_solver import (
    solve_from_matrix_result,
    solve_from_output_dir,
    solve_from_payloads,
)

__all__ = [
    "payloads_from_matrix_build_result",
    "payloads_from_output_dir",
    "AxisDiffusionResult",
    "AxisMatrixPayload",
    "PICKLE_CONTRACT_VERSION",
    "compute_axis_diffusion_coefficient",
    "read_axis_payloads",
    "solve_from_matrix_result",
    "solve_from_output_dir",
    "solve_from_payloads",
    "solve_all_axes",
    "write_axis_payloads",
    "write_diffusion_coefficients",
]
