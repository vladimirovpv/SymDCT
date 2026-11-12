"""High-level Stage-3 orchestration for symbolic diffusion coefficients."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Mapping, Sequence

import sympy as sp

from src.common_utils.utils import configure_cli_logging
from .adapters import payloads_from_matrix_build_result, payloads_from_output_dir
from .contracts import AxisMatrixPayload
from .io import write_diffusion_coefficients
from .solver import solve_all_axes


def solve_from_payloads(
    payloads: Mapping[tuple[int, int, int], AxisMatrixPayload]
    | Sequence[AxisMatrixPayload],
    *,
    min_jump_length: sp.Expr | int | float | None = None,
    apply_w_substitutions: bool = True,
) -> dict[tuple[int, int, int], sp.Expr]:
    """Solve symbolic diffusion coefficients from Stage-3 payloads."""
    return solve_all_axes(
        payloads,
        min_jump_length=min_jump_length,
        apply_w_substitutions=apply_w_substitutions,
    )


def solve_from_output_dir(
    output_dir: str | Path,
    *,
    min_jump_length: sp.Expr | int | float | None = None,
    apply_w_substitutions: bool = True,
    export: bool = True,
) -> dict[tuple[int, int, int], sp.Expr]:
    """Load payloads from disk and solve symbolic diffusion coefficients."""
    payloads = payloads_from_output_dir(output_dir)
    coeffs = solve_from_payloads(
        payloads,
        min_jump_length=min_jump_length,
        apply_w_substitutions=apply_w_substitutions,
    )
    if export:
        write_diffusion_coefficients(coeffs, output_dir)
    return coeffs


def solve_from_matrix_result(
    matrix_result,
    *,
    min_jump_length: sp.Expr | int | float | None = None,
    apply_w_substitutions: bool = True,
) -> dict[tuple[int, int, int], sp.Expr]:
    """Solve symbolic diffusion coefficients directly from Stage-2 in-memory output."""
    payloads = payloads_from_matrix_build_result(matrix_result)
    return solve_from_payloads(
        payloads,
        min_jump_length=min_jump_length,
        apply_w_substitutions=apply_w_substitutions,
    )


def _parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Solve symbolic diffusion coefficients from Stage-2 axis payload artifacts."
    )
    parser.add_argument(
        "path",
        help="Path to output dir (containing tc_axis_payload_*.pickle) or input YAML file.",
    )
    parser.add_argument(
        "--no-substitutions",
        action="store_true",
        help="Do not apply w_substitutions from payload metadata.",
    )
    parser.add_argument(
        "--no-export",
        action="store_true",
        help="Do not write dc_symbolic_<h>_<k>_<l>.txt/.m files.",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    configure_cli_logging()
    args = _parse_args(argv)
    path = Path(args.path)
    output_dir = path.parent.parent / "output" if path.suffix.lower() in {".yaml", ".yml"} else path

    coeffs = solve_from_output_dir(
        output_dir=output_dir,
        apply_w_substitutions=not args.no_substitutions,
        export=not args.no_export,
    )
    print("Solved symbolic diffusion coefficients for axes:")
    for axis, expr in sorted(coeffs.items()):
        print(f"  {axis}: {expr}")
    return 0
