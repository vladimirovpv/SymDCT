from __future__ import annotations

import argparse
import logging
import os
from dataclasses import dataclass
from pathlib import Path

from sympy import simplify
from sympy.parsing.sympy_parser import parse_expr

from src.common_utils.energetics_io import (
    ensure_energetics_in_sync,
    get_jump_info,
    get_lattice_parameters,
    load_energetics_yaml,
)
from src.common_utils.label_utils import (
    map_jump_label_to_id,
    mathematica_subscript_to_python,
)
from src.common_utils.utils import configure_cli_logging, get_example_name, get_hkl_from_filename
from src.transport_coefficient_evaluator.dc_numeric import DiffCoeff, plot_arrhenius_with_regression

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

configure_cli_logging()


@dataclass
class EvaluationContext:
    """Shared context derived from the production energetics file.

    Loaded once and reused for every DC-expression file (axis) to avoid
    redundant I/O and computation when evaluating multiple axes.
    """
    energetics: dict
    jump_info: list
    em: tuple
    fq: tuple
    map_label2index: dict
    lattice_parameters: tuple  # (a, b, c)


def load_evaluation_context(
    production_path: Path,
    template_path: Path | None = None,
    length_tolerance: float | None = None,
) -> EvaluationContext:
    """Load energetics and derive all axis-independent quantities."""
    production_path = ensure_energetics_in_sync(
        template_path, production_path, length_tolerance=length_tolerance
    )
    energetics = load_energetics_yaml(str(production_path))
    jump_info = get_jump_info(energetics)
    _, _, _, _, em, fq = zip(*sorted(jump_info))
    alias_map = energetics.get("symbol_alias_map", {})
    map_label2index = map_jump_label_to_id(jump_info, alias_map)
    lattice_parameters = get_lattice_parameters(energetics)
    return EvaluationContext(
        energetics=energetics,
        jump_info=jump_info,
        em=em,
        fq=fq,
        map_label2index=map_label2index,
        lattice_parameters=lattice_parameters,
    )


def evaluate_dc_file(dc_file: Path, ctx: EvaluationContext) -> int:
    """Evaluate one symbolic DC expression file and produce an Arrhenius plot."""
    if not dc_file.exists():
        logger.error("dc expression file not found: %s", dc_file)
        raise SystemExit(2)

    with open(dc_file) as fh:
        dc_expr = fh.read().replace("^", "**")

    print(f"Read diffusion coefficient expression from file: {dc_file}")
    print(f"Expression:\n   {dc_expr}")

    dc_expr = mathematica_subscript_to_python(dc_expr)
    s_expr = parse_expr(dc_expr, evaluate=False)
    reduced_expr = s_expr.subs(ctx.map_label2index)
    dc_expr = str(simplify(reduced_expr))

    folder_name = get_example_name(dc_file)
    axis = get_hkl_from_filename(dc_file.name)

    dc = DiffCoeff(dc_expr, ctx.em, ctx.fq, *ctx.lattice_parameters)

    model = dc.get_model()
    title = (
        f"Diffusion Coefficient for {folder_name} along axis {axis}"
        + f"\nwithin a {model}-frequency model"
    )

    T_min = 100
    T_max = 2000
    # to show the plots instead of saving them to files, set figure_filename=None
    plot_name = dc_file.parent / f"arrhenius_plot_{folder_name}_{axis[0]}_{axis[1]}_{axis[2]}_freq{model}.png"
    plot_arrhenius_with_regression(dc, T_min, T_max, title=title, figure_filename=plot_name)

    return plot_name


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Evaluator for transport coefficients using a production energetics YAML."
    )
    p.add_argument("dc_expr", help="Path to file containing the symbolic DC expression (text file)")
    p.add_argument("production", help="Path to production energetics YAML (required)")
    p.add_argument(
        "template",
        nargs="?",
        help=(
            "Path to energetics_template.yaml (optional). "
            "If omitted the tool will look next to the production file in "
            "../output/energetics_template.yaml"
        ),
    )
    p.add_argument(
        "--length-tolerance",
        type=float,
        help="Override length tolerance used when comparing jumps",
    )
    p.add_argument(
        "--auto-sync",
        action="store_true",
        help="Attempt automatic synchronization (same as setting DC_AUTO_SYNC=1)",
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = _parse_args(argv)

    if args.auto_sync:
        os.environ["DC_AUTO_SYNC"] = "1"

    ctx = load_evaluation_context(
        production_path=Path(args.production),
        template_path=Path(args.template) if args.template else None,
        length_tolerance=args.length_tolerance,
    )
    return evaluate_dc_file(Path(args.dc_expr), ctx)


if __name__ == "__main__":
    raise SystemExit(main())
