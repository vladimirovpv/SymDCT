"""Utility helpers shared across public modules.

This is the common home for the trimmed ``utils`` functions, so that
they can be imported as ``from src.common_utils.utils import ...``.
"""

from __future__ import annotations

import logging
import re
import sys
from pathlib import Path
from typing import Any, Dict, Mapping, Tuple

import sympy
from sympy import Rational, nsimplify


logger = logging.getLogger(__name__)
# set logging level to INFO to avoid debug spam from SymPy
logger.setLevel(logging.INFO)


def configure_cli_logging(level: int = logging.INFO) -> None:
    """Install a simple root logger formatter for command-line entry points."""

    logging.basicConfig(
        level=level,
        format="%(levelname)s: %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
        force=True,
    )


def to_rational(value: float, tolerance: float = 1e-8, max_denominator: int = 100) -> Rational:
    """Convert a float value to a SymPy :class:`Rational`."""

    try:
        simplified_value = nsimplify(value, tolerance=tolerance)
        if not simplified_value.is_number:
            raise ValueError(f"Simplified value is not a number: {simplified_value}")
        rational_value = (
            Rational(simplified_value).limit_denominator(max_denominator)
            if abs(value) > tolerance
            else Rational(0, 1)
        )
    except Exception as exc:  # pragma: no cover - defensive
        raise TypeError(f"Failed to convert value to rational: {exc}") from exc
    return rational_value


def get_args(prompt: str = "Enter name of file with parameters (yaml format expected):") -> str:
    """Get name of input YAML file from command line or prompt interactively."""

    if len(sys.argv) >= 2:
        inp_file = sys.argv[1]
    else:  # pragma: no cover - interactive path
        inp_file = input(prompt + " ")
        print()
    return inp_file


def get_input_parameters(inp_file: str) -> Dict[str, Any]:
    """Load input parameters from YAML file and print a compact summary."""

    separation_line_length = 80
    print(f"\nLoaded input parameters from: {inp_file}")
    print("-" * separation_line_length)
    # Import here to avoid a circular import between ``utils`` and
    # ``jump_network_builder.geometry``.  The geometry module needs
    # :func:`to_rational` from this file, while the YAML loader lives in
    # :mod:`src.common_utils.yaml_with_comments`.  Delaying the import keeps
    # ``utils`` free of hard dependencies at module-import time.
    from src.common_utils.yaml_with_comments import load_data

    params: Dict[str, Any] = load_data(inp_file)
    if not params:
        return {}

    # max length of the parameter names for formatting
    max_length = max(len(key) for key in params.keys())
    print(f"{{:<{max_length}}} | {{}}".format("Parameter", "Value"))
    print("-" * separation_line_length)
    for key, value in params.items():
        if isinstance(value, list):
            for i, v in enumerate(value):
                if i == 0:
                    print(f"{key:<{max_length}} | {v}")
                else:
                    print(f"{'':<{max_length}} | {v}")
        else:
            print(f"{key:<{max_length}} | {value}")

        logger.debug("%s = %r", key, value)
    print("-" * separation_line_length)
    print()
    return params


def replace_W_with_subscript(input_string: str) -> str:
    """Replace tokens like ``W1`` or ``W_ab`` by Mathematica-style subscripts."""

    # W followed by digits -> Subscript[W,n]
    def replacement_digits(match: re.Match[str]) -> str:
        return f"Subscript[W,{match.group(1)}]"

    # W followed by two letters -> Subscript[W,a,b]
    def replacement_letters(match: re.Match[str]) -> str:
        return f"Subscript[W,{match.group(1)},{match.group(2)}]"

    # W_{ab} -> Subscript[W,a,b]
    def replacement_braces(match: re.Match[str]) -> str:
        return f"Subscript[W,{match.group(1)},{match.group(2)}]"

    # W_{ab}^{n} -> Supersubscript[W,ab,n]
    def replacement_superscript(match: re.Match[str]) -> str:
        return f"Supersubscript[W,{match.group(1)}{match.group(2)}, {match.group(3)}]"

    result = re.sub(r"W(\d+)", replacement_digits, input_string)
    result = re.sub(r"W([a-zA-Z])([a-zA-Z])", replacement_letters, result)
    result = re.sub(r"W_\{([a-zA-Z])([a-zA-Z])\}\^\{(\d+)\}", replacement_superscript, result)
    result = re.sub(r"W_\{([a-zA-Z])([a-zA-Z])\}", replacement_braces, result)
    return result


def export_to_mathematica(matrix: sympy.Matrix, filename_wo_extension: str, is_replace_W_with_subscript: bool = False) -> None:
    """Export a SymPy matrix to a Mathematica ``.m`` file."""

    math_filename = filename_wo_extension + ".m"
    logger.info("Saving results to %s", math_filename)
    with open(math_filename, "w", encoding="utf-8") as f:
        n_rows, n_cols = matrix.shape
        f.write("{")
        for i in range(n_rows):
            f.write("{")
            row = []
            for j in range(n_cols):
                m_ij = f"{sympy.printing.mathematica.mathematica_code(matrix[i, j])}"
                m_ij_subscript = replace_W_with_subscript(m_ij) if is_replace_W_with_subscript else m_ij
                row.append(m_ij_subscript.replace(" ", ""))
            f.write(", ".join(row))
            if i != n_rows - 1:
                f.write("},\n")
        f.write("}}\n")

def get_example_name(path: Path) -> str:
    """Extract the example name from the file path."""
    # Folder name following "examples/"
    folder_match = re.search(r'examples/([^/]+)/', str(path))
    folder_name = folder_match.group(1) if folder_match else None
    return folder_name

def format_axis_label(hkl_conventional, hkl_structure) -> str:
    """Format an axis label using the conventional-cell indices -- as typed
    in the input YAML ``axes`` entry -- as the primary, canonical value,
    noting the primitive-cell indices actually used internally for the
    geometric calculation whenever they differ (which happens whenever the
    POSCAR is a primitive cell)."""
    if hkl_structure is not None and tuple(hkl_structure) != tuple(hkl_conventional):
        return f"{hkl_conventional} [structure-cell basis {hkl_structure}]"
    return str(hkl_conventional)


def load_cutoffs_config(
    yaml_file: str | Path,
    params: Mapping[str, Any],
    default_max_jump_length: float = 2.5,
) -> Tuple[dict, float]:
    """Auto-discover a ``cutoffs.yaml`` sibling of the input ``yaml_file``
    and derive the neighbour-search radius from it.

    If ``cutoffs.yaml`` exists in the same input folder as ``yaml_file``,
    it is loaded and the maximum bond length across its ``bonds`` entries
    is used as the search radius (neighbours are later filtered per-pair
    against the individual cutoffs by the geometry stage). If it is
    missing, or has no usable ``bonds`` entries, this falls back to the
    single global ``max_jump_length`` from ``params`` (or
    ``default_max_jump_length`` when that key is also absent), and returns
    an empty cutoffs config -- ``cutoffs = {}`` is only the pre-file-read
    default, not a final value once a cutoffs.yaml is found.
    """
    from src.common_utils.yaml_with_comments import load_data

    cutoffs_path = Path(yaml_file).parent / "cutoffs.yaml"
    if cutoffs_path.exists():
        cutoffs_config = load_data(str(cutoffs_path)) or {}
        bonds = cutoffs_config.get("bonds", {}) if isinstance(cutoffs_config, dict) else {}
        if bonds:
            max_jump_length = max(float(v) for v in bonds.values())
            return cutoffs_config, max_jump_length

    return {}, params.get("max_jump_length", default_max_jump_length)


def get_hkl_from_filename(filename: str) -> tuple[int, int, int]:
    """Extract h, k, l from the filename."""
    match = re.search(r'dc_symbolic_(-?\d+)_(-?\d+)_(-?\d+)\.txt$', filename)
    if not match:
        raise ValueError(f"Filename does not match expected pattern: {filename}")
    return tuple(int(x) for x in match.groups())
