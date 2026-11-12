"""Shared initialization helpers for high-level orchestration scripts.

This module provides :func:`init_from_yaml`, which encapsulates the
common pattern used by both the jump-network builder and the
matrix-builder stage in this subproject:

* obtain the YAML path from the command line via :func:`src.utils.get_args`;
* parse it using :func:`src.utils.get_input_parameters`;
* derive an ``output/`` directory next to the YAML ``input/`` folder.

By keeping this in ``src/common_utils/`` we avoid polluting the parent
project's top-level ``src`` namespace.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Dict, Tuple

from src.common_utils.utils import get_args, get_input_parameters


def init_from_yaml() -> Tuple[str, Dict[str, Any], Path]:
    """Initialise a calculation from a YAML input file.

    Returns
    -------
    inp_file
        Path to the YAML file as obtained from :func:`get_args`.
    params
        Parsed parameter dictionary from :func:`get_input_parameters`.
    output_dir
        Directory where results should be written.  This is constructed
        as ``yaml_dir / ".." / "output"`` and created if necessary.
    """

    inp_file = get_args()
    params: Dict[str, Any] = get_input_parameters(inp_file)
    yaml_path = Path(inp_file)
    if "fname" in params:
        fname_path = Path(str(params["fname"]))
        if not fname_path.is_absolute():
            # First try to resolve relative to the YAML file location (common case)
            resolved = yaml_path.parent / fname_path
            # If that does not exist, fall back to a project‑root relative path.
            if not resolved.exists():
                fallback = Path.cwd() / fname_path
                if fallback.exists():
                    resolved = fallback
                else:
                    # As a last resort, keep the original relative path – downstream
                    # code will raise a clear error if the file truly does not exist.
                    resolved = yaml_path.parent / fname_path
            if not resolved.exists():
                raise FileNotFoundError(
                    f"POSCAR file '{resolved}' not found. "
                    f"Check the 'fname' parameter in {yaml_path}."
                )
            params["fname"] = str(resolved.resolve())
            output_dir = yaml_path.parent.parent / "output"
            output_dir.mkdir(parents=True, exist_ok=True)

    return inp_file, params, output_dir
