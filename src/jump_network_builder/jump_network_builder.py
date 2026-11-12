"""Orchestration script for the geometry stage (jump_network_builder).

This script reads a POSCAR and cutoffs.yaml, constructs the jump network,
registers jumps, dumps an adjacency list and an energetics template YAML.
"""

from __future__ import annotations

import pathlib as pl
import sys

# When this module is executed directly as a script (e.g. via the IDE),
# Python may not have the package root on sys.path which causes
# `ModuleNotFoundError: No module named 'src'`.  Insert a small
# bootstrap that adds the repository `public` directory (two levels up)
# to sys.path so `from src...` imports work the same as when running
# with `python -m src...`.
if __package__ is None or __package__ == "":
    import os
    _this_dir = os.path.dirname(os.path.abspath(__file__))
    # src/<package>/... -> go two levels up to reach the repository root
    _public_root = os.path.abspath(os.path.join(_this_dir, "..", ".."))
    if _public_root not in sys.path:
        sys.path.insert(0, _public_root)
    # Make the current working directory the `public` root so that
    # relative paths (like examples/...) resolve the same as when the
    # module is executed with `python -m` from the `public` directory.
    try:
        os.chdir(_public_root)
    except Exception:
        # Non-fatal: if chdir fails, proceed without changing cwd.
        pass

# from src.rw_yaml import load_data
from src.common_utils.yaml_with_comments import load_data
from src.jump_network_builder.geometry import build_jump_network_from_poscar
from src.common_utils.common_init import init_from_yaml
from src.common_utils.utils import configure_cli_logging, load_cutoffs_config


def main(argv=None):
    if argv is None:
        argv = sys.argv

    configure_cli_logging()
    inp_file, params, output_dir = init_from_yaml()

    poscar_fname = params["fname"]
    crystal_matrix = params.get("matrix", None)

    if len(argv) > 2:
        # Explicit override: use the given cutoffs.yaml regardless of what
        # (if anything) sits next to the input yaml.
        cutoffs = load_data(argv[2])
        if cutoffs and "bonds" in cutoffs:
            max_jump_length = max(bond_length for bond_length in cutoffs["bonds"].values())
        else:
            max_jump_length = params.get("max_jump_length", 2.5)
    else:
        # Auto-discover cutoffs.yaml next to the input yaml; falls back to
        # the global max_jump_length from params (cutoffs = {}) when absent.
        cutoffs, max_jump_length = load_cutoffs_config(inp_file, params)

    G, jr, jump_sites = build_jump_network_from_poscar(
        poscar_fname,
        cutoffs or {},
        str(output_dir),
        crystal_matrix,
        max_jump_length,
    )
    print(
        f"Registered {jr.get_number_of_jumps()} jump types; "
        f"graph JSON at: {output_dir / pl.Path('jump_structure.node_link.json')}"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
