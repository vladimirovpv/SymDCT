"""Run Stage 1 and 2 for an example and test lumpability of every axis.

Usage: uv run python scripts/check_lumpability_examples.py examples/hcp_vac/input/hcp_vac.yaml
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import src.matrix_builder.matrix_builder as mb
from src.common_utils.common_init import init_from_yaml
from src.common_utils.utils import load_cutoffs_config
from src.jump_network_builder.geometry import build_jump_network_from_poscar
from src.matrix_builder.lumpability import check_lumpability
from src.symmetry_context.jump_symmetry_context import JumpNetworkSymmetryContext

_original = mb.find_equiv_and_reduce_matrix
reports = []


def _checked(equiv_sites, full_tp_matrix, hkl_structure, k_symbols, t_symbol):
    rep = check_lumpability(equiv_sites, full_tp_matrix, hkl_structure, k_symbols, t_symbol)
    reports.append((tuple(hkl_structure), len(set(equiv_sites)), len(equiv_sites), rep))
    return _original(equiv_sites, full_tp_matrix, hkl_structure, k_symbols, t_symbol)


mb.find_equiv_and_reduce_matrix = _checked


def main() -> int:
    yaml_file, params, output_dir = init_from_yaml()
    cutoffs, max_jump_length = load_cutoffs_config(yaml_file, params)
    graph, _, _ = build_jump_network_from_poscar(
        params["fname"], cutoffs, str(output_dir), params.get("matrix", None), max_jump_length
    )
    ctx = JumpNetworkSymmetryContext(params["fname"], params.get("matrix", None))
    mb.build_matrices(graph=graph, symmetry_ctx=ctx, params=params)
    print("\n=== LUMPABILITY CHECK ===")
    for hkl, n_cls, n_sites, rep in reports:
        print(f"axis {hkl}: {n_sites} sites -> {n_cls} classes: {rep.summary()}")
    return 0 if all(r.ok for *_, r in reports) else 1


if __name__ == "__main__":
    raise SystemExit(main())
