"""Adapters between Stage-2 matrix outputs and Stage-3 solver contracts."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

import sympy as sp

from .contracts import AxisMatrixPayload
from .io import read_axis_payloads

if TYPE_CHECKING:
    from src.matrix_builder.matrix_builder import MatrixBuildResult


def _build_w_substitutions_from_jump_network(jump_network) -> dict[sp.Symbol, sp.Expr]:
    if not hasattr(jump_network, "jump_probabilities"):
        return {}
    jump_probabilities = jump_network.jump_probabilities
    if not hasattr(jump_probabilities, "get_jump_probabilities"):
        return {}

    values = jump_probabilities.get_jump_probabilities()
    substitutions: dict[sp.Symbol, sp.Expr] = {}
    for idx, expr in values.items():
        substitutions[sp.Symbol(f"W{int(idx)}", real=True)] = sp.nsimplify(expr)
    return substitutions


def _build_jump_label_map_from_graph(graph) -> dict[int, str]:
    """Build global_jump_index -> 'A->B' mapping from graph edges."""
    label_map: dict[int, str] = {}
    
    # Handle both networkx graphs and JumpNetworkGraph wrappers
    nx_graph = graph
    if hasattr(graph, "graph") and hasattr(graph.graph, "edges"):
        # It's a JumpNetworkGraph, extract the networkx graph
        nx_graph = graph.graph
    elif not hasattr(graph, "edges"):
        # Not a valid graph
        return label_map
    
    for start, end, _key, attrs in nx_graph.edges(keys=True, data=True):
        glob_idx = int(attrs["glob_jump_cnt"])
        # Try both 'label' and 'lbl' attributes for node labels
        start_label = nx_graph.nodes[start].get("label") or nx_graph.nodes[start].get("lbl", str(start))
        end_label = nx_graph.nodes[end].get("label") or nx_graph.nodes[end].get("lbl", str(end))
        label = f"{start_label}->{end_label}"
        label_map[glob_idx] = label
    return label_map


def payloads_from_matrix_build_result(
    result: "MatrixBuildResult",
) -> dict[tuple[int, int, int], AxisMatrixPayload]:
    """Adapt Stage-2 `MatrixBuildResult` into Stage-3 axis payloads."""
    jump_network = result.jump_network
    t_symbol = jump_network.get_t() if hasattr(jump_network, "get_t") else None
    w_substitutions = _build_w_substitutions_from_jump_network(jump_network)
    
    # Build jump label map from the graph
    graph = jump_network.graph if hasattr(jump_network, "graph") else None
    jump_label_map = _build_jump_label_map_from_graph(graph) if graph else {}
    jump_label_display_map: dict[str, str] = {}
    alias_map = {}
    if hasattr(result, "symmetry_ctx") and hasattr(result.symmetry_ctx, "get_symbol_alias_map"):
        alias_map = dict(result.symmetry_ctx.get_symbol_alias_map())
    canonical_labels = list(getattr(result, "labels", ()))
    if not canonical_labels and result.axis_results:
        canonical_labels = list(next(iter(result.axis_results.values())).labels)
    for canon in canonical_labels:
        jump_label_display_map.setdefault(str(canon), str(alias_map.get(str(canon), canon)))

    # result.axis_results is keyed by hkl_conventional (see
    # MatrixBuildResult.axis_results), so payloads inherit the same,
    # always-available canonical identifier -- no fallback needed.
    payloads: dict[tuple[int, int, int], AxisMatrixPayload] = {}
    for hkl_conventional, axis_result in result.axis_results.items():
        payloads[hkl_conventional] = AxisMatrixPayload(
            axis_hkl_structure=tuple(int(v) for v in axis_result.hkl_structure),
            axis_hkl_conventional=tuple(int(v) for v in hkl_conventional),
            transition_matrix=sp.Matrix(axis_result.reordered_tp_matrix),
            labels=tuple(str(lbl) for lbl in axis_result.labels),
            min_jump_length=axis_result.min_jump_length,
            interplanar_distance_squared=axis_result.interplanar_distance_squared,
            t_symbol=t_symbol,
            w_substitutions=w_substitutions,
            jump_label_map=jump_label_map if jump_label_map else None,
            jump_label_display_map=jump_label_display_map if jump_label_display_map else None,
        )
    return payloads


def payloads_from_output_dir(output_dir: str | Path) -> dict[tuple[int, int, int], AxisMatrixPayload]:
    """Load Stage-3 payloads from Stage-2 output artifacts."""
    payloads = read_axis_payloads(output_dir)
    if not payloads:
        raise FileNotFoundError(
            f"No Stage-3 payload artifacts found in {output_dir!r}. "
            "Run matrix_builder to generate tc_axis_payload_<h>_<k>_<l>.pickle files."
        )
    return payloads
