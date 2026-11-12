from __future__ import annotations

import json
import pickle
import pathlib as pl
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

import networkx as nx
import sympy as sp
from pymatgen.core import Structure

from src.common_utils.common_init import init_from_yaml
from src.common_utils.utils import configure_cli_logging, export_to_mathematica
from src.matrix_builder.matrix_reduction_helpers import (
    find_equiv_and_reduce_matrix,
    find_min_jump_length,
    prepare_label_data_for_reordered_matrix,
    prepare_reordered_tp_matrix,
)
from src.symmetry_context.axis_equivalence import (
    JumpAxisEquivalenceAnalyzer,
)
from src.symmetry_context.jump_symmetry_context import JumpNetworkSymmetryContext
from src.matrix_builder.lumpability import (
    check_lumpability,
    refine_to_lumpable,
    suggest_rate_relations,
)
from src.jump_network_builder.geometry import _node_link_graph



def _coerce_integer_hkl(axis: Any) -> tuple[sp.Integer, sp.Integer, sp.Integer]:
    """Validate and convert axis definition to integer Miller components."""

    if not isinstance(axis, (list, tuple)) or len(axis) != 3:
        raise ValueError(
            f"Axis must be a 3-component [h, k, l] vector, got: {axis!r}"
        )

    hkl: list[sp.Integer] = []
    for idx, component in enumerate(axis):
        expr = sp.nsimplify(component)
        try:
            as_int = int(expr)
        except Exception as exc:
            raise ValueError(
                f"Axis component at index {idx} is not integer-like: {component!r}"
            ) from exc

        if sp.simplify(expr - as_int) != 0:
            raise ValueError(
                f"Axis component at index {idx} must be integer-valued, got: {component!r}"
            )
        hkl.append(sp.Integer(as_int))

    return hkl[0], hkl[1], hkl[2]


# The data classes and the ``MatrixStageJumpNetwork`` implementation have been
# moved to ``matrix_builder/types.py`` to keep this module focused on orchestration.
from .types import AxisMatrixResult, MatrixBuildResult, MatrixStageJumpNetwork

def _load_jump_graph(graph_path: str | pl.Path) -> nx.MultiDiGraph:
    path = pl.Path(graph_path)
    if not path.exists():
        raise FileNotFoundError(
            f"Graph file not found: {path}. "
            "Run Stage 1 first: python -m src.jump_network_builder <yaml_config>"
        )

    if path.suffix.lower() == ".adjlist":
        raise ValueError(
            "Stage 2 requires jump_structure.node_link.json. "
            "jump_structure.adjlist does not preserve the node/edge metadata required for matrix construction."
        )

    if path.suffix.lower() == ".json":
        with path.open("r", encoding="utf-8") as fh:
            data = json.load(fh)
        graph = _node_link_graph(data)
    else:
        raise ValueError(
            f"Unsupported graph format '{path.suffix}'. "
            "Use jump_structure.node_link.json from Stage 1."
        )

    if not isinstance(graph, nx.MultiDiGraph):
        graph = nx.MultiDiGraph(graph)

    required_edge_attrs = {"h", "nn_dist", "image", "glob_jump_cnt"}
    for start, end, key, attrs in graph.edges(keys=True, data=True):
        missing = required_edge_attrs.difference(attrs.keys())
        if missing:
            raise ValueError(
                f"Serialized graph edge ({start}, {end}, {key}) misses attributes: {sorted(missing)}"
            )

    return graph


def _extract_node_labels(
    graph: nx.MultiDiGraph,
    symmetry_ctx: JumpNetworkSymmetryContext,
) -> list[str]:
    labels_from_graph = {
        int(node): str(attrs["lbl"])
        for node, attrs in graph.nodes(data=True)
        if "lbl" in attrs
    }
    fallback_labels = symmetry_ctx.get_node_lbls()
    nsites = symmetry_ctx.get_refined_jump_structure().num_sites
    labels: list[str] = []
    for idx in range(nsites):
        labels.append(labels_from_graph.get(idx, fallback_labels[idx]))
    return labels


def _number_of_jumps_from_graph(graph: nx.MultiDiGraph) -> int:
    if graph.number_of_edges() == 0:
        return 0
    max_jump_index = max(int(attrs["glob_jump_cnt"]) for *_edge, attrs in graph.edges(data=True))
    return max_jump_index + 1


def build_full_transition_matrix(jump_net: Any, subst_w: Mapping[str, Any]) -> sp.Matrix:
    """Return the full transition–probability matrix for a given network."""
    return jump_net.get_full_tp_matrix(dict(subst_w))


def export_axis_matrices_to_mathematica(
    result: MatrixBuildResult,
    params: Mapping[str, Any],
    output_dir: pl.Path,
) -> None:
    """Export axis-dependent reduced matrices and Stage-3 payload pickles.

    Filenames (``tp_matrix_axis_<h>_<k>_<l>.m``,
    ``tc_axis_meta_<h>_<k>_<l>.m`` and
    ``tc_axis_payload_<h>_<k>_<l>.pickle``) use the conventional-cell plane
    indices -- ``result.axis_results`` is keyed by ``hkl_conventional``,
    since ``JumpNetworkSymmetryContext`` always refines to a conventional
    structure regardless of the input POSCAR, so the conventional indices
    are always available and unambiguous.
    """
    t_symbol = result.jump_network.get_t() if hasattr(result.jump_network, "get_t") else None
    w_substitutions: dict[sp.Symbol, sp.Expr] = {}
    if hasattr(result.jump_network, "jump_probabilities"):
        jump_probabilities = result.jump_network.jump_probabilities
        if hasattr(jump_probabilities, "get_jump_probabilities"):
            for idx, expr in jump_probabilities.get_jump_probabilities().items():
                w_substitutions[sp.Symbol(f"W{int(idx)}", real=True)] = sp.nsimplify(expr)

    # Build jump label map from graph edges
    jump_label_map: dict[int, str] = {}
    jump_label_display_map: dict[str, str] = {}
    if hasattr(result.jump_network, "graph") and hasattr(result.jump_network.graph, "graph"):
        nx_graph = result.jump_network.graph.graph
        alias_map = {}
        if hasattr(result.symmetry_ctx, "get_symbol_alias_map"):
            alias_map = dict(result.symmetry_ctx.get_symbol_alias_map())
        for start, end, _key, attrs in nx_graph.edges(keys=True, data=True):
            glob_idx = int(attrs["glob_jump_cnt"])
            # Try both 'label' and 'lbl' attributes for node labels
            start_label = nx_graph.nodes[start].get("label") or nx_graph.nodes[start].get("lbl", str(start))
            end_label = nx_graph.nodes[end].get("label") or nx_graph.nodes[end].get("lbl", str(end))
            label = f"{start_label}->{end_label}"
            jump_label_map[glob_idx] = label
            start_disp = alias_map.get(str(start_label), str(start_label))
            end_disp = alias_map.get(str(end_label), str(end_label))
            jump_label_display_map.setdefault(str(start_label), str(start_disp))
            jump_label_display_map.setdefault(str(end_label), str(end_disp))

    for hkl_conventional, ax_res in result.axis_results.items():
        hkl_label = "_".join(str(v) for v in hkl_conventional)
        fname_wo_extension = output_dir / f"tp_matrix_axis_{hkl_label}"
        export_to_mathematica(ax_res.reordered_tp_matrix, str(fname_wo_extension))
        _write_axis_meta_file(
            output_dir / f"tc_axis_meta_{hkl_label}.m",
            hkl_conventional=hkl_conventional,
            ax_res=ax_res,
            w_substitutions=w_substitutions,
            jump_label_map=jump_label_map,
        )

        payload = {
            "contract": "transport-coefficient-axis-payload",
            "contract_version": 1,
            "axis_hkl_structure": tuple(int(v) for v in ax_res.hkl_structure),
            "axis_hkl_conventional": tuple(int(v) for v in hkl_conventional),
            "transition_matrix": ax_res.reordered_tp_matrix,
            "labels": list(ax_res.labels),
            "min_jump_length": ax_res.min_jump_length,
            "interplanar_distance_squared": ax_res.interplanar_distance_squared,
            "t_symbol": t_symbol,
            "w_substitutions": w_substitutions,
            "jump_label_map": jump_label_map,
            "jump_label_display_map": jump_label_display_map,
        }
        payload_file = output_dir / f"tc_axis_payload_{hkl_label}.pickle"
        with payload_file.open("wb") as fh:
            pickle.dump(payload, fh)


def _write_axis_meta_file(
    meta_file: pl.Path,
    *,
    hkl_conventional: Sequence[Any],
    ax_res: Any,
    w_substitutions: Mapping[sp.Symbol, sp.Expr],
    jump_label_map: Mapping[int, str],
) -> None:
    """Write the scale factors of one axis as a Mathematica association.

    The file ``tc_axis_meta_<h>_<k>_<l>.m`` accompanies
    ``tp_matrix_axis_<h>_<k>_<l>.m``.  It holds the quantities that
    ``scripts/solve_dc.wls`` needs in addition to the matrix: the squared
    interplanar distance and the minimal jump length.
    """
    from sympy.printing.mathematica import mathematica_code

    def code(expr: Any) -> str:
        return mathematica_code(sp.nsimplify(expr))

    min_jump_length = 1 if ax_res.min_jump_length is None else ax_res.min_jump_length
    axis = ", ".join(str(int(v)) for v in hkl_conventional)
    labels = ", ".join(f'"{label}"' for label in ax_res.labels)
    w_rules = ", ".join(
        f"{code(lhs)} -> {code(rhs)}" for lhs, rhs in w_substitutions.items() if lhs != rhs
    )
    jump_labels = ", ".join(
        f'"W{int(idx)}" -> "{label}"' for idx, label in sorted(jump_label_map.items())
    )
    lines = [
        "<|",
        f'  "axis" -> {{{axis}}},',
        f'  "interplanarDistanceSquared" -> {code(ax_res.interplanar_distance_squared)},',
        f'  "minJumpLength" -> {code(min_jump_length)},',
        f'  "labels" -> {{{labels}}},',
        f'  "wSubstitutions" -> {{{w_rules}}},',
        f'  "jumpLabels" -> {{{jump_labels}}}',
        "|>",
    ]
    meta_file.write_text("\n".join(lines) + "\n", encoding="utf-8")


MERGE_MODES = ("check", "force")


def _apply_merge_policy(
    equiv_sites: list[int],
    full_tp: sp.Matrix,
    hkl_structure: Sequence[Any],
    k_symbols: Sequence[sp.Symbol],
    t_symbol: sp.Symbol,
    mode: str,
    axis_label: str,
) -> list[int]:
    """Decide whether the proposed equivalence classes may be merged.

    Merging with ``A' = R A R^+`` is exact only if the partition is lumpable
    (see :mod:`src.matrix_builder.lumpability`).  ``mode`` is the value of the
    ``merge_classes`` input parameter:

    * ``"check"`` (default): if the partition is not lumpable, classes are split
      to the coarsest lumpable refinement and the rate relations that would
      allow the merge are printed.
    * ``"force"``: the proposed classes are merged anyway.  The result is valid
      only if the printed rate relations hold; apply them afterwards.
    """

    report = check_lumpability(equiv_sites, full_tp, hkl_structure, k_symbols, t_symbol)
    if report.ok:
        n_sites = len(equiv_sites)
        n_classes = len(set(equiv_sites))
        if n_classes < n_sites:
            print(
                f"\naxis {axis_label}: merging {n_sites} sites into {n_classes} classes "
                "is lumpable, the reduction is exact."
            )
        return equiv_sites

    relations = suggest_rate_relations(report)
    relation_text = (
        ", ".join(f"{r.lhs} = {r.rhs}" for r in relations)
        if relations
        else "(could not be expressed as linear relations between rates)"
    )
    print(f"\nWARNING: axis {axis_label}: proposed merging of equivalent sites is not lumpable.")
    print(f"  {report.summary()}")
    print(f"  The merge is exact only if: {relation_text}")
    if mode == "force":
        print("  merge_classes=force: classes are merged anyway; impose the relations above in the result.")
        return equiv_sites

    refined = refine_to_lumpable(equiv_sites, full_tp, hkl_structure, k_symbols, t_symbol)
    print(
        f"  merge_classes=check: classes split ({len(set(equiv_sites))} -> {len(set(refined))}); "
        "set merge_classes: force in the input YAML to merge them anyway."
    )
    return refined


def build_matrices(
    *,
    graph_path: str | pl.Path | None = None,
    jump_sites_poscar: str | pl.Path | None = None,
    graph: nx.MultiDiGraph | None = None,
    symmetry_ctx: JumpNetworkSymmetryContext | None = None,
    crystal_matrix: Sequence[str] | None = None,
    params: Mapping[str, Any],
) -> MatrixBuildResult:
    """Build full and axis-reduced matrices from files or in-memory objects."""

    if graph is None:
        if graph_path is None:
            raise ValueError("graph_path must be provided when graph is not supplied")
        graph = _load_jump_graph(graph_path)

    if symmetry_ctx is None:
        if jump_sites_poscar is None:
            raise ValueError("jump_sites_poscar must be provided when symmetry_ctx is not supplied")
        jump_sites_poscar_path = pl.Path(jump_sites_poscar).resolve()
        matrix_species: Sequence[str] | None = crystal_matrix
        if matrix_species is None and jump_sites_poscar_path.name != "jump_sites.poscar":
            matrix_species = params.get("matrix")
        symmetry_ctx = JumpNetworkSymmetryContext(
            poscar_file=str(jump_sites_poscar_path),
            crystal_matrix=matrix_species,
        )

    jump_sites = symmetry_ctx.get_refined_jump_structure()
    node_labels = _extract_node_labels(graph, symmetry_ctx)
    # Obtain the symbolic lattice‑vector matrix (conventional basis).
    lat_const = symmetry_ctx.get_symbolic_lattice_vectors()
    number_of_jumps = _number_of_jumps_from_graph(graph)

    jump_net = MatrixStageJumpNetwork(
        graph=graph,
        jump_sites=jump_sites,
        node_labels=node_labels,
        lat_const=lat_const,
        number_of_jumps=number_of_jumps,
    )

    full_tp = build_full_transition_matrix(
        jump_net,
        params.get("tweak_transition_rates", {}),
    )

    axis_analyzer = JumpAxisEquivalenceAnalyzer(
        symmetry_ctx, label_source=str(params.get("site_equivalence", "axis_orbit")).lower()
    )
    merge_mode = str(params.get("merge_classes", "check")).lower()
    if merge_mode not in MERGE_MODES:
        raise ValueError(f"merge_classes must be one of {MERGE_MODES}, got {merge_mode!r}")
    axis_results: Dict[Tuple[int, int, int], AxisMatrixResult] = {}

    axes = params.get("axes", []) or []
    for axis in axes:
        if (
            isinstance(axis, (list, tuple))
            and len(axis) == 2
            and isinstance(axis[0], (list, tuple))
        ):
            raise ValueError(
                "Axis definitions must be plain [h, k, l] vectors; "
                "the screw parameter format [[h, k, l], n] is not supported."
            )
        hkl = _coerce_integer_hkl(axis)
        hkl_conventional = tuple(int(v) for v in hkl)
        hkl_structure = symmetry_ctx.convert_hkl(hkl)

        # Determine equivalence classes of sites for the current transport axis.
        equiv_sites = axis_analyzer.find_equiv_classes_for_axis(hkl_structure)
        site_projections = axis_analyzer.last_site_projections

        k = jump_net.get_k()
        t = jump_net.get_t()
        equiv_sites = _apply_merge_policy(
            equiv_sites, full_tp, hkl_structure, k, t, merge_mode, str(hkl_conventional)
        )

        hkl_vec = sp.Matrix(hkl_structure)
        interplanar_distance_squared = jump_net.compute_squared_interplanar_distance(hkl_vec)
        print_equiv_coords(
            equiv_sites, hkl_vec, interplanar_distance_squared, jump_sites, node_labels,
            hkl_conventional=hkl_conventional,
        )

        ncols, nrows, reduced_tp_cos = find_equiv_and_reduce_matrix(
            equiv_sites,
            full_tp,
            hkl_structure,
            k,
            t,
        )

        # Apply minimal-argument normalization to obtain a transition matrix
        # expressed in terms of the smallest trigonometric argument scale.
        #
        # When this succeeds, ``find_min_jump_length`` returns a matrix where
        # the wave-vector symbol has effectively been rescaled via the
        # substitution constructed in ``find_min_arg``. The returned
        # ``min_jump_length`` then encodes the corresponding jump-length scale,
        # and Stage-3 is allowed to multiply the diffusion coefficient by
        # ``min_jump_length**2``.
        #
        # If normalization fails for any reason, we fall back to the
        # unnormalized ``reduced_tp_cos`` and do not report ``min_jump_length``;
        # in that case, Stage-3 must *not* apply the additional scaling.
        try:
            normalized_tp, min_jump_length = find_min_jump_length(
                ncols,
                nrows,
                reduced_tp_cos,
                t,
                axis_label=str(hkl_structure),
            )
        except Exception:
            normalized_tp = reduced_tp_cos
            min_jump_length = None

        labels, pretty_site_coords = prepare_label_data_for_reordered_matrix(
            equiv_sites=equiv_sites,
            node_labels=jump_net.node_labels,
            frac_coords=jump_net.jump_sites.frac_coords,
            site_projections=site_projections,
        )
        reordered_tp, labels = prepare_reordered_tp_matrix(
            normalized_tp,
            labels,
            nrows=nrows,
            ncols=ncols,
            verbose=False,
            pretty_site_coords=pretty_site_coords,
        )

        axis_results[hkl_conventional] = AxisMatrixResult(
            hkl_structure=hkl_structure,
            hkl_conventional=hkl_conventional,
            equiv_sites=equiv_sites,
            reordered_tp_matrix=reordered_tp,
            labels=labels,
            min_jump_length=min_jump_length,
            interplanar_distance_squared=interplanar_distance_squared,
            diffusion_coefficient=None,
        )

    return MatrixBuildResult(
        jump_network=jump_net,
        symmetry_ctx=symmetry_ctx,
        labels=tuple(jump_net.node_labels),
        full_tp_matrix=full_tp,
        axis_results=axis_results,
    )


def print_equiv_coords(equiv_sites: list[int], hkl_structure: tuple[int, ...] | tuple[int, int, int],
                    hkl_norm_sq: sp.Expr, jump_sites: Structure, node_labels: list[str],
                    hkl_conventional: tuple[int, int, int] | None = None):
    # ------------------------------------------------------------------
    # Print direct (fractional) coordinates of equivalent sites per class,
    # so the grouping of sites that contribute identically to transport
    # along a given axis can be inspected. The output is concise: for each
    # class we list the site indices and their fractional coordinates.
    # ------------------------------------------------------------------
    from collections import defaultdict
    class_map: dict[int, list[int]] = defaultdict(list)
    for site_idx, class_id in enumerate(equiv_sites):
        class_map[class_id].append(site_idx)

    hkl_vec = sp.Matrix(hkl_structure)

    # Display-only rounding: floating-point round-off (e.g. from
    # convert_hkl's basis change) can leave a projection that is physically
    # zero as a tiny nonzero float, which nsimplify then renders as a long,
    # meaningless rational below. This only affects what gets printed; the
    # underlying `proj`/`equiv_sites` values used elsewhere are untouched.
    PROJ_DISPLAY_EPSILON = 1e-5

    structure_label = ",".join(str(hkl) for hkl in hkl_structure)
    if hkl_conventional is not None and tuple(hkl_conventional) != tuple(hkl_structure):
        conventional_label = ",".join(str(hkl) for hkl in hkl_conventional)
        print(
            f"\n--- Equivalent sites for axis ({structure_label}) "
            f"[input as conventional ({conventional_label}); the POSCAR is a "
            f"primitive cell, so plane indices were converted to the "
            f"primitive-cell basis for this calculation] ---"
        )
    else:
        print(f"\n--- Equivalent sites for axis ({structure_label}) ---")
    for cls_id, sites in sorted(class_map.items()):
        lines = []
        for site_idx in sites:
            label = node_labels[site_idx]
            # Wrap to the canonical [0, 1) representative before projecting,
            # matching JumpAxisEquivalenceAnalyzer.find_equiv_classes_for_axis:
            # a symmetrized structure can leave a site's representative
            # sitting exactly on the periodic boundary (a frac component of
            # 1.0 rather than 0.0), which would otherwise print a projection
            # a whole interplanar spacing away from other members of the very
            # same class this site was just grouped into above.
            fcoords = jump_sites.frac_coords[site_idx] % 1.0
            fvec = sp.Matrix(list(fcoords))
            proj = (hkl_vec.dot(fvec))
            if abs(float(proj)) < PROJ_DISPLAY_EPSILON:
                frac = sp.Integer(0)
            else:
                # Normalised fraction; guard against zero-length hkl (should not occur)
                frac = sp.nsimplify(proj, constants=[sp.sqrt(3), sp.sqrt(2)], tolerance=1e-6) / sp.sqrt(
                    hkl_norm_sq) if hkl_norm_sq != 0 else sp.Integer(0)
            lines.append(
                f"\n     {label} idx={site_idx}"
                f" coord={[sp.nsimplify(fc) for fc in fcoords.tolist()]} proj={frac}"
            )
        class_str = ", ".join(lines)
        print(f"  Class {cls_id}: {class_str}")


def build_matrices_from_yaml(
    yaml_file: str | pl.Path,
    params: Mapping[str, Any],
    output_dir: pl.Path,
) -> MatrixBuildResult:
    """Build matrices using Stage-1 artifacts in ``output_dir``."""

    graph_path = output_dir / "jump_structure.node_link.json"
    jump_sites_poscar = output_dir / "jump_sites.poscar"
    if not graph_path.exists():
        raise FileNotFoundError(
            f"Missing Stage-1 graph file: {graph_path}. "
            "Run jump_network_builder first."
        )
    if not jump_sites_poscar.exists():
        raise FileNotFoundError(
            f"Missing Stage-1 jump-site structure: {jump_sites_poscar}. "
            "Run jump_network_builder first."
        )

    return build_matrices(
        graph_path=graph_path,
        jump_sites_poscar=jump_sites_poscar,
        params=params,
    )


def print_axis_build_summary(result: MatrixBuildResult) -> None:
    """Print a compact summary for each transport direction."""
    if not result.axis_results:
        print("No axes specified in YAML; only full matrix was built.")
        return

    print("Built reduced matrices for axes:")
    for hkl_conventional, ax_res in sorted(result.axis_results.items()):
        hkl_structure = getattr(ax_res, "hkl_structure", None)
        effective_distance_squared = None
        if ax_res.interplanar_distance_squared is not None:
            effective_distance_squared = ax_res.interplanar_distance_squared
            if ax_res.min_jump_length is not None:
                effective_distance_squared = sp.simplify(
                    effective_distance_squared * sp.nsimplify(ax_res.min_jump_length) ** 2
                )
        print(
            f"  hkl={hkl_conventional}"
            + (
                f" (structure-cell basis {hkl_structure})"
                if hkl_structure is not None and tuple(hkl_structure) != tuple(hkl_conventional)
                else ""
            )
            + f": size={ax_res.reordered_tp_matrix.shape}, "
            f"min_jump_length={ax_res.min_jump_length}, "
            f"interplanar_distance_squared={ax_res.interplanar_distance_squared}, "
            f"effective_transport_distance_squared={effective_distance_squared}"
        )


def main(argv: Optional[Sequence[str]] = None) -> int:
    """CLI entry point for the matrix-builder stage."""

    configure_cli_logging()
    _yaml_file, params, output_dir = init_from_yaml()

    graph_path = output_dir / pl.Path("jump_structure.node_link.json")
    jump_sites_poscar = output_dir / pl.Path("jump_sites.poscar")

    result = build_matrices(
        graph_path=graph_path,
        jump_sites_poscar=jump_sites_poscar,
        params=params,
    )

    print("Full transition matrix shape:", result.full_tp_matrix.shape)
    print_axis_build_summary(result)
    export_axis_matrices_to_mathematica(result, params, output_dir)
    return 0


if __name__ == "__main__":
    main()
