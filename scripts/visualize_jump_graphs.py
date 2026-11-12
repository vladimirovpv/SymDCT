#!/usr/bin/env python3
"""Render full and axis-quotient jump graphs from Stage-1/Stage-2 artifacts.

For each requested transport direction [hkl] this produces an axis-oriented
full jump graph (every physical periodic jump, carrying the signed transport
displacement Delta_s along that axis) and an explicit quotient graph (one
node per axis-dependent equivalence class, self-loops preserved), each as
PNG (quick look) and GraphML (the authoritative output, for polishing in
yEd).

The quotient graph is built explicitly -- NOT via ``nx.quotient_graph()``.
That function can never represent a self-loop: its block-pair enumeration is
built from ``permutations``/``combinations``, which structurally exclude a
block paired with itself, for any partition or edge_relation. That would
silently drop same-type jumps (e.g. A-A) whenever a class is a singleton,
even when they matter physically for that axis's diffusion coefficient (a
self-loop C->C does not mean Delta_s == 0 -- an A-A jump can map A1->A1 in
the quotient while the actual periodic jump has Delta_s != 0).

Module structure (each stage is a separate, independently testable step):

    axis-oriented geometry (xi_hat, d_hkl, oriented unit cell)
        -> full graph + physical jump geometry (PhysicalJump per edge)
        -> axis-specific equivalence classes (JumpAxisEquivalenceAnalyzer)
        -> explicit quotient construction (build_quotient_graph_explicit)
        -> graph labels / display attributes
        -> GraphML export
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pymatgen as pmg

logging.getLogger("matplotlib").setLevel(logging.WARNING)

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.common_utils.jump_vector_codec import decode_jump_vector
from src.common_utils.yaml_io import load_data
from src.jump_network_builder.geometry import _node_link_graph
from src.matrix_builder.matrix_reduction_helpers import (
    representative_site_indices_from_equiv,
)
from src.symmetry_context.axis_equivalence import JumpAxisEquivalenceAnalyzer
from src.symmetry_context.jump_symmetry_context import JumpNetworkSymmetryContext
from src.symmetry_context.node_labeler import build_non_equivalent_site_label_data

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Path / artifact resolution and Stage-1 graph loading (unchanged)
# ---------------------------------------------------------------------------


def _resolve_paths(path_arg: Path) -> tuple[Path, Path]:
    if path_arg.suffix.lower() in {".yaml", ".yml"}:
        yaml_path = path_arg.resolve()
        output_dir = yaml_path.parent.parent / "output"
        return yaml_path, output_dir

    output_dir = path_arg.resolve()
    input_dir = output_dir.parent / "input"
    candidates = sorted(input_dir.glob("*.yaml")) + sorted(input_dir.glob("*.yml"))
    if not candidates:
        raise FileNotFoundError(f"No YAML input found in: {input_dir}")
    return candidates[0], output_dir


def _load_graph(node_link_path: Path) -> nx.MultiDiGraph:
    if not node_link_path.exists():
        raise FileNotFoundError(
            f"Missing graph artifact: {node_link_path}. Run jump_network_builder first."
        )
    with node_link_path.open("r", encoding="utf-8") as fh:
        data = json.load(fh)
    # Reuse Stage 1's own loader so the edge-record key ("edges", not
    # NetworkX's "links" default) matches what jump_network_builder wrote.
    graph = _node_link_graph(data)
    return graph if isinstance(graph, nx.MultiDiGraph) else nx.MultiDiGraph(graph)


def _extract_node_labels(graph: nx.MultiDiGraph) -> dict[int, str]:
    return {
        int(node): str(attrs.get("lbl", node))
        for node, attrs in graph.nodes(data=True)
    }


def _safe_axes(raw_axes: Iterable[Any]) -> list[tuple[int, int, int]]:
    axes: list[tuple[int, int, int]] = []
    for axis in raw_axes:
        if isinstance(axis, (list, tuple)) and len(axis) == 3:
            axes.append((int(axis[0]), int(axis[1]), int(axis[2])))
    return axes


# ---------------------------------------------------------------------------
# Topology-only overview of the full network (all axes combined) -- kept for
# a quick "what jumps exist at all" picture, independent of any one
# transport direction. Unrelated to the axis-specific graphs below.
# ---------------------------------------------------------------------------


def _collapse_multidigraph(graph: nx.MultiDiGraph) -> nx.DiGraph:
    collapsed = nx.DiGraph()
    for node in graph.nodes:
        collapsed.add_node(int(node))
    for start, end, attrs in graph.edges(data=True):
        u = int(start)
        v = int(end)
        nn_dist = attrs.get("nn_dist")
        jump_idx = attrs.get("glob_jump_cnt")
        if not collapsed.has_edge(u, v):
            collapsed.add_edge(u, v, count=0, lengths=[], jump_indices=[])
        collapsed[u][v]["count"] += 1
        if nn_dist is not None:
            collapsed[u][v]["lengths"].append(float(nn_dist))
        if jump_idx is not None:
            collapsed[u][v]["jump_indices"].append(int(jump_idx))
    return collapsed


def _topology_edge_label(attrs: dict[str, Any]) -> str:
    count = int(attrs.get("count", 1))
    lengths = attrs.get("lengths", [])
    if not lengths:
        return f"n={count}"
    lo = min(lengths)
    hi = max(lengths)
    if abs(lo - hi) < 1e-9:
        return f"n={count}, d={lo:.3f}"
    return f"n={count}, d=[{lo:.3f},{hi:.3f}]"


def _draw_graph(
    graph: nx.DiGraph,
    node_labels: dict[Any, str],
    title: str,
    output_file: Path,
    pos: dict[Any, tuple[float, float]] | None = None,
    edge_labels: dict[tuple[Any, Any], str] | None = None,
    figsize: tuple[float, float] = (12, 9),
) -> None:
    """Generic node-link renderer, reused for the topology overview and the
    axis-specific full/quotient graphs (the latter pass an explicit ``pos``
    so the layout reflects plane order and in-plane position rather than a
    force-directed spring layout).

    NOTE: this PNG is a quick-look picture only, not the authoritative
    output -- GraphML (written separately) carries the complete, exact
    per-jump data and is meant for manual layout/polishing in yEd. In
    particular, self-loops (e.g. a same-type A-A jump) can render as very
    thin arcs when a transport axis places every site at the same in-plane
    coordinate, since matplotlib scales the two screen axes independently
    to fill the figure; this does not affect the underlying data.
    """

    if pos is None:
        pos = nx.spring_layout(graph, seed=42)
    plt.figure(figsize=figsize)
    nx.draw_networkx_nodes(graph, pos, node_color="#b9dcff", node_size=900)
    nx.draw_networkx_edges(
        graph, pos, arrowstyle="->", arrowsize=14, width=1.5, connectionstyle="arc3,rad=0.08"
    )
    nx.draw_networkx_labels(graph, pos, labels=node_labels, font_size=10)

    if edge_labels is None:
        edge_labels = {(u, v): _topology_edge_label(data) for u, v, data in graph.edges(data=True)}

    # Self-loops and "through" edges need different label treatment: a
    # through edge's label can be rotated to align with that edge's own
    # direction (rotate=True), which keeps it from sitting at an arbitrary
    # angle across other lines; a self-loop has no such direction, and
    # draw_networkx_edge_labels' automatic curve-point placement puts its
    # label right on top of the loop arc (which itself hugs the node), so
    # self-loop labels are instead placed by hand, offset above the node.
    self_loop_labels = {(u, v): text for (u, v), text in edge_labels.items() if u == v}
    through_labels = {(u, v): text for (u, v), text in edge_labels.items() if u != v}

    if through_labels:
        # connectionstyle must match draw_networkx_edges' above: this
        # function has its own independent default ('arc3', no rad), and a
        # mismatch between the two can make some edge/self-loop
        # combinations crash inside networkx's label-placement geometry.
        nx.draw_networkx_edge_labels(
            graph, pos, edge_labels=through_labels, font_size=8, rotate=True, label_pos=0.45,
            connectionstyle="arc3,rad=0.08",
        )

    if self_loop_labels:
        # Offset by a fixed on-screen distance (matplotlib "points", via the
        # current display transform) rather than a fraction of the overall
        # inter-node data span: the self-loop arc's own size is tied to
        # node_size in points, not to how far apart *other* nodes happen to
        # be, so a data-span-relative offset is either invisible (small
        # span) or wildly oversized (large span, as hcp_vac's two-plane
        # axes are) -- both observed when this used data_span.
        ax = plt.gca()
        node_radius_points = math.sqrt(900) / 2
        offset_points = node_radius_points + 12
        for (node, _node2), text in self_loop_labels.items():
            x_disp, y_disp = ax.transData.transform(pos[node])
            x_text, y_text = ax.transData.inverted().transform((x_disp, y_disp + offset_points))
            plt.text(x_text, y_text, text, fontsize=8, ha="center", va="bottom")

    plt.title(title)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(output_file, dpi=200)
    plt.close()


# ---------------------------------------------------------------------------
# 1. Axis-oriented geometry
# ---------------------------------------------------------------------------


def _extended_gcd(a: int, b: int) -> tuple[int, int, int]:
    """Return (g, x, y) with a*x + b*y = g = gcd(a, b) (g >= 0)."""
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r != 0:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    if old_r < 0:
        old_r, old_s, old_t = -old_r, -old_s, -old_t
    return old_r, old_s, old_t


def _unimodular_plane_basis(hkl: tuple[int, int, int]) -> tuple[
    tuple[int, int, int], tuple[int, int, int], tuple[int, int, int]
]:
    """Return integer vectors (w0, u0, v0) via exact-integer extended-GCD:
    h.u0 = h.v0 = 0 (u0, v0 span the (hkl) plane through the origin) and
    h.w0 = 1, with |det([w0, u0, v0])| == 1 -- i.e. (w0, u0, v0) is a
    unimodular integer basis change of the whole lattice, not just of the
    2-D sublattice. This is the "exact-integer extended-GCD construction"
    step: two elementary unimodular eliminations (first zeroing k against
    h, then l against the result) via the extended Euclidean algorithm,
    composed as a single 3x3 unimodular matrix R with (h,k,l) @ R = (1,0,0).
    """
    h, k, l = (int(v) for v in hkl)
    g_hkl = math.gcd(math.gcd(abs(h), abs(k)), abs(l))
    if g_hkl == 0:
        raise ValueError("hkl must not be (0, 0, 0)")
    h, k, l = h // g_hkl, k // g_hkl, l // g_hkl  # reduce to coprime Miller indices

    g1, x1, y1 = _extended_gcd(h, k)
    r1 = np.eye(3, dtype=object)
    if g1 == 0:
        # h == k == 0 (so l == +-1 after reduction): first elimination is trivial.
        r1[:] = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
    else:
        r1[0, 0], r1[0, 1] = x1, -k // g1
        r1[1, 0], r1[1, 1] = y1, h // g1
    step1 = np.array([h, k, l], dtype=object) @ r1  # == (g1, 0, l)
    assert int(step1[1]) == 0

    g2, x2, y2 = _extended_gcd(int(step1[0]), int(step1[2]))
    r2 = np.eye(3, dtype=object)
    if g2 == 0:
        r2[:] = [[0, 0, 1], [0, 1, 0], [1, 0, 0]]
    else:
        r2[0, 0], r2[0, 2] = x2, -int(step1[2]) // g2
        r2[2, 0], r2[2, 2] = y2, int(step1[0]) // g2

    r_full = r1 @ r2
    result = np.array([h, k, l], dtype=object) @ r_full
    if int(result[0]) != 1 or int(result[1]) != 0 or int(result[2]) != 0:
        raise AssertionError(
            f"Extended-GCD plane-basis construction failed for hkl={hkl}: "
            f"(h,k,l) @ R = {tuple(int(v) for v in result)}, expected (1, 0, 0)."
        )
    det = int(round(np.linalg.det(r_full.astype(float))))
    if abs(det) != 1:
        raise AssertionError(f"Plane-basis matrix R is not unimodular for hkl={hkl}: det={det}")

    w0 = tuple(int(v) for v in r_full[:, 0])
    u0 = tuple(int(v) for v in r_full[:, 1])
    v0 = tuple(int(v) for v in r_full[:, 2])
    return w0, u0, v0


def _gauss_reduce_2d(
    u: tuple[int, int, int], v: tuple[int, int, int], metric: np.ndarray
) -> tuple[tuple[int, int, int], tuple[int, int, int]]:
    """2-D Lagrange-Gauss lattice reduction of integer vectors (u, v) using
    the structure's actual (generally anisotropic) metric tensor -- *not*
    an isotropic/Euclidean-in-integer-coordinates assumption, which would
    be wrong whenever the lattice is not cubic (e.g. tetragonal a != c).
    Returns the shortest possible basis spanning the same 2-D sublattice.
    """
    u_arr = np.array(u, dtype=object)
    v_arr = np.array(v, dtype=object)

    def _dot(a: np.ndarray, b: np.ndarray) -> float:
        af = np.array(a, dtype=float)
        bf = np.array(b, dtype=float)
        return float(af @ metric @ bf)

    for _ in range(1000):  # generous bound; converges in a handful of steps in practice
        if _dot(v_arr, v_arr) < _dot(u_arr, u_arr):
            u_arr, v_arr = v_arr, u_arr
        mu = round(_dot(u_arr, v_arr) / _dot(v_arr, v_arr))
        if mu == 0:
            break
        u_arr = u_arr - mu * v_arr
    else:
        raise AssertionError("2-D Gauss lattice reduction did not converge")

    return tuple(int(x) for x in u_arr), tuple(int(x) for x in v_arr)


def _minimal_transverse_vector(
    w0: tuple[int, int, int],
    u: tuple[int, int, int],
    v: tuple[int, int, int],
    metric: np.ndarray,
    search_range: int = 3,
) -> tuple[int, int, int]:
    """Brute-force search, over integer shifts w0 + m*u + n*v (which all
    satisfy h.w == 1, since u and v lie in the plane h.x == 0), for the
    shortest such completion of the primitive cell."""
    w0_arr = np.array(w0, dtype=object)
    u_arr = np.array(u, dtype=object)
    v_arr = np.array(v, dtype=object)
    best: tuple[int, int, int] | None = None
    best_len2: float | None = None
    for m in range(-search_range, search_range + 1):
        for n in range(-search_range, search_range + 1):
            w = w0_arr + m * u_arr + n * v_arr
            w_float = np.array(w, dtype=float)
            len2 = float(w_float @ metric @ w_float)
            if best_len2 is None or len2 < best_len2 - 1e-9:
                best_len2 = len2
                best = tuple(int(x) for x in w)
    assert best is not None
    return best


@dataclass(frozen=True)
class OrientedCell:
    """A provably primitive oriented cell for one transport direction:
    ``u``, ``v`` (integer + Cartesian) span the (hkl) plane; ``w`` is the
    shortest transverse completion. ``|det([w, u, v])| == 1`` exactly."""

    u: tuple[int, int, int]
    v: tuple[int, int, int]
    w: tuple[int, int, int]
    u_cart: np.ndarray
    v_cart: np.ndarray
    w_cart: np.ndarray


def oriented_unit_cell_for_axis(
    structure: pmg.core.Structure, hkl: tuple[int, int, int]
) -> OrientedCell:
    """Return the oriented unit cell for transport direction ``hkl``.

    Replaces pymatgen's ``SlabGenerator``, which was found to pick a
    non-minimal -- and for some {hkl} not even confined to the plane --
    oriented unit cell for body-centered tetragonal structures (e.g. Be12Ti,
    I4/mmm, a != c), especially {211}-type planes: for hkl=(2,1,1) on that
    structure, SlabGenerator's cell had area*height = 4x the true primitive
    cell volume. Reimplemented from scratch instead of patching
    SlabGenerator's basis search:

    1. ``structure.lattice.matrix`` is already the correctly constructed
       primitive cell (from ``JumpStructureSymmetryRefiner``/spglib, not
       rebuilt here); no separate conventional-to-primitive step is needed
       in this function.
    2. Find the minimal 2-D in-plane basis for {hkl} via exact-integer
       extended-GCD construction (``_unimodular_plane_basis``), then reduce
       it with a metric-aware 2-D Gauss/Lagrange lattice reduction
       (``_gauss_reduce_2d``) using the structure's actual (generally
       anisotropic) metric tensor -- required whenever a != c.
    3. Brute-force search (``_minimal_transverse_vector``) for the minimal
       third (out-of-plane) vector completing a provably primitive cell.
    4. Validate with the area * height == volume identity (and
       |det([w, u, v])| == 1) every time -- raises AssertionError if either
       fails, rather than silently returning a wrong cell.
    """
    lattice_matrix = structure.lattice.matrix
    metric = lattice_matrix @ lattice_matrix.T  # anisotropic metric tensor, not assumed isotropic

    w0, u0, v0 = _unimodular_plane_basis(hkl)
    u, v = _gauss_reduce_2d(u0, v0, metric)
    w = _minimal_transverse_vector(w0, u, v, metric)

    det = round(np.linalg.det(np.array([w, u, v], dtype=float)))
    if abs(det) != 1:
        raise AssertionError(
            f"Oriented cell for hkl={hkl} is not primitive: |det([w,u,v])|={abs(det)} != 1."
        )

    u_cart = np.array(u, dtype=float) @ lattice_matrix
    v_cart = np.array(v, dtype=float) @ lattice_matrix
    w_cart = np.array(w, dtype=float) @ lattice_matrix

    area = float(np.linalg.norm(np.cross(u_cart, v_cart)))
    normal = np.cross(u_cart, v_cart)
    normal_hat = normal / np.linalg.norm(normal)
    height = float(abs(np.dot(w_cart, normal_hat)))
    volume = float(abs(np.linalg.det(lattice_matrix)))
    if abs(area * height - volume) > 1e-6 * max(volume, 1e-12):
        raise AssertionError(
            f"Oriented cell for hkl={hkl} fails area*height == volume: "
            f"{area:.6f} * {height:.6f} = {area * height:.6f} != {volume:.6f}."
        )

    return OrientedCell(u=u, v=v, w=w, u_cart=u_cart, v_cart=v_cart, w_cart=w_cart)


def axis_unit_vector_and_spacing(
    structure: pmg.core.Structure, hkl: tuple[int, int, int]
) -> tuple[np.ndarray, float]:
    """Cartesian unit vector xi_hat along transport direction ``hkl``, and d_hkl.

    The transport direction is the (hkl) plane normal: the reciprocal
    lattice vector ``h*b1* + k*b2* + l*b3*``, which pymatgen's ``Lattice``
    already provides in Cartesian form via ``reciprocal_lattice`` -- this is
    the physically correct way to turn Miller indices into a real-space
    Cartesian direction for a general (non-cubic) lattice, and is *not*
    something to compute by hand from the direct lattice matrix.
    """
    h, k, l = (int(v) for v in hkl)
    recip = structure.lattice.reciprocal_lattice.matrix  # rows = b1*, b2*, b3* (Cartesian)
    g_cart = np.array([h, k, l], dtype=float) @ recip
    norm = float(np.linalg.norm(g_cart))
    if norm < 1e-12:
        raise ValueError(f"Degenerate transport direction for hkl={hkl}")
    xi_hat = g_cart / norm
    d_hkl = float(structure.lattice.d_hkl((h, k, l)))
    return xi_hat, d_hkl


def axis_inplane_basis(
    structure: pmg.core.Structure, hkl: tuple[int, int, int], xi_hat: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return an in-plane orthonormal basis (e2, e3) perpendicular to xi_hat,
    used *only* for the 2-D drawing's in-plane projection -- never for
    computing Delta_s, which depends purely on xi_hat and the physical jump
    vector.

    Uses the oriented unit cell's own in-plane lattice vectors (``u``, ``v``
    -- in-plane *by construction*, not by picking the two least aligned
    with xi_hat), so the projection reflects real, minimal crystallographic
    in-plane directions. Falls back to an arbitrary deterministic orthogonal
    basis only if that construction itself raises (e.g. a genuinely
    degenerate lattice), which is unexpected in normal operation.
    """
    try:
        oriented = oriented_unit_cell_for_axis(structure, hkl)
        e2 = oriented.u_cart - np.dot(oriented.u_cart, xi_hat) * xi_hat
        e2 = e2 / np.linalg.norm(e2)
        e3_raw = oriented.v_cart - np.dot(oriented.v_cart, xi_hat) * xi_hat
        e3_raw = e3_raw - np.dot(e3_raw, e2) * e2  # Gram-Schmidt for a clean 2-D plot
        if np.linalg.norm(e3_raw) < 1e-8:
            raise ValueError("In-plane lattice vectors are collinear after orthogonalization")
        e3 = e3_raw / np.linalg.norm(e3_raw)
        return e2, e3
    except Exception as exc:  # noqa: BLE001 -- deliberate broad fallback, logged
        logger.warning(
            "Falling back to an arbitrary in-plane basis for hkl=%s "
            "(oriented-unit-cell construction failed: %s)",
            hkl, exc,
        )
        ref = np.array([1.0, 0.0, 0.0])
        if abs(np.dot(ref, xi_hat)) > 0.9:
            ref = np.array([0.0, 1.0, 0.0])
        e2 = np.cross(xi_hat, ref)
        e2 = e2 / np.linalg.norm(e2)
        e3 = np.cross(xi_hat, e2)
        e3 = e3 / np.linalg.norm(e3)
        return e2, e3


# ---------------------------------------------------------------------------
# 2. Full jump graph with physical jump geometry
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class PhysicalJump:
    """One directed, physically real periodic jump u -> v for a given axis."""

    u: int
    v: int
    glob_jump_cnt: int
    pair_key: str  # canonical "A->B" label, from the existing node labels
    rate_symbol: str  # "W<glob_jump_cnt>", the existing per-jump rate identifier
    length: float  # ordinary Cartesian jump length (nn_dist)
    image: tuple[int, int, int]  # periodic lattice translation R
    delta_s: float  # signed transport displacement (r_v - r_u + R) . xi_hat
    delta_r_cart: tuple[float, float, float]  # full Cartesian displacement r_v - r_u + R


def extract_physical_jumps(
    graph: nx.MultiDiGraph,
    structure: pmg.core.Structure,
    node_labels: dict[int, str],
    xi_hat: np.ndarray,
) -> list[PhysicalJump]:
    """Enumerate every physical periodic jump with its signed transport
    displacement along ``xi_hat``.

    ``h`` (the encoded jump vector already on every edge) is
    ``r_v - r_u + R`` in *fractional* coordinates -- pymatgen's neighbor
    list already folds the periodic image ``R`` into the neighbor's
    fractional coordinates when the edge was built (geometry.py), so no
    separate reconstruction from the ``image`` attribute is needed here;
    ``image`` is retained purely for provenance/export.
    """
    jumps: list[PhysicalJump] = []
    for u, v, attrs in graph.edges(data=True):
        u = int(u)
        v = int(v)
        edge_id = f"({u},{v},glob_jump_cnt={attrs.get('glob_jump_cnt')})"
        h_frac = [float(c) for c in decode_jump_vector(attrs["h"], edge_id=edge_id)]
        delta_r_cart = np.asarray(structure.lattice.get_cartesian_coords(h_frac), dtype=float)
        delta_s = float(np.dot(delta_r_cart, xi_hat))

        length = float(attrs.get("nn_dist", np.linalg.norm(delta_r_cart)))
        recon_len = float(np.linalg.norm(delta_r_cart))
        if abs(recon_len - length) > 1e-3:
            raise AssertionError(
                f"Jump {edge_id}: stored nn_dist={length:.6f} does not match "
                f"the Cartesian norm of the encoded jump vector 'h' ({recon_len:.6f}); "
                "the jump-vector/nn_dist edge data is inconsistent."
            )

        image = tuple(int(c) for c in attrs.get("image", (0, 0, 0)))
        glob_jump_cnt = int(attrs["glob_jump_cnt"])
        pair_key = f"{node_labels[u]}->{node_labels[v]}"
        jumps.append(
            PhysicalJump(
                u=u,
                v=v,
                glob_jump_cnt=glob_jump_cnt,
                pair_key=pair_key,
                rate_symbol=f"W{glob_jump_cnt}",
                length=length,
                image=image,
                delta_s=delta_s,
                delta_r_cart=(float(delta_r_cart[0]), float(delta_r_cart[1]), float(delta_r_cart[2])),
            )
        )
    return jumps


# ---------------------------------------------------------------------------
# 3. Axis-specific equivalence classes and site labels
# ---------------------------------------------------------------------------


def axis_site_labels(
    structure: pmg.core.Structure,
    analyzer: JumpAxisEquivalenceAnalyzer,
    hkl_struct: tuple[int, int, int],
    node_labels: dict[int, str],
) -> tuple[list[int], dict[int, str], list[float], list[int]]:
    """Return (equiv_sites, site_label, projections, plane_of_site) for this axis.

    Reuses ``JumpAxisEquivalenceAnalyzer`` for the equivalence classes and
    ``build_non_equivalent_site_label_data`` for the site-*character* part
    of the "Wyckoff/site character + plane number" label convention already
    used throughout the project (e.g. A0, B1) -- but that function numbers
    same-character sites by rank *within their own character* (A's own 0,
    1, 2, ...; B's own 0, 1, 2, ... independently), which is exactly right
    for the matrix-table printouts it was built for, and exactly wrong here:
    on this diagram, the trailing number must mean "this plane", the same
    quantity the x-axis encodes, or two different characters at two
    different planes can end up both labelled "0" (e.g. bcc_int's O0 and T0
    on axis (1,1,1): O0 is the lowest-projection O site, T0 the lowest-
    projection T site, but they are *not* at the same projection -- nothing
    here implies they should be). So the character comes from the reused
    label, but its trailing number is replaced with the site's actual
    global plane index (``group_into_planes``, shared across all
    characters), which is also returned so the caller doesn't need to
    recompute it separately -- guaranteeing the label and the x-position it
    is plotted at always agree.
    """
    equiv_sites = analyzer.find_equiv_classes_for_axis(hkl_struct)
    projections = list(analyzer.last_site_projections)
    plane_of_site = group_into_planes(projections, tol=analyzer.dist_tol)

    node_label_list = [node_labels[i] for i in range(len(equiv_sites))]
    rep_indices = representative_site_indices_from_equiv(equiv_sites)
    rep_coords = [structure.frac_coords[i] for i in rep_indices]
    rep_projections = [projections[i] for i in rep_indices]

    label_data = build_non_equivalent_site_label_data(
        node_labels=node_label_list,
        representative_site_indices=rep_indices,
        representative_site_coords=rep_coords,
        representative_projections=rep_projections,
    )
    class_id_by_rep = [equiv_sites[i] for i in rep_indices]
    class_character = {
        cls: label.rstrip("0123456789") for cls, label in zip(class_id_by_rep, label_data.labels)
    }

    site_label = {
        i: f"{class_character[equiv_sites[i]]}{plane_of_site[i]}" for i in range(len(equiv_sites))
    }
    return equiv_sites, site_label, projections, plane_of_site


def group_into_planes(projections: Sequence[float], tol: float) -> list[int]:
    """Assign each site a sequential plane index (0-based), ordered by
    increasing transport-axis projection.

    This is a purely geometric grouping by projection value alone -- unlike
    the axis equivalence classes (which additionally require equal Wyckoff
    label), several different site labels can and do share the same plane.
    Used only for the 2-D layout's left-to-right plane ordering, never for
    the site *label* text (which keeps the existing per-label convention).
    """
    order = sorted(range(len(projections)), key=lambda i: projections[i])
    plane_of_site = [-1] * len(projections)
    current_plane = -1
    reference = None
    for idx in order:
        proj = projections[idx]
        if reference is None or abs(proj - reference) > tol:
            current_plane += 1
            reference = proj
        plane_of_site[idx] = current_plane
    return plane_of_site


# ---------------------------------------------------------------------------
# 4. Explicit quotient construction -- NOT nx.quotient_graph()
# ---------------------------------------------------------------------------


@dataclass
class QuotientEdge:
    """One axis-dependent quotient edge C_u -> C_v, aggregating every
    physical jump mapped to that class pair without discarding any of
    them (retained in full in ``jumps``)."""

    c_u: int
    c_v: int
    jumps: list[PhysicalJump] = field(default_factory=list)

    @property
    def count(self) -> int:
        return len(self.jumps)


def build_quotient_graph_explicit(
    physical_jumps: Sequence[PhysicalJump],
    equiv_sites: Sequence[int],
) -> dict[tuple[int, int], QuotientEdge]:
    """Explicit axis-dependent quotient construction.

    1. full_node -> equivalence_class is ``equiv_sites`` (already computed).
    2. One quotient node per equivalence class (implicit: any class id that
       appears in ``equiv_sites``).
    3.-4. For every physical directed edge u->v, compute C_u=class(u),
       C_v=class(v) and add the jump to quotient edge C_u->C_v.
    5. Several physical jumps sharing (C_u, C_v) aggregate into the same
       quotient edge, keeping every individual jump (QuotientEdge.jumps).
    6. C_u == C_v produces a genuine quotient self-loop -- never discarded
       (unlike ``nx.quotient_graph``, whose block-pair enumeration
       structurally excludes a block paired with itself).
    """
    edges: dict[tuple[int, int], QuotientEdge] = {}
    for jump in physical_jumps:
        c_u = int(equiv_sites[jump.u])
        c_v = int(equiv_sites[jump.v])
        key = (c_u, c_v)
        if key not in edges:
            edges[key] = QuotientEdge(c_u=c_u, c_v=c_v)
        edges[key].jumps.append(jump)
    return edges


def validate_quotient_graph(
    physical_jumps: Sequence[PhysicalJump],
    equiv_sites: Sequence[int],
    quotient_edges: dict[tuple[int, int], QuotientEdge],
    n_sites: int,
) -> None:
    """Explicit consistency checks for one axis's quotient graph. Fails
    loudly (AssertionError with a specific message) on any violation."""

    # (1) Every full-graph node belongs to exactly one equivalence class.
    assert len(equiv_sites) == n_sites, (
        f"equiv_sites has {len(equiv_sites)} entries, expected {n_sites} (one per site)."
    )
    for i, cls in enumerate(equiv_sites):
        assert isinstance(cls, (int, np.integer)) and cls >= 0, (
            f"Site {i} has no valid equivalence class assigned: {cls!r}."
        )

    # (2) & (5) Every full-graph edge assigned to exactly one quotient edge;
    # total physical-jump multiplicity is conserved.
    total_full = len(physical_jumps)
    total_quotient = sum(qe.count for qe in quotient_edges.values())
    assert total_quotient == total_full, (
        f"Quotient graph represents {total_quotient} physical jumps, "
        f"but the full graph has {total_full}; jumps were lost or duplicated."
    )
    # Identity, not (u, v, glob_jump_cnt): several geometrically distinct
    # edges (different image/h) can share one glob_jump_cnt -- JumpRegister
    # groups jumps by (species pair, length), not by individual geometry --
    # so that tuple is not a unique key for a single physical-jump instance.
    def _describe(jump: PhysicalJump) -> str:
        return f"(u={jump.u}, v={jump.v}, glob_jump_cnt={jump.glob_jump_cnt}, image={jump.image})"

    seen_jump_ids: set[int] = set()
    for qe in quotient_edges.values():
        for jump in qe.jumps:
            assert id(jump) not in seen_jump_ids, (
                f"Physical jump {_describe(jump)} assigned to more than one quotient edge."
            )
            seen_jump_ids.add(id(jump))
    for jump in physical_jumps:
        assert id(jump) in seen_jump_ids, (
            f"Physical jump {_describe(jump)} from the full graph is missing from the quotient graph."
        )

    # (3) & (4) No same-class edge lost; every required self-loop exists.
    same_class_pairs = {
        (int(equiv_sites[j.u]), int(equiv_sites[j.v]))
        for j in physical_jumps
        if equiv_sites[j.u] == equiv_sites[j.v]
    }
    for pair in same_class_pairs:
        assert pair in quotient_edges, (
            f"Same-class jumps map to quotient self-loop {pair}, but no such "
            f"quotient edge exists -- a self-loop was dropped."
        )

    # (6) & (7) Every physical jump has a well-defined signed Delta_s,
    # computed from the periodic jump geometry, not from quotient-node
    # positions -- structurally guaranteed, since Delta_s is computed by
    # extract_physical_jumps() before any quotient graph exists.
    for jump in physical_jumps:
        assert jump.delta_s is not None and np.isfinite(jump.delta_s), (
            f"Physical jump {(jump.u, jump.v, jump.glob_jump_cnt)} has no valid Delta_s: {jump.delta_s!r}."
        )


def validate_axis_specificity(seen_axes: list[dict[str, Any]]) -> None:
    """Cross-axis check (8): distinct requested transport directions must
    never collapse onto the same transport geometry, even when their
    quotient *topology* happens to coincide (e.g. [0,0,1] vs [-2,1,1])."""
    for i in range(len(seen_axes)):
        for j in range(i + 1, len(seen_axes)):
            a, b = seen_axes[i], seen_axes[j]
            if a["hkl_conventional"] == b["hkl_conventional"]:
                continue
            same_xi = np.allclose(a["xi_hat"], b["xi_hat"], atol=1e-9)
            same_d = abs(a["d_hkl"] - b["d_hkl"]) < 1e-9
            assert not (same_xi and same_d), (
                f"Distinct transport directions {a['hkl_conventional']} and "
                f"{b['hkl_conventional']} produced identical transport geometry "
                f"(xi_hat, d_hkl) -- axis-specific data is not being kept distinct."
            )


# ---------------------------------------------------------------------------
# 5. Display labels / GraphML serialization-safe attribute helpers
# ---------------------------------------------------------------------------


def _quotient_node_label(members: Sequence[str]) -> str:
    return "/".join(sorted(set(members)))


def _jump_edge_display_label(jumps: Sequence[PhysicalJump]) -> str:
    """Full edge label for GraphML, e.g. '3 W_AB\\nDelta s = [+2.212, +2.212, -2.212]'
    -- every individual contribution listed, never averaged away."""
    if not jumps:
        return ""
    rate = jumps[0].rate_symbol if len(jumps) == 1 else jumps[0].pair_key.replace("->", "")
    prefix = jumps[0].rate_symbol if len(jumps) == 1 else f"{len(jumps)}x W_{rate}"
    deltas = ", ".join(f"{j.delta_s:+.3f}" for j in sorted(jumps, key=lambda j: j.glob_jump_cnt))
    return f"{prefix}\nΔs = [{deltas}]"


def _compact_png_label(jumps: Sequence[PhysicalJump], d_hkl: float) -> str:
    """Short PNG-only label (min/max Delta_s rather than every value) -- the
    full per-jump list stays in GraphML, in Angstrom, unaffected by this.

    Displayed here as a fraction of d_hkl (the spacing between transport-
    axis unit-cell boundaries) rather than an absolute Angstrom distance:
    on a plane-indexed schematic, "traversed 1/2 a plane spacing" is the
    legible quantity, not an absolute length in a unit the diagram has no
    other use for.
    """
    if not jumps:
        return ""
    rate = jumps[0].rate_symbol if len(jumps) == 1 else jumps[0].pair_key.replace("->", "")
    prefix = jumps[0].rate_symbol if len(jumps) == 1 else f"{len(jumps)}x W_{rate}"
    fractions = sorted(j.delta_s / d_hkl for j in jumps)
    if abs(fractions[0] - fractions[-1]) < 1e-6:
        span = f"{fractions[0]:+.3f}"
    else:
        span = f"[{fractions[0]:+.3f} .. {fractions[-1]:+.3f}]"
    return f"{prefix}\nΔs/d = {span}"


def _json_str(value: Any) -> str:
    """Deterministic, serialization-safe string encoding for a composite
    attribute value (GraphML attributes must be scalar str/int/float/bool)."""
    return json.dumps(value, sort_keys=True)


def _flatten_quotient_edge_attrs(qe: QuotientEdge) -> dict[str, Any]:
    ordered = sorted(qe.jumps, key=lambda j: j.glob_jump_cnt)
    return {
        "count": qe.count,
        "jump_indices": _json_str([j.glob_jump_cnt for j in ordered]),
        "rate_symbols": _json_str([j.rate_symbol for j in ordered]),
        "pair_keys": _json_str(sorted({j.pair_key for j in ordered})),
        "delta_s": _json_str([round(j.delta_s, 9) for j in ordered]),
        "lengths": _json_str([round(j.length, 9) for j in ordered]),
        "images": _json_str([list(j.image) for j in ordered]),
        "endpoints": _json_str([[j.u, j.v] for j in ordered]),
        "label": _jump_edge_display_label(ordered),
    }


def _flatten_full_edge_attrs(jump: PhysicalJump) -> dict[str, Any]:
    return {
        "glob_jump_cnt": jump.glob_jump_cnt,
        "pair_key": jump.pair_key,
        "rate_symbol": jump.rate_symbol,
        "length": jump.length,
        "delta_s": jump.delta_s,
        "image": _json_str(list(jump.image)),
        "delta_r_cart": _json_str([round(c, 9) for c in jump.delta_r_cart]),
        "label": f"{jump.rate_symbol}\nΔs = {jump.delta_s:+.3f}",
    }


# ---------------------------------------------------------------------------
# 6. GraphML export
# ---------------------------------------------------------------------------


def _axis_graph_metadata(
    hkl_conventional: tuple[int, int, int],
    hkl_struct: tuple[int, int, int],
    xi_hat: np.ndarray,
    d_hkl: float,
) -> dict[str, Any]:
    return {
        "transport_hkl_conventional": _json_str(list(hkl_conventional)),
        "transport_hkl_structure": _json_str(list(hkl_struct)),
        "transport_xi_hat_cartesian": _json_str([round(float(c), 9) for c in xi_hat]),
        "interplanar_spacing_d_hkl": float(d_hkl),
        "equivalence_relation": (
            "Two sites are equivalent iff they share the same Wyckoff label and "
            "have equal projection (hkl . frac_coords) onto the transport axis "
            "within tolerance (JumpAxisEquivalenceAnalyzer.find_equiv_classes_for_axis)."
        ),
        "oriented_uc_transform": (
            "oriented_unit_cell_for_axis(structure, hkl_structure): exact-integer "
            "extended-GCD plane basis + metric-aware 2-D Gauss/Lagrange reduction + "
            "brute-force minimal transverse vector, validated against "
            "area(u,v)*height(w) == primitive_cell_volume and |det([w,u,v])| == 1 "
            "(not pymatgen.core.surface.SlabGenerator, which was found to produce a "
            "non-minimal, non-primitive oriented cell for some {hkl} in "
            "body-centered tetragonal structures)."
        ),
    }


def write_full_graph_graphml(
    physical_jumps: Sequence[PhysicalJump],
    site_label: dict[int, str],
    axis_metadata: dict[str, Any],
    output_file: Path,
) -> nx.MultiDiGraph:
    G = nx.MultiDiGraph()
    for key, value in axis_metadata.items():
        G.graph[key] = value
    for node, label in site_label.items():
        G.add_node(int(node), label=label)
    for jump in physical_jumps:
        G.add_edge(jump.u, jump.v, **_flatten_full_edge_attrs(jump))
    nx.write_graphml(G, output_file)
    return G


def write_quotient_graph_graphml(
    quotient_edges: dict[tuple[int, int], QuotientEdge],
    equiv_sites: Sequence[int],
    site_label: dict[int, str],
    axis_metadata: dict[str, Any],
    output_file: Path,
) -> nx.MultiDiGraph:
    members_by_class: dict[int, list[int]] = defaultdict(list)
    for site_idx, cls in enumerate(equiv_sites):
        members_by_class[int(cls)].append(site_idx)

    G = nx.MultiDiGraph()
    for key, value in axis_metadata.items():
        G.graph[key] = value
    for cls_id, members in sorted(members_by_class.items()):
        member_labels = [site_label[m] for m in sorted(members)]
        G.add_node(
            cls_id,
            label=_quotient_node_label(member_labels),
            member_site_ids=_json_str(sorted(members)),
            member_site_labels=_json_str(sorted(set(member_labels))),
        )
    for (c_u, c_v), qe in sorted(quotient_edges.items()):
        G.add_edge(c_u, c_v, **_flatten_quotient_edge_attrs(qe))
    nx.write_graphml(G, output_file)
    return G


# ---------------------------------------------------------------------------
# Drawing the axis-specific graphs (plane-ordered 2-D layout)
# ---------------------------------------------------------------------------


def _axis_positions(
    site_indices: Iterable[int],
    structure: pmg.core.Structure,
    plane_of_site: Sequence[int],
    e2: np.ndarray,
    e3: np.ndarray,
) -> dict[int, tuple[float, float]]:
    """2-D layout: x = plane index (left to right), y = in-plane position
    (projection onto e2; e3 is available but a single in-plane axis keeps
    same-plane sites visibly separated without a 3rd screen dimension).
    Sites are never collapsed just because they share a plane index.

    Both axes are independently rescaled to a fixed target span. This is a
    schematic-diagram concern only (matplotlib node markers have a fixed
    on-screen size irrespective of the data coordinate range, so a genuinely
    tiny span -- e.g. two planes 1 unit apart with near-zero in-plane spread,
    common for high-symmetry axes -- makes markers and self-loops overlap
    into illegibility); it does not change plane order, in-plane ordering,
    or which sites coincide, and has no effect on any physical quantity
    (Delta_s, labels, GraphML data), which are computed independently of
    this layout.
    """
    raw: dict[int, tuple[float, float]] = {}
    for i in site_indices:
        cart = np.asarray(structure.lattice.get_cartesian_coords(structure.frac_coords[i]))
        y = float(np.dot(cart, e2))
        raw[i] = (float(plane_of_site[i]), y)

    target_span = 6.0
    xs = [p[0] for p in raw.values()]
    ys = [p[1] for p in raw.values()]
    x_span = (max(xs) - min(xs)) if xs else 0.0
    y_span = (max(ys) - min(ys)) if ys else 0.0
    x_scale = target_span / x_span if x_span > 1e-9 else 1.0
    y_scale = target_span / y_span if y_span > 1e-9 else 1.0
    return {i: (x * x_scale, y * y_scale) for i, (x, y) in raw.items()}


def draw_axis_full_graph(
    physical_jumps: Sequence[PhysicalJump],
    site_label: dict[int, str],
    structure: pmg.core.Structure,
    plane_of_site: Sequence[int],
    e2: np.ndarray,
    e3: np.ndarray,
    title: str,
    output_file: Path,
    d_hkl: float,
) -> None:
    G = nx.MultiDiGraph()
    for i in site_label:
        G.add_node(i)
    for jump in physical_jumps:
        G.add_edge(jump.u, jump.v)
    pos = _axis_positions(site_label.keys(), structure, plane_of_site, e2, e3)

    # Several physical jumps can share the same (u, v) pair (e.g. 6
    # symmetry-equivalent A-A jumps); a plain {(u,v): label} dict can only
    # hold one label per pair, so group them into one compact label per
    # pair for this quick-look PNG -- every individual jump still has its
    # own edge/attributes in G (and in the GraphML export) unaffected by
    # this display-only grouping.
    jumps_by_pair: dict[tuple[int, int], list[PhysicalJump]] = defaultdict(list)
    for jump in physical_jumps:
        jumps_by_pair[(jump.u, jump.v)].append(jump)
    edge_labels = {pair: _compact_png_label(jumps, d_hkl) for pair, jumps in jumps_by_pair.items()}
    _draw_graph(
        nx.DiGraph(G), site_label, title, output_file, pos=pos, edge_labels=edge_labels,
        figsize=(14, 10),
    )


def draw_axis_quotient_graph(
    quotient_edges: dict[tuple[int, int], QuotientEdge],
    equiv_sites: Sequence[int],
    site_label: dict[int, str],
    structure: pmg.core.Structure,
    plane_of_site: Sequence[int],
    e2: np.ndarray,
    e3: np.ndarray,
    title: str,
    output_file: Path,
    d_hkl: float,
) -> None:
    members_by_class: dict[int, list[int]] = defaultdict(list)
    for site_idx, cls in enumerate(equiv_sites):
        members_by_class[int(cls)].append(site_idx)

    G = nx.DiGraph()
    qlabels: dict[int, str] = {}
    for cls_id, members in members_by_class.items():
        G.add_node(cls_id)
        qlabels[cls_id] = _quotient_node_label([site_label[m] for m in members])
    for (c_u, c_v), qe in quotient_edges.items():
        G.add_edge(c_u, c_v, label=_compact_png_label(qe.jumps, d_hkl))

    # Position each class node at the (mean plane, mean in-plane coordinate)
    # of its members -- for display only; never used to (re)compute Delta_s.
    per_site_pos = _axis_positions(site_label.keys(), structure, plane_of_site, e2, e3)
    pos = {}
    for cls_id, members in members_by_class.items():
        xs = [per_site_pos[m][0] for m in members]
        ys = [per_site_pos[m][1] for m in members]
        pos[cls_id] = (float(np.mean(xs)), float(np.mean(ys)))

    edge_labels = {(u, v): data["label"] for u, v, data in G.edges(data=True)}
    _draw_graph(G, qlabels, title, output_file, pos=pos, edge_labels=edge_labels, figsize=(14, 10))


# ---------------------------------------------------------------------------
# Orchestration
# ---------------------------------------------------------------------------


def process_axis(
    axis: tuple[int, int, int],
    graph: nx.MultiDiGraph,
    structure: pmg.core.Structure,
    symmetry_ctx: JumpNetworkSymmetryContext,
    analyzer: JumpAxisEquivalenceAnalyzer,
    node_labels: dict[int, str],
    save_dir: Path,
) -> dict[str, Any]:
    hkl_conventional = tuple(int(v) for v in axis)
    hkl_struct = tuple(int(v) for v in symmetry_ctx.convert_hkl(axis))

    xi_hat, d_hkl = axis_unit_vector_and_spacing(structure, hkl_struct)
    e2, e3 = axis_inplane_basis(structure, hkl_struct, xi_hat)

    physical_jumps = extract_physical_jumps(graph, structure, node_labels, xi_hat)
    equiv_sites, site_label, projections, plane_of_site = axis_site_labels(
        structure, analyzer, hkl_struct, node_labels
    )

    quotient_edges = build_quotient_graph_explicit(physical_jumps, equiv_sites)
    validate_quotient_graph(physical_jumps, equiv_sites, quotient_edges, n_sites=structure.num_sites)

    axis_metadata = _axis_graph_metadata(hkl_conventional, hkl_struct, xi_hat, d_hkl)
    label = "_".join(str(v) for v in hkl_conventional)

    full_graphml = save_dir / f"full_axis_{label}.graphml"
    write_full_graph_graphml(physical_jumps, site_label, axis_metadata, full_graphml)
    print(f"Saved: {full_graphml}")

    quotient_graphml = save_dir / f"quotient_axis_{label}.graphml"
    write_quotient_graph_graphml(quotient_edges, equiv_sites, site_label, axis_metadata, quotient_graphml)
    print(f"Saved: {quotient_graphml}")

    # PNG rendering is a secondary quick-look (GraphML above is the
    # authoritative output); a matplotlib/networkx edge-label-placement
    # quirk for some self-loop-heavy graphs must never prevent GraphML
    # export for this or any subsequent axis, so failures here are logged
    # and skipped rather than propagated.
    full_png = save_dir / f"full_axis_{label}.png"
    try:
        draw_axis_full_graph(
            physical_jumps, site_label, structure, plane_of_site, e2, e3,
            f"Full jump graph for axis {hkl_conventional} (d={d_hkl:.4f})", full_png,
            d_hkl,
        )
        print(f"Saved: {full_png}")
    except Exception:
        logger.warning("Failed to render %s (GraphML export above is unaffected)", full_png, exc_info=True)

    quotient_png = save_dir / f"quotient_axis_{label}.png"
    try:
        draw_axis_quotient_graph(
            quotient_edges, equiv_sites, site_label, structure, plane_of_site, e2, e3,
            f"Quotient graph for axis {hkl_conventional} (d={d_hkl:.4f})", quotient_png,
            d_hkl,
        )
        print(f"Saved: {quotient_png}")
    except Exception:
        logger.warning("Failed to render %s (GraphML export above is unaffected)", quotient_png, exc_info=True)

    return {"hkl_conventional": hkl_conventional, "xi_hat": xi_hat, "d_hkl": d_hkl}


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Render full and quotient jump-network graphs from example output artifacts."
    )
    parser.add_argument(
        "path",
        help="Path to input YAML or output directory (examples/<name>/output).",
    )
    parser.add_argument(
        "--save-dir",
        default=None,
        help="Directory to write images/GraphML (default: <output>/graph_visualization).",
    )
    args = parser.parse_args()

    yaml_path, output_dir = _resolve_paths(Path(args.path))
    save_dir = Path(args.save_dir).resolve() if args.save_dir else (output_dir / "graph_visualization")
    save_dir.mkdir(parents=True, exist_ok=True)

    params = load_data(str(yaml_path))
    graph = _load_graph(output_dir / "jump_structure.node_link.json")
    node_labels = _extract_node_labels(graph)

    full_png = save_dir / "full_jump_network.png"
    _draw_graph(_collapse_multidigraph(graph), node_labels, "Full jump network", full_png)
    print(f"Saved: {full_png}")

    jump_sites_poscar = output_dir / "jump_sites.poscar"
    if not jump_sites_poscar.exists():
        raise FileNotFoundError(
            f"Missing jump_sites artifact: {jump_sites_poscar}. Run jump_network_builder first."
        )

    # jump_sites.poscar is Stage 1's *matrix-already-removed* output, so
    # there is no matrix species left in it to separate out -- passing
    # params["matrix"] here (pre-existing bug, inherited from before this
    # rewrite) makes PoscarJumpStructureLoader assert on any example whose
    # YAML declares a non-empty matrix (e.g. bcc_int's ["Al"]).
    symmetry_ctx = JumpNetworkSymmetryContext(
        poscar_file=str(jump_sites_poscar),
        crystal_matrix=None,
    )
    structure = symmetry_ctx.get_refined_jump_structure()
    analyzer = JumpAxisEquivalenceAnalyzer(symmetry_ctx)

    axes = _safe_axes(params.get("axes", []))
    if not axes:
        print("No valid [h,k,l] axes in YAML; only the full network overview was rendered.")
        return 0

    seen_axes: list[dict[str, Any]] = []
    for axis in axes:
        seen_axes.append(
            process_axis(axis, graph, structure, symmetry_ctx, analyzer, node_labels, save_dir)
        )

    validate_axis_specificity(seen_axes)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
