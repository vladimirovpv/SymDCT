"""Classification of directed jumps by space-group orbits.

Two directed jumps share a jump class (and therefore one rate symbol) when a
symmetry operation of the crystal maps one onto the other.  The operations are
taken from spglib for the refined jump structure, with node labels acting as
atom types, so an operation never exchanges sites that carry different labels.

The older classification by the triplet (label of the start site, label of the
end site, jump length) is retained as ``jump_classification: triplet``.  The
triplet cannot separate symmetry-inequivalent jumps of equal length between
sites with equal labels, whereas orbits can.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import numpy as np
import spglib

# A directed edge: (start site index, end site index, jump vector in fractional
# coordinates of the structure, including the lattice translation of the end site).
Edge = tuple[int, int, Sequence[float]]

_KEY_DECIMALS = 4


@dataclass(frozen=True)
class JumpOrbits:
    """Result of :func:`classify_jump_orbits`.

    Attributes
    ----------
    ids
        Orbit index of every edge, numbered by first appearance.
    reverse_ids
        Orbit index of the reversed jump for every orbit index.
    n_ops
        Number of symmetry operations reported by spglib.
    n_skipped_ops
        Operations that could not be applied because they do not map the
        sites onto themselves within ``site_tol``.
    n_unmapped_edges
        Number of (operation, edge) pairs whose image is not an edge of the
        network.  A non-zero value signals that the neighbour selection is
        not symmetric.
    """

    ids: list[int]
    reverse_ids: dict[int, int]
    n_ops: int
    n_skipped_ops: int
    n_unmapped_edges: int


def _edge_key(i: int, j: int, h) -> tuple:
    return (int(i), int(j), tuple(np.round(np.asarray(h, dtype=float), _KEY_DECIMALS) + 0.0))


def _site_map(positions: np.ndarray, lattice: np.ndarray, type_ids: Sequence[int],
              rot: np.ndarray, trans: np.ndarray, site_tol: float) -> list[int] | None:
    """Map every site onto the site that (rot, trans) carries it to, or None."""

    images = positions @ rot.T + trans
    mapping: list[int] = []
    for i, x in enumerate(images):
        found = None
        for k, xk in enumerate(positions):
            if type_ids[k] != type_ids[i]:
                continue
            d = x - xk
            d -= np.round(d)
            if np.linalg.norm(d @ lattice) < site_tol:
                found = k
                break
        if found is None:
            return None
        mapping.append(found)
    return mapping


def classify_jump_orbits(
    lattice: np.ndarray,
    frac_coords: np.ndarray,
    labels: Sequence[str],
    edges: Sequence[Edge],
    symprec: float = 0.1,
    site_tol: float | None = None,
) -> JumpOrbits:
    """Partition directed edges into orbits of the space group.

    Parameters
    ----------
    lattice
        Lattice vectors as rows, in Angstrom.
    frac_coords
        Fractional coordinates of the sites.
    labels
        Node label of every site.  Labels play the role of atom types.
    edges
        Directed edges ``(i, j, h)`` with ``h`` the jump vector in fractional
        coordinates.
    symprec
        Symmetry tolerance in Angstrom, passed to spglib.
    site_tol
        Tolerance in Angstrom for matching a transformed site with a site of
        the structure.  Defaults to ``symprec / 2``.
    """

    lattice = np.asarray(lattice, dtype=float)
    positions = np.asarray(frac_coords, dtype=float)
    site_tol = symprec / 2 if site_tol is None else site_tol

    label_ids = {lbl: k for k, lbl in enumerate(sorted(set(labels)))}
    type_ids = [label_ids[lbl] for lbl in labels]
    dataset = spglib.get_symmetry((lattice, positions, type_ids), symprec=symprec)
    if dataset is None:
        raise RuntimeError("spglib could not determine the symmetry of the jump structure")

    keys: dict[tuple, int] = {}
    for e, (i, j, h) in enumerate(edges):
        keys[_edge_key(i, j, h)] = e

    parent = list(range(len(edges)))

    def find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    n_skipped = 0
    n_unmapped = 0
    for rot, trans in zip(dataset["rotations"], dataset["translations"]):
        mapping = _site_map(positions, lattice, type_ids, rot, trans, site_tol)
        if mapping is None:
            n_skipped += 1
            continue
        for e, (i, j, h) in enumerate(edges):
            image = keys.get(_edge_key(mapping[i], mapping[j], rot @ np.asarray(h, dtype=float)))
            if image is None:
                n_unmapped += 1
                continue
            a, b = find(e), find(image)
            if a != b:
                parent[a] = b

    numbering: dict[int, int] = {}
    ids: list[int] = []
    for e in range(len(edges)):
        ids.append(numbering.setdefault(find(e), len(numbering)))

    reverse_ids: dict[int, int] = {}
    for e, (i, j, h) in enumerate(edges):
        reverse = keys.get(_edge_key(j, i, -np.asarray(h, dtype=float)))
        if reverse is not None:
            reverse_ids.setdefault(ids[e], ids[reverse])

    return JumpOrbits(
        ids=ids,
        reverse_ids=reverse_ids,
        n_ops=len(dataset["rotations"]),
        n_skipped_ops=n_skipped,
        n_unmapped_edges=n_unmapped,
    )
