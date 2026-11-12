"""Assign labels to symmetrically equivalent sites.

This module provides :class:`SymmetryNodeLabeler`, which converts the
``equivalent_sites`` information stored on a symmetrized structure into
per‑site labels, implemented here as a small, independent utility.
"""

from __future__ import annotations

from dataclasses import dataclass
import string
from typing import Sequence

import pymatgen as pmg

from src.common_utils.utils import to_rational


class SymmetryNodeLabeler:
    """Generate node labels based on spacegroup equivalent sites.

    Parameters
    ----------
    symmetrized_structure
        A pymatgen structure returned by
        ``SpacegroupAnalyzer.get_symmetrized_structure()`` with an
        additional integer site property ``"site_idx"`` defining a
        unique index for each site.
    """

    def __init__(
        self,
        symmetrized_structure: pmg.core.Structure,
        site_labels: Sequence[str] | None = None,
    ) -> None:
        self._sym_struct = symmetrized_structure

    # ------------------------------------------------------------------
    # public API
    # ------------------------------------------------------------------
    def compute_labels(self) -> list[str]:
        """Return a list of labels for each site in the structure.

        The default label for a site is the chemical symbol of its
        species. If the number of symmetry‑equivalent site groups exceeds
        the number of distinct symbols, the groups are relabelled with
        capital letters ``"A"``, ``"B"``, ... based on their
        equivalence class index, in the same way as in the original
        implementation.
        """

        sym_struct = self._sym_struct
        equivalent_sites = getattr(sym_struct, "equivalent_sites", [])
        node_labels = [specie.symbol for specie in sym_struct.species]
        sublat_set = set(node_labels)

        # If there are more equivalence classes than unique species,
        # assign explicit sublattice labels A, B, C, ...
        if len(sublat_set) < len(equivalent_sites):
            for eq_class, equiv_sites in enumerate(equivalent_sites):
                label_char = string.ascii_uppercase[eq_class]
                for site in equiv_sites:
                    idx = site.properties["site_idx"]
                    node_labels[idx] = label_char

        return node_labels


@dataclass(frozen=True)
class NonEquivalentSiteLabelData:
    """Labels and pretty coordinate strings for representative sites."""

    labels: list[str]
    pretty_site_coords: dict[str, str]


def build_non_equivalent_site_label_data(
    node_labels: Sequence[str],
    representative_site_indices: Sequence[int],
    representative_site_coords: Sequence[Sequence[float]],
    representative_projections: Sequence[float] | None = None,
) -> NonEquivalentSiteLabelData:
    """Construct non-equivalent site labels from explicit inputs.

    The trailing integer in each label numbers sites through in order of
    increasing projection onto the transport axis (i.e. their plane
    position along the axis), separately for each Wyckoff-derived label,
    when ``representative_projections`` is supplied. Without projections
    (e.g. callers not tied to a specific axis), classes are numbered in
    the order they are given, matching the historical behaviour.
    """

    def _site_coords_value(site_index: int, vector: Sequence[float]) -> str:
        values = " ".join(f"{str(to_rational(float(ev), tolerance=1e-6)): >6}" for ev in vector)
        return f"{values} ({site_index:>2})"

    n = len(representative_site_indices)
    if representative_projections is None:
        plane_numbers = list(range(n))
    else:
        # Rank sites that share the same Wyckoff-derived label by increasing
        # projection along the axis, so the number reflects plane position.
        by_label: dict[str, list[int]] = {}
        for i, site_idx in enumerate(representative_site_indices):
            by_label.setdefault(node_labels[site_idx], []).append(i)
        plane_numbers = [0] * n
        for positions in by_label.values():
            ordered = sorted(positions, key=lambda i: representative_projections[i])
            for rank, i in enumerate(ordered):
                plane_numbers[i] = rank

    labels: list[str] = []
    pretty_site_coords: dict[str, str] = {}
    for i, (site_idx, coords) in enumerate(zip(representative_site_indices, representative_site_coords)):
        lbl = f"{node_labels[site_idx]}{plane_numbers[i]}"
        labels.append(lbl)
        pretty_site_coords[lbl] = _site_coords_value(site_idx, coords)

    return NonEquivalentSiteLabelData(labels=labels, pretty_site_coords=pretty_site_coords)


def print_non_equivalent_site_coordinates(
    sorted_labels: Sequence[str],
    pretty_site_coords: dict[str, str],
) -> None:
    """Print representative-site coordinates in a compact, label-indexed table."""

    print(78 * "#")
    print("Fractional coordinates of nonequivalent jump sites ")
    print(78 * "#")
    for lbl in sorted_labels:
        print(f"{lbl}: {pretty_site_coords[lbl]}")
    print(78 * "#")
