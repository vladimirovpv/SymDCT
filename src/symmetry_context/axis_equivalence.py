"""Axis-preserving symmetry and site-equivalence helpers.

This module adds an additional *layer* on top of
``JumpStructureSymmetryRefiner`` that is responsible for partitioning jump
sites into equivalence classes for a given transport axis (specified by
Miller indices).

Equivalence follows Ishioka and Koiwa (Philos. Mag. A 52, 267, 1985, Sec. 3).
Two sites can be merged for a transport direction ``xi`` (perpendicular to
the planes ``hkl``) if a space-group operation ``x' = W x + w`` (fractional
coordinates) maps one onto the other and keeps the transport direction:

* ``hkl W = hkl``: the normal of the planes is preserved, not reversed.  This
  excludes inversion and every operation that turns the direction around.

Such an operation leaves the xi component of every displacement unchanged.
Jump rates and the Fourier phases ``exp(i k . h)`` of the jumps are therefore
mapped onto each other, and the sum of the matrix elements of ``A(k)`` into a
class is the same from every site of the class, for every ``theta``.  The
merge is exact.

Two cases are distinguished, as in the original paper:

* (a) ``hkl . w`` is an integer: the operation keeps every site in its plane
  (rotations about xi, mirror planes containing xi, glide planes whose glide
  vector lies in the plane).  These are the "equivalent sublattices" of the
  paper.
* (b) ``hkl . w`` is not an integer: a screw axis or a glide plane with a
  component along xi.  It carries a site into another plane, and the period of
  the reduced one-dimensional lattice becomes ``a/n``.  Sites related this
  way also belong to one class.

Values of ``site_equivalence`` (``label_source``):

* ``axis_orbit`` (default): orbit under the operations with ``hkl W = hkl``,
  cases (a) and (b).  A class may contain sites from several planes.
* ``plane_orbit``: only case (a).  A class lies in one plane.
  This is the conservative choice.
* ``orbit``: classical spglib orbit (``equivalent_indices``) and the same
  plane.  It is weaker than the two above, because two sites of one orbit may
  be related only by an operation that reverses the transport direction
  (for example inversion about a point in the plane).
* ``wyckoff``: the bare Wyckoff letter and the same plane.  It also merges
  sites with different labels at the same position (orientation
  pseudo-species) and is meant for use together with ``merge_classes: force``.

The merge is checked for lumpability in :mod:`src.matrix_builder.lumpability`
as a safeguard.  The class numbering follows the smallest plane coordinate
``s = (hkl . frac) mod 1`` among the members.

The dot product is meaningful without any change of basis because
fractional coordinates are components in the *direct* lattice basis while
Miller indices are components in the *reciprocal* lattice basis by
construction; their pairing (``hkl . r``) is exactly the (unnormalised)
projection of ``r`` onto the plane normal defined by ``hkl`` and must not be
computed in Cartesian coordinates, where this direct/reciprocal duality no
longer holds.

Classes are numbered in order of increasing plane coordinate, ties broken by
orbit index.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence, Tuple

import numpy as np
import pymatgen as pmg

from src.symmetry_context.jump_symmetry_context import JumpNetworkSymmetryContext


# ``Vector3`` can be either a sequence of floats (list/tuple) or a NumPy array.
# The ``|`` union operator is only available from Python 3.10 onward, but the
# project runs on Python 3.9.  Use ``typing.Union`` for compatibility.
from typing import Union
from fractions import Fraction
from functools import reduce
from math import gcd
Vector3 = Union[Sequence[float], np.ndarray]

def _reduce_hkl_vector(vec: Vector3, tol: float = 1e-6, max_denominator: int = 12) -> Tuple[int, int, int]:
    """Reduce a Miller-index vector to the smallest integer triplet.

    Floating-point noise (e.g. 1.999999999) is cleaned up via limit_denominator,
    but genuine fractional components (e.g. 0.5) are handled by scaling the
    whole vector up to clear denominators *before* rounding — not by rounding
    each component independently, which would destroy fractional information
    (np.rint's round-half-to-even turns 0.5 into 0).
    """
    arr = np.asarray(vec, dtype=float)

    if not np.any(np.abs(arr) > tol):
        return 0, 0, 0

    fracs = [Fraction(v).limit_denominator(max_denominator) for v in arr]
    denoms = [f.denominator for f in fracs]
    lcm = reduce(lambda a, b: a * b // gcd(a, b), denoms, 1)

    scaled = [int(round(v * lcm)) for v in arr]

    div = reduce(gcd, (abs(v) for v in scaled if v != 0))
    reduced = tuple(v // div for v in scaled)
    return reduced  # type: ignore[return-value]


@dataclass
class JumpAxisEquivalenceAnalyzer:
    """Analyse axis-preserving symmetry and equivalent sites.

    This class is designed to sit on top of
    :class:`JumpStructureSymmetryRefiner` and provide a small API for
    finding site equivalence classes for a given axis (specified by
    Miller indices in the *structure* coordinates).
    """

    symmetry_context: JumpNetworkSymmetryContext
    dist_tol: float = 1e-3
    # "axis_orbit" (default): orbit under the operations that preserve the
    # transport direction (hkl W = hkl); screw axes and glide planes along xi
    # may carry a site into another plane.
    # "plane_orbit": the same, but only operations that keep each site in its
    # plane (hkl . w integer).
    # "orbit": classical spglib orbit (``equivalent_indices``).
    # "wyckoff": the bare Wyckoff letter is used instead; this also merges
    # sites that carry different labels at the same position (for example
    # orientation pseudo-species) and is meant for use together with
    # ``merge_classes: force``.
    label_source: str = "axis_orbit"

    def __post_init__(self) -> None:
        self._sym_struct: pmg.core.Structure = self.symmetry_context.get_refined_jump_structure()
        # Populated by ``find_equiv_classes_for_axis`` with the per-site
        # projection (hkl . frac_coords) onto the most recently requested
        # transport axis, index-aligned with the refined jump structure.
        # Exposed so callers can number equivalence classes by their plane
        # position along the axis (e.g. for matrix-printing labels) without
        # recomputing the dot products.
        self.last_site_projections: list[float] = []
        # Number of space-group operations that preserved the plane family
        # in the last plane-orbit analysis (diagnostics).
        self.last_plane_ops_used: int = 0

    # ------------------------------------------------------------------
    # public API
    # ------------------------------------------------------------------
    def find_equiv_classes_for_axis(self, hkl_struct: Vector3) -> List[int]:
        """Return equivalence-class index for each site for a given axis.

        Two sites are equivalent iff they share the same label (plane orbit
        by default, see the module docstring) *and* the plane coordinate
        ``(hkl_struct . frac) mod 1`` is equal within ``dist_tol``. Classes
        are numbered in order of increasing plane coordinate.

        Parameters
        ----------
        hkl_struct
            Miller indices in the *structure* lattice coordinates.

        Returns
        -------
        list of int
            For each site index ``i`` in the symmetrized jump structure,
            the value is the equivalence-class index for the axis
            specified by ``hkl_struct``.
        """

        hkl = np.asarray(hkl_struct, dtype=float)
        jump = self._sym_struct

        # ``wyckoff_letters`` (set by SpacegroupAnalyzer.get_symmetrized_structure)
        # is index-aligned with the structure's own site order, exactly like
        # ``frac_coords`` -- no grouping/reverse-lookup needed to read it.
        # Use the symmetry-orbit index instead of the bare Wyckoff letter as the
        # site label: the letter alone does not distinguish sites that carry
        # different internal states (e.g. pseudo-species for orientations),
        # whereas the orbit index is the space-group equivalence class.
        if self.label_source == "wyckoff":
            wyckoff_letters = jump.wyckoff_letters
        elif self.label_source == "orbit":
            orbit_of_site = [-1] * len(jump)
            for orbit_idx, group in enumerate(jump.equivalent_indices):
                for site_idx in group:
                    orbit_of_site[site_idx] = orbit_idx
            wyckoff_letters = orbit_of_site
        elif self.label_source in ("axis_orbit", "plane_orbit"):
            wyckoff_letters = self._plane_orbit_labels(
                hkl, same_plane=(self.label_source == "plane_orbit")
            )
        else:
            raise ValueError(
                "label_source must be 'axis_orbit', 'plane_orbit', 'orbit' or 'wyckoff', "
                f"got {self.label_source!r}"
            )

        # Wrap each site's fractional coordinates to their canonical [0, 1)
        # representative before projecting. hkl . frac is only translation-
        # invariant up to an integer (shifting frac by any lattice vector n
        # shifts hkl . frac by hkl . n, an integer for reduced Miller indices)
        # -- it is emphatically NOT confined to [0, 1) for a general hkl with
        # more than one nonzero component (e.g. (1,1,1): each frac component
        # is <1, but their sum ranges over [0, 3)). Symmetrized structures can
        # also leave a site's representative sitting exactly on the periodic
        # boundary (frac component == 1.0 instead of 0.0) rather than always
        # wrapped; skipping this step let such a site "look" one whole
        # interplanar spacing away from Wyckoff-equivalent sites that are
        # actually on the very same plane, splitting one true equivalence
        # class into two.
        wrapped_frac_coords = [np.mod(site.frac_coords, 1.0) for site in jump]
        # Plane coordinate s_alpha = (hkl . frac) mod 1 in [0, 1): sites related
        # by a lattice translation along the transport axis share the same
        # plane family.  Values within ``dist_tol`` of 1 are mapped to 0.
        projections = []
        for fc in wrapped_frac_coords:
            phi = float(np.dot(hkl, fc))
            plane = phi - np.floor(phi + self.dist_tol)
            projections.append(0.0 if abs(plane) < self.dist_tol else float(plane))
        self.last_site_projections = projections

        # Group sites by (equal Wyckoff label, equal projection) without
        # relying on any external equivalence-class utility: sort by
        # (wyckoff label, projection), then walk the sorted order, starting a
        # new class whenever the Wyckoff label changes or the projection
        # departs from the current class's reference projection by more than
        # ``dist_tol``.
        # For ``axis_orbit`` the plane coordinate does not split a class:
        # sites related by a screw axis or glide plane along xi lie in
        # different planes and still form one class.
        split_by_plane = self.label_source != "axis_orbit"
        order = sorted(
            range(len(projections)),
            key=lambda site_idx: (
                wyckoff_letters[site_idx],
                projections[site_idx] if split_by_plane else 0.0,
            ),
        )

        raw_classes: list[tuple[float, list[int]]] = []
        current_members: list[int] = []
        current_letter = None
        reference_projection = None
        for idx in order:
            letter = wyckoff_letters[idx]
            proj = projections[idx]
            if (
                current_letter is None
                or letter != current_letter
                or (split_by_plane and abs(proj - reference_projection) > self.dist_tol)
            ):
                if current_members:
                    raw_classes.append((reference_projection, current_members))
                current_members = []
                current_letter = letter
                reference_projection = proj
            current_members.append(idx)
        if current_members:
            raw_classes.append((reference_projection, current_members))
        # A class spanning several planes is placed at its smallest plane
        # coordinate.
        raw_classes = [
            (min(projections[m] for m in members), members)
            for _, members in raw_classes
        ]

        # Number classes in order of increasing projection (global ordering
        # across Wyckoff labels).
        raw_classes.sort(key=lambda item: item[0])

        equiv_sites = [-1] * len(projections)
        for cls_id, (_, members) in enumerate(raw_classes):
            for idx in members:
                equiv_sites[idx] = cls_id

        return equiv_sites

    # ------------------------------------------------------------------
    # internal helpers
    # ------------------------------------------------------------------
    def _plane_orbit_labels(self, hkl: np.ndarray, same_plane: bool) -> List[int]:
        """Return an equivalence-class id for every site of the refined structure.

        Uses the space-group operations of the primitive cell that satisfy
        ``hkl W = hkl`` (and, if ``same_plane`` is true, also ``hkl . w``
        integer).  Sites that such an operation maps onto each other (same
        species, same position modulo the lattice, within ``dist_tol``) are
        put into one equivalence class.  An operation that does not map the
        structure onto itself within tolerance is skipped, which can only
        make the classes smaller.
        """

        struct = self._sym_struct
        n_sites = len(struct)
        frac = [np.mod(site.frac_coords, 1.0) for site in struct]
        parent = list(range(n_sites))

        def find(i: int) -> int:
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        self.last_plane_ops_used = 0
        for op in self.symmetry_context.get_symmetry_operations():
            rot = np.asarray(op.rotation_matrix, dtype=float)
            trans = np.asarray(op.translation_vector, dtype=float)
            if not np.allclose(hkl @ rot, hkl, atol=1e-6):
                continue
            shift = float(hkl @ trans)
            if same_plane and abs(shift - round(shift)) > self.dist_tol:
                continue
            images = []
            for i in range(n_sites):
                x = rot @ frac[i] + trans
                match = -1
                for j in range(n_sites):
                    if struct[j].species != struct[i].species:
                        continue
                    diff = x - frac[j]
                    diff -= np.round(diff)
                    if np.all(np.abs(diff) < self.dist_tol):
                        match = j
                        break
                if match < 0:
                    images = None
                    break
                images.append(match)
            if images is None:
                continue
            self.last_plane_ops_used += 1
            for i, j in enumerate(images):
                parent[find(i)] = find(j)

        roots = sorted({find(i) for i in range(n_sites)})
        root_id = {r: k for k, r in enumerate(roots)}
        return [root_id[find(i)] for i in range(n_sites)]
