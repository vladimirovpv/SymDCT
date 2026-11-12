"""Symmetry analysis and primitive/symmetrized jump structure.

This module provides :class:`JumpStructureSymmetryRefiner`, a lightweight
utility to:

* analyse the symmetry of a *jump* structure (and optionally a *matrix*
  structure for consistency),
* construct a refined conventional cell,
* reduce it to a primitive cell, and
* obtain a *symmetrized* primitive structure with ``equivalent_sites``
  information.

The implementation covers just the essentials required by the
jump‑network builder.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np
import pymatgen as pmg
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


@dataclass
class JumpStructureSymmetryRefiner:
    """Analyse symmetry and build symmetrized primitive jump structure.

    Parameters
    ----------
    matrix_structure
        Structure containing matrix atoms.
    jump_structure
        Structure containing only jump / interstitial sites.
    symprec
        Symmetry tolerance passed to :class:`SpacegroupAnalyzer`.
    check_symmetry_consistency
        If ``True`` (default), checks that matrix and jump structures
        have the same number of symmetry operations and raises
        :class:`ValueError` otherwise.
    """

    matrix_structure: pmg.core.Structure
    jump_structure: pmg.core.Structure
    symprec: float = 0.1
    check_symmetry_consistency: bool = True

    def __post_init__(self) -> None:
        self._lattice_type: str | None = None
        self._s_matrix: np.ndarray | None = None
        self._symmetrized: pmg.core.Structure | None = None
        self._spg_analyzer: Optional[SpacegroupAnalyzer] = None
        self._analyse()

    # ------------------------------------------------------------------
    # internal helpers
    # ------------------------------------------------------------------
    def _analyse(self) -> None:
        # Spacegroup analyzers for jump and matrix structures
        spg_jump = SpacegroupAnalyzer(self.jump_structure, symprec=self.symprec)
        spg_matrix = SpacegroupAnalyzer(self.matrix_structure, symprec=self.symprec)

        if self.check_symmetry_consistency:
            n_ops_matrix = len(spg_matrix.get_symmetry_operations(cartesian=False))
            n_ops_jump = len(spg_jump.get_symmetry_operations(cartesian=False))
            if n_ops_matrix != n_ops_jump:
                raise ValueError(
                    "Matrix and jump structures have different numbers of symmetry "
                    f"operations: matrix={n_ops_matrix}, jump={n_ops_jump}. "
                    "Check input coordinates or adjust symprec."
                )

        # Basic lattice information from the jump structure
        self._lattice_type = spg_jump.get_lattice_type()

        # Refine atomic positions according to symmetry (conventional cell)
        refined_conv = spg_jump.get_refined_structure()
        spg_conv = SpacegroupAnalyzer(refined_conv, symprec=self.symprec)

        # Reduce conventional to primitive cell
        primitive = spg_conv.find_primitive()
        self._s_matrix = np.array(
            spg_conv.get_conventional_to_primitive_transformation_matrix()
        )

        # Symmetrized structure for primitive cell
        spg_prim = SpacegroupAnalyzer(primitive, symprec=self.symprec)
        self._spg_analyzer = spg_prim
        sym_struct = spg_prim.get_symmetrized_structure()

        # Add sequential site index property (used for labelling and
        # graph construction in the jump‑network builder).
        sym_struct.add_site_property("site_idx", list(range(sym_struct.num_sites)))
        self._symmetrized = sym_struct

    # ------------------------------------------------------------------
    # public API
    # ------------------------------------------------------------------
    @property
    def lattice_type(self) -> str:
        """Lattice type of the jump structure (e.g. "cubic")."""

        assert self._lattice_type is not None
        return self._lattice_type

    @property
    def s_matrix(self) -> np.ndarray:
        """Conventional → primitive transformation matrix as a 3×3 array."""

        assert self._s_matrix is not None
        return self._s_matrix

    @property
    def symmetrized_structure(self) -> pmg.core.Structure:
        """Symmetrized primitive structure used for jump‑network building."""

        assert self._symmetrized is not None
        return self._symmetrized

    @property
    def spacegroup_analyzer(self) -> SpacegroupAnalyzer:
        """Spacegroup analyzer associated with the primitive cell."""

        assert self._spg_analyzer is not None
        return self._spg_analyzer
