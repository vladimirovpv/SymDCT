"""High‑level symmetry context for jump‑network construction.

This module exposes :class:`JumpNetworkSymmetryContext`, a small façade
tailored to the needs of ``build_jump_network_from_poscar``:

* load a POSCAR file and separate matrix/jump structures,
* perform symmetry analysis to build a symmetrized primitive structure,
* provide per‑site labels based on spacegroup equivalent sites.

The implementation is intentionally self‑contained and lightweight,
covering just the capabilities required for jump‑network building.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional
import numpy as np

import pymatgen as pmg
import sympy as sp  # Imported here for type hints and matrix utilities

from src.jump_network_builder.loader import PoscarJumpStructureLoader
from .symmetry_refiner import JumpStructureSymmetryRefiner
from .node_labeler import SymmetryNodeLabeler


@dataclass
class JumpNetworkSymmetryContext:
    """Symmetry context used by the jump‑network builder.

    Parameters
    ----------
    poscar_file
        Path to a POSCAR/CONTCAR file.
    crystal_matrix
        Optional list of species symbols that define the matrix. If
        provided, they are removed from the *jump* structure and all
        other species are removed from the *matrix* structure.
    symprec
        Symmetry tolerance passed to :class:`SpacegroupAnalyzer`.
    check_symmetry_consistency
        Whether to check that matrix and jump structures share the same
        number of symmetry operations. Enabled by default as a
        consistency guard.
    """

    poscar_file: str
    crystal_matrix: Optional[List[str]] = None
    symprec: float = 0.1
    check_symmetry_consistency: bool = True

    def __post_init__(self) -> None:
        # 1) Load POSCAR and separate matrix/jump structures
        self._loader = PoscarJumpStructureLoader(self.poscar_file, self.crystal_matrix)

        # 2) Analyse symmetry and build symmetrized primitive structure
        self._refiner = JumpStructureSymmetryRefiner(
            self._loader.matrix_structure,
            self._loader.jump_structure,
            symprec=self.symprec,
            check_symmetry_consistency=self.check_symmetry_consistency,
        )
        self._symbol_alias_map = self._loader.raw_symbol_alias_map

        # 3) Compute node labels based on equivalent sites
        self._labeler = SymmetryNodeLabeler(self._refiner.symmetrized_structure)
        self._node_labels: list[str] = self._labeler.compute_labels()

        # Initialise centralised S‑matrix handler for later use.
        from .smatrix_handler import SMatrixHandler
        self._smatrix_handler = SMatrixHandler(self._refiner)

    # ------------------------------------------------------------------
    # public API used by build_jump_network_from_poscar (conceptually)
    # ------------------------------------------------------------------
    def get_refined_jump_structure(self) -> pmg.core.Structure:
        """Return the symmetrized primitive jump structure."""

        return self._refiner.symmetrized_structure

    def get_node_lbls(self) -> list[str]:
        """Return per‑site labels derived from symmetry equivalence."""

        return self._node_labels

    def get_symbol_alias_map(self) -> dict[str, str]:
        """Return canonical species -> original POSCAR alias mapping."""
        return dict(self._symbol_alias_map)

    # New accessor returning the handler instance.
    def get_smatrix_handler(self):
        """Return the :class:`SMatrixHandler` instance for centralized operations."""
        return self._smatrix_handler

    # ------------------------------------------------------------------
    # Miller‑index conversion helper
    # ------------------------------------------------------------------
    def convert_hkl(self, hkl):
        """Convert Miller indices from the conventional cell to the primitive‑cell basis if required.
        Delegates to :meth:`SMatrixHandler.convert_hkl_to_structure`.
        """
        if _is_conventional := self.is_conventional_cell():
            # Conventional cell – use Miller indices directly.
            hkl_converted= tuple(int(v) for v in hkl)
        else:
            # Primitive cell – convert to primitive basis using the symmetry context.
            hkl_converted = self.to_structure_hkl(hkl)

        # print(f"Simulation cell is {'conventional' if _is_conventional else 'primitive'}.")

        return hkl_converted

    def to_structure_hkl(self, hkl):
        """Convert Miller indices given for the conventional cell to the basis of the structure.

        The result is expressed in the basis of the simulation cell (the primitive
        cell if the structure is primitive).  Delegates to
        :meth:`SMatrixHandler.convert_hkl_to_structure`.
        """
        return self._smatrix_handler.convert_hkl_to_structure(hkl)

    # ------------------------------------------------------------------
    # New helper: determine if the simulation cell is a conventional cell.
    # ------------------------------------------------------------------
    def is_conventional_cell(self) -> bool:
        """Return ``True`` if the cell is conventional.

        The original implementation used a two‑pass ``SpacegroupAnalyzer``
        approach.  A cheaper and equivalent check is to examine the determinant
        of the conventional→primitive transformation matrix (the S‑matrix).  If
        ``|det(S)| > 1`` the cell contains a conventional supercell; otherwise it
        is already primitive.
        """
        # The S‑matrix handler stores the matrix as a NumPy array.
        smatrix = self._smatrix_handler.matrix()
        det = np.linalg.det(smatrix)
        return abs(det) > 1

    def get_lattice_type(self) -> str:
        """Return lattice type string ("cubic", "hexagonal", ...)."""

        return self._refiner.lattice_type

    def get_symbolic_lattice_vectors(self) -> sp.Matrix:
        """Return the symbolic lattice‑vector matrix.

        The matrix is expressed in the *conventional* cell basis.  If the
        underlying structure is primitive, the matrix is automatically converted
        to the primitive basis using the centralised :class:`SMatrixHandler`.
        This provides a single source of truth for the conventional→primitive
        transformation and hides the conversion from callers.
        """

        a, b, c = sp.symbols("a b c", real=True, positive=True)
        lattice_type = self.get_lattice_type()
        if lattice_type in {"cubic"}:
            conventional = sp.matrices.diag(a, a, a)
        elif lattice_type in {"hexagonal"}:
            conventional = sp.Matrix(
                [
                    [a, 0, 0],
                    [-a / 2, sp.sqrt(3) * a / 2, 0],
                    [0, 0, c],
                ]
            )
        elif lattice_type in {"tetragonal"}:
            conventional = sp.matrices.diag(a, a, c)
        else:
            conventional = sp.matrices.diag(a, b, c)

        # If the cell is primitive, convert to the primitive basis.
        if not self.is_conventional_cell():
            # Use the handler to obtain the conversion matrix as a SymPy object.
            conv2prim = self.get_smatrix_handler().as_sympy()
            # Columns of ``conventional`` are the basis vectors; multiply on the left.
            return conv2prim * conventional
        return conventional

    # ------------------------------------------------------------------
    # advanced symmetry accessors
    # ------------------------------------------------------------------
    def get_symmetry_operations(self, cartesian: bool = False):
        """Return the list of symmetry operations for the primitive cell.

        Parameters
        ----------
        cartesian
            Passed through to ``SpacegroupAnalyzer.get_symmetry_operations``.

        Returns
        -------
        list
            A list of ``SymmOp`` objects representing the symmetry
            operations of the primitive cell.
        """

        return self._refiner.spacegroup_analyzer.get_symmetry_operations(
            cartesian=cartesian
        )

    # def get_refiner(self) -> JumpStructureSymmetryRefiner:
    #     """Return the underlying symmetry refiner used by this context."""
    #
    #     return self._refiner
