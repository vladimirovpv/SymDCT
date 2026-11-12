"""Centralised handling of the conventional→primitive transformation matrix (S‑matrix).

The :class:`SMatrixHandler` encapsulates the transformation matrix and provides
clear, consistent accessors for the rest of the code base.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np

from .symmetry_refiner import JumpStructureSymmetryRefiner


@dataclass
class SMatrixHandler:
    """Utility class for the conventional→primitive transformation matrix.

    Parameters
    ----------
    refiner:
        Instance of :class:`JumpStructureSymmetryRefiner` from which the raw
        NumPy matrix is obtained via ``refiner.s_matrix``.
    """

    refiner: JumpStructureSymmetryRefiner
    _convertion_matrix: np.ndarray | None = None

    def __post_init__(self) -> None:
        # Pull the matrix from the refiner once – it is a 3×3 NumPy array.
        self._convertion_matrix = np.array(self.refiner.s_matrix, dtype=float)

    # ---------------------------------------------------------------------
    # Public accessors
    # ---------------------------------------------------------------------
    def matrix(self) -> np.ndarray:
        """Return the raw NumPy 3×3 transformation matrix.

        The returned array is a *copy* to prevent accidental in‑place
        modification of the internal state.
        """
        assert self._convertion_matrix is not None
        return self._convertion_matrix.copy()

    def as_sympy(self):
        """Return the matrix as a :class:`sympy.Matrix` for symbolic work.

        The import is performed lazily to avoid a hard dependency at module import
        time, which resolves the IDE's unresolved reference warning when the
        optional ``sympy`` package is not present in the environment.
        """
        import sympy as sp
        return sp.Matrix([[sp.nsimplify(el, rational=True) for el in row] for row in self.matrix()])

    # ---------------------------------------------------------------------
    # Convenience helpers used by the existing ``axis_equivalence`` module.
    # ---------------------------------------------------------------------
    def convert_hkl_to_structure(self, hkl: Any) -> tuple[int, int, int]:
        """Convert Miller indices from conventional to *structure* coordinates.

        This mirrors the logic previously implemented in
        ``axis_equivalence.convert_hkl_conventional_to_structure`` but routes
        through the centralised handler.  The conversion uses the inverse‑
        transpose of the conventional→primitive matrix because Miller indices
        are reciprocal‑space vectors.
        """
        # Import locally to avoid circular imports.
        from .axis_equivalence import _reduce_hkl_vector

        conv2prim = self.matrix()
        hkl_arr = np.asarray(hkl, dtype=float)

        # Identity or singular matrix – fall back to reduction only.
        if np.allclose(conv2prim, np.eye(3)) or abs(np.linalg.det(conv2prim)) < 1e-12:
            return _reduce_hkl_vector(hkl_arr)

        hkl_struct = conv2prim.T @ hkl_arr
        return _reduce_hkl_vector(hkl_struct)

    def pretty_print(self, tol: float = 1e-6) -> None:
        """Pretty‑print the S‑matrix using ``sympy.pprint``.

        This method mirrors ``jump_network_builder.geometry.pretty_s_matrix``
        but is kept here for convenience when the handler is used directly.
        """
        sp.pprint(self.as_sympy().applyfunc(lambda e: sp.nsimplify(e, tolerance=tol)))
