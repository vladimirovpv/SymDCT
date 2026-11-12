"""Small crystallographic rotation helpers used across the public code.

Currently this module only provides :func:`rotation_matches_G`, a utility
for checking whether a symmetry operation preserves a given Miller axis.

This is the light‑weight home for this functionality. Keeping it under
``src.common_utils`` avoids polluting the top‑level namespace.
"""

from __future__ import annotations

from typing import Iterable

import numpy as np


def rotation_matches_G(
    op_rotation_matrix: Iterable[Iterable[float]],
    G: Iterable[float],
    tol: float = 1e-8,
) -> bool:
    """Return ``True`` if a symmetry rotation preserves a Miller vector.

    Parameters
    ----------
    op_rotation_matrix
        3x3 matrix acting on *direct‑space fractional* coordinates (as
        provided by ``pymatgen`` symmetry operations).
    G
        Miller index vector ``(h, k, l)`` in reciprocal fractional coordinates.
    tol
        Absolute tolerance used in the ``np.allclose`` comparison.

    Notes
    -----
    All current call sites *only* treat rotations that map ``G``
    to ``+G`` as axis‑preserving – a rotation that flips the direction
    is **not** considered equivalent for our purposes.

    To make that intent explicit, this helper returns a single
    boolean and only reports ``True`` when the induced action on
    reciprocal space maps ``G`` onto itself (within ``tol``).
    """

    G_arr = np.asarray(G, dtype=float)
    W_direct = np.asarray(op_rotation_matrix, dtype=float)

    try:
        W_rec = np.linalg.inv(W_direct).T
    except np.linalg.LinAlgError as exc:  # pragma: no cover - defensive
        raise ValueError("Rotation matrix is singular and cannot be inverted") from exc

    G_rot = W_rec @ G_arr

    # Only accept +G, not -G, as preserving the axis direction.
    return bool(np.allclose(G_rot, G_arr, atol=tol))


__all__ = ["rotation_matches_G"]
