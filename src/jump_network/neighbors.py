"""Neighbor discovery helpers for jump-network construction."""

from __future__ import annotations

from typing import Any


def get_all_neighbors(jump_sites: Any, lat_const: Any, params: dict) -> Any:
    """Return all neighbors within configured cutoff distance.

    ``lat_const`` is accepted for call-site compatibility and is not used.
    """

    max_jump_length = params.get("max_jump_length") or params.get("max_jump_length_angstrom")
    if max_jump_length is None:
        raise KeyError("'max_jump_length' (or 'max_jump_length_angstrom') must be set in params")

    all_neighbors = jump_sites.get_all_neighbors(max_jump_length)
    if not any(all_neighbors):
        raise ValueError(
            f"No neighbors found within cutoff distance {max_jump_length}. "
            "Try to increase the cutoff."
        )
    return all_neighbors

