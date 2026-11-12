"""Encoding/decoding helpers for serialized jump vectors.

Stage 1 writes jump vectors into node-link JSON edge attributes. To keep the
contract stable and avoid ad-hoc conversions in each stage, this module
provides one canonical schema and parser:

- Canonical JSON form per component: ``{"num": <int>, "den": <int>}``
  representing an exact rational ``num/den``.
- Decoder also accepts plain scalar/string forms for backward compatibility
  with older serialized data.
"""

from __future__ import annotations

from typing import Any, Sequence

import sympy as sp


def encode_jump_vector(jump_vector: Sequence[Any]) -> list[dict[str, int]]:
    """Encode jump-vector components as exact rational pairs."""

    encoded: list[dict[str, int]] = []
    for idx, component in enumerate(jump_vector):
        expr = sp.nsimplify(component, rational=True)
        if expr.is_real is False:
            raise ValueError(f"Jump vector component must be real at index {idx}: {component!r}")
        if not expr.is_rational:
            raise ValueError(
                f"Jump vector component must be rational at index {idx}: {component!r}"
            )
        rat = sp.Rational(expr)
        encoded.append({"num": int(rat.p), "den": int(rat.q)})
    return encoded


def decode_jump_vector(jump_vector: Sequence[Any], *, edge_id: str = "<edge>") -> list[sp.Expr]:
    """Decode serialized jump-vector components to numeric SymPy expressions."""

    decoded: list[sp.Expr] = []
    for idx, component in enumerate(jump_vector):
        if isinstance(component, dict):
            if "num" not in component or "den" not in component:
                raise ValueError(
                    f"Invalid jump-vector dict at {edge_id}[{idx}]: expected keys 'num' and 'den'"
                )
            value = sp.Rational(int(component["num"]), int(component["den"]))
        else:
            try:
                value = sp.nsimplify(component, rational=True)
            except Exception as exc:
                raise ValueError(
                    f"Invalid jump vector component at {edge_id}[{idx}]: {component!r}"
                ) from exc

        if value.is_real is False:
            raise ValueError(
                f"Jump vector component must be real at {edge_id}[{idx}]: {component!r}"
            )
        decoded.append(value)
    return decoded

