"""Centralized utilities for label manipulation and pair-key normalization."""

from __future__ import annotations

import re
from typing import Any


def sort_pair_key(a: str, b: str) -> str:
    """Return sorted key for bond pair lookup, order-independent.
    
    Example: sort_pair_key("B", "A") -> "A-B"
    """
    return "-".join(sorted([a, b]))


def map_pair_key_to_jump_label(
    pair_key: str,
    alias_map: dict[str, str],
) -> str:
    """Map a pair key like 'A-B' or 'H-T' to a compact jump label like 'AB' or 'TT'.
    
    Uses symbol_alias_map to resolve display names (e.g., "H" -> "T").
    
    Parameters
    ----------
    pair_key
        A pair key like "A-B", "H-T", "O-T" from the energetics YAML.
    alias_map
        Mapping from canonical site labels to display names, e.g. 
        {"H": "T", "O": "O", "A": "A", "B": "B"}.
    
    Returns
    -------
    str
        Compact label like "AB", "TO", "TT" (no separator).
    """
    parts = pair_key.split("-")
    if len(parts) != 2:
        raise ValueError(f"Invalid pair key: {pair_key!r}")
    
    resolved = [alias_map.get(part.strip(), part.strip()) for part in parts]
    return "".join(resolved)


def remap_pair_keys_in_template(
    energetics_template: dict[str, Any],
    alias_map: dict[str, str] | None = None,
) -> None:
    """Remap pair keys in energetics template jumps section using symbol alias map.
    
    This is called ONCE during template creation (geometry.py) before writing to disk.
    The evaluator reads the already-processed YAML and does NOT call this again.
    
    Modifies the template in-place.
    """
    if alias_map is None:
        alias_map = energetics_template.get("symbol_alias_map", {})
    
    if not alias_map:
        return
    
    template_jumps = energetics_template.get("jumps", {})
    mapped_tmpl_jumps: dict[str, list[dict[str, Any]]] = {}
    
    for pair_key, entries in template_jumps.items():
        mapped_labels = [alias_map.get(label.strip(), label.strip()) 
                        for label in pair_key.split("-")]
        mapped_pair_key = "-".join(mapped_labels)
        mapped_tmpl_jumps.setdefault(mapped_pair_key, []).extend(entries)
    
    energetics_template["jumps"] = mapped_tmpl_jumps


def max_W_index_from_symbols(symset) -> int | None:
    """Return the highest integer suffix among ``W<n>`` symbols in *symset*.

    Parameters
    ----------
    symset:
        Iterable of sympy.Symbol objects (or anything whose ``str()`` yields
        a name like ``"W0"``, ``"W12"``, etc.).

    Returns
    -------
    int or None
        Maximum index found, or ``None`` if no matching symbol is present.
    """
    idxs = [
        int(m.group(1))
        for s in symset
        if (m := re.search(r'(\d+)$', str(s)))
    ]
    return max(idxs) if idxs else None


def map_jump_label_to_id(
    jump_info: list[tuple],
    alias_map: dict[str, str],
) -> dict[str, str]:
    """Build a mapping from symbolic jump labels to numeric jump-ID symbols.

    For each entry in *jump_info* the pair key is resolved via *alias_map* to
    produce a compact label (e.g. ``"AB"``), which is then used to construct
    the symbolic name ``"W_AB"``.  The corresponding jump-ID symbol ``"W0"``
    is the value.

    Parameters
    ----------
    jump_info:
        Sequence of tuples ``(jump_id, pair_key, direction, length, em, fq)``
        as returned by :func:`~src.common_utils.utils.get_jump_info`.
    alias_map:
        Symbol alias map, e.g. ``{"H": "T"}``, used to rename site labels
        in pair keys before building the compact label.

    Returns
    -------
    dict
        Mapping ``{"W_AB": "W0", "W_BA": "W1", ...}``.
    """
    result: dict[str, str] = {}
    for jump_id, pair_key, *_rest in jump_info:
        compact_label = map_pair_key_to_jump_label(pair_key, alias_map)
        result[f'W_{compact_label}'] = f'W{jump_id}'
    return result


def mathematica_subscript_to_python(expr_str: str) -> str:
    """Convert Mathematica-style subscript notation to plain Python identifiers.

    Turns ``W_{OT}`` → ``W_OT``, ``W_{TO}`` → ``W_TO``, etc.

    Parameters
    ----------
    expr_str:
        Expression string that may contain ``W_{…}`` patterns.

    Returns
    -------
    str
        Expression string with curly-brace subscripts removed.
    """
    return re.sub(r'_\{([^{}]+)\}', r'_\1', expr_str)
