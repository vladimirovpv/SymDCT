"""Utility functions for normalising ``tweak_transition_rates``.

The original implementation lived as a large private method inside
``MatrixStageJumpNetwork``.  It performed three tasks:

1. Build a map from canonical jump labels (e.g. ``"A->B"``) to a list of
   records containing the global jump index and the jump length.
2. Parse user‑provided keys and values which may be integers, label strings,
   or label strings with a ``:length=`` qualifier.  The ``:all`` suffix forces
   the operation to apply to *all* matching jumps.
3. Emit warnings for ambiguous or missing specifications.

In the refactor we extracted this logic to a dedicated module so that it can
be unit‑tested independently and reused elsewhere.  The behaviour is unchanged
except that warning messages are now emitted via the :mod:`logging` module
instead of ``print``.
"""

from __future__ import annotations

import logging
from typing import Dict, List, Tuple, Union

# The record format used throughout the original code.
_Record = Dict[str, Union[int, float]]


def _parse_label_spec(spec: str) -> Tuple[str, Union[float, None], bool]:
    """Parse a label specification.

    Supported forms are ``"A->B"``, ``"A->B:length=1.23"``, ``"A->B:1.23"``
    (a bare number is shorthand for ``length=``) and ``"A->B:all"``. Returns
    a tuple ``(label, length, all_flag)`` where ``length`` is ``None`` unless
    a length qualifier is present, and ``all_flag`` is ``True`` when the
    ``:all`` suffix is used.
    """
    if ":" not in spec:
        return spec, None, False

    label, qualifier = spec.split(":", 1)
    label = label.strip()
    qualifier = qualifier.strip().lower()

    if qualifier == "all":
        return label, None, True

    if qualifier.startswith("length=") or qualifier.startswith("len=") or qualifier.startswith("l="):
        try:
            length = float(qualifier.split("=", 1)[1])
            return label, length, False
        except Exception:
            # Fall through to treat as an unrecognised spec – the caller will
            # handle the warning.
            return label, None, False

    # Bare number shorthand for length=<value>, e.g. "A->B:2.213".
    try:
        return label, float(qualifier), False
    except ValueError:
        return label, None, False


def normalize_tweak_rates(
    tweak_rates: Dict[object, object],
    records: Dict[str, List[_Record]],
    tolerance: float,
) -> Dict[int, Union[int, float, str]]:
    """Normalize ``tweak_transition_rates`` to a mapping of global indices.

    Parameters
    ----------
    tweak_rates:
        The user‑provided mapping where keys and values may be integers (global
        jump indices), label strings, or label strings with a length qualifier.
    records:
        Mapping from canonical label to a list of records ``{'glob': int,
        'length': float}`` as built by ``MatrixStageJumpNetwork._build_label_records_map``.
    tolerance:
        Numerical tolerance used when matching lengths.

    Returns
    -------
    dict[int, Union[int, float, str]]
        Normalised mapping where each key is a global jump index and the value
        is the target (global index, numeric value, or original string if it
        could not be resolved).
    """
    if not tweak_rates:
        return {}

    normalized: Dict[int, Union[int, float, str]] = {}

    for raw_k, raw_v in dict(tweak_rates).items():
        # Resolve source indices
        src_indices: List[int] = []
        try:
            src_indices = [int(raw_k)]
        except Exception:
            k_label, k_length, k_all = _parse_label_spec(str(raw_k))
            matches = records.get(k_label, [])
            if not matches:
                logging.warning("tweak_transition_rates source label '%s' did not match any jumps.", raw_k)
                continue
            if k_length is not None:
                matched = [r["glob"] for r in matches if abs(r["length"] - k_length) <= tolerance]
                if not matched:
                    logging.warning("no jump with length ~%s for label '%s'", k_length, k_label)
                    continue
                src_indices = matched
            else:
                # No length qualifier – either a single match or apply to all
                if len(matches) == 1 and not k_all:
                    src_indices = [matches[0]["glob"]]
                else:
                    if len(matches) > 1 and not k_all:
                        logging.warning(
                            "label '%s' matches %d jumps; applying tweak to ALL matching jumps.",
                            k_label,
                            len(matches),
                        )
                    src_indices = [r["glob"] for r in matches]

        # Resolve target value
        resolved_val: Union[int, float, str, None] = None
        try:
            resolved_val = int(raw_v)
        except Exception:
            if isinstance(raw_v, str):
                v_label, v_length, v_all = _parse_label_spec(raw_v)
                v_matches = records.get(v_label, [])
                if not v_matches:
                    # Try to interpret as a raw float literal
                    try:
                        resolved_val = float(raw_v)
                    except Exception:
                        logging.warning(
                            "tweak_transition_rates target '%s' did not match jumps and is not numeric.",
                            raw_v,
                        )
                        resolved_val = raw_v
                else:
                    if v_length is not None:
                        found = [r["glob"] for r in v_matches if abs(r["length"] - v_length) <= tolerance]
                        if found:
                            resolved_val = found[0]
                        else:
                            logging.warning(
                                "target label '%s' matched no jump by length.", raw_v,
                            )
                            resolved_val = raw_v
                    else:
                        if len(v_matches) > 1:
                            logging.warning(
                                "target label '%s' matches %d jumps; using first match.", v_label, len(v_matches)
                            )
                        resolved_val = v_matches[0]["glob"]
            else:
                # Pass through numeric or other types unchanged
                resolved_val = raw_v

        for sidx in src_indices:
            normalized[sidx] = resolved_val  # type: ignore[assignment]

    return normalized
