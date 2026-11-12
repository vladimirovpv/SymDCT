"""Shared YAML load/dump helpers for the public codebase.

This module centralises YAML I/O so that both the jump‑network builder
and the matrix‑builder stage can depend on a single implementation.  It
wraps :mod:`ruamel.yaml` and preserves the extended comment‑handling
behaviour provided by :mod:`src.common_utils.yaml_with_comments`.

Public entry points
-------------------

* :func:`load_data` – load a YAML file into a Python ``dict``.
* :func:`dump_data` – dump a mapping to YAML.

Both helpers expose a small, stable ``rw_yaml``-style signature so
call-sites can depend on it without touching :mod:`ruamel.yaml` directly.
"""

from __future__ import annotations

from typing import Any, Dict

from src.common_utils.yaml_with_comments import load_data as _jn_load_data
from src.common_utils.yaml_with_comments import dump_data as _jn_dump_data


def load_data(fname: str, dbg_print: bool = False) -> Dict[str, Any]:
    """Load YAML data from *fname*.

    Parameters
    ----------
    fname
        Path to the YAML file.
    dbg_print
        If ``True``, print a short debug representation after loading
        (primarily for interactive use).
    """

    data: Dict[str, Any] = _jn_load_data(fname)
    if dbg_print:
        # Lazy import to avoid circular logging setup.
        import pprint

        pprint.pp(data)
    return data


def dump_data(fname: str, data: Dict[str, Any]) -> None:
    """Dump *data* as YAML into *fname*.

    This is a thin wrapper around the Stage‑1 helper so that all public
    code can import YAML I/O from :mod:`src.common_utils.yaml_io` rather
    than reaching into the geometry package directly.
    """

    _jn_dump_data(fname, data)


__all__ = ["load_data", "dump_data"]
