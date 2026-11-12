"""Centralized I/O operations for energetics YAML files.

This module provides all operations for reading, writing, comparing,
and synchronizing energetics YAML files used by the diffusion coefficient
evaluator pipeline.
"""

from __future__ import annotations

import logging
import math
import os
from pathlib import Path
from typing import Any, Dict, List, Tuple

from src.common_utils.yaml_io import load_data, dump_data
from src.common_utils.yaml_with_comments import CommentedData


logger = logging.getLogger(__name__)


def load_energetics_yaml(fname: str) -> Dict[str, Any]:
    """Load and perform light validation of an energetics template YAML file.

    Parameters
    ----------
    fname:
        Path to the YAML file.

    Returns
    -------
    dict
        Parsed YAML content as a dictionary.

    Raises
    ------
    ValueError
        If required keys are missing or structure is invalid.
    """
    data = load_data(fname)
    if not isinstance(data, dict):
        raise ValueError(f"Energetics template {fname!r} did not yield a mapping")

    required_keys = {"settings", "site frequencies", "site energies", "jumps"}
    missing = required_keys - set(data.keys())
    if missing:
        raise ValueError(f"Energetics template {fname!r} is missing keys: {sorted(missing)}")

    # Basic structure checks for jumps: expect a mapping of pair_key -> list
    jumps = data.get("jumps", {})
    if not isinstance(jumps, dict):
        raise ValueError("'jumps' section must be a mapping of pair_key -> list of jump entries")

    for pair_key, entries in jumps.items():
        if not isinstance(entries, list):
            raise ValueError(f"Jump entry for {pair_key!r} must be a list")
        for e in entries:
            if not isinstance(e, dict):
                raise ValueError(f"Each jump entry under {pair_key!r} must be a mapping")
            if "jump_id" not in e or "length" not in e:
                raise ValueError(f"Jump entry under {pair_key!r} missing 'jump_id' or 'length': {e}")

    return data


def get_lattice_parameters(energetics: Dict[str, Any]) -> Tuple[Any, Any, Any]:
    """Extract lattice parameters (a, b, c) from an energetics mapping.

    Parameters
    ----------
    energetics:
        Dict loaded by load_energetics_yaml.

    Returns
    -------
    tuple
        (a, b, c) values from energetics["lattice constants"].
    """
    return tuple(energetics["lattice constants"][axis] for axis in ("a", "b", "c"))


def get_jump_info(energetics: Dict[str, Any]) -> List[Tuple]:
    """Extract per-jump migration energy and frequency data from an energetics dict.

    For every jump entry in the ``jumps`` section of *energetics* the function
    computes the migration energy E_m = E_sp - E_site for both forward and
    backward directions and collects the result as a flat list.

    Parameters
    ----------
    energetics:
        Dict loaded by load_energetics_yaml.

    Returns
    -------
    list of tuple
        Each element is a tuple (jump_id, pair_key, direction, length, e_m, frequency)
        where direction is "forward" or "backward".
    """
    site_energies = energetics["site energies"]
    site_labels = site_energies.keys()
    site_frequencies = energetics["site frequencies"]
    fixed_frequencies = {
        site: site_frequencies.get(site, site_frequencies["default_frequency"])
        for site in site_labels
    }
    jumps = energetics["jumps"]
    jump_info: list = []
    for pair_key, jump_entries in jumps.items():
        start, end = pair_key.split("-")
        for e in jump_entries:
            e_m_fwd = e["energy_sp"] - site_energies[start]
            e_m_bkw = e["energy_sp"] - site_energies[end]
            j_id = e["jump_id"]
            jump_info.append((
                j_id["forward"], pair_key, "forward",
                e["length"], e_m_fwd, fixed_frequencies[start],
            ))
            if start != end:
                jump_info.append((
                    j_id["backward"], f"{end}-{start}", "backward",
                    e["length"], e_m_bkw, fixed_frequencies[end],
                ))

    print("Migration energies (sorted by jump_id):")
    for jump_id, pk, _, jump_length, em, fq in sorted(jump_info):
        print(f"  {jump_id}: {pk} L={jump_length} Å Em={em:.3f} eV freq={fq} THz")

    return jump_info


def _collect_global_jumps(jumps_mapping: Dict[str, List[Dict[str, Any]]]) -> Dict[int, Dict[str, Any]]:
    """Return mapping global_idx -> {'pair': pair_key, 'length': float, 'entry': entry}

    jumps_mapping is YAML 'jumps' section: pair_key -> list of entries
    """
    result: Dict[int, Dict[str, Any]] = {}
    for pair_key, entries in jumps_mapping.items():
        for e in entries:
            jid = e.get("jump_id")
            if not isinstance(jid, dict):
                continue
            # prefer forward index as representative if available
            f = jid.get("forward")
            b = jid.get("backward")
            # take both indices — jumps may store different indices for directions
            for idx in (f, b):
                if idx is None:
                    continue
                result[int(idx)] = {"pair": pair_key, "length": float(e.get("length", math.nan)), "entry": e}
    return result


def compare_energetics(template_path: str, production_path: str, length_tolerance: float | None = None) -> Tuple[bool, Dict[str, Any]]:
    """Compare template and production energetics YAML.

    Checks:
      - same number of unique global jumps
      - same global indices map to same pair keys and similar lengths within tolerance

    Returns: (ok: bool, details: dict)
    """
    tmpl = load_data(template_path)
    prod = load_data(production_path)

    # determine tolerance: production settings override template if present
    tol = None
    if length_tolerance is not None:
        tol = float(length_tolerance)
    else:
        # look for settings.length_tolerance in either file
        for d in (prod, tmpl):
            try:
                st = d.get("settings", {})
                if st and st.get("length_tolerance") is not None:
                    tol = float(st.get("length_tolerance"))
                    break
            except Exception:
                pass
    if tol is None:
        tol = 1e-6

    # New comparison strategy: compare by (pair_key, length) rather than global indices.
    tmpl_jumps = tmpl.get("jumps", {})
    prod_jumps = prod.get("jumps", {})

    # Flatten into lists of (pair_key, length, entry)
    def _flatten(jumps_dict: Dict[str, List[Dict[str, Any]]]) -> List[Tuple[str, float, Dict[str, Any]]]:
        items = []
        for pair_key, entries in jumps_dict.items():
            for e in entries:
                try:
                    length = float(e.get("length", math.nan))
                except Exception:
                    length = math.nan
                items.append((pair_key, length, e))
        return items

    tmpl_list = _flatten(tmpl_jumps)
    prod_list = _flatten(prod_jumps)

    details: Dict[str, Any] = {
        "template_n_jumps": len(tmpl_list),
        "production_n_jumps": len(prod_list),
        "missing_template_entries": [],
        "extra_production_entries": [],
        "mismatches": [],
        "length_tolerance": tol,
    }

    # Match production entries to template entries by pair_key and length within tolerance.
    used_prod = [False] * len(prod_list)
    for i, (tpair, tlen, tentry) in enumerate(tmpl_list):
        matched = False
        for j, (ppair, plen, pentry) in enumerate(prod_list):
            if used_prod[j]:
                continue
            if tpair == ppair and abs(tlen - plen) <= tol:
                matched = True
                used_prod[j] = True
                break
        if not matched:
            details["missing_template_entries"].append({"pair": tpair, "length": tlen})

    # Any production entries not matched are extra
    for j, (ppair, plen, pentry) in enumerate(prod_list):
        if not used_prod[j]:
            details["extra_production_entries"].append({"pair": ppair, "length": plen})

    ok = len(details["missing_template_entries"]) == 0 and len(details["extra_production_entries"]) == 0
    return ok, details


def synchronize_energetics(
    template_path: str,
    production_path: str,
    out_dir: str | None = None,
    length_tolerance: float | None = None,
) -> Tuple[str, Dict[str, Any]]:
    """Produce a synchronized production YAML aligned to template.

    The new file will contain all jumps from the template with placeholder
    energies (null). Existing production values are preserved where possible.

    Returns: (new_path, report_dict)
    """
    tmpl = load_data(template_path)
    prod = load_data(production_path)

    # determine tolerance
    tol = None
    if length_tolerance is not None:
        tol = float(length_tolerance)
    else:
        for d in (prod, tmpl):
            try:
                st = d.get("settings", {})
                if st and st.get("length_tolerance") is not None:
                    tol = float(st.get("length_tolerance"))
                    break
            except Exception:
                pass
    if tol is None:
        tol = 1e-6

    # Build lookup for production jumps by (pair_key, length)
    prod_jumps = prod.get("jumps", {})
    prod_lookup: Dict[Tuple[str, float], Dict[str, Any]] = {}
    for pair_key, entries in prod_jumps.items():
        for e in entries:
            try:
                length = round(float(e.get("length", math.nan)), 6)
            except Exception:
                continue
            prod_lookup[(pair_key, length)] = e

    # Build new jumps section from template, copying matching production values
    tmpl_jumps = tmpl.get("jumps", {})
    new_jumps: Dict[str, List[Dict[str, Any]]] = {}
    copied = 0
    new_added = 0

    for pair_key, tmpl_entries in tmpl_jumps.items():
        new_entries: List[Dict[str, Any]] = []
        for te in tmpl_entries:
            try:
                tlen = round(float(te.get("length", math.nan)), 6)
            except Exception:
                tlen = math.nan

            key = (pair_key, tlen)
            if key in prod_lookup:
                # Copy production entry (including user-filled energies)
                new_entries.append(dict(prod_lookup[key]))
                copied += 1
            else:
                # Create placeholder entry
                jid = te.get("jump_id", {})
                new_entry = {
                    "jump_id": dict(jid) if isinstance(jid, dict) else {"forward": 0, "backward": 0},
                    "length": te.get("length"),
                    "energy_sp": CommentedData(None, commented=True),
                }
                new_entries.append(new_entry)
                new_added += 1
        new_jumps[pair_key] = new_entries

    # Assemble new production dict
    new_prod = dict(prod)
    new_prod["jumps"] = new_jumps

    # Determine output path
    if out_dir is None:
        prod_dir = os.path.dirname(production_path)
    else:
        prod_dir = out_dir

    # Count unique jump indices for filename
    unique_indices = set()
    for entries in new_jumps.values():
        for e in entries:
            jid = e.get("jump_id", {})
            for k in ("forward", "backward"):
                try:
                    unique_indices.add(int(jid.get(k)))
                except (TypeError, ValueError):
                    pass

    n_jumps = len(unique_indices)
    base = os.path.splitext(os.path.basename(production_path))[0]
    new_name = f"{base}_{n_jumps}freq.YAML"
    new_path = os.path.join(prod_dir, new_name)

    dump_data(new_path, new_prod)

    report = {
        "copied_entries": copied,
        "new_placeholder_entries": new_added,
        "n_jumps_new_file": n_jumps,
        "new_file": new_path,
    }
    return new_path, report


def ensure_energetics_in_sync(
    template: Path | None,
    production: Path,
    auto_sync_env: str = "DC_AUTO_SYNC",
    length_tolerance: float | None = None,
) -> Path:
    """Ensure template and production energetics are compatible.

    If *template* is ``None`` the function will attempt to locate it as
    ``production.parent.parent / "output" / "energetics_template.yaml"``.

    If a mismatch is detected and the environment variable named by
    *auto_sync_env* is set to ``"1"``, the function will call
    :func:`synchronize_energetics` and return the path to the newly written
    synchronized production file.  Otherwise it logs a warning and returns
    the original *production* path unchanged so the caller can decide how to
    proceed.

    Parameters
    ----------
    template:
        Path to the geometry-stage template YAML, or ``None`` to auto-locate.
    production:
        Path to the filled-in production energetics YAML.
    auto_sync_env:
        Name of the environment variable that enables automatic sync
        (default: ``"DC_AUTO_SYNC"``).
    length_tolerance:
        Override the length tolerance used when comparing jumps.

    Returns
    -------
    Path
        *production* (possibly replaced by the newly synchronized file).
    """
    if template is None:
        template = production.parent.parent / "output" / "energetics_template.yaml"

    if not template.exists():
        logger.info("Template not found at %s; skipping template/production check.", template)
        return production

    ok, details = compare_energetics(str(template), str(production), length_tolerance=length_tolerance)
    if not ok:
        if os.getenv(auto_sync_env) == "1":
            logger.info("Template/production mismatch detected. Attempting automatic synchronization (%s=1)...", auto_sync_env)
            try:
                new_path, report = synchronize_energetics(
                    str(template), str(production),
                    out_dir=str(production.parent),
                    length_tolerance=length_tolerance,
                )
                logger.info("Synchronization written to: %s", new_path)
                logger.debug("Sync report: %s", report)
                return Path(new_path)
            except Exception as exc:
                logger.error("Automatic synchronization failed: %s", exc)
                logger.error("Please run the synchronizer manually and fix placeholders. See src/common_utils/README.md for details.")
                raise SystemExit(1) from exc
        else:
            if not length_tolerance:
                tol = 1e-6
            missing = details.get("missing_template_entries", [])
            extra = details.get("extra_production_entries", [])
            logger.warning(
                f"Energetics template/production mismatch (tolerance {tol:.4e} A)"
                f"\nmissing from template: {[' @ '.join((jump['pair'],f"{jump['length']:.4f}")) for jump in missing]}"
                f"\nextra in production: {[' @ '.join((jump['pair'],f"{jump['length']:.4f}")) for jump in extra]}"
                "\nSet %s=1 to auto-fix.\n",
                auto_sync_env,
            )
    return production
