"""Read a crystal structure from any file format that pymatgen can parse.

:func:`read_structure` is the single entry point used by the loader. The format
is chosen from the file name:

1. Names that contain ``POSCAR`` or ``CONTCAR`` (any case) or end in ``.vasp``
   are read as VASP files. This path keeps the species symbols exactly as
   written, so that aliases such as ``He1`` survive (see ``poscar_aliases``).
2. All other names go to :meth:`pymatgen.core.Structure.from_file`. It covers
   CIF, mCIF, CSSR, XSF, JSON/YAML (pymatgen serialisation), vasprun.xml,
   CHGCAR/LOCPOT, ATAT (rndstr.in, lat.in, bestsqs), FHI-aims ``geometry.in``,
   Fleur, exciting, LMTO, SHELX ``.res``, PWmat ``.config`` and netCDF.
3. If pymatgen does not know the extension, ASE is tried when it is installed
   (this adds extxyz, xyz, cell, gen and the other ASE formats). Last, the file
   is read as a POSCAR, which was the only format of earlier releases.

The key ``structure_format`` of the input YAML overrides the choice by name. It
accepts ``poscar``, ``ase`` or any format string of
:meth:`pymatgen.core.Structure.from_str`.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from pymatgen.core import DummySpecies, Species, Structure
from pymatgen.io.vasp import Poscar


def _is_poscar_name(name: str) -> bool:
    low = name.lower()
    return "poscar" in low or "contcar" in low or low.endswith(".vasp")


def _read_poscar(filename: str) -> tuple[Structure, dict[str, str]]:
    poscar = Poscar.from_file(
        filename, check_for_potcar=False, read_velocities=False, true_names=False
    )
    return poscar.structure, _poscar_aliases(filename, list(poscar.site_symbols))


def _poscar_aliases(filename: str, canonical: list[str]) -> dict[str, str]:
    """Map canonical symbols back to the species tokens written in the file."""
    lines = [ln.strip() for ln in Path(filename).read_text(encoding="utf-8").splitlines() if ln.strip()]
    if len(lines) < 7:
        return {}
    tokens = lines[5].split()
    if not tokens or len(tokens) != len(canonical):
        return {}
    return {str(c): str(r) for r, c in zip(tokens, canonical)}


def _read_ase(filename: str) -> Structure:
    try:
        from ase.io import read as ase_read
        from pymatgen.io.ase import AseAtomsAdaptor
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise ImportError(
            "Reading this file format needs the optional package 'ase' (pip install ase)."
        ) from exc
    return AseAtomsAdaptor.get_structure(ase_read(filename))


def read_structure(filename: str, fmt: Optional[str] = None) -> tuple[Structure, dict[str, str]]:
    """Return ``(structure, aliases)`` for *filename*.

    ``aliases`` maps canonical element symbols to the original species tokens
    of a POSCAR. It is empty for all other formats.
    """
    path = Path(filename)
    if not path.is_file():
        raise FileNotFoundError(f"Structure file '{filename}' not found.")

    aliases: dict[str, str] = {}
    fmt_low = fmt.lower() if fmt else None

    if fmt_low == "poscar" or (fmt_low is None and _is_poscar_name(path.name)):
        struct, aliases = _read_poscar(str(path))
    elif fmt_low == "ase":
        struct = _read_ase(str(path))
    elif fmt_low is not None:
        struct = Structure.from_str(path.read_text(encoding="utf-8", errors="replace"), fmt=fmt_low)
    else:
        try:
            struct = Structure.from_file(str(path))
        except ValueError as exc:
            if "Unrecognized extension" not in str(exc):
                raise
            try:
                struct = _read_ase(str(path))
            except Exception:
                try:
                    struct, aliases = _read_poscar(str(path))
                except Exception:
                    raise ValueError(
                        f"Cannot read '{filename}': the extension is unknown to pymatgen, "
                        "and the file is neither readable by ASE nor a POSCAR. "
                        "Set 'structure_format' in the input YAML."
                    ) from exc

    if not isinstance(struct, Structure):
        raise ValueError(f"'{filename}' does not contain a periodic crystal structure.")
    if not struct.is_ordered:
        raise ValueError(
            f"'{filename}' contains partially occupied sites. "
            "Provide an ordered structure, with one species per site."
        )
    charged = {
        sp: sp.element
        for sp in struct.types_of_species
        if isinstance(sp, Species) and not isinstance(sp, DummySpecies)
    }
    if charged:
        struct.replace_species(charged)
    return struct, aliases
