"""Lightweight POSCAR loader for jump-network construction.

This module provides :class:`PoscarJumpStructureLoader`, which encapsulates
just the responsibilities needed by the jump‑network builder:

- Read a VASP POSCAR/CONTCAR file using pymatgen.
- Optionally separate *matrix* atoms from *jump* (interstitial) sites
  using a list of species symbols that define the matrix.
- Expose ``matrix_structure``, ``jump_structure`` and the number of
  jump sites ``nsites``.

Implemented here in isolation so that this package can be used as a
standalone, lightweight component.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional

from pymatgen.core import Structure
from pymatgen.io.vasp import Poscar


@dataclass
class PoscarJumpStructureLoader:
    """Load matrix and jump structures from a POSCAR file.

    Parameters
    ----------
    filename
        Path to a VASP POSCAR/CONTCAR file.
    crystal_matrix
        Optional list of chemical symbols (e.g. ["Ti", "Be"]) that
        define the *matrix* species. If given, atoms of these species
        are removed from the *jump* structure, and all other species are
        removed from the *matrix* structure. If omitted, the matrix and
        jump structures are identical copies of the input.
    """

    filename: str
    crystal_matrix: Optional[List[str]] = None

    def __post_init__(self) -> None:
        self._matrix_structure: Structure | None = None
        self._jump_structure: Structure | None = None
        self._raw_symbol_alias_map: dict[str, str] = {}
        self._nsites: int = 0
        self._load()

    # ------------------------------------------------------------------
    # internal helpers
    # ------------------------------------------------------------------
    def _read_raw_symbol_alias_map(self, canonical_symbols: list[str]) -> dict[str, str]:
        """Read original POSCAR species symbols and map canonical symbols back."""

        path = Path(self.filename)
        lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
        if len(lines) < 7:
            return {}

        symbol_tokens = lines[5].split()
        if not symbol_tokens or len(symbol_tokens) != len(canonical_symbols):
            return {}

        return {
            str(canonical): str(raw)
            for raw, canonical in zip(symbol_tokens, canonical_symbols)
        }

    def _load(self) -> None:
        """Read POSCAR and build matrix / jump structures.

        This is intentionally minimal: no logging side effects, no
        process termination – errors are surfaced as Python exceptions.
        """

        poscar = Poscar.from_file(
            self.filename,
            check_for_potcar=False,
            read_velocities=False,
            true_names=False,
        )

        base_struct = poscar.structure
        self._raw_symbol_alias_map = self._read_raw_symbol_alias_map(list(poscar.site_symbols))

        matrix_struct = base_struct.copy()
        jump_struct = base_struct.copy()

        if self.crystal_matrix:
            # Remove all matrix species from the jump structure so that
            # it contains only interstitial / jump sites.
            jump_struct.remove_species(self.crystal_matrix)

            # In the matrix structure keep *only* matrix species.
            non_matrix_species = [
                sp.symbol for sp in matrix_struct.species if sp.symbol not in self.crystal_matrix
            ]
            if non_matrix_species:
                matrix_struct.remove_species(non_matrix_species)

        self._matrix_structure = matrix_struct if matrix_struct.num_sites > 0 else None
        self._jump_structure = jump_struct if jump_struct.num_sites > 0 else None
        # intentionally removed: use structure.num_sites instead
        # self._nsites = len(jump_struct)

    # ------------------------------------------------------------------
    # public API
    # ------------------------------------------------------------------
    @property
    def matrix_structure(self) -> Structure:
        """Structure containing only matrix atoms (if ``crystal_matrix``
        was given) or a copy of the full structure otherwise."""

        assert self._matrix_structure is not None
        return self._matrix_structure

    @property
    def jump_structure(self) -> Structure:
        """Structure containing only jump / interstitial sites.

        If ``crystal_matrix`` was omitted, this is a copy of the full
        input structure.
        """

        assert self._jump_structure is not None
        return self._jump_structure

    @property
    def raw_symbol_alias_map(self) -> dict[str, str]:
        """Map canonical species symbols back to the original POSCAR aliases."""

        return dict(self._raw_symbol_alias_map)

    # this property is removed intentionally as both _matrix_structure and _jump_structure
    # have their own attribute .num_sites. These should be used instead of the .nsites.
    # @property
    # def nsites(self) -> int:
    #     """Number of sites in :attr:`jump_structure`."""
    #
    #     return self._nsites
