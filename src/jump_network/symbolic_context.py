"""Symbol registry for jump-network symbolic calculations."""

from __future__ import annotations

from typing import List

import sympy as sp


class SymbolicContext:
    """Provide k-vectors, time symbol, substitutions, and tolerance."""

    def __init__(self, tolerance: float = 1e-6):
        self.tolerance = tolerance
        self._k: List[sp.Symbol] = [sp.symbols(f"k{i+1}", real=True) for i in range(3)]
        self._t: sp.Symbol = sp.Symbol("t", real=True, positive=True)
        self._subst_w: dict = {}

    @property
    def k(self) -> List[sp.Symbol]:
        return self._k

    @property
    def t(self) -> sp.Symbol:
        return self._t

    def get_tolerance(self) -> float:
        return self.tolerance

    def set_substitutions(self, subst_w: dict) -> None:
        self._subst_w = dict(subst_w) if subst_w else {}

    def get_substitutions(self, replace: bool = True) -> dict:
        return self._subst_w if replace else {}

