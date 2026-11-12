"""Transition-matrix builder for jump networks."""

from __future__ import annotations

import sys
from typing import Any
import time

import sympy as sp
from sympy.matrices import zeros
from tqdm import tqdm

from .core import JumpNetworkGraph
from .symbolic_context import SymbolicContext

THRESHOLD_ROWS = 20   # below this, just print a one-liner instead of a progress bar


class TransitionMatrixBuilder:
    """Fill symbolic transition-probability matrix from graph edges."""

    def __init__(
        self,
        graph: JumpNetworkGraph,
        symbolic_ctx: SymbolicContext,
        jump_probabilities: Any,
    ):
        self.graph = graph
        self.symbolic_ctx = symbolic_ctx
        self.jump_probabilities = jump_probabilities
        self._matrix = None

    def fill_matrix(self, subst_w: dict | None = None) -> sp.Matrix:
        if self._matrix is not None:
            return self._matrix

        subst_w = subst_w or {}
        self.symbolic_ctx.set_substitutions(subst_w)

        num_sites = self.graph.num_sites
        k = self.symbolic_ctx.k
        w = self.jump_probabilities.get_jump_probabilities()
        tol = self.symbolic_ctx.get_tolerance()

        matrix = zeros(num_sites, num_sites)

        show_bar = num_sites >= THRESHOLD_ROWS
        t0 = time.perf_counter()
        row_iter = tqdm(range(num_sites),
                        desc="\nFilling transition matrix",
                        disable= not show_bar,
                        file=sys.stdout)

        if not show_bar:
            print(f"\nFilling transition matrix ({num_sites}x{num_sites})...", end="")
        else:
            print("\n")

        for irow in row_iter:
            outgoing = 0
            for start, end, key, _ in self.graph.get_out_edges(irow):
                edge_data = self.graph.get_edge_data(start, end, key)
                outgoing += w[edge_data["glob_jump_cnt"]]

            matrix[irow, irow] -= outgoing

            for start_j, end_i, key, h in self.graph.get_in_edges(irow):
                h_sym = sp.Matrix([sp.nsimplify(h_el, tolerance=tol, rational=True) for h_el in h])
                pot = h_sym.T.dot(sp.Matrix(k))
                edge_data = self.graph.get_edge_data(start_j, end_i, key)
                global_idx = edge_data["glob_jump_cnt"]
                matrix[irow, start_j] += w[global_idx] * sp.exp(sp.I * pot)

        if not show_bar:
            print(f" done in {time.perf_counter() - t0:.2f} s")

        self._matrix = matrix
        return self._matrix

    def get_matrix(self) -> sp.Matrix | None:
        return self._matrix

    def reset(self) -> None:
        self._matrix = None

