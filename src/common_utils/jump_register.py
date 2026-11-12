"""Common JumpRegister shared between jump_network_builder and matrix_builder.

This module hosts the lightweight ``JumpRegister`` implementation so it
can be imported as ``from src.common_utils.jump_register import JumpRegister``
without polluting the top-level :mod:`src` namespace.

The implementation is based on ``src/JumpRegister.py`` and kept
API‑compatible for its current usages.
"""

from __future__ import annotations

from dataclasses import dataclass
import logging
from typing import Iterable, List

import numpy as np
from sympy import Matrix, Rational, fraction, nsimplify


log = logging.getLogger(__name__)


def custom_isclose(a, b, rtol: float = 1e-5, atol: float = 1e-8, equal_nan: bool = False) -> bool:
    """Numerical closeness helper that also works for SymPy Rationals."""

    if isinstance(a, Rational) and isinstance(b, Rational):
        return abs(a - b) <= max(rtol * max(abs(a), abs(b)), atol)
    else:
        return np.isclose(a, b, rtol=rtol, atol=atol, equal_nan=equal_nan)


def custom_allclose(a: Iterable, b: Iterable, rtol: float = 1e-5, atol: float = 1e-8, equal_nan: bool = False) -> bool:
    return all(custom_isclose(x, y, rtol=rtol, atol=atol, equal_nan=equal_nan) for x, y in zip(a, b))


@dataclass(frozen=True)
class Jump:
    """Descriptor for a jump type between two sublattices."""

    sub_a: str = ""
    sub_b: str = ""
    jump_length: float = 0.0


@dataclass(frozen=False)
class JumpRegInfo:
    """Bookkeeping information for a registered jump type."""

    global_idx: int = 0
    local_idx: int = 0


class JumpRegister:
    """Register and group jumps between sublattices.

    The register answers two main questions:

    * Given a jump (A -> B, length h, vector, image), what is its global
      jump index?  (:meth:`register`).
    * How many distinct jump types have been seen so far?
      (:meth:`get_number_of_jumps`).

    It also stores the list of jump vectors associated with each jump
    type which is useful for diagnostic printing.
    """

    def __init__(self, tol: float = 0.008, rounding_precision: int = 3) -> None:
        # tolerance for distance comparison
        self._tol = tol
        self._rounding_precision = rounding_precision
        # mapping Jump(sub_a, sub_b, length) -> (JumpRegInfo, [vectors], [images])
        self._jumps: dict[Jump, tuple[JumpRegInfo, List[Matrix], List[List[int]]]] = {}
        # current maximal jump count, increased when a new jump type is found
        self._jump_cnt: int = 0
        # list of global indices of jumps added to a specific network
        self._added: list[int] = []
        # flag reserved for future support of symmetric transition matrices
        self._symmetric_matrix: bool = False

    # ------------------------------------------------------------------
    # Basic queries
    # ------------------------------------------------------------------
    def get_symmetric_matrix(self) -> bool:
        return self._symmetric_matrix

    def get_number_of_jumps(self) -> int:
        """Return the number of distinct jump types registered so far."""

        return self._jump_cnt

    # ------------------------------------------------------------------
    # Registration and bookkeeping
    # ------------------------------------------------------------------
    def _normalise_length(self, jump_length):
        """Normalise jump length so nearby floats collapse to one value."""

        if isinstance(jump_length, float):
            # check if a close-by length was already registered
            for k in self._jumps.keys():
                if abs(k.jump_length - jump_length) < self._tol:
                    jump_length = k.jump_length
                    break
            jump_length = round(jump_length, self._rounding_precision)
        elif isinstance(jump_length, Rational):
            nomi, denom = fraction(jump_length)
            log.debug("Rational jump length: numerator=%s, denominator=%s", nomi, denom)
        else:
            log.warning("jump_length is of an unexpected type: %r", type(jump_length))
        return jump_length

    def register(
        self,
        sub_a: str,
        sub_b: str,
        jump_length: float,
        jump_vector,
        image: List[int],
    ) -> int:
        """Register a jump between ``sub_a`` and ``sub_b`` and return its global index."""

        jump_length = self._normalise_length(jump_length)

        # Convert jump_vector to a SymPy Matrix for consistent handling,
        # but keep the original iterable shape for comparisons.
        if not isinstance(jump_vector, Matrix):
            try:
                jump_vector = Matrix(jump_vector)
            except Exception:  # pragma: no cover - defensive
                log.error("Could not convert jump_vector %r to Matrix", jump_vector)
                jump_vector = Matrix(jump_vector)

        one_jump = Jump(sub_a, sub_b, jump_length)

        if one_jump in self._jumps:
            jump_info, jump_vectors, images = self._jumps[one_jump]
            if not any(custom_allclose(jump_vector, v, rtol=0.0, atol=self._tol) for v in jump_vectors):
                jump_vectors.append(jump_vector)
                images.append(list(map(int, image)))
            jump_ind = jump_info.global_idx
            local_cnt = jump_info.local_idx
        else:
            jump_ind = self._jump_cnt
            # count how many jumps between (sub_a, sub_b) exist already
            local_cnt = self.num_jumps(one_jump.sub_a, one_jump.sub_b)
            self._jumps[one_jump] = (
                JumpRegInfo(jump_ind, local_cnt),
                [jump_vector],
                [list(map(int, image))],
            )
            if self._symmetric_matrix:
                reversed_jump = Jump(sub_b, sub_a, jump_length)
                self._jumps[reversed_jump] = (
                    JumpRegInfo(jump_ind, local_cnt),
                    [(-jump_vector)],
                    [(-1 * Matrix(image))],
                )
            self._jump_cnt += 1

            log.debug(
                "Jump length %9.6f %s->%s global=%d local=%d image=%s",
                float(one_jump.jump_length),
                sub_a,
                sub_b,
                jump_ind,
                local_cnt,
                image,
            )

        return jump_ind

    # ------------------------------------------------------------------
    # Convenience helpers used by the builder
    # ------------------------------------------------------------------
    def num_jumps(self, sub_a: str, sub_b: str) -> int:
        """Return how many jump types connect ``sub_a`` and ``sub_b``."""

        return sum(1 for a_jump in self._jumps if (a_jump.sub_a, a_jump.sub_b) == (sub_a, sub_b))

    def add_to_netwk(self, jump_count: int) -> None:
        """Mark a given global jump index as present in the network."""

        self._added.append(jump_count)

    def added_to_netwk(self, jump_count: int) -> bool:
        """Return True if the given jump index has been added to the network."""

        return jump_count in self._added

    def print_jumps(self, compact: bool = True) -> None:
        """Print registered jump statistics (compatibility method)."""

        print(f"Total jump types: {self._jump_cnt}")
        if not compact:
            for jump, (info, vectors, images) in self._jumps.items():
                print(f"  {jump.sub_a}->{jump.sub_b} (length={jump.jump_length}): "
                      f"global={info.global_idx}, local={info.local_idx}, "
                      f"vectors={len(vectors)}")

    def map_tprob(self, w_dict: dict) -> dict:
        """Map jump probabilities using the register's global indices."""
        return dict(w_dict)

