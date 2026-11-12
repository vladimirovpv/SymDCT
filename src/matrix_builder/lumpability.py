"""Lumpability check for merging equivalent sites in the Fourier-space matrix.

Merging sites into classes with ``A' = R A R^+`` is exact only if the
class partition is *lumpable* (equitable) for the wave vector ``k = hkl * t``
used for the transport direction, i.e. if

    R A(k) = A'(k) R ,

where ``R`` is the class indicator matrix.  In components, for every site
``j`` of a class ``d`` the sum ``sum_{i in c} A_ij(k)`` over the rows of any
class ``c`` must be the same for all ``j`` in ``d``.  If this holds, the
probabilities summed over classes obey a closed master equation and the
reduced matrix keeps the eigenvalue ``lambda_1(k)`` exactly.  If it does not
hold, the reduction assumes equal occupation inside a class, which is
correct only at ``k = 0``.

The check is symbolic: it forms ``R A`` (a product with rows of ones that
adds up the elements of each column within a class) and tests that the
difference of any two columns of one class is identically zero.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Sequence

import numpy as np
import sympy as sp


@dataclass
class LumpabilityReport:
    """Result of :func:`check_lumpability`."""

    ok: bool
    # (source class d, target class c, site j0, site j, residual)
    # ``residual`` is a symbolic expression (method="symbolic") or a
    # float (method="numeric").
    violations: list[tuple] = field(default_factory=list)

    def __bool__(self) -> bool:
        return self.ok

    def summary(self) -> str:
        if self.ok:
            return "lumpable"
        d, c, j0, j, res = self.violations[0]
        return (
            f"NOT lumpable: sites {j0} and {j} of class {d} differ in their total "
            f"rate into class {c} by {res}"
        )


def _is_zero(expr: sp.Expr) -> tuple[bool, sp.Expr]:
    """Return (is identically zero, residual).

    The residual is returned in expanded exponential form, which is the form
    that :func:`suggest_rate_relations` needs.
    """

    if expr == 0:
        return True, sp.Integer(0)
    # Cheap canonical test: exponentials with identical arguments cancel
    # automatically after expansion.
    cheap = sp.expand(expr.rewrite(sp.exp))
    if cheap == 0:
        return True, sp.Integer(0)
    # Fall back to the full simplifier before declaring a violation.
    if sp.simplify(cheap) == 0:
        return True, sp.Integer(0)
    return False, cheap


def check_lumpability(
    equiv_sites: Sequence[int],
    full_matrix: sp.Matrix,
    hkl_structure: Sequence[float],
    k_symbols: Sequence[sp.Symbol],
    t_symbol: sp.Symbol,
    method: str = "symbolic",
    n_samples: int = 3,
    tol: float = 1e-9,
    seed: int = 0,
) -> LumpabilityReport:
    """Check that merging ``equiv_sites`` is exact along ``k = hkl * t``.

    ``method="symbolic"`` (default) multiplies the matrix from the left by
    the class indicator rows and requires the column differences within a
    class to vanish identically.  ``method="numeric"`` evaluates the same
    condition at random points and is kept as an independent cross-check.
    """

    equiv = [int(v) for v in equiv_sites]
    n = len(equiv)
    if full_matrix.shape != (n, n):
        raise ValueError(f"full_matrix has shape {full_matrix.shape}, expected {(n, n)}")
    if method not in ("symbolic", "numeric"):
        raise ValueError("method must be 'symbolic' or 'numeric'")

    subst = {kk: sp.nsimplify(hkl_structure[i]) * t_symbol for i, kk in enumerate(k_symbols)}
    matrix = full_matrix.subs(subst)

    classes = sorted(set(equiv))
    members = {c: [i for i, e in enumerate(equiv) if e == c] for c in classes}
    violations: list[tuple] = []

    if method == "symbolic":
        indicator = sp.Matrix(
            [[1 if equiv[j] == c else 0 for j in range(n)] for c in classes]
        )
        ra = indicator * matrix  # entry (c, j): total rate of site j into class c
        for d in classes:
            cols = members[d]
            j0 = cols[0]
            for j in cols[1:]:
                for ci, c in enumerate(classes):
                    zero, residual = _is_zero(ra[ci, j] - ra[ci, j0])
                    if not zero:
                        violations.append((d, c, j0, j, residual))
        return LumpabilityReport(ok=not violations, violations=violations)

    symbols = sorted(matrix.free_symbols, key=lambda s: s.name)
    func = sp.lambdify(symbols, matrix, modules="numpy")
    rng = np.random.default_rng(seed)
    for _ in range(n_samples):
        values = rng.uniform(0.5, 2.0, size=len(symbols))
        a_num = np.array(func(*values), dtype=complex)
        for d in classes:
            cols = members[d]
            j0 = cols[0]
            for c in classes:
                rows = np.array(members[c])
                ref = a_num[rows, j0].sum()
                for j in cols[1:]:
                    dev = abs(a_num[rows, j].sum() - ref)
                    if dev > tol:
                        violations.append((d, c, j0, j, float(dev)))
    return LumpabilityReport(ok=not violations, violations=violations)


def _substituted(full_matrix, hkl_structure, k_symbols, t_symbol) -> sp.Matrix:
    subst = {kk: sp.nsimplify(hkl_structure[i]) * t_symbol for i, kk in enumerate(k_symbols)}
    return full_matrix.subs(subst)


def refine_to_lumpable(
    equiv_sites: Sequence[int],
    full_matrix: sp.Matrix,
    hkl_structure: Sequence[float],
    k_symbols: Sequence[sp.Symbol],
    t_symbol: sp.Symbol,
) -> list[int]:
    """Split classes until the partition is lumpable (coarsest refinement).

    Two sites of one class stay together only if their columns of ``R A`` are
    identically equal for the current partition.  The split is repeated until
    nothing changes.  The returned class labels are renumbered by the position
    of the first member (input order of ``equiv_sites`` is preserved as far as
    possible).
    """

    matrix = _substituted(full_matrix, hkl_structure, k_symbols, t_symbol)
    n = matrix.shape[0]
    labels = [int(v) for v in equiv_sites]

    while True:
        classes = sorted(set(labels))
        indicator = sp.Matrix([[1 if labels[j] == c else 0 for j in range(n)] for c in classes])
        ra = indicator * matrix
        new_labels = list(labels)
        next_label = max(labels) + 1
        changed = False
        for d in classes:
            members = [j for j in range(n) if labels[j] == d]
            groups: list[list[int]] = []
            for j in members:
                for grp in groups:
                    j0 = grp[0]
                    if all(_is_zero(ra[ci, j] - ra[ci, j0])[0] for ci in range(len(classes))):
                        grp.append(j)
                        break
                else:
                    groups.append([j])
            for grp in groups[1:]:
                changed = True
                for j in grp:
                    new_labels[j] = next_label
                next_label += 1
        labels = new_labels
        if not changed:
            break

    # Renumber classes in the order in which their first member appears
    # when the input classes are walked in increasing class id.
    order = {}
    for _, j in sorted((equiv_sites[j], j) for j in range(n)):
        order.setdefault(labels[j], len(order))
    return [order[l] for l in labels]


def suggest_rate_relations(
    report: LumpabilityReport,
    rate_prefix: str = "W",
) -> list[sp.Eq]:
    """Linear relations between rates that would make the partition lumpable.

    Every symbolic residual of :func:`check_lumpability` is a sum of
    exponentials in ``t`` whose coefficients are linear in the jump rates.
    All coefficients must vanish; this gives linear equations, which are
    solved for the rate symbols whose name starts with ``rate_prefix``
    (symbols with larger index are expressed through those with smaller
    index).  An empty list is returned if the partition is already lumpable
    or if the relations cannot be expressed through rate symbols.
    """

    equations: list[sp.Expr] = []
    for *_, residual in report.violations:
        if not isinstance(residual, sp.Basic):
            continue
        expr = sp.expand(residual)
        exps = sorted(expr.atoms(sp.exp), key=str)
        dummies = sp.symbols(f"_e0:{len(exps)}") if exps else ()
        if len(exps) == 1:
            dummies = (dummies,) if not isinstance(dummies, tuple) else dummies
        replaced = expr.xreplace(dict(zip(exps, dummies)))
        if dummies:
            poly = sp.Poly(replaced, *dummies)
            equations.extend(sp.expand(c) for c in poly.coeffs())
        else:
            equations.append(replaced)

    unknowns = sorted(
        {s for e in equations for s in e.free_symbols if s.name.startswith(rate_prefix)},
        key=lambda s: _natural_key(s.name),
        reverse=True,
    )
    if not equations or not unknowns:
        return []
    solutions = sp.solve(equations, unknowns, dict=True)
    if not solutions:
        return []
    return [sp.Eq(lhs, rhs) for lhs, rhs in sorted(solutions[0].items(), key=lambda kv: _natural_key(kv[0].name))]


def _natural_key(name: str):
    import re

    return [int(p) if p.isdigit() else p for p in re.split(r"(\d+)", name)]
