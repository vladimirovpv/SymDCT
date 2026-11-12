"""Helpers for axis-dependent transition-matrix reduction."""

from __future__ import annotations

import logging
from typing import Any, Sequence

import sympy as sp
from sympy import MutableDenseMatrix

from src.symmetry_context.node_labeler import (
    build_non_equivalent_site_label_data,
    print_non_equivalent_site_coordinates,
)


def _ordered_unique(values: Sequence[int]) -> list[int]:
    """Return first-seen unique values, preserving order."""

    return list(dict.fromkeys(int(v) for v in values))


def representative_site_indices_from_equiv(equiv_sites: Sequence[int]) -> list[int]:
    """Return one representative site index for each equivalence class."""

    class_to_site: dict[int, int] = {}
    for site_index, cls in enumerate(equiv_sites):
        class_to_site.setdefault(int(cls), site_index)
    return [class_to_site[cls] for cls in _ordered_unique(equiv_sites)]


def prepare_label_data_for_reordered_matrix(
    *,
    equiv_sites: Sequence[int],
    node_labels: Sequence[str],
    frac_coords: Sequence[Sequence[float]],
    site_projections: Sequence[float] | None = None,
) -> tuple[list[str], dict[str, str]]:
    """Build labels and pretty coordinate strings from explicit inputs.

    ``site_projections``, when given, is the per-site projection onto the
    transport axis (index-aligned with ``frac_coords``/``node_labels``, as
    produced by ``JumpAxisEquivalenceAnalyzer.find_equiv_classes_for_axis``
    via ``last_site_projections``). It is used to number same-label classes
    by their plane position along the axis rather than an arbitrary order.
    A class that spans several planes (screw axis or glide plane along the
    axis) is placed at the smallest projection of its sites, not at the
    projection of its first site.
    """

    representative_indices = representative_site_indices_from_equiv(equiv_sites)
    representative_coords = [frac_coords[idx] for idx in representative_indices]
    representative_projections = None
    if site_projections is not None:
        # Position of a class: smallest projection over all of its sites.
        class_projection: dict[int, float] = {}
        for site_index, cls in enumerate(equiv_sites):
            proj = float(site_projections[site_index])
            class_projection[int(cls)] = min(proj, class_projection.get(int(cls), proj))
        representative_projections = [
            class_projection[int(equiv_sites[idx])] for idx in representative_indices
        ]
    label_data = build_non_equivalent_site_label_data(
        node_labels=node_labels,
        representative_site_indices=representative_indices,
        representative_site_coords=representative_coords,
        representative_projections=representative_projections,
    )
    return label_data.labels, label_data.pretty_site_coords


def make_real(expr: sp.Expr) -> sp.Expr:
    if expr.is_real:
        return expr
    if (expr / sp.I).is_real:
        return expr / sp.I
    return expr


def find_min_arg(
    my_transition_matrix: sp.Matrix,
    t_: sp.Symbol,
    axis_label: str | None = None,
) -> tuple[dict[Any, Any], Any]:
    """Find substitution to normalize trigonometric arguments to a minimal scale."""

    func_atoms = list(my_transition_matrix.atoms(sp.Function))
    if not func_atoms:
        return {}, 1

    my_arg_list = []
    for elem in func_atoms:
        try:
            inner = elem.args[0]
            if getattr(inner, "args", None) and len(inner.args) > 0:
                if inner.is_Mul:
                    non_t_factors = [a for a in inner.args if a != t_]
                    factor = sp.Abs(sp.Mul(*non_t_factors)) if non_t_factors else sp.Integer(1)
                else:
                    factor = sp.Abs(inner.args[0])
            else:
                factor = 1
            my_arg_list.append((factor, inner))
        except Exception:
            my_arg_list.append((sp.Abs(elem), elem))

    try:
        min_factor, expr = min(my_arg_list, key=lambda x: x[0])
    except Exception:
        expr = func_atoms[0]
        min_factor = abs(expr)

    expr_real = make_real(expr)
    prefix = f"** Axis {axis_label}; " if axis_label else "** "
    logging.debug(f"{prefix}Min_factor = {sp.nsimplify(min_factor)}; Expr = {sp.nsimplify(abs(expr_real))}")

    subst = {e: sp.nsimplify(e / expr_real, tolerance=0.02) * t_ for _, e in my_arg_list}
    return subst, sp.nsimplify(abs(expr_real / t_), tolerance=0.02)


def print_transition_matrix(
    nrows_: int,
    ncols_: int,
    my_transition_matrix: sp.Matrix,
    table_lbls_: Sequence[str],
    title: str = "Transition matrix",
) -> None:
    """Pretty-print a transition matrix with aligned columns."""

    field_length = [
        max(len(str(sp.nsimplify(my_transition_matrix[r_, c_]))) for r_ in range(nrows_))
        for c_ in range(ncols_)
    ]

    print(f"\n{title}\n")
    print("  " + ", ".join(f"{table_lbls_[l]: >{field_length[l]}}" for l in range(len(table_lbls_))))
    print("--" + "".join("-" * field_length[l] for l in range(len(table_lbls_))))
    for row in range(nrows_):
        print(
            table_lbls_[row],
            ", ".join(
                f"{str(sp.simplify(my_transition_matrix[row, col])): >{field_length[col]}}"
                for col in range(ncols_)
            ),
        )
    print("\n")


def reorder_matrix(
    nrows_: int,
    ncols_: int,
    my_transition_matrix: sp.Matrix,
    table_lbls_: Sequence[str],
) -> tuple[sp.Matrix, list[str]]:
    """Reorder matrix rows and columns according to lexicographically sorted labels."""

    labels = [str(label) for label in table_lbls_]
    label_indices = {label: index for index, label in enumerate(labels)}
    sorted_labels = sorted(labels)
    reordered_matrix = sp.Matrix(
        nrows_,
        ncols_,
        lambda i, j: my_transition_matrix[label_indices[sorted_labels[i]], label_indices[sorted_labels[j]]],
    )
    return reordered_matrix, sorted_labels


def prepare_reordered_tp_matrix(
    final_tp_matrix: sp.Matrix,
    labels: Sequence[str],
    *,
    nrows: int | None = None,
    ncols: int | None = None,
    verbose: bool = False,
    pretty_site_coords: dict[str, str] | None = None,
) -> tuple[sp.Matrix, list[str]]:
    """Reorder a reduced transition matrix using precomputed explicit labels."""

    if nrows is None or ncols is None:
        nrows, ncols = final_tp_matrix.shape

    reordered_tp_matrix, sorted_labels = reorder_matrix(nrows, ncols, final_tp_matrix, labels)

    if verbose and pretty_site_coords is not None:
        print_non_equivalent_site_coordinates(sorted_labels, pretty_site_coords)
    if verbose:
        print_transition_matrix(nrows, ncols, reordered_tp_matrix, sorted_labels)

    return reordered_tp_matrix, sorted_labels


def combine_equiv_sites(equiv: Sequence[int], full_matrix: sp.Matrix) -> sp.Matrix:
    """Reduce a matrix by summing rows/cols belonging to equivalent classes."""

    classes = _ordered_unique(equiv)
    reducer = sp.Matrix(
        [sp.Matrix([1 if eqv == cls else 0 for eqv in equiv]).T for cls in classes]
    )
    return reducer * full_matrix * reducer.pinv()


def find_equiv_and_reduce_matrix(
    equiv_sites: Sequence[int],
    full_tp_matrix: sp.Matrix,
    hkl_structure: Sequence[int],
    k_symbols: Sequence[sp.Symbol],
    t_symbol: sp.Symbol,
) -> tuple[int, int, sp.Matrix]:
    """Reduce by equivalence classes and substitute ``k -> hkl * t``."""

    reduced_tp_matrix = combine_equiv_sites(equiv_sites, full_tp_matrix)
    nrows, ncols = reduced_tp_matrix.shape

    subst_dict = {kk: hkl_structure[i] * t_symbol for i, kk in enumerate(k_symbols)}
    reduced_tp_matrix_cos = sp.Matrix(
        [
            [sp.simplify(reduced_tp_matrix[row, col].subs(subst_dict).rewrite(sp.cos)) for col in range(ncols)]
            for row in range(nrows)
        ]
    )
    return ncols, nrows, reduced_tp_matrix_cos


def find_min_jump_length(
    ncols: int,
    nrows: int,
    reduced_tp_matrix_cos: sp.Matrix,
    t: sp.Symbol,
    axis_label: str | None = None,
) -> tuple[MutableDenseMatrix, Any]:
    """Apply minimal-argument normalization and return jump-length scale."""

    subst_dict2, min_jump_length = find_min_arg(reduced_tp_matrix_cos, t, axis_label=axis_label)
    if sp.im(min_jump_length) != 0 and sp.re(min_jump_length) == 0:
        min_jump_length = -sp.I * min_jump_length

    # simulaneous=True is important to avoid issues with nested substitutions
    # otherwise substitutions can cascade and produce incorrect results
    final_tp_matrix = sp.Matrix(
        [
            [sp.simplify(reduced_tp_matrix_cos[row, col].subs(subst_dict2, simultaneous=True))
             for col in range(ncols)]
            for row in range(nrows)
        ]
    )
    return final_tp_matrix, min_jump_length
