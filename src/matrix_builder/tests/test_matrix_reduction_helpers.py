from __future__ import annotations

import sympy as sp

from src.matrix_builder.matrix_reduction_helpers import (
    find_equiv_and_reduce_matrix,
    find_min_jump_length,
    prepare_label_data_for_reordered_matrix,
    prepare_reordered_tp_matrix,
)


def test_prepare_reordered_tp_matrix_reorders_by_label():
    matrix = sp.Matrix([[1, 2], [3, 4]])

    reordered, sorted_labels = prepare_reordered_tp_matrix(
        matrix,
        labels=["B0", "A0"],
    )

    assert sorted_labels == ["A0", "B0"]
    assert reordered == sp.Matrix([[4, 3], [2, 1]])


def test_prepare_label_data_from_explicit_inputs():
    labels, pretty = prepare_label_data_for_reordered_matrix(
        equiv_sites=[1, 1, 0, 0],
        node_labels=["A", "A", "B", "B"],
        frac_coords=[
            [0.0, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            [0.0, 0.5, 0.0],
            [0.5, 0.5, 0.0],
        ],
    )

    assert labels == ["A0", "B1"]
    assert set(pretty.keys()) == {"A0", "B1"}


def test_prepare_label_data_numbers_by_axis_projection():
    # Structure order (0,1,2,3) intentionally does not match projection
    # order, to verify numbering follows plane position along the axis
    # rather than order of appearance.
    labels, pretty = prepare_label_data_for_reordered_matrix(
        equiv_sites=[0, 1, 2, 3],
        node_labels=["A", "A", "B", "B"],
        frac_coords=[
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.25],
            [0.0, 0.0, 0.5],
            [0.0, 0.0, 0.75],
        ],
        site_projections=[5.0, 2.0, 7.0, 1.0],
    )

    # Site 1 (proj=2) is the lower plane among the "A" sites -> A0;
    # site 0 (proj=5) is the higher plane -> A1. Likewise for "B".
    assert labels == ["A1", "A0", "B1", "B0"]
    assert set(pretty.keys()) == {"A0", "A1", "B0", "B1"}


def test_prepare_label_data_places_multi_plane_class_at_smallest_projection():
    # Class 0 spans two planes (screw axis): projections 0.6 (first site) and
    # 0.1. Class 1 has projections 0.3 (first site) and 0.8. By the smallest
    # projection class 0 (0.1) lies below class 1 (0.3). Using the first site
    # of each class (0.6 and 0.3) would reverse the numbering.
    labels, _ = prepare_label_data_for_reordered_matrix(
        equiv_sites=[0, 1, 0, 1],
        node_labels=["O", "O", "O", "O"],
        frac_coords=[[0.0, 0.0, 0.0]] * 4,
        site_projections=[0.6, 0.3, 0.1, 0.8],
    )

    assert labels == ["O0", "O1"]


def test_prepare_label_data_without_projections_falls_back_to_order():
    labels, _ = prepare_label_data_for_reordered_matrix(
        equiv_sites=[1, 1, 0, 0],
        node_labels=["A", "A", "B", "B"],
        frac_coords=[
            [0.0, 0.0, 0.0],
            [0.5, 0.0, 0.0],
            [0.0, 0.5, 0.0],
            [0.5, 0.5, 0.0],
        ],
    )

    assert labels == ["A0", "B1"]


def test_find_equiv_and_reduce_matrix_shape_and_substitution():
    k1, k2, k3, t = sp.symbols("k1 k2 k3 t", real=True)
    full = sp.Matrix(
        [
            [sp.cos(k1), 0, 1, 0],
            [0, sp.cos(k1), 0, 1],
            [1, 0, sp.cos(k1), 0],
            [0, 1, 0, sp.cos(k1)],
        ]
    )

    ncols, nrows, reduced = find_equiv_and_reduce_matrix(
        equiv_sites=[0, 0, 1, 1],
        full_tp_matrix=full,
        hkl_structure=(1, 0, 0),
        k_symbols=[k1, k2, k3],
        t_symbol=t,
    )

    assert (nrows, ncols) == (2, 2)
    assert reduced.shape == (2, 2)
    assert reduced[0, 0].has(sp.cos(t))


def test_find_min_jump_length_normalizes_arguments():
    t = sp.Symbol("t", real=True)
    reduced = sp.Matrix([[sp.cos(2 * t), 0], [0, sp.cos(4 * t)]])

    final, min_jump_length = find_min_jump_length(2, 2, reduced, t)

    assert sp.simplify(min_jump_length - 2) == 0
    assert sp.simplify(final[0, 0] - sp.cos(t)) == 0
    assert sp.simplify(final[1, 1] - sp.cos(2 * t)) == 0

