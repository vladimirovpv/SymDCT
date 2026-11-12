"""Tests for the lumpability check of the equivalent-site reduction."""

from __future__ import annotations

import time

import pytest
import sympy as sp

from src.matrix_builder.lumpability import check_lumpability

k1, k2, k3, t = sp.symbols("k1 k2 k3 t", real=True)
W0, W1 = sp.symbols("W0 W1", positive=True)
KS = [k1, k2, k3]
METHODS = ["symbolic", "numeric"]


@pytest.mark.parametrize("method", METHODS)
def test_symmetric_pair_is_lumpable(method):
    m = -2 * W1
    a12 = W1 * (sp.exp(sp.I * k1 / 2) + sp.exp(-sp.I * k1 / 2))
    full = sp.Matrix([[m, a12], [a12, m]])
    assert check_lumpability([0, 0], full, (1, 0, 0), KS, t, method=method).ok


@pytest.mark.parametrize("method", METHODS)
def test_one_sided_pair_is_not_lumpable(method):
    # Chain 1(l) -W0- 2(l) -W1- 1(l+1): A12 = W0 + W1 exp(-i k1), A21 = conj(A12).
    m = -(W0 + W1)
    a12 = W0 + W1 * sp.exp(-sp.I * k1)
    a21 = W0 + W1 * sp.exp(sp.I * k1)
    full = sp.Matrix([[m, a12], [a21, m]])
    rep = check_lumpability([0, 0], full, (1, 0, 0), KS, t, method=method)
    assert not rep.ok


@pytest.mark.parametrize("method", METHODS)
def test_one_sided_pair_is_lumpable_perpendicular_to_chain(method):
    m = -(W0 + W1)
    a12 = W0 + W1 * sp.exp(-sp.I * k1)
    a21 = W0 + W1 * sp.exp(sp.I * k1)
    full = sp.Matrix([[m, a12], [a21, m]])
    assert check_lumpability([0, 0], full, (0, 1, 0), KS, t, method=method).ok


def test_different_classes_are_not_compared():
    full = sp.Matrix([[-1, W0], [W0, -1]])
    assert check_lumpability([0, 1], full, (1, 0, 0), KS, t).ok


def _hcp_vacancy_matrix():
    a, c, Wa, Wb = sp.symbols("a c W_a W_b", positive=True)
    kx, ky, kz = KS
    a11 = -6 * Wa + 2 * Wb * (
        sp.cos(kx * a) + 2 * sp.cos(kx * a / 2) * sp.cos(ky * a * sp.sqrt(3) / 2) - 3
    )
    a12 = 2 * Wa * sp.cos(kz * c / 2) * (
        2 * sp.exp(-sp.I * ky * a * sp.sqrt(3) / 6) * sp.cos(kx * a / 2)
        + sp.exp(sp.I * ky * a * sp.sqrt(3) / 3)
    )
    a21 = sp.conjugate(a12).subs({sp.conjugate(s): s for s in KS})
    a21 = 2 * Wa * sp.cos(kz * c / 2) * (
        2 * sp.exp(sp.I * ky * a * sp.sqrt(3) / 6) * sp.cos(kx * a / 2)
        + sp.exp(-sp.I * ky * a * sp.sqrt(3) / 3)
    )
    return sp.Matrix([[a11, a12], [a21, a11]])


@pytest.mark.parametrize("method", METHODS)
def test_hcp_vacancy_ab_merge_along_a_and_c(method):
    # Paper example: merging the two basis sites A and B of the hcp vacancy is
    # exact along the transport directions used (A12 is real there).
    full = _hcp_vacancy_matrix()
    assert check_lumpability([0, 0], full, (1, 0, 0), KS, t, method=method).ok
    assert check_lumpability([0, 0], full, (0, 0, 1), KS, t, method=method).ok


@pytest.mark.parametrize("method", METHODS)
def test_hcp_vacancy_ab_merge_fails_in_general_direction(method):
    # In a general direction (ky != 0) Im A12 != 0 and the merge is only approximate.
    full = _hcp_vacancy_matrix()
    assert not check_lumpability([0, 0], full, (1, 1, 0), KS, t, method=method).ok


def test_large_lumpable_matrix_is_fast():
    # 12 sites in 4 classes of 3, lumpable by construction:
    # A[(c,a),(d,b)] = q_cd(t) * w_ab with sum_a w_ab = 1 for every b.
    rates = sp.symbols("r0:16", positive=True)
    q = sp.Matrix(
        4, 4,
        lambda c, d: rates[4 * c + d] * (sp.exp(sp.I * k1 * (c - d) / 2) + sp.exp(-sp.I * k1 * (c - d) / 2)),
    )
    u, v = sp.symbols("u v", positive=True)
    w = sp.Matrix([[u, v, 1 / sp.Integer(3)], [1 - u - sp.Rational(1, 3), v, 1 / sp.Integer(3)],
                   [sp.Rational(1, 3), 1 - 2 * v, sp.Rational(1, 3)]])
    # make column sums of w equal to one
    w = sp.Matrix([[u, v, sp.Rational(1, 3)], [sp.Rational(2, 3) - u, sp.Rational(1, 3), sp.Rational(1, 3)],
                   [sp.Rational(1, 3), sp.Rational(2, 3) - v, sp.Rational(1, 3)]])
    big = sp.Matrix(12, 12, lambda i, j: q[i // 3, j // 3] * w[i % 3, j % 3])
    equiv = [i // 3 for i in range(12)]
    start = time.perf_counter()
    rep = check_lumpability(equiv, big, (1, 0, 0), KS, t)
    elapsed = time.perf_counter() - start
    assert rep.ok, rep.summary()
    assert elapsed < 10


# ---------------------------------------------------------------------------
# refinement and rate relations
# ---------------------------------------------------------------------------
from src.matrix_builder.lumpability import refine_to_lumpable, suggest_rate_relations


def _diplon_like(independent: bool):
    w = sp.symbols("W0:9", positive=True)

    def diag(wd, r1, r2):
        return 3 * wd * sp.exp(sp.I * k1) + 3 * wd * sp.exp(-sp.I * k1) - 6 * wd - r1 - r2

    a = sp.Matrix([
        [diag(w[0], w[1], w[2]), w[5], w[8]],
        [w[2], diag(w[3], w[4], w[5]), w[7]],
        [w[1], w[4], diag(w[6], w[7], w[8])],
    ])
    if independent:
        return a, w
    T, R = sp.symbols("wT wR", positive=True)
    sub = {w[0]: T, w[3]: T, w[6]: T}
    sub.update({w[i]: R for i in (1, 2, 4, 5, 7, 8)})
    return a.subs(sub), w


def test_refine_splits_nonlumpable_class():
    full, _ = _diplon_like(independent=True)
    labels = refine_to_lumpable([0, 0, 0], full, (1, 0, 0), KS, t)
    assert sorted(set(labels)) == [0, 1, 2]


def test_refine_keeps_lumpable_class():
    full, _ = _diplon_like(independent=False)
    assert refine_to_lumpable([0, 0, 0], full, (1, 0, 0), KS, t) == [0, 0, 0]


def test_refine_keeps_other_classes_and_order():
    full = sp.Matrix([[-1, 0, 0], [0, -1, 0], [0, 0, -1]])
    assert refine_to_lumpable([0, 1, 1], full, (1, 0, 0), KS, t) == [0, 1, 1]


def test_suggested_relations_for_diplon():
    full, w = _diplon_like(independent=True)
    rep = check_lumpability([0, 0, 0], full, (1, 0, 0), KS, t)
    rels = suggest_rate_relations(rep)
    assert not rep.ok and rels
    solved = {r.lhs: r.rhs for r in rels}
    # translation rates of the three orientations must be identified
    assert solved[w[3]] == w[0] and solved[w[6]] == w[0]
    # applying the relations makes the partition lumpable
    fixed = full.subs(solved)
    assert check_lumpability([0, 0, 0], fixed, (1, 0, 0), KS, t).ok


def test_suggested_relations_empty_when_lumpable():
    full, _ = _diplon_like(independent=False)
    rep = check_lumpability([0, 0, 0], full, (1, 0, 0), KS, t)
    assert suggest_rate_relations(rep) == []


@pytest.mark.parametrize("direction", [(1, 0, 0), (0, 1, 0), (1, 1, 0)])
def test_hcp_vacancy_merge_keeps_D_and_isotropy(direction):
    # Even where the merge is not lumpable (Im A12 != 0), the quadratic part of
    # the merged matrix gives the same D as the full 2x2 matrix, because Im A12
    # starts at order theta^3 (no dipole moment of the sites).
    full = _hcp_vacancy_matrix()
    sub = {KS[i]: direction[i] * t for i in range(3)}
    A = full.subs(sub)
    B = A.subs(t, 0)
    B[1, :] = sp.Matrix([[1, 1]])
    d_full = -sp.diff(A.det(), t, 2).subs(t, 0) / (2 * B.det())
    merged = A[0, 0] + (A[0, 1] + A[1, 0]) / 2
    d_merged = -sp.diff(merged, t, 2).subs(t, 0) / 2
    assert sp.simplify(sp.expand_complex(sp.simplify(d_full - d_merged))) == 0


# ---------------------------------------------------------------------------
# message printed by the merge policy
# ---------------------------------------------------------------------------
def test_merge_policy_reports_lumpable_merge(capsys):
    from src.matrix_builder.matrix_builder import _apply_merge_policy

    full, _ = _diplon_like(independent=False)
    labels = _apply_merge_policy([0, 0, 0], full, (1, 0, 0), KS, t, "check", "(1, 0, 0)")
    out = capsys.readouterr().out
    assert labels == [0, 0, 0]
    assert "lumpable" in out and "WARNING" not in out


def test_merge_policy_is_silent_without_merge(capsys):
    from src.matrix_builder.matrix_builder import _apply_merge_policy

    full, _ = _diplon_like(independent=False)
    _apply_merge_policy([0, 1, 2], full, (1, 0, 0), KS, t, "check", "(1, 0, 0)")
    assert capsys.readouterr().out == ""
