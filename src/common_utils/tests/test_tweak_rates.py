"""Tests for :mod:`src.common_utils.tweak_rates`.

Covers ``_parse_label_spec``'s qualifier grammar -- in particular the bare
number shorthand for ``length=<value>``, added after ``hcp_vac.yaml``'s
``"A->B: 2.213"`` key was found to be silently ignored (the bare number was
not recognised, so the tweak fell back to matching every same-label jump
instead of just the one at that length) -- and ``normalize_tweak_rates``'s
exact-match-within-tolerance semantics for length-qualified specs.

Length values in params YAML files are written to 3 decimal places by
convention (see ``hcp_vac.yaml``), so the tolerance actually used at the
call site, ``src.matrix_builder.types.TWEAK_LENGTH_TOLERANCE`` (0.005 --
10x the max rounding error at 3 decimals), is deliberately looser than
the generic 1e-6 default exercised directly below.
"""

from __future__ import annotations

from src.common_utils.tweak_rates import _parse_label_spec, normalize_tweak_rates
from src.matrix_builder.types import TWEAK_LENGTH_TOLERANCE


def test_parse_label_spec_bare_number_is_length_shorthand():
    assert _parse_label_spec("A->B: 2.213") == ("A->B", 2.213, False)


def test_parse_label_spec_explicit_length_prefixes():
    assert _parse_label_spec("A->B: length=2.213") == ("A->B", 2.213, False)
    assert _parse_label_spec("A->B: len=2.213") == ("A->B", 2.213, False)
    assert _parse_label_spec("A->B: l=2.213") == ("A->B", 2.213, False)


def test_parse_label_spec_all_and_unqualified():
    assert _parse_label_spec("A->A: all") == ("A->A", None, True)
    assert _parse_label_spec("A->A") == ("A->A", None, False)


def test_parse_label_spec_unparseable_qualifier_falls_back_to_no_filter():
    assert _parse_label_spec("A->B: not_a_number") == ("A->B", None, False)


def test_normalize_tweak_rates_length_match_is_exact_within_tolerance():
    records = {
        "A->B": [{"glob": 1, "length": 2.211735}],
        "B->A": [{"glob": 3, "length": 2.211735}],
    }
    # Exact match (within tolerance) resolves and applies the tweak.
    result = normalize_tweak_rates({"A->B: 2.211735": "B->A"}, records, tolerance=1e-6)
    assert result == {1: 3}

    # An approximate length that misses a tight tolerance window matches
    # nothing -- this is what broke when hcp_vac.yaml's key said "2.213"
    # (off by ~1.3e-3 from the true length) while the call site still used
    # the tight 1e-6 tolerance meant for exact symbolic work.
    result = normalize_tweak_rates({"A->B: 2.213": "B->A"}, records, tolerance=1e-6)
    assert result == {}


def test_normalize_tweak_rates_project_tolerance_covers_3_decimal_rounding():
    # At the project's actual TWEAK_LENGTH_TOLERANCE (0.005), a length
    # written to 3 decimal places -- the params YAML convention -- matches
    # the true computed length even though they differ past the 3rd decimal.
    records = {
        "A->B": [{"glob": 1, "length": 2.211735}],
        "B->A": [{"glob": 3, "length": 2.211735}],
    }
    result = normalize_tweak_rates({"A->B: 2.212": "B->A"}, records, tolerance=TWEAK_LENGTH_TOLERANCE)
    assert result == {1: 3}
