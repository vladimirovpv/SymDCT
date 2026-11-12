"""Tests for the plane-orbit site equivalence in ``axis_equivalence``."""

import pathlib as pl

import pytest
import pymatgen as pmg
from pymatgen.core import Lattice, Structure

from src.symmetry_context.axis_equivalence import JumpAxisEquivalenceAnalyzer
from src.symmetry_context.jump_symmetry_context import JumpNetworkSymmetryContext


def _context(tmp_path: pl.Path, lattice: Lattice, species, coords) -> JumpNetworkSymmetryContext:
    struct = Structure(lattice, species, coords)
    path = tmp_path / "POSCAR"
    pmg.io.vasp.Poscar(struct).write_file(path)
    return JumpNetworkSymmetryContext(str(path), [])


def _classes(ctx, label_source, axis=(0, 0, 1)):
    analyzer = JumpAxisEquivalenceAnalyzer(ctx, label_source=label_source)
    hkl = ctx.convert_hkl(axis)
    eq = analyzer.find_equiv_classes_for_axis(hkl)
    return eq, analyzer


def _triclinic():
    # Generic triclinic cell: only inversion (P-1) is a symmetry.
    return Lattice([[3.0, 0.0, 0.0], [0.4, 3.5, 0.0], [0.3, 0.5, 4.0]])


PAIR = [(0.2, 0.1, 0.0), (0.8, 0.9, 0.0)]


def test_inversion_only_pair_is_not_merged_by_plane_orbit(tmp_path):
    """Sites related only by inversion lie in one plane but are not merged."""
    ctx = _context(tmp_path, _triclinic(), ["H", "H"], PAIR)
    eq_orbit, _ = _classes(ctx, "orbit")
    eq_plane, analyzer = _classes(ctx, "plane_orbit")
    assert len(set(eq_orbit)) == 1
    assert len(set(eq_plane)) == 2
    assert analyzer.last_plane_ops_used == 1  # identity only


def test_rotation_about_normal_merges_pair(tmp_path):
    """In an orthorhombic cell a 2-fold axis along z maps the pair onto itself."""
    ctx = _context(tmp_path, Lattice.orthorhombic(3.0, 3.5, 4.0), ["H", "H"], PAIR)
    eq_plane, analyzer = _classes(ctx, "plane_orbit")
    assert len(set(eq_plane)) == 1
    assert analyzer.last_plane_ops_used >= 2


def test_pair_along_in_plane_axis_uses_other_operations(tmp_path):
    """For the x axis the 2-fold about z reverses the normal, so no merge."""
    ctx = _context(tmp_path, Lattice.orthorhombic(3.0, 3.5, 4.0), ["H", "H"], PAIR)
    eq_plane, _ = _classes(ctx, "plane_orbit", axis=(1, 0, 0))
    assert len(set(eq_plane)) == 2


def test_different_species_are_never_merged(tmp_path):
    ctx = _context(tmp_path, Lattice.orthorhombic(3.0, 3.5, 4.0), ["H", "Li"], PAIR)
    eq_plane, _ = _classes(ctx, "plane_orbit")
    assert len(set(eq_plane)) == 2


def test_classes_are_sorted_by_plane_coordinate(tmp_path):
    coords = [(0.0, 0.0, 0.7), (0.0, 0.0, 0.2)]
    ctx = _context(tmp_path, Lattice.orthorhombic(3.0, 3.5, 4.0), ["H", "Li"], coords)
    eq_plane, analyzer = _classes(ctx, "plane_orbit")
    proj = analyzer.last_site_projections
    order = sorted(range(len(proj)), key=lambda i: eq_plane[i])
    assert [proj[i] for i in order] == sorted(proj)


def test_invalid_label_source_raises(tmp_path):
    ctx = _context(tmp_path, Lattice.cubic(3.0), ["H"], [(0.0, 0.0, 0.0)])
    analyzer = JumpAxisEquivalenceAnalyzer(ctx, label_source="nonsense")
    with pytest.raises(ValueError):
        analyzer.find_equiv_classes_for_axis(ctx.convert_hkl((0, 0, 1)))


def _hcp_pair_context(tmp_path):
    lattice = Lattice.hexagonal(2.95, 4.82)
    coords = [(1 / 3, 2 / 3, 0.25), (2 / 3, 1 / 3, 0.75)]
    return _context(tmp_path, lattice, ["H", "H"], coords)


def test_screw_axis_along_c_merges_planes_in_axis_orbit(tmp_path):
    """The 6_3 screw axis maps the two hcp sites onto each other along c."""
    ctx = _hcp_pair_context(tmp_path)
    eq_axis, _ = _classes(ctx, "axis_orbit", axis=(0, 0, 1))
    eq_plane, _ = _classes(ctx, "plane_orbit", axis=(0, 0, 1))
    assert len(set(eq_axis)) == 1
    assert len(set(eq_plane)) == 2


def test_hcp_pair_not_merged_perpendicular_to_a(tmp_path):
    """No operation that keeps this direction relates the two sites."""
    ctx = _hcp_pair_context(tmp_path)
    eq_axis, _ = _classes(ctx, "axis_orbit", axis=(1, 0, 0))
    assert len(set(eq_axis)) == 2


def test_inversion_only_pair_not_merged_by_axis_orbit(tmp_path):
    ctx = _context(tmp_path, _triclinic(), ["H", "H"], PAIR)
    eq_axis, _ = _classes(ctx, "axis_orbit")
    assert len(set(eq_axis)) == 2
