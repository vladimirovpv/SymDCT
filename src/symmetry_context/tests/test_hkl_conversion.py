"""Regression tests for the conventional -> structure (primitive) hkl conversion.

The conversion is checked against an independent geometric criterion. The plane
normal built from the conventional indices in the refined conventional cell must
be parallel to the plane normal built from the converted indices in the primitive
cell used by the pipeline.  Fixed expected indices are not used on purpose,
because spglib's choice of primitive vectors is not part of the contract.
"""

import numpy as np
import pytest
from pymatgen.core import Lattice, Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from src.symmetry_context.smatrix_handler import SMatrixHandler
from src.symmetry_context.symmetry_refiner import JumpStructureSymmetryRefiner

CASES = {
    "fcc": ("Fm-3m", Lattice.cubic(4.0), [[0, 0, 0]]),
    "bcc": ("Im-3m", Lattice.cubic(3.0), [[0, 0, 0]]),
    "hcp": ("P6_3/mmc", Lattice.hexagonal(2.9, 4.7), [[1 / 3, 2 / 3, 0.25]]),
    "C-orthorhombic": ("Cmcm", Lattice.orthorhombic(3.2, 9.9, 4.2), [[0, 0.18, 0.25]]),
    "R-hexagonal": ("R-3m", Lattice.hexagonal(4.1, 25.0), [[0, 0, 0.23]]),
    "I-tetragonal": ("I4/mmm", Lattice.tetragonal(3.1, 8.0), [[0, 0, 0]]),
}

HKLS = [(1, 0, 0), (0, 0, 1), (1, 1, 0), (1, 1, 1), (2, -1, 0), (1, 2, 3)]


def _normal(lattice: Lattice, hkl) -> np.ndarray:
    return lattice.reciprocal_lattice_crystallographic.matrix.T @ np.asarray(hkl, float)


@pytest.mark.parametrize("name", list(CASES))
@pytest.mark.parametrize("hkl", HKLS)
def test_converted_hkl_keeps_the_plane_normal(name, hkl):
    spacegroup, lattice, coords = CASES[name]
    species = ["X"] * len(coords)
    struct = Structure.from_spacegroup(spacegroup, lattice, species, coords)

    refiner = JumpStructureSymmetryRefiner(struct, struct, symprec=1e-3)
    handler = SMatrixHandler(refiner)

    conv_lattice = SpacegroupAnalyzer(struct, symprec=1e-3).get_refined_structure().lattice
    prim_lattice = refiner.symmetrized_structure.lattice

    h_struct = handler.convert_hkl_to_structure(hkl)

    n_conv = _normal(conv_lattice, hkl)
    n_prim = _normal(prim_lattice, h_struct)
    cross = np.cross(n_conv, n_prim)
    assert np.linalg.norm(cross) < 1e-6 * np.linalg.norm(n_conv) * np.linalg.norm(n_prim)
    assert np.dot(n_conv, n_prim) > 0


@pytest.mark.parametrize("name", list(CASES))
def test_s_matrix_rows_are_primitive_vectors_in_conventional_basis(name):
    spacegroup, lattice, coords = CASES[name]
    struct = Structure.from_spacegroup(spacegroup, lattice, ["X"] * len(coords), coords)
    refiner = JumpStructureSymmetryRefiner(struct, struct, symprec=1e-3)

    conv_lattice = SpacegroupAnalyzer(struct, symprec=1e-3).get_refined_structure().lattice
    prim_lattice = refiner.symmetrized_structure.lattice

    s = refiner.s_matrix
    np.testing.assert_allclose(s @ conv_lattice.matrix, prim_lattice.matrix, atol=1e-8)
    assert abs(np.linalg.det(s)) == pytest.approx(prim_lattice.volume / conv_lattice.volume)


@pytest.mark.parametrize("name", list(CASES))
def test_s_matrix_is_a_unimodular_change_of_the_tabulated_matrix(name):
    """Frame check independent of the construction of S.

    The tabulated matrix and the matrix of the pipeline describe the same primitive
    lattice, so they differ by an integer matrix of determinant +-1. A rotation of the
    Cartesian frame between the two cells would break this.
    """
    spacegroup, lattice, coords = CASES[name]
    struct = Structure.from_spacegroup(spacegroup, lattice, ["X"] * len(coords), coords)
    refiner = JumpStructureSymmetryRefiner(struct, struct, symprec=1e-3)

    spg_conv = SpacegroupAnalyzer(
        SpacegroupAnalyzer(struct, symprec=1e-3).get_refined_structure(), symprec=1e-3
    )
    s_tab = np.array(spg_conv.get_conventional_to_primitive_transformation_matrix())

    u = refiner.s_matrix @ np.linalg.inv(s_tab)
    np.testing.assert_allclose(u, np.round(u), atol=1e-6)
    assert abs(round(np.linalg.det(u))) == 1


# ---------------------------------------------------------------------------
# Setting of the user's indices.  hkl refer to the input cell if it is a
# conventional cell (also in a non-standard setting) and to the standardised
# conventional cell if the input cell is primitive.  The check is frame
# independent, because the standardisation may rotate the Cartesian frame.
# The Gram matrix of the plane normals must be the same in both descriptions.
# ---------------------------------------------------------------------------


def _cyclic(struct: Structure) -> Structure:
    """Same crystal with the axes permuted cyclically (a, b, c) -> (b, c, a)."""
    perm = [1, 2, 0]
    lattice = Lattice(struct.lattice.matrix[perm])
    return Structure(lattice, struct.species, struct.frac_coords[:, perm])


def _input_cases():
    pair = Structure(
        Lattice.orthorhombic(3.0, 3.5, 4.0), ["H", "H"], [(0.2, 0.1, 0.0), (0.8, 0.9, 0.0)]
    )
    cmcm = Structure.from_spacegroup(*CASES["C-orthorhombic"][:2], ["X"], CASES["C-orthorhombic"][2])
    fcc = Structure.from_spacegroup(*CASES["fcc"][:2], ["X"], CASES["fcc"][2])
    fcc_prim = SpacegroupAnalyzer(fcc, symprec=1e-3).get_primitive_standard_structure()
    return {
        "P-orthorhombic, non-ascending refined axes": (pair, "input"),
        "Cmcm, cyclic axis permutation": (_cyclic(cmcm), "input"),
        "fcc primitive input": (fcc_prim, "standardised"),
    }


@pytest.mark.parametrize("label", list(_input_cases()))
def test_hkl_setting_follows_input_cell_if_conventional(label):
    struct, reference = _input_cases()[label]
    refiner = JumpStructureSymmetryRefiner(struct, struct, symprec=1e-3)
    handler = SMatrixHandler(refiner)

    if reference == "input":
        ref_lattice = struct.lattice
    else:
        ref_lattice = SpacegroupAnalyzer(struct, symprec=1e-3).get_refined_structure().lattice
    prim_lattice = refiner.symmetrized_structure.lattice

    # Scaling by the reduction to coprime integers does not change the plane,
    # so compare normalised normals through the Gram matrix of unit vectors.
    converted = [handler.convert_hkl_to_structure(h) for h in HKLS]
    n_ref = np.array([_normal(ref_lattice, h) for h in HKLS])
    n_prim = np.array([_normal(prim_lattice, h) for h in converted])
    n_ref /= np.linalg.norm(n_ref, axis=1)[:, None]
    n_prim /= np.linalg.norm(n_prim, axis=1)[:, None]
    np.testing.assert_allclose(n_ref @ n_ref.T, n_prim @ n_prim.T, atol=1e-8)
