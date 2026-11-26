"""Consistency of the symbolic lattice vectors with the primitive cell used by the pipeline.

The symbols a, b and c are replaced by the lengths of the standardised conventional
vectors. The Gram matrix of the resulting primitive vectors must equal the Gram matrix
of the primitive cell of the pipeline. The Gram matrix is independent of the Cartesian
frame, which the standardisation may rotate.
"""

import numpy as np
import pytest
import sympy as sp
from pymatgen.core import Lattice, Structure

from src.symmetry_context.jump_symmetry_context import JumpNetworkSymmetryContext
from src.symmetry_context.tests.test_hkl_conversion import CASES, _cyclic


def _context(struct: Structure, tmp_path) -> JumpNetworkSymmetryContext:
    path = tmp_path / "POSCAR"
    struct.to(filename=str(path), fmt="poscar")
    return JumpNetworkSymmetryContext(poscar_file=str(path), crystal_matrix=None)


def _symbolic_gram(ctx: JumpNetworkSymmetryContext) -> np.ndarray:
    a, b, c = sp.symbols("a b c", real=True, positive=True)
    conv = ctx._refiner.conventional_lattice
    lat = np.array(ctx.get_symbolic_lattice_vectors().subs({a: conv.a, b: conv.b, c: conv.c}).evalf(), float)
    return lat @ lat.T


def _structures():
    out = {name: Structure.from_spacegroup(sg, lat, ["X"] * len(co), co) for name, (sg, lat, co) in CASES.items()}
    cmcm = out["C-orthorhombic"]
    out["Cmcm, cyclic axis permutation"] = _cyclic(cmcm)
    out["P-orthorhombic pair"] = Structure(
        Lattice.orthorhombic(3.0, 3.5, 4.0), ["H", "H"], [(0.2, 0.1, 0.0), (0.8, 0.9, 0.0)]
    )
    return out


@pytest.mark.parametrize("label", list(_structures()))
def test_symbolic_lattice_vectors_reproduce_the_primitive_cell(label, tmp_path):
    ctx = _context(_structures()[label], tmp_path)
    prim = ctx.get_refined_jump_structure().lattice.matrix
    np.testing.assert_allclose(_symbolic_gram(ctx), prim @ prim.T, atol=1e-6)


def test_monoclinic_cell_with_oblique_angle_is_rejected(tmp_path):
    struct = Structure(Lattice.monoclinic(3.0, 4.0, 5.0, 100.0), ["H", "H"], [(0.1, 0.2, 0.3), (0.7, 0.1, 0.6)])
    ctx = _context(struct, tmp_path)
    with pytest.raises(NotImplementedError):
        ctx.get_symbolic_lattice_vectors()
