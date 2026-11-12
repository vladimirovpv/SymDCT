import pathlib as pl

import pytest
import pymatgen as pmg
from pymatgen.core import Lattice, Structure

from ..loader import PoscarJumpStructureLoader
from src.symmetry_context.symmetry_refiner import JumpStructureSymmetryRefiner
from src.symmetry_context.node_labeler import SymmetryNodeLabeler
from src.symmetry_context.jump_symmetry_context import JumpNetworkSymmetryContext


def _make_simple_bcc_structure(tmp_path: pl.Path) -> str:
    """Create a minimal bcc-like POSCAR on disk and return its path.

    The structure is deliberately simple: a cubic cell with two atoms of
    the same species at (0, 0, 0) and (1/2, 1/2, 1/2).
    """

    lattice = Lattice.cubic(3.0)
    species = ["X", "X"]
    coords = [
        (0.0, 0.0, 0.0),
        (0.5, 0.5, 0.5),
    ]
    struct = Structure(lattice, species, coords)
    poscar = pmg.io.vasp.Poscar(struct)
    poscar_path = tmp_path / "POSCAR"
    poscar.write_file(poscar_path)
    return str(poscar_path)


def _make_ot_poscar(tmp_path: pl.Path) -> str:
    poscar_path = tmp_path / "POSCAR_OT"
    poscar_path.write_text(
        "\n".join(
            [
                "OT test",
                "1.0",
                "3 0 0",
                "0 3 0",
                "0 0 3",
                "O T",
                "1 1",
                "direct",
                "0 0 0",
                "0.5 0.5 0.5",
                "",
            ]
        ),
        encoding="utf-8",
    )
    return str(poscar_path)


def test_poscar_loader_matrix_and_jump_structures(tmp_path):
    poscar_path = _make_simple_bcc_structure(tmp_path)

    # With no crystal_matrix provided, matrix and jump structures should
    # both be identical copies of the input structure.
    loader = PoscarJumpStructureLoader(poscar_path)
    assert loader.matrix_structure.num_sites == 2
    assert loader.jump_structure.num_sites == 2

    # With a matrix species specified, the jump structure should lose
    # those sites. Since no interstitial sites remain, the jump_structure
    # property asserts rather than returning an empty Structure.
    loader2 = PoscarJumpStructureLoader(poscar_path, crystal_matrix=["X"])
    assert loader2.matrix_structure.num_sites == 2
    with pytest.raises(AssertionError):
        loader2.jump_structure


def test_poscar_loader_preserves_raw_t_label(tmp_path):
    poscar_path = _make_ot_poscar(tmp_path)
    loader = PoscarJumpStructureLoader(poscar_path)

    assert loader.raw_symbol_alias_map == {"O": "O", "H": "T"}


def test_symmetry_refiner_builds_symmetrized_structure(tmp_path):
    poscar_path = _make_simple_bcc_structure(tmp_path)
    loader = PoscarJumpStructureLoader(poscar_path)

    refiner = JumpStructureSymmetryRefiner(
        loader.matrix_structure,
        loader.jump_structure,
        symprec=1e-3,
        check_symmetry_consistency=True,
    )

    sym_struct = refiner.symmetrized_structure
    assert sym_struct.num_sites > 0

    # Each site in the symmetrized structure should carry a "site_idx"
    # property (used by the node labeler and graph builder).
    assert "site_idx" in sym_struct.site_properties
    assert len(sym_struct.site_properties["site_idx"]) == sym_struct.num_sites

    # A spacegroup analyzer must be available and expose symmetry ops.
    spg = refiner.spacegroup_analyzer
    ops = spg.get_symmetry_operations(cartesian=False)
    assert len(ops) > 0


def test_symmetry_node_labeler_produces_labels(tmp_path):
    poscar_path = _make_simple_bcc_structure(tmp_path)
    loader = PoscarJumpStructureLoader(poscar_path)
    refiner = JumpStructureSymmetryRefiner(
        loader.matrix_structure,
        loader.jump_structure,
        symprec=1e-3,
        check_symmetry_consistency=True,
    )

    sym_struct = refiner.symmetrized_structure
    labeler = SymmetryNodeLabeler(sym_struct)
    labels = labeler.compute_labels()

    assert isinstance(labels, list)
    assert len(labels) == sym_struct.num_sites
    # For a simple two-site bcc-like structure of the same species, we
    # expect at least one label to be present.
    assert all(isinstance(lbl, str) for lbl in labels)

def test_jump_network_symmetry_context_symmetry_operations(tmp_path):
    poscar_path = _make_simple_bcc_structure(tmp_path)

    ctx = JumpNetworkSymmetryContext(poscar_path, crystal_matrix=None, symprec=1e-3)

    ops_frac = ctx.get_symmetry_operations(cartesian=False)

    assert len(ops_frac) > 0


def test_jump_network_symmetry_context_preserves_t_label(tmp_path):
    poscar_path = _make_ot_poscar(tmp_path)

    ctx = JumpNetworkSymmetryContext(poscar_path, crystal_matrix=None, symprec=1e-3)
    labels = ctx.get_node_lbls()
    alias_map = ctx.get_symbol_alias_map()

    assert "H" in labels
    assert alias_map.get("H") == "T"
