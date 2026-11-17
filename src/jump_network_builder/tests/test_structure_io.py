"""Tests for reading crystal structures from the formats that pymatgen supports."""

import json
import pathlib as pl

import pytest
from pymatgen.core import Lattice, Structure
from pymatgen.io.vasp import Poscar

from ..geometry import build_jump_network_from_poscar
from ..loader import PoscarJumpStructureLoader
from ..structure_io import read_structure

EXAMPLES = pl.Path(__file__).resolve().parents[3] / "examples"
HCP = EXAMPLES / "hcp_vac" / "input" / "hcp_Be.poscar"


def _reference() -> Structure:
    """hcp host with one real-element jump site, readable from every format."""
    dummy = read_structure(str(HCP))[0]
    return Structure(dummy.lattice, ["Be", "He"], dummy.frac_coords)


@pytest.fixture()
def hcp_real(tmp_path):
    path = tmp_path / "hcp_real.poscar"
    Poscar(_reference()).write_file(path)
    return str(path)


@pytest.mark.parametrize("name,fmt", [
    ("s.cif", None), ("s.json", None), ("s.xsf", None), ("s.cssr", None),
    ("s.yaml", None), ("s.res", None), ("s_cif.txt", "cif"),
])
def test_formats_give_same_structure(tmp_path, name, fmt):
    ref = _reference()
    path = tmp_path / name
    if name.endswith(".yaml"):
        import io
        from ruamel.yaml import YAML
        buf = io.StringIO()
        YAML().dump(ref.as_dict(), buf)
        path.write_text(buf.getvalue())
    elif name.endswith(".json"):
        path.write_text(json.dumps(ref.as_dict()))
    else:
        ref.to(filename=str(path), **({"fmt": fmt} if fmt else {}))
    struct, aliases = read_structure(str(path), fmt)
    assert aliases == {}
    assert len(struct) == len(ref)
    assert struct.lattice.abc == pytest.approx(ref.lattice.abc, rel=1e-4)
    assert struct.matches(ref)


def test_cif_gives_same_jump_network(tmp_path, hcp_real):
    ref = _reference()
    cif = tmp_path / "hcp.cif"
    ref.to(filename=str(cif))
    out_a, out_b = tmp_path / "a", tmp_path / "b"
    Ga, jra, _ = build_jump_network_from_poscar(hcp_real, {}, str(out_a), None, 3.0)
    Gb, jrb, _ = build_jump_network_from_poscar(str(cif), {}, str(out_b), None, 3.0)
    assert jra.get_number_of_jumps() == jrb.get_number_of_jumps()
    assert Ga.number_of_nodes() == Gb.number_of_nodes()
    assert Ga.number_of_edges() == Gb.number_of_edges()


def test_poscar_with_lowercase_extension_keeps_aliases(tmp_path):
    s = Structure(Lattice.cubic(3.0), ["He", "He"], [(0, 0, 0), (0.5, 0.5, 0.5)])
    path = tmp_path / "x.poscar"
    Poscar(s).write_file(path)
    lines = path.read_text().splitlines()
    lines[5] = lines[5].replace("He", "He1")
    path.write_text("\n".join(lines) + "\n")
    _, aliases = read_structure(str(path))
    assert aliases.get("He") == "He1"


def test_unknown_extension_falls_back_to_poscar(tmp_path):
    path = tmp_path / "structure.dat"
    path.write_text(HCP.read_text())
    ref = read_structure(str(HCP))[0]
    struct, _ = read_structure(str(path))
    assert len(struct) == len(ref)


def test_unreadable_file_raises_clear_error(tmp_path):
    path = tmp_path / "nonsense.dat"
    path.write_text("this is not a structure\n")
    with pytest.raises(ValueError, match="structure_format"):
        read_structure(str(path))


def test_missing_file(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_structure(str(tmp_path / "nope.cif"))


def test_partial_occupancy_rejected(tmp_path):
    s = Structure(Lattice.cubic(3.0), [{"Fe": 0.5, "Ni": 0.5}], [(0, 0, 0)])
    path = tmp_path / "dis.json"
    path.write_text(json.dumps(s.as_dict()))
    with pytest.raises(ValueError, match="partially occupied"):
        read_structure(str(path))


def test_oxidation_states_removed_and_matrix_selection_works(tmp_path):
    s = Structure(Lattice.cubic(4.0), ["Ti4+", "O2-"], [(0, 0, 0), (0.5, 0.5, 0.5)])
    path = tmp_path / "ox.cif"
    s.to(filename=str(path))
    loader = PoscarJumpStructureLoader(str(path), crystal_matrix=["Ti"])
    assert [sp.symbol for sp in loader.jump_structure.species] == ["O"]
    assert [sp.symbol for sp in loader.matrix_structure.species] == ["Ti"]


def test_ase_formats(tmp_path):
    ase_io = pytest.importorskip("ase.io")
    pytest.importorskip("pymatgen.io.ase")
    from pymatgen.io.ase import AseAtomsAdaptor
    ref = _reference()
    path = tmp_path / "s.extxyz"
    ase_io.write(str(path), AseAtomsAdaptor.get_atoms(ref))
    struct, _ = read_structure(str(path))
    assert struct.matches(ref)
