"""Tests for the classification of jumps by space-group orbits."""

import json
import logging
import pathlib as pl

import numpy as np
import pytest
from pymatgen.core import Lattice, Structure
from pymatgen.io.vasp import Poscar

from ..geometry import build_jump_network_from_poscar
from ..jump_orbits import classify_jump_orbits
from src.common_utils.yaml_io import load_data

EXAMPLES = pl.Path(__file__).resolve().parents[3] / "examples"
EXAMPLE_NAMES = ["Be12Ti_intHe", "bcc_int", "hcp_int_OT", "hcp_vac", "host-host_diplon"]


def _hexagonal_poscar(tmp_path: pl.Path, c_over_a: float) -> str:
    a = 2.0
    struct = Structure(Lattice.hexagonal(a, c_over_a * a), ["X"], [(0.0, 0.0, 0.0)])
    path = tmp_path / "POSCAR"
    Poscar(struct).write_file(path)
    return str(path)


def _partition(ids):
    classes = {}
    for k, i in enumerate(ids):
        classes.setdefault(i, []).append(k)
    return sorted(classes.values())


@pytest.mark.parametrize("name", EXAMPLE_NAMES)
def test_orbits_reproduce_triplet_classes_of_shipped_examples(name):
    """On the reference examples both classifications give the same partition."""

    assets = EXAMPLES / name / "assets"
    if not (assets / "jump_structure.node_link.json").is_file():
        pytest.skip(f"no reference assets for {name}")
    graph = json.loads((assets / "jump_structure.node_link.json").read_text())
    structure = Poscar.from_file(assets / "jump_sites.poscar").structure
    labels = [n["lbl"] for n in sorted(graph["nodes"], key=lambda n: n["id"])]
    edges = [
        (e["source"], e["target"], [h["num"] / h["den"] for h in e["h"]])
        for e in graph["edges"]
    ]
    orbits = classify_jump_orbits(
        structure.lattice.matrix, structure.frac_coords, labels, edges, symprec=0.1
    )
    assert orbits.n_unmapped_edges == 0
    assert orbits.n_skipped_ops == 0
    assert _partition(orbits.ids) == _partition([e["glob_jump_cnt"] for e in graph["edges"]])


def test_tetragonal_lattice_gives_two_orbits():
    lattice = np.diag([2.0, 2.0, 3.0])
    edges = [(0, 0, h) for h in
             ([1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1])]
    orbits = classify_jump_orbits(lattice, np.zeros((1, 3)), ["A"], edges)
    assert _partition(orbits.ids) == [[0, 1, 2, 3], [4, 5]]
    assert orbits.reverse_ids == {0: 0, 1: 1}


def test_equal_length_inequivalent_jumps_are_separated_by_orbits(tmp_path, caplog):
    """Simple hexagonal lattice with c = a: six in-plane and two axial jumps
    have the same triplet but are not related by symmetry."""

    poscar = _hexagonal_poscar(tmp_path, c_over_a=1.0)
    cutoffs = {"bonds": {"X-X": 2.1}}

    with caplog.at_level(logging.WARNING):
        G, jr, _ = build_jump_network_from_poscar(
            poscar, cutoffs, str(tmp_path / "orbit"), max_jump_length=2.1,
            jump_classification="orbit")
    classes = sorted({d["glob_jump_cnt"] for *_e, d in G.edges(data=True)})
    assert len(classes) == 2
    assert jr.get_number_of_jumps() == 2
    sizes = sorted(sum(1 for *_e, d in G.edges(data=True) if d["glob_jump_cnt"] == c)
                   for c in classes)
    assert sizes == [2, 6]
    assert "inequivalent orbits" in caplog.text

    template = load_data(str(tmp_path / "orbit" / "energetics_template.yaml"))
    entries = template["jumps"]["X-X"]
    assert len(entries) == 2
    assert {e["jump_id"]["forward"] for e in entries} == {0, 1}

    G_t, jr_t, _ = build_jump_network_from_poscar(
        poscar, cutoffs, str(tmp_path / "triplet"), max_jump_length=2.1,
        jump_classification="triplet")
    assert {d["glob_jump_cnt"] for *_e, d in G_t.edges(data=True)} == {0}
    assert jr_t.get_number_of_jumps() == 1


def test_unknown_classification_is_rejected(tmp_path):
    poscar = _hexagonal_poscar(tmp_path, c_over_a=1.0)
    with pytest.raises(ValueError):
        build_jump_network_from_poscar(poscar, {}, str(tmp_path / "o"),
                                       jump_classification="hash")
