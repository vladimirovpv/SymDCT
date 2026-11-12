import os
import pathlib as pl

import networkx as nx

from ..geometry import build_jump_network_from_poscar
from src.common_utils.yaml_io import load_data


def _write_minimal_cutoffs_yaml(tmp_path: pl.Path) -> str:
    """Write a minimal cutoffs.yaml compatible with geometry.

    Uses a single species "X" with a modest cutoff so the neighbour
    search in the bcc-like test structure finds at least one jump.
    """

    content = """\
units:
  length: angstrom

bonds:
  "X-X": 2.6
"""
    path = tmp_path / "cutoffs.yaml"
    path.write_text(content)
    return str(path)


def _make_simple_bcc_poscar(tmp_path: pl.Path) -> str:
    # Reuse the helper from the other test module if tests are run as a
    # package; fall back to a local copy otherwise to keep this file
    # self-contained.
    try:
        from .test_loader_and_symmetry import _make_simple_bcc_structure

        return _make_simple_bcc_structure(tmp_path)
    except Exception:  # pragma: no cover - fallback path
        import pymatgen as pmg
        from pymatgen.core import Lattice, Structure

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


def pathlib_name(path):
    import os
    return os.path.basename(str(path))


def test_build_jump_network_from_poscar_creates_graph_and_template(tmp_path, monkeypatch):
    """End-to-end smoke test for geometry + energetics template pipeline.

    This uses a minimal synthetic POSCAR and a tiny cutoffs.yaml to
    exercise:

    * JumpNetworkSymmetryContext (via geometry.build_jump_network_from_poscar)
    * neighbour search and cutoffs filtering
    * JumpRegister integration (glob_jump_cnt on edges)
    * energetics template generation and YAML dumping
    """

    poscar_path = _make_simple_bcc_poscar(tmp_path)
    cutoffs_path = _write_minimal_cutoffs_yaml(tmp_path)

    cutoffs = load_data(cutoffs_path)

    output_dir = tmp_path / "output"
    G, jr, jump_sites = build_jump_network_from_poscar(
        poscar_path,
        cutoffs,
        str(output_dir),
        crystal_matrix=None,
        max_jump_length=2.6,
    )

    # Basic graph sanity: we should have nodes and at least one edge.
    assert isinstance(G, nx.MultiDiGraph)
    assert G.number_of_nodes() == jump_sites.num_sites
    assert G.number_of_edges() >= 1

    # Every edge must carry a global jump index and a length.
    for _, _, data in G.edges(data=True):
        assert "glob_jump_cnt" in data
        assert "nn_dist" in data

    # Graph exports and energetics template must have been written.
    node_link_path = output_dir / "jump_structure.node_link.json"
    adjlist_path = output_dir / "jump_structure.adjlist"
    jump_sites_poscar_path = output_dir / "jump_sites.poscar"
    energetics_path = output_dir / "energetics_template.yaml"
    assert node_link_path.is_file()
    assert adjlist_path.is_file()
    assert jump_sites_poscar_path.is_file()
    assert energetics_path.is_file()

    # The header of jump_sites.poscar names the input file without its directory.
    header = jump_sites_poscar_path.read_text(encoding="utf-8").splitlines()[0]
    assert header == f"Refined jump sites from {pathlib_name(poscar_path)}"
    assert str(tmp_path) not in header

    # Energetics YAML should be parseable and contain the expected keys.
    energetics = load_data(str(energetics_path))
    assert "settings" in energetics
    assert "site frequencies" in energetics
    assert "site energies" in energetics
    assert "jumps" in energetics
