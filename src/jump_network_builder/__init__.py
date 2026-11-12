"""jump_network_builder package - geometry stage only

This package provides helpers to build a jump network (connectivity graph)
from a POSCAR and a cutoffs YAML. It registers every discovered jump using
JumpRegister and writes out a jump template (energetics) file without
energies so that energetics can be filled later.
"""
from .geometry import build_jump_network_from_poscar, pretty_s_matrix
from .jump_network_builder import main as run_builder

__all__ = ["build_jump_network_from_poscar", "pretty_s_matrix", "run_builder"]

