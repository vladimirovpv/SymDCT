"""High–level helpers for constructing and reducing transition matrices.

The :mod:`matrix_builder` package encapsulates the *matrix stage* of the
workflow: given a YAML input file describing a jump network problem, it

* loads Stage-1 artifacts (`jump_structure.node_link.json` and
  `jump_sites.poscar`) and reconstructs the jump graph plus symmetry context;
* constructs the full symbolic transition–probability matrix; and
* optionally reduces this matrix along selected diffusion axes using the
  symmetry tools from :mod:`src.symmetry_context` and
  :mod:`src.matrix_builder.matrix_reduction_helpers`.

The public API is deliberately small:

* :func:`build_matrices_from_yaml` – one–shot helper for Python callers.
* :func:`main` – CLI entry point that mirrors the initialization pattern
  used in :mod:`src.jump_network_builder.jump_network_builder`.

This module does *not* perform any heavy I/O or plotting; it focuses on
matrix construction and optimization and can be reused from scripts,
tests or notebooks.
"""

from .matrix_builder import MatrixBuildResult, build_matrices_from_yaml, main

__all__ = [
    "MatrixBuildResult",
    "build_matrices_from_yaml",
    "main",
]
