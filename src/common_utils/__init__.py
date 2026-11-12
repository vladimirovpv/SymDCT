"""Utility package for the jump-network and matrix-builder tooling.

Currently exposes :mod:`common_init`, which holds helpers shared between
the orchestration scripts of this standalone project, without adding
top-level modules to the parent project's ``src`` directory.
"""

from .common_init import init_from_yaml
from .jump_vector_codec import decode_jump_vector, encode_jump_vector

__all__ = ["init_from_yaml", "encode_jump_vector", "decode_jump_vector"]
