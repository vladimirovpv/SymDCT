"""Namespace package for the SymDCT pipeline.

Having this file ensures that ``src`` is recognised as a package when the
test suite is executed via ``python -m pytest`` (or ``uv run pytest``).

The subpackages are:

* ``src.common_utils``
* ``src.jump_network_builder``
* ``src.matrix_builder``

which are all used by the tests.
"""

# Allow the `public` package to resolve imports originating from the
# top-level `src` directory in the repository. When modules inside
# `public` import `src.*`, Python sometimes resolves `src` relative to
# the `public` package.  Ensure the repository's top-level `src`
# directory is on the package search path so imports like
# `from src.common_utils...` work consistently.
import os

# Ensure spglib uses modern error handling to avoid deprecation warnings
# This must be set before any submodule imports that may pull in spglib.
os.environ.setdefault("OLD_ERROR_HANDLING", "false")

_repo_root_src = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if os.path.isdir(_repo_root_src) and _repo_root_src not in __path__:
	__path__.append(_repo_root_src)

