"""Pytest configuration for the SymDCT tests.

The spglib legacy error-handling flag is switched off to avoid the
deprecation message emitted at import time. Two approaches are used:

- Set the environment variable `SPGLIB_OLD_ERROR_HANDLING=0` before any
  imports that may trigger spglib to configure its error handling.
- If spglib is already importable, attempt to flip an in-module flag
  to the modern behaviour.

This addresses the root cause (spglib's legacy error-handling mode)
instead of merely filtering the warning text.
"""
import os

# Ensure the environment variable is set before any test imports occur.
os.environ.setdefault("SPGLIB_OLD_ERROR_HANDLING", "0")


def pytest_configure(config):
    """Attempt to set spglib's runtime flag if the module is available."""
    try:
        import spglib

        # Some builds expose a module-level OLD_ERROR_HANDLING attribute;
        # try to set it to False so spglib stops using legacy paths.
        import spglib.error as spglib_error

        spglib_error.OLD_ERROR_HANDLING = False
    except Exception:
        # spglib not available at configure time; env var will still be used
        # when the module is later imported.
        pass

