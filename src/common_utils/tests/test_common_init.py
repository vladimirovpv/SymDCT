"""Tests for :mod:`src.common_utils.common_init`.

These tests focus on the behaviour of :func:`init_from_yaml`, ensuring
that it derives the ``output/`` directory relative to the YAML input
path as ``yaml_dir / ".." / "output"`` and that the directory is
created.

We monkeypatch :mod:`src.utils` functions so no real command-line
parsing or YAML loading is performed.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from src.common_utils.common_init import init_from_yaml


def test_init_from_yaml_uses_expected_output_dir(monkeypatch, tmp_path):
    """``init_from_yaml`` should create ``output/`` next to the YAML's parent.

    Given a YAML path ``<tmp>/some/input/config.yaml``, the helper is
    expected to return an ``output_dir`` equal to
    ``<tmp>/some/output`` and to ensure that this directory exists.
    """

    # Arrange: craft a fake YAML path under a temporary directory.
    yaml_path = tmp_path / "some" / "input" / "config.yaml"
    yaml_path.parent.mkdir(parents=True, exist_ok=True)
    # init_from_yaml now requires the referenced POSCAR to actually exist.
    (yaml_path.parent / "dummy.poscar").write_text("POSCAR", encoding="utf-8")

    # Monkeypatch the functions that ``init_from_yaml`` imports from
    # :mod:`src.utils` so that no real I/O or argument parsing occurs.
    def fake_get_args():
        return str(yaml_path)

    def fake_get_input_parameters(path):  # noqa: ARG001 - exercised indirectly
        # Return a minimal, but valid, parameter dict.
        return {"fname": "dummy.poscar"}

    monkeypatch.setattr(
        "src.common_utils.common_init.get_args",
        fake_get_args,
        raising=True,
    )
    monkeypatch.setattr(
        "src.common_utils.common_init.get_input_parameters",
        fake_get_input_parameters,
        raising=True,
    )

    # Act
    inp_file, params, output_dir = init_from_yaml()

    # Assert
    assert Path(inp_file) == yaml_path
    assert params["fname"] == str((yaml_path.parent / "dummy.poscar").resolve())

    expected_output_dir = yaml_path.parent.parent / "output"
    assert output_dir == expected_output_dir
    assert output_dir.is_dir()


def test_init_from_yaml_raises_when_fname_missing(monkeypatch, tmp_path):
    """``fname`` is used exactly as written in the YAML — no prefix fallback.

    If the referenced POSCAR does not exist next to the YAML (nor relative
    to the current working directory), init_from_yaml raises
    FileNotFoundError rather than guessing at alternate (e.g. "jump_"
    prefixed) filenames.
    """

    yaml_path = tmp_path / "examples" / "input" / "config.yaml"
    yaml_path.parent.mkdir(parents=True, exist_ok=True)
    # Only a differently-named file exists; "dummy.poscar" itself is absent.
    (yaml_path.parent / "jump_dummy.poscar").write_text("POSCAR", encoding="utf-8")

    monkeypatch.setattr("src.common_utils.common_init.get_args", lambda: str(yaml_path), raising=True)
    monkeypatch.setattr(
        "src.common_utils.common_init.get_input_parameters",
        lambda path: {"fname": "dummy.poscar"},
        raising=True,
    )

    with pytest.raises(FileNotFoundError):
        init_from_yaml()
