from __future__ import annotations

from pathlib import Path

import sympy as sp

import main as public_main
from unittest.mock import patch
from src.common_utils.yaml_with_comments import load_data


class _FakeJumpNetwork:
    pass


class _FakeSymmetryContext:
    """Minimal stand-in exposing the methods main() actually calls on it."""

    def is_conventional_cell(self) -> bool:
        return True

    def get_lattice_type(self) -> str:
        return "cubic"

    def get_symbol_alias_map(self) -> dict:
        return {}


class _FakeStage2Result:
    def __init__(self) -> None:
        self.full_tp_matrix = sp.eye(2)
        self.axis_results = {
            (0, 0, 1): type(
                "_AxisResult",
                (),
                {
                    "reordered_tp_matrix": sp.eye(2),
                    "min_jump_length": sp.Rational(1, 2),
                    "interplanar_distance_squared": sp.Symbol("d_sq"),
                    "labels": [("A", "0"), ("B", "1")],
                    "hkl_structure": (0, 0, 1),
                    "hkl_conventional": (0, 0, 1),
                },
            )()
        }


_ENERGETICS_YAML = (
    "settings:\n  length_tolerance: 0.01\n  units:\n    energy: eV\n    length: angstrom\n    frequency: THz\n"
    "lattice constants:\n  a: 4.0\n  b: 4.0\n  c: 4.0\nsite frequencies:\n  default_frequency: 10.0\n  A: 5.0\n"
    "  B: 5.0\nsite energies:\n  A: 0.0\n  B: 0.1\njumps:\n  A-B:\n    - jump_id: {forward: 0, backward: 1}\n"
    "      length: 2.5\n      energy_sp: 0.5\nsymbol_alias_map: {}\n"
)


def test_main_runs_stage3(monkeypatch, tmp_path, capsys) -> None:
    # Use tmp_path to keep test artifacts isolated from repo root
    out_dir = tmp_path / "out"
    input_dir = tmp_path / "input"
    out_dir.mkdir(parents=True, exist_ok=True)
    input_dir.mkdir(parents=True, exist_ok=True)

    yaml_path = input_dir / "input.yaml"
    yaml_path.write_text("fname: POSCAR\nmatrix: []\n")

    monkeypatch.setattr(public_main, "verbose", True)
    monkeypatch.setattr(public_main, "init_from_yaml", lambda: (str(yaml_path), {"fname": "POSCAR"}, out_dir))
    monkeypatch.setattr(public_main, "build_jump_network_from_poscar", lambda *args, **kwargs: (None, None, None))
    monkeypatch.setattr(public_main, "JumpNetworkSymmetryContext", lambda *args, **kwargs: _FakeSymmetryContext())
    monkeypatch.setattr(public_main, "build_matrices", lambda *args, **kwargs: _FakeStage2Result())
    monkeypatch.setattr(
        public_main,
        "solve_from_matrix_result",
        lambda result: {(0, 0, 1): sp.Symbol("W0")},
    )

    # energetics_template.yaml and dc_symbolic files are Stage-1/2 output
    # artifacts, so they live in out_dir; the production energetics.yaml is
    # hand-authored and lives next to the input YAML (main.py looks for it
    # at yaml_path.parent / "energetics.yaml").
    (out_dir / "energetics_template.yaml").write_text(_ENERGETICS_YAML)
    (out_dir / "dc_symbolic_0_0_1.txt").write_text("W0\n")
    (input_dir / "energetics.yaml").write_text(_ENERGETICS_YAML)

    # Mock the evaluator to avoid numeric evaluation of symbolic expressions.
    # main.py imports these two names directly from
    # src.transport_coefficient_evaluator.evaluator (Stage 4), so they are
    # patched as attributes of the main module itself.
    with patch.object(public_main, "load_evaluation_context", lambda *a, **k: object()), \
         patch.object(public_main, "evaluate_dc_file", lambda dc_file, ctx: Path(str(dc_file) + ".png")):
        assert public_main.main() == 0
    out = capsys.readouterr().out
    assert "interplanar_distance_squared=d_sq" in out
    assert "effective_transport_distance_squared=d_sq/4" in out
    assert "Solved symbolic diffusion coefficients for axes:" in out
    assert "(0, 0, 1): W0" in out


def test_bcc_int_main_accounts_for_min_jump_length(monkeypatch, tmp_path, capsys) -> None:
    example_dir = Path(__file__).resolve().parents[1] / "examples" / "bcc_int" / "input"
    input_dir = tmp_path / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    yaml_path = input_dir / "bcc_int.yaml"
    # bcc_int.yaml's `fname` points at prim_int_bcc.poscar (the actual input
    # structure); jump_prim_int_bcc.poscar is a Stage-1 *output* artifact
    # (the refined jump-sites structure), not an input file.
    poscar_path = input_dir / "prim_int_bcc.poscar"
    yaml_path.write_text((example_dir / "bcc_int.yaml").read_text(encoding="utf-8"), encoding="utf-8")
    poscar_path.write_text((example_dir / "prim_int_bcc.poscar").read_text(encoding="utf-8"), encoding="utf-8")

    params = load_data(str(yaml_path))
    monkeypatch.setattr(public_main, "verbose", True)
    monkeypatch.setattr("src.common_utils.common_init.get_args", lambda: str(yaml_path), raising=True)
    monkeypatch.setattr(
        "src.common_utils.common_init.get_input_parameters",
        lambda path: dict(params),
        raising=True,
    )

    # main()'s Stage 4 requires a hand-authored production energetics.yaml to
    # exist next to the input YAML *before* it runs (main.py looks for it at
    # yaml_path.parent / "energetics.yaml"); pre-seed it here from the
    # example's own fixture.
    (input_dir / "energetics.yaml").write_text(
        (example_dir / "energetics.yaml").read_text(encoding="utf-8"),
        encoding="utf-8",
    )

    # Mock the evaluator to avoid numeric evaluation of symbolic expressions.
    with patch.object(public_main, "load_evaluation_context", lambda *a, **k: object()), \
         patch.object(public_main, "evaluate_dc_file", lambda dc_file, ctx: Path(str(dc_file) + ".png")):
        assert public_main.main() == 0
    out = capsys.readouterr().out
    # axis (1,1,1): interplanar_distance_squared=a**2/12, min_jump_length=1/2,
    # so effective_transport_distance_squared = a**2/12 * (1/2)**2 = a**2/48
    assert "interplanar_distance_squared=a**2/12" in out
    assert "effective_transport_distance_squared=a**2/48" in out
    assert "W_{TO}" in out
    assert "W_{TT}" in out
    assert "W_{H" not in out


def test_main_skips_stage4_and_warns_when_energetics_missing(monkeypatch, tmp_path, capsys) -> None:
    out_dir = tmp_path / "out"
    input_dir = tmp_path / "input"
    out_dir.mkdir(parents=True, exist_ok=True)
    input_dir.mkdir(parents=True, exist_ok=True)

    yaml_path = input_dir / "input.yaml"
    yaml_path.write_text("fname: POSCAR\nmatrix: []\n")

    monkeypatch.setattr(public_main, "verbose", True)
    monkeypatch.setattr(public_main, "init_from_yaml", lambda: (str(yaml_path), {"fname": "POSCAR"}, out_dir))
    monkeypatch.setattr(public_main, "build_jump_network_from_poscar", lambda *args, **kwargs: (None, None, None))
    monkeypatch.setattr(public_main, "JumpNetworkSymmetryContext", lambda *args, **kwargs: _FakeSymmetryContext())
    monkeypatch.setattr(public_main, "build_matrices", lambda *args, **kwargs: _FakeStage2Result())
    monkeypatch.setattr(
        public_main,
        "solve_from_matrix_result",
        lambda result: {(0, 0, 1): sp.Symbol("W0")},
    )

    # No energetics.yaml written to input_dir: main() should warn and return
    # 0 without attempting Stage 4, rather than raising.
    with patch.object(public_main, "load_evaluation_context") as load_ctx, \
         patch.object(public_main, "evaluate_dc_file") as evaluate_dc:
        assert public_main.main() == 0
    load_ctx.assert_not_called()
    evaluate_dc.assert_not_called()
    out = capsys.readouterr().out
    assert "energetics" in out
    assert str(input_dir / "energetics.yaml") in out
