"""Disk I/O for Stage-3 transport-coefficient payloads and results."""

from __future__ import annotations

import pickle
from pathlib import Path
from typing import Mapping

import sympy as sp

from .contracts import AxisMatrixPayload


AXIS_PAYLOAD_PREFIX = "tc_axis_payload"


def axis_payload_filename(axis_hkl_structure: tuple[int, int, int]) -> str:
    h, k, l = (int(v) for v in axis_hkl_structure)
    return f"{AXIS_PAYLOAD_PREFIX}_{h}_{k}_{l}.pickle"


def write_axis_payload(payload: AxisMatrixPayload, output_file: str | Path) -> Path:
    path = Path(output_file)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as fh:
        pickle.dump(payload.to_pickle_payload(), fh)
    return path


def read_axis_payload(payload_file: str | Path) -> AxisMatrixPayload:
    path = Path(payload_file)
    with path.open("rb") as fh:
        data = pickle.load(fh)
    return AxisMatrixPayload.from_pickle_payload(data)


def write_axis_payloads(
    payloads: Mapping[tuple[int, int, int], AxisMatrixPayload],
    output_dir: str | Path,
) -> list[Path]:
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    files: list[Path] = []
    for axis, payload in sorted(payloads.items()):
        files.append(write_axis_payload(payload, output_path / axis_payload_filename(axis)))
    return files


def read_axis_payloads(output_dir: str | Path) -> dict[tuple[int, int, int], AxisMatrixPayload]:
    """Load Stage-2 payload pickles, keyed by ``axis_hkl_conventional``.

    Conventional is the canonical, always-available identifier (see
    ``AxisMatrixPayload``), consistent with how in-memory payloads are
    keyed by ``payloads_from_matrix_build_result``.
    """
    output_path = Path(output_dir)
    payload_files = sorted(output_path.glob(f"{AXIS_PAYLOAD_PREFIX}_*.pickle"))
    payloads: dict[tuple[int, int, int], AxisMatrixPayload] = {}
    for payload_file in payload_files:
        payload = read_axis_payload(payload_file)
        payloads[payload.axis_hkl_conventional] = payload
    return payloads


def write_diffusion_coefficients(
    diffusion_coefficients: Mapping[tuple[int, int, int], sp.Expr],
    output_dir: str | Path,
    prefix: str = "dc_symbolic",
) -> list[Path]:
    """Write dc_symbolic_<h>_<k>_<l>.txt/.m files, one pair per axis.

    ``diffusion_coefficients`` is keyed by ``axis_hkl_conventional`` (see
    ``solve_all_axes``) -- the plane indices exactly as given in the input
    YAML, and always available regardless of whether the POSCAR is
    primitive or conventional -- so the dict key is used directly for the
    filename.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for axis, expr in sorted(diffusion_coefficients.items()):
        axis_label = "_".join(str(int(v)) for v in axis)
        txt_file = output_path / f"{prefix}_{axis_label}.txt"
        txt_file.write_text(str(expr), encoding="utf-8")
        m_file = output_path / f"{prefix}_{axis_label}.m"
        m_file.write_text(str(expr).replace("**", "^"), encoding="utf-8")
        written.extend([txt_file, m_file])
    return written

