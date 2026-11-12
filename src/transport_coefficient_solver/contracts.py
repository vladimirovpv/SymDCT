"""Stage-3 data contracts for symbolic transport-coefficient solving."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Optional, Sequence, Tuple

import sympy as sp


PICKLE_CONTRACT_VERSION = 1


@dataclass(frozen=True)
class AxisMatrixPayload:
    """Per-axis Stage-2 payload consumed by the transport-coefficient solver.

    ``axis_hkl_conventional`` -- the plane indices exactly as given in the
    input YAML ``axes`` entry -- is the canonical, always-available
    identifier for this axis; it is never ``None``, since
    ``JumpNetworkSymmetryContext`` always refines to a conventional
    structure regardless of the input POSCAR. It is what output filenames
    and dict keys should use. ``axis_hkl_structure`` is the same axis in
    whatever basis the working structure actually used internally
    (identical to ``axis_hkl_conventional`` unless that structure is
    primitive), needed only for reproducing the Stage-2 geometry, not for
    labelling.
    """

    axis_hkl_structure: Tuple[int, int, int]
    axis_hkl_conventional: Tuple[int, int, int]
    transition_matrix: sp.Matrix
    labels: Tuple[str, ...]
    min_jump_length: Optional[sp.Expr] = None
    interplanar_distance_squared: Optional[sp.Expr] = None
    t_symbol: Optional[sp.Symbol] = None
    w_substitutions: Optional[Mapping[sp.Symbol, sp.Expr]] = None
    jump_label_map: Optional[Mapping[int, str]] = None
    jump_label_display_map: Optional[Mapping[str, str]] = None

    def to_pickle_payload(self) -> dict[str, Any]:
        """Return a validated pickle-ready payload dictionary."""
        return {
            "contract": "transport-coefficient-axis-payload",
            "contract_version": PICKLE_CONTRACT_VERSION,
            "axis_hkl_structure": tuple(int(v) for v in self.axis_hkl_structure),
            "axis_hkl_conventional": tuple(int(v) for v in self.axis_hkl_conventional),
            "transition_matrix": self.transition_matrix,
            "labels": list(self.labels),
            "min_jump_length": self.min_jump_length,
            "interplanar_distance_squared": self.interplanar_distance_squared,
            "t_symbol": self.t_symbol,
            "w_substitutions": dict(self.w_substitutions or {}),
            "jump_label_map": dict(self.jump_label_map or {}),
            "jump_label_display_map": dict(self.jump_label_display_map or {}),
        }

    @classmethod
    def from_pickle_payload(cls, payload: Mapping[str, Any]) -> "AxisMatrixPayload":
        """Parse and validate payload dictionary loaded from pickle."""
        if payload.get("contract") != "transport-coefficient-axis-payload":
            raise ValueError("Unsupported payload contract type")
        version = int(payload.get("contract_version", 0))
        if version != PICKLE_CONTRACT_VERSION:
            raise ValueError(
                f"Unsupported payload version {version}; expected {PICKLE_CONTRACT_VERSION}"
            )

        axis = payload.get("axis_hkl_structure")
        if not isinstance(axis, (list, tuple)) or len(axis) != 3:
            raise ValueError("axis_hkl_structure must be a 3-component sequence")
        axis_tuple = (int(axis[0]), int(axis[1]), int(axis[2]))

        matrix = payload.get("transition_matrix")
        if not isinstance(matrix, sp.MatrixBase):
            raise ValueError("transition_matrix must be a sympy Matrix")
        matrix = sp.Matrix(matrix)

        labels_raw = payload.get("labels")
        if not isinstance(labels_raw, Sequence):
            raise ValueError("labels must be a sequence")
        labels = tuple(str(lbl) for lbl in labels_raw)
        if matrix.shape[0] != len(labels) or matrix.shape[1] != len(labels):
            raise ValueError(
                "labels length must match transition_matrix dimensions (square matrix expected)"
            )

        interplanar = payload.get("interplanar_distance_squared")
        if interplanar is not None:
            interplanar = sp.nsimplify(interplanar)

        min_jump_length = payload.get("min_jump_length")
        if min_jump_length is not None:
            min_jump_length = sp.nsimplify(min_jump_length)

        t_symbol = payload.get("t_symbol")
        if t_symbol is not None and not isinstance(t_symbol, sp.Symbol):
            raise ValueError("t_symbol must be a sympy.Symbol when provided")

        w_substitutions_raw = payload.get("w_substitutions") or {}
        if not isinstance(w_substitutions_raw, Mapping):
            raise ValueError("w_substitutions must be a mapping")
        w_substitutions = {
            sp.Symbol(str(k)) if not isinstance(k, sp.Symbol) else k: sp.nsimplify(v)
            for k, v in w_substitutions_raw.items()
        }

        jump_label_map_raw = payload.get("jump_label_map") or {}
        if not isinstance(jump_label_map_raw, Mapping):
            raise ValueError("jump_label_map must be a mapping")
        jump_label_map = {
            int(k): str(v) for k, v in jump_label_map_raw.items()
        }

        jump_label_display_map_raw = payload.get("jump_label_display_map") or {}
        if not isinstance(jump_label_display_map_raw, Mapping):
            raise ValueError("jump_label_display_map must be a mapping")
        jump_label_display_map = {
            str(k): str(v) for k, v in jump_label_display_map_raw.items()
        }

        conventional = payload.get("axis_hkl_conventional")
        if not isinstance(conventional, (list, tuple)) or len(conventional) != 3:
            raise ValueError("axis_hkl_conventional must be a 3-component sequence")
        conventional_tuple = (int(conventional[0]), int(conventional[1]), int(conventional[2]))

        return cls(
            axis_hkl_structure=axis_tuple,
            axis_hkl_conventional=conventional_tuple,
            transition_matrix=matrix,
            labels=labels,
            min_jump_length=min_jump_length,
            interplanar_distance_squared=interplanar,
            t_symbol=t_symbol,
            w_substitutions=w_substitutions,
            jump_label_map=jump_label_map if jump_label_map else None,
            jump_label_display_map=jump_label_display_map if jump_label_display_map else None,
        )


@dataclass(frozen=True)
class AxisDiffusionResult:
    """Stage-3 symbolic result for a specific diffusion axis.

    See ``AxisMatrixPayload`` for the ``axis_hkl_conventional`` vs.
    ``axis_hkl_structure`` distinction; ``axis_hkl_conventional`` is the
    canonical identifier ``solve_all_axes`` keys its result dict by.
    """

    axis_hkl_structure: Tuple[int, int, int]
    axis_hkl_conventional: Tuple[int, int, int]
    diffusion_coefficient: sp.Expr
