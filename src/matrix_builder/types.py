"""Data classes used by the matrix builder.

The original implementation defined ``AxisMatrixResult``, ``MatrixBuildResult``
and ``MatrixStageJumpNetwork`` directly in ``matrix_builder.py``.  Moving them to
their own module reduces the size of the main file and makes the public API
clearer.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence, Tuple, List, Optional

import sympy as sp

# Local imports required for the ``MatrixStageJumpNetwork`` implementation
from src.jump_network.core import JumpNetworkGraph
from src.common_utils.jump_vector_codec import decode_jump_vector
from src.common_utils.tweak_rates import normalize_tweak_rates
from src.common_utils.jump_probabilities import JumpProbabilities
from src.jump_network.transition_matrix import TransitionMatrixBuilder
from src.matrix_builder.matrix_reduction_helpers import (
    find_equiv_and_reduce_matrix,
    find_min_jump_length,
    prepare_label_data_for_reordered_matrix,
    prepare_reordered_tp_matrix,
)
from src.jump_network.symbolic_context import SymbolicContext
from src.symmetry_context.symmetry_refiner import JumpStructureSymmetryRefiner

# Tolerance (in the POSCAR's length unit, typically Angstrom) used when
# matching a tweak_transition_rates length qualifier (e.g. "A->B:2.213")
# against actual jump lengths. Kept separate from SymbolicContext.tolerance,
# which is the much tighter tolerance nsimplify uses to rationalize jump
# *vector* components -- these are unrelated concerns despite both being
# "tolerances". 0.005 comfortably covers jump lengths written to 3 decimal
# places in params YAML files (max rounding error 0.0005) with margin to
# spare.
TWEAK_LENGTH_TOLERANCE = 0.005


@dataclass
class AxisMatrixResult:
    """Per-axis Stage-2 result.

    ``hkl_conventional`` -- the plane indices exactly as given in the input
    YAML ``axes`` entry -- is the canonical, always-available identifier for
    this axis (see ``MatrixBuildResult.axis_results``); it is never
    ``None``, since ``JumpNetworkSymmetryContext`` always refines to a
    conventional structure regardless of the input POSCAR. ``hkl_structure``
    is the same axis converted into whatever basis the working structure
    actually uses internally (identical to ``hkl_conventional`` unless that
    structure is primitive), needed for the geometric calculations in this
    module but not for output labelling/filenames.
    """

    hkl_structure: Tuple[int, int, int]
    hkl_conventional: Tuple[int, int, int]
    equiv_sites: List[int]
    reordered_tp_matrix: sp.Matrix
    labels: Sequence[str]
    min_jump_length: Optional[sp.Expr] = None
    interplanar_distance_squared: Optional[sp.Expr] = None
    diffusion_coefficient: Optional[sp.Expr] = None


@dataclass
class MatrixBuildResult:
    jump_network: Any
    symmetry_ctx: Any
    labels: Sequence[str]
    full_tp_matrix: sp.Matrix
    # Keyed by hkl_conventional (see AxisMatrixResult), not hkl_structure --
    # this is the canonical axis identifier used for all output labelling.
    axis_results: Mapping[Tuple[int, int, int], AxisMatrixResult]


class MatrixStageJumpNetwork:
    """Minimal JumpNetwork‑like adapter used by Stage 2.

    This is a direct copy of the original implementation from
    ``matrix_builder.py``.  It provides the same public API required by the
    rest of the pipeline while living in a dedicated module.
    """

    def __init__(
        self,
        graph: Any,
        jump_sites: Any,
        node_labels: Sequence[str],
        lat_const: sp.Matrix,
        number_of_jumps: int,
    ) -> None:
        # The conventional→primitive matrix is now accessed via the centralised
        # ``SMatrixHandler`` on the symmetry context; it is no longer stored here.
        self.jump_sites = jump_sites
        self.node_labels = list(node_labels)
        self.lat_const = lat_const
        self.symbolic_ctx = SymbolicContext(tolerance=1e-6)
        self.jump_probabilities = JumpProbabilities(number_of_jumps)
        self._full_tp_matrix: sp.Matrix | None = None

        self.graph = JumpNetworkGraph(jump_sites.num_sites, self.node_labels)
        for start, end, _key, attrs in graph.edges(keys=True, data=True):
            edge_id = f"({start}, {end}, {_key})"
            self.graph.add_edge(
                int(start),
                int(end),
                jump_vector=decode_jump_vector(attrs["h"], edge_id=edge_id),
                distance=float(attrs["nn_dist"]),
                image=list(attrs["image"]),
                global_jump_index=int(attrs["glob_jump_cnt"]),
            )

    def _build_label_records_map(self) -> dict:
        """Build map from canonical label -> list of records.
        
        Returns dict: 'A->B' -> [{'glob': idx, 'length': distance}, ...]
        """
        label_map = {}
        for start, end, _key, attrs in self.graph.graph.edges(keys=True, data=True):
            start_label = self.node_labels[int(start)]
            end_label = self.node_labels[int(end)]
            canonical = f"{start_label}->{end_label}"
            rec = {
                'glob': int(attrs['glob_jump_cnt']),
                'length': float(attrs['nn_dist'])
            }
            label_map.setdefault(canonical, []).append(rec)
        for k in label_map:
            label_map[k] = sorted(label_map[k], key=lambda r: r['glob'])
        return label_map

    @staticmethod
    def _parse_label_spec(spec: str) -> Tuple[str, Optional[float], bool]:
        """Parse label spec like 'A->B', 'A->B:length=1.234' or 'A->B: all'."""
        if ':' not in spec:
            return spec, None, False
        parts = spec.split(':', 1)
        label = parts[0].strip()
        qualifier = parts[1].strip().lower()
        if qualifier == 'all':
            return label, None, True
        if qualifier.startswith('length=') or qualifier.startswith('len=') or qualifier.startswith('l='):
            try:
                length_str = qualifier.split('=', 1)[1]
                length = float(length_str)
                return label, length, False
            except (ValueError, IndexError):
                pass
        return label, None, False

    def _normalize_tweak_rates(self, tweak_rates: dict) -> dict:
        """Instance-bound wrapper forwarding to ``normalize_tweak_rates``."""
        return normalize_tweak_rates(
            tweak_rates,
            self._build_label_records_map(),
            TWEAK_LENGTH_TOLERANCE,
        )

    def get_k(self) -> List[sp.Symbol]:
        return self.symbolic_ctx.k

    def get_t(self) -> sp.Symbol:
        return self.symbolic_ctx.t

    def get_full_tp_matrix(self, subst_w: dict | None = None) -> sp.Matrix:
        if self._full_tp_matrix is not None:
            return self._full_tp_matrix
        builder = TransitionMatrixBuilder(
            self.graph,
            self.symbolic_ctx,
            self.jump_probabilities,
        )
        if subst_w:
            normalized_w = self._normalize_tweak_rates(dict(subst_w))
            self.jump_probabilities.customize_jump_probabilities(normalized_w)
            self._full_tp_matrix = builder.fill_matrix(normalized_w)
        else:
            self._full_tp_matrix = builder.fill_matrix({})
        return self._full_tp_matrix

    def compute_squared_interplanar_distance(self, hkl_structure: sp.Matrix) -> sp.Expr:
        """Return the squared interplanar distance for conventional Miller indices."""
        lattice_matrix = sp.Matrix(self.lat_const)
        g_tensor = lattice_matrix * lattice_matrix.T
        interplanar_distance_squared = 1 / (hkl_structure.T * g_tensor.inv() * hkl_structure)[0, 0]
        return sp.simplify(interplanar_distance_squared)
