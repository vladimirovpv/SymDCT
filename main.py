"""Top-level orchestration script for the SymDCT pipeline.

For now this module serves primarily as an *executable example* of how to
wire together the geometry (jump_network_builder) and matrix-builder
stages in a single in-memory pipeline.

The idea is:

1. Use a YAML file compatible with ``jump_network_builder`` to describe
   the structure (POSCAR), matrix species and axes of interest.
2. Run the geometry stage to build the jump network and adjacency list.
3. Immediately feed the resulting graph and symmetry context into the
   matrix builder.
4. Pass the stage-2 matrix result into the transport-coefficient solver.
5. Feed the symbolic diffusion coefficient(s) into the transport-coefficient
   evaluator to generate Arrhenius plots.

This mirrors the two entry points supported by
``src.matrix_builder.matrix_builder.build_matrices``:

* ``from files`` – via the dedicated CLI in
  :mod:`src.matrix_builder.matrix_builder`;
* ``from memory`` – demonstrated by :func:`main` below.
"""

from __future__ import annotations

import logging
from pathlib import Path

from src.common_utils.common_init import init_from_yaml
from src.common_utils.utils import configure_cli_logging, format_axis_label, load_cutoffs_config
from src.jump_network_builder.geometry import build_jump_network_from_poscar
from src.matrix_builder.matrix_builder import build_matrices, print_axis_build_summary
from src.symmetry_context.jump_symmetry_context import JumpNetworkSymmetryContext
from src.transport_coefficient_solver.transport_coefficient_solver import solve_from_matrix_result
from scripts.print_transition_matrices import print_transition_matrix_table
from src.transport_coefficient_evaluator.evaluator import (
    evaluate_dc_file,
    load_evaluation_context,
)

export_dc = True  # Set to False to skip writing diffusion coefficients to disk
verbose = False  # suppress verbose output for now; can be toggled via CLI args in future

logger = logging.getLogger(__name__)

def print_header(header: str):
    """Print a header for SymDCT or a single module."""
    print()
    print("=" * 80)
    print(header)
    print("=" * 80)

def main() -> int:
    """Run geometry and matrix stages in a single in-memory pipeline.

    The configuration YAML is interpreted using
    :func:`src.common_utils.common_init.init_from_yaml`, which mirrors the
    behaviour of both the geometry and matrix-stage CLIs.  This function
    intentionally prints only a short summary of the shapes of the
    constructed matrices; it does not write any side-effect files.
    """


    print_header("SymDCT:\n  Jump Network + Transport Matrix + DiffCoeff Solver & Evaluator")

    configure_cli_logging()
    yaml_file, params, output_dir = init_from_yaml()

    poscar_fname = params["fname"]
    crystal_matrix = params.get("matrix", None)

    # Geometry stage: build jump network and refined jump structure.
    # Here we reuse the same logic as the geometry CLI but we keep all
    # data in memory.
    # Auto-discover cutoffs.yaml next to the input yaml; falls back to the
    # global max_jump_length from params (cutoffs = {}) when absent.
    cutoffs, max_jump_length = load_cutoffs_config(yaml_file, params)

    print_header("Stage 1: Building jump network")
    G, jr, jump_sites = build_jump_network_from_poscar(
        poscar_fname,
        cutoffs,
        str(output_dir),
        crystal_matrix,
        max_jump_length,
    )

    # Build a symmetry context compatible with the geometry stage and
    # pass both graph and context into the matrix builder.  The current
    # implementation of ``build_matrices`` still reconstructs its own
    # internal JumpNetwork from ``params``, but the call signature makes
    # the intended "from memory" usage explicit.
    symmetry_ctx = JumpNetworkSymmetryContext(poscar_fname, crystal_matrix)
    print(f"Simulation cell is {'conventional' if symmetry_ctx.is_conventional_cell() else 'primitive'}")
    lattice_type = symmetry_ctx.get_lattice_type()
    print(f"Crystal lattice type: {lattice_type}")
    if lattice_type == "cubic":
        print("  Note: For cubic lattices, the diffusion coefficient is isotropic and the solver will return a single value.")
    elif lattice_type == "hexagonal":
        print("  Note: For hexagonal lattices, the diffusion coefficient is anisotropic.\n"
              "        The solver will return two values: in-plane and out-of-plane.")
    elif lattice_type == "tetragonal":
        print("  Note: For tetragonal lattices, the diffusion coefficient is anisotropic.\n"
              "        The solver will return three values.")
    else:
        print("  Note: For other lattice types, the diffusion coefficient is anisotropic.\n"
              "        The solver will return three values.")

    print_header("Stage 2: Constructing transport matrix")
    result = build_matrices(
        graph=G,
        symmetry_ctx=symmetry_ctx,
        params=params,
    )


    print(f"\nFull transition matrix (from memory): {result.full_tp_matrix.shape}")
    if verbose:
        print_axis_build_summary(result)
    for axis, axis_result in result.axis_results.items():
        shape = axis_result.reordered_tp_matrix.shape
        # remap labels before printing matrix tables
        symbol_alias_map = symmetry_ctx.get_symbol_alias_map()
        remapped_labels = []
        for lbl in axis_result.labels:
            label = str(lbl)
            symbol = label.rstrip("0123456789")
            suffix = label[len(symbol):]
            remapped_labels.append(f"{symbol_alias_map.get(symbol, symbol)}{suffix}")
        axis_label = format_axis_label(axis, axis_result.hkl_structure)
        print_transition_matrix_table(axis_result.reordered_tp_matrix, labels=remapped_labels,
                                      title=f"Transition Matrix ({shape[0]}x{shape[1]}) for axis {axis_label}")

    print_header("Stage 3: Solving for diffusion coefficients (symbolic)")
    # coeffs is keyed by hkl_conventional -- see MatrixBuildResult.axis_results
    # and solve_all_axes -- so it can be handed straight to
    # write_diffusion_coefficients below with no separate conventional lookup.
    coeffs = solve_from_matrix_result(result)
    print("Solved symbolic diffusion coefficients for axes:")
    for axis, expr in sorted(coeffs.items()):
        hkl_structure = result.axis_results[axis].hkl_structure
        axis_label = format_axis_label(axis, hkl_structure)
        print(f"  {axis_label}: {expr}")

    if export_dc:
        from src.transport_coefficient_solver.transport_coefficient_solver import write_diffusion_coefficients
        write_diffusion_coefficients(coeffs, output_dir)

    print_header("Stage 4: Evaluating diffusion coefficients (numeric)")
    yaml_path = Path(yaml_file)
    # energetics.yaml is hand-authored input data, so -- like the config YAML
    # and the POSCAR it references -- it lives in input/ next to yaml_path,
    # not in output/ (which holds only Stage 1-3 generated artifacts,
    # including the auto-generated energetics_template.yaml below).
    production_path = yaml_path.parent / "energetics.yaml"
    if not production_path.exists():
        logger.warning(
            "Production energetics file not found: %s; skipping Stage 4 (numeric evaluation).",
            production_path,
        )
        return 0

    template_path = output_dir / "energetics_template.yaml"
    if not template_path.exists():
        raise FileNotFoundError(f"Geometry template not found: {template_path}")

    dc_files = sorted(output_dir.glob("dc_symbolic_*.txt"))
    if not dc_files:
        raise FileNotFoundError(f"No symbolic diffusion coefficient files found in {output_dir}")

    ctx = load_evaluation_context(production_path, template_path)
    plot_files = []
    for dc_file in dc_files:
        plot_file = evaluate_dc_file(dc_file, ctx)
        plot_files.append(plot_file)

    print_header(f"Done - results are written to {output_dir}"
                 "\n  Diffusion coefficients:\n    "
                 f"{"\n    ".join([str(file) for file in dc_files])}"
                 "\n  Arrhenius plots:\n    "
                 f"{"\n    ".join([str(file) for file in plot_files])}"
                 )

    return 0
if __name__ == "__main__":
    raise SystemExit(main())
