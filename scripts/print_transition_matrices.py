#!/usr/bin/env python3
"""Pretty-print transition matrices from the Stage-2 payload files.

Usage:
    python -m src.jump_network_builder examples/hcp_vac/input/hcp_vac.yaml
    python -m src.matrix_builder examples/hcp_vac/input/hcp_vac.yaml
    python scripts/print_transition_matrices.py examples/hcp_vac/input/hcp_vac.yaml

The script reads the files tc_axis_payload_<h>_<k>_<l>.pickle written by Stage 2
and displays the matrices as formatted tables with symbolic expressions.
"""

import sys
import pickle
from pathlib import Path
from typing import Sequence, Union

import sympy as sp


def print_transition_matrix_table(
    matrix: sp.Matrix,
    labels: Sequence[str],
    title: str = "Transition Matrix",
) -> None:
    """Pretty-print a transition matrix as an aligned table."""
    
    nrows, ncols = matrix.shape
    
    # Calculate field widths for alignment
    field_width = {}
    for col in range(ncols):
        max_len = len(str(labels[col]))
        for row in range(nrows):
            expr_str = str(sp.simplify(matrix[row, col]))
            max_len = max(max_len, len(expr_str))
        field_width[col] = max_len + 2
    
    row_label_width = max(len(str(lbl)) for lbl in labels) + 2
    
    # Print header
    print(f"\n{title}")
    print("─" * (row_label_width + sum(field_width.values()) + 5))
    
    # Column headers
    header = " " * row_label_width + "|"
    for col in range(ncols):
        header += f" {str(labels[col]):>{field_width[col]-2}} |"
    print(header)
    print("─" * (row_label_width + sum(field_width.values()) + 5))
    
    # Data rows
    for row in range(nrows):
        row_str = f" {str(labels[row]):<{row_label_width-2}} |"
        for col in range(ncols):
            expr_str = str(sp.simplify(matrix[row, col]))
            row_str += f" {expr_str:>{field_width[col]-2}} |"
        print(row_str)
    
    print("─" * (row_label_width + sum(field_width.values()) + 5))
    print()


def print_matrices_from_output_dir(output_dir: Union[Path, str]) -> None:
    """Load and display all transition matrices from pickle files."""
    
    output_dir = Path(output_dir)
    
    # Find the results filename prefix from YAML
    yaml_files = list(output_dir.parent.glob("input/*.yaml"))
    if not yaml_files:
        print("Error: Could not find YAML configuration file")
        return
    
    yaml_file = yaml_files[0]
    
    # Stage 2 writes one payload file per axis
    pickle_files = sorted(output_dir.glob("tc_axis_payload_*.pickle"))
    if not pickle_files:
        print(f"Error: No tc_axis_payload_*.pickle files found in {output_dir}")
        print("Run Stages 1 and 2 first:")
        print("  python -m src.jump_network_builder <yaml_config>")
        print("  python -m src.matrix_builder <yaml_config>")
        return

    print("\n" + "=" * 100)
    print(f"TRANSITION MATRICES - {yaml_file.stem.upper()}")
    print("=" * 100)
    print(f"Configuration: {yaml_file}")
    print(f"Output directory: {output_dir}")

    results = {}
    for pkl_file in pickle_files:
        with open(pkl_file, 'rb') as f:
            data = pickle.load(f)
        axis = data.get("axis_hkl_conventional") if isinstance(data, dict) else None
        axis_name = str(tuple(axis)) if axis is not None else pkl_file.stem.replace("tc_axis_payload_", "")
        results[axis_name] = data

    # Display each matrix
    for axis_name in sorted(results.keys()):
        data = results[axis_name]
        print(f"\n{'─' * 100}")
        print(f"Diffusion Axis: {axis_name}")
        print(f"{'─' * 100}")
        
        if isinstance(data, dict) and 'transition_matrix' in data:
            matrix = data['transition_matrix']
            labels = data.get('labels', [f'S{i}' for i in range(matrix.shape[0])])
            jump_length = data.get('min_jump_length', 1.0)
            
            print(f"Jump length factor: {jump_length}")
            print(f"Sites: {labels}")
            
            print_transition_matrix_table(matrix, labels, f"Transition Matrix for axis {axis_name}")
        else:
            print(f"Unexpected data format: {type(data)}")
    
    print("=" * 100)
    print("Legend:")
    print("  W_AA = Jump probability for A->A jumps (W0)")
    print("  W_AB = Jump probability for A->B jumps (W1, W2, ...)")
    print("  cos(t), exp(I*t), etc. = Wave vector components")
    print("=" * 100 + "\n")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python scripts/print_transition_matrices.py <yaml_config_or_output_dir>")
        sys.exit(1)
    
    arg = Path(sys.argv[1])
    
    # If argument is a YAML file, use its output directory
    if arg.suffix == '.yaml':
        output_dir = arg.parent.parent / 'output'
    else:
        output_dir = arg
    
    print_matrices_from_output_dir(output_dir)
