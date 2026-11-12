"""Symbolic Stage-3 transport-coefficient solver."""

from __future__ import annotations

from typing import Mapping, Sequence

import sympy as sp

from .contracts import AxisDiffusionResult, AxisMatrixPayload


def infer_t_symbol(matrix: sp.Matrix, preferred_name: str = "t") -> sp.Symbol:
    """Infer the wave-vector symbol used in transition matrix expressions."""
    free_symbols = sorted(matrix.free_symbols, key=lambda s: str(s))
    preferred = [s for s in free_symbols if str(s) == preferred_name]
    if preferred:
        return preferred[0]
    if len(free_symbols) == 1:
        return free_symbols[0]
    raise ValueError(
        f"Cannot infer t symbol from matrix free symbols: {[str(s) for s in free_symbols]}"
    )


def apply_jump_label_substitutions(
    dc_expr: sp.Expr, 
    jump_label_map: Mapping[int, str] | None,
    jump_label_display_map: Mapping[str, str] | None = None,
    w_substitutions: Mapping[sp.Symbol, sp.Expr] | None = None,
) -> sp.Expr:
    """Apply human-readable jump label substitutions after DC is computed.
    
    Replaces W_i symbols with W_{label} notation (e.g., W1 -> W_AB).
    This is done POST-computation to avoid slowing down SymPy during solving.
    """
    if not jump_label_map or not dc_expr.free_symbols:
        return dc_expr

    # Find which final W symbols appear in dc_expr
    w_symbols_in_expr = {s for s in dc_expr.free_symbols if str(s).startswith('W')}
    
    if not w_symbols_in_expr:
        return dc_expr

    # Build mapping from W symbol string to (global_idx, label) pairs
    w_str_to_labels = {}  # "W2" -> [(0, 'A->A'), (2, 'B->B')]
    
    # Convert w_substitutions dict to string keys for reliable lookup
    w_subs_by_str = {}
    if w_substitutions:
        w_subs_by_str = {str(k): str(v) for k, v in w_substitutions.items()}
    
    for global_idx, label in jump_label_map.items():
        # Start with the original symbol W_i (where i=global_idx)
        current_w_str = f"W{global_idx}"
        
        # Follow the substitution chain using string-based lookup
        if w_subs_by_str:
            for _ in range(len(jump_label_map) + 1):  # Safety limit
                if current_w_str in w_subs_by_str:
                    current_w_str = w_subs_by_str[current_w_str]
                else:
                    break
        
        # current_w_str is now the final W symbol string
        if current_w_str not in w_str_to_labels:
            w_str_to_labels[current_w_str] = []
        w_str_to_labels[current_w_str].append((global_idx, label))

    # Now build substitutions for W symbols that appear in the expression
    substitutions = {}
    for w_symbol in w_symbols_in_expr:
        w_str = str(w_symbol)
        if w_str in w_str_to_labels:
            labels_list = w_str_to_labels[w_str]
            
            # Separate identity and non-identity jumps
            non_identity_labels = []
            identity_labels = []
            
            for global_idx, label in labels_list:
                if "->" in label:
                    a, b = label.split("->")
                    a, b = a.strip(), b.strip()
                    if a != b:
                        non_identity_labels.append(label)
                    else:
                        identity_labels.append(label)
            
            # Prefer non-identity labels (A->B over A->A)
            # Sort alphabetically within each category for deterministic output
            candidates = sorted(non_identity_labels) if non_identity_labels else sorted(identity_labels)
            selected_label = candidates[0] if candidates else None
            
            if selected_label:
                a, b = selected_label.split("->")
                a_disp = jump_label_display_map.get(a.strip(), a.strip()) if jump_label_display_map else a.strip()
                b_disp = jump_label_display_map.get(b.strip(), b.strip()) if jump_label_display_map else b.strip()
                new_symbol = sp.Symbol(f"W_{{{a_disp}{b_disp}}}")
                substitutions[w_symbol] = new_symbol

    if substitutions:
        dc_expr = dc_expr.subs(substitutions)
    return dc_expr


def get_denom(matrix: sp.Matrix, t_symbol: sp.Symbol) -> sp.Expr:
    """Return denominator used in symbolic diffusion-coefficient evaluation."""
    irow = 0
    augmented = matrix.row_insert(irow, sp.ones(1, matrix.shape[1]))
    augmented.row_del(irow + 1)
    return sp.factor(sp.det(augmented.subs(t_symbol, 0)))


def get_num(matrix: sp.Matrix, t_symbol: sp.Symbol) -> sp.Expr:
    """Return numerator used in symbolic diffusion-coefficient evaluation."""
    return sp.factor(sp.diff(sp.det(matrix), t_symbol, 2).subs(t_symbol, 0))


def compute_axis_diffusion_coefficient(
    payload: AxisMatrixPayload,
    *,
    min_jump_length: sp.Expr | int | float | None = None,
    apply_w_substitutions: bool = True,
) -> AxisDiffusionResult:
    """Compute symbolic diffusion coefficient for one axis payload."""
    t_symbol = payload.t_symbol or infer_t_symbol(payload.transition_matrix)
    interplanar_distance_squared = payload.interplanar_distance_squared
    if interplanar_distance_squared is None:
        raise ValueError("interplanar_distance_squared is required for solver input")
    interplanar_distance_squared = sp.nsimplify(interplanar_distance_squared)
    effective_min_jump_length = (
        payload.min_jump_length if min_jump_length is None else min_jump_length
    )
    if effective_min_jump_length is None:
        effective_min_jump_length = 1

    denom = get_denom(payload.transition_matrix, t_symbol)
    if sp.simplify(denom) == 0:
        raise ValueError("Denominator is zero; transition matrix is singular for t=0")

    num = get_num(payload.transition_matrix, t_symbol)
    dc_expr = -sp.Rational(1, 2) * interplanar_distance_squared * (
        sp.nsimplify(effective_min_jump_length) ** 2
    ) * sp.nsimplify(num / denom)
    
    # Simplify first (without expensive label substitutions)
    dc_expr = sp.factor(sp.nsimplify(dc_expr))
    
    # Apply stored w_substitutions if provided (these are special values, not labels)
    if apply_w_substitutions and payload.w_substitutions:
        dc_expr = dc_expr.subs(dict(payload.w_substitutions))
    
    # Apply human-readable label substitutions AFTER solving (post-processing only)
    if apply_w_substitutions and payload.jump_label_map:
        dc_expr = apply_jump_label_substitutions(
            dc_expr,
            payload.jump_label_map,
            payload.jump_label_display_map,
            payload.w_substitutions,
        )

    return AxisDiffusionResult(
        axis_hkl_structure=payload.axis_hkl_structure,
        axis_hkl_conventional=payload.axis_hkl_conventional,
        diffusion_coefficient=dc_expr,
    )


def solve_all_axes(
    payloads: Sequence[AxisMatrixPayload] | Mapping[tuple[int, int, int], AxisMatrixPayload],
    *,
    min_jump_length: sp.Expr | int | float | None = None,
    apply_w_substitutions: bool = True,
) -> dict[tuple[int, int, int], sp.Expr]:
    """Solve symbolic diffusion coefficients for all provided axes.

    The returned dict is keyed by ``axis_hkl_conventional`` (see
    ``AxisMatrixPayload``), the canonical, always-available axis identifier,
    so callers can hand it straight to ``write_diffusion_coefficients``
    without any separate conventional-hkl lookup.
    """
    if isinstance(payloads, Mapping):
        iterable = payloads.values()
    else:
        iterable = payloads

    results: dict[tuple[int, int, int], sp.Expr] = {}
    for payload in iterable:
        axis_result = compute_axis_diffusion_coefficient(
            payload,
            min_jump_length=min_jump_length,
            apply_w_substitutions=apply_w_substitutions,
        )
        results[axis_result.axis_hkl_conventional] = axis_result.diffusion_coefficient
    return results
