from __future__ import annotations

import json
import logging
import os
from typing import List, Any
import networkx as nx
import pymatgen as pmg
import sympy as sp
from src.common_utils.jump_vector_codec import encode_jump_vector
from src.common_utils.jump_register import JumpRegister
from src.common_utils.utils import to_rational, configure_cli_logging
from src.common_utils.yaml_with_comments import dump_data
from src.common_utils.label_utils import sort_pair_key, remap_pair_keys_in_template
from .jump_orbits import classify_jump_orbits


configure_cli_logging()
logger = logging.getLogger(__name__)

# Key under which nx.node_link_data/node_link_graph store edge records.
# Kept explicit (rather than relying on the library default) so Stage 1's
# writer and Stage 2's reader agree regardless of the installed NetworkX
# version. Mirrors the convention used by src.jump_network.core.LINKS.
LINKS = "edges"


def _node_link_data(graph: nx.Graph) -> dict:
    """Return node-link data with a NetworkX-version-compatible edge key."""

    try:
        return nx.node_link_data(graph, edges=LINKS)
    except TypeError:
        return nx.node_link_data(graph)


def _node_link_graph(data: dict) -> nx.Graph:
    """Load node-link data with a NetworkX-version-compatible edge key."""

    try:
        return nx.node_link_graph(data, edges=LINKS)
    except TypeError:
        return nx.node_link_graph(data)

def pretty_s_matrix(ss_matrix, tol: float = 1e-6) -> sp.Matrix:
    """Return and print S-matrix (conventional -> primitive)."""
    print("S-matrix (used for conversion from primitive to conventional unit cell):")
    s_matrix = sp.Matrix([sp.nsimplify(elem, tolerance=tol) for vec in ss_matrix for elem in vec]).reshape(3, 3)
    sp.pprint(s_matrix)
    return s_matrix


def get_all_neighbors(jump_sites: pmg.core.Structure, max_jump_length: float):
    _all_nbs = jump_sites.get_all_neighbors(max_jump_length)
    if any(_all_nbs) == 0:
        import sys
        sys.exit(f"ERROR: No neighbors found within the cutoff distance {max_jump_length}. Try to increase the cutoff.")
    else:
        return _all_nbs





def _select_neighbors_by_cutoffs(jump_sites: pmg.core.Structure, all_neighbors: list, cutoffs: dict,
                                 symbol_alias_map: dict, default_unit='angstrom'):
    """Filter neighbor entries using cutoffs mapping. Returns filtered neighbor list.

    cutoffs is expected to be a dict with keys like "A-B" and float values.
    """
    # normalize cutoffs keys
    normalized = {"-".join(sorted(k.split("-"))): v for k, v in cutoffs.items()}
    filtered = []
    no_cutoff_warn = {}
    for i, neighbors in enumerate(all_neighbors):
        kept = []
        for nb in neighbors:
            type_i = symbol_alias_map.get(jump_sites.sites[i].specie.symbol, jump_sites.sites[i].specie.symbol)
            type_j = symbol_alias_map.get(nb.specie.symbol, nb.specie.symbol)
            key = sort_pair_key(type_i, type_j)
            maxd = normalized.get(key, None)
            if maxd is None:
                no_cutoff_warn[key] = no_cutoff_warn[key] + 1 if key in no_cutoff_warn else 1
                kept.append(nb)
                continue
            # if no cutoff provided for this pair keep the neighbor (user omitted pair)
            if nb.nn_distance <= maxd:
                kept.append(nb)
        filtered.append(kept)

    if no_cutoff_warn:
        for key, count in no_cutoff_warn.items():
            logger.warning(f"{count} neighbors found for bond {key[0]}-{key[2]} without a specified cutoff.")

    return filtered


def _report_orbit_diagnostics(orbits, triplet_ids, records, node_lbls) -> None:
    """Log where the orbit partition differs from the triplet partition."""

    if orbits.n_skipped_ops:
        logger.warning(
            "%d of %d symmetry operations do not map the jump sites onto themselves "
            "within the site tolerance and were not used.", orbits.n_skipped_ops, orbits.n_ops)
    if orbits.n_unmapped_edges:
        logger.warning(
            "The jump network is not invariant under the space group: %d images of edges "
            "are missing. Check the cutoffs and tweak_transition_rates; jumps that the "
            "symmetry relates may be classified separately.", orbits.n_unmapped_edges)

    members: dict[int, set] = {}
    for tid, oid in zip(triplet_ids, orbits.ids):
        members.setdefault(tid, set()).add(oid)
    for tid, oids in members.items():
        if len(oids) > 1:
            i, j, _h, nb, _ = next(rec for rec in records if rec[4] == tid)
            logger.warning(
                "Jumps %s->%s of length %.4f A form one class by the triplet rule but %d "
                "inequivalent orbits by symmetry. The orbit classification separates them; "
                "use jump_classification: triplet to reproduce version 0.1.0.",
                node_lbls[i], node_lbls[j], float(nb.nn_distance), len(oids))
    orbit_members: dict[int, set] = {}
    for tid, oid in zip(triplet_ids, orbits.ids):
        orbit_members.setdefault(oid, set()).add(tid)
    for oid, tids in orbit_members.items():
        if len(tids) > 1:
            logger.warning(
                "Symmetry-equivalent jumps fall into %d triplet classes (orbit %d). "
                "The orbit classification merges them.", len(tids), oid)


def build_jump_network_from_poscar(
        poscar_file: str, cutoffs_config: dict, output_dir: str,
        crystal_matrix: List[str] | None = None, max_jump_length: float = 2.5, tol: float = 1e-6,
        jump_classification: str = "orbit",
        structure_format: str | None = None,
):
    """Build jump network from a structure file (POSCAR, CIF, ...) and cutoffs and dump graph/template outputs.

    ``jump_classification`` selects how directed jumps are grouped into jump
    classes, i.e. into rate symbols.  ``"orbit"`` (default) groups jumps that
    are related by a space-group operation.  ``"triplet"`` groups jumps with
    equal start label, end label and length (the classification of version
    0.1.0).  Both partitions are computed, and a warning is logged when they
    differ.

    Returns: graph (MultiDiGraph), JumpRegister instance, pymatgen refined jump structure
    """
    if jump_classification not in ("orbit", "triplet"):
        raise ValueError(
            f"jump_classification must be 'orbit' or 'triplet', got {jump_classification!r}"
        )
    from src.symmetry_context.jump_symmetry_context import JumpNetworkSymmetryContext
    from .energetics_template import prepare_energetics_template

    # prepare output dir
    os.makedirs(output_dir, exist_ok=True)

    # build refined jump structure
    ctx = JumpNetworkSymmetryContext(
        poscar_file,
        crystal_matrix,
        symprec=0.1,  # or tie to tol if desired
        structure_format=structure_format,
    )

    jump_sites = ctx.get_refined_jump_structure()
    symbol_alias_map = ctx.get_symbol_alias_map()
    node_lbls = [symbol_alias_map.get(lbl,lbl) for lbl in ctx.get_node_lbls()]

    # neighbors within global cutoff
    all_nbs = get_all_neighbors(jump_sites, max_jump_length)
    # apply per-pair cutoffs from config
    cutoffs = cutoffs_config.get('bonds', {}) if isinstance(cutoffs_config, dict) else {}
    all_nbs_filtered = _select_neighbors_by_cutoffs(jump_sites, all_nbs, cutoffs, symbol_alias_map)

    # create JumpNetwork-like graph and register jumps
    G = nx.MultiDiGraph()
    G.add_nodes_from(range(jump_sites.num_sites))
    jr = JumpRegister()

    records = []
    for inb, site in enumerate(jump_sites.sites):
        # store site labels on nodes (used later by energetics template)
        G.nodes[inb]["lbl"] = node_lbls[inb]
        for nb in all_nbs_filtered[inb]:
            jnb = nb.properties.get('site_idx', getattr(nb, 'index', None))
            h = [
                to_rational(fc1, tolerance=tol) - to_rational(fc2, tolerance=tol)
                for fc1, fc2 in zip(nb.frac_coords, site.frac_coords)
            ]
            # classification by the triplet (start label, end label, length)
            triplet_id = jr.register(
                node_lbls[inb],
                nb.specie.symbol,
                float(nb.nn_distance),
                jump_vector=list(h),
                image=list(map(int, nb.image)),
            )
            jr.add_to_netwk(triplet_id)
            records.append((inb, jnb, h, nb, triplet_id))

    jump_ids = [rec[4] for rec in records]
    reverse_ids = None
    n_classes = jr.get_number_of_jumps()
    if records:
        orbits = classify_jump_orbits(
            jump_sites.lattice.matrix,
            jump_sites.frac_coords,
            node_lbls,
            [(i, j, [float(x) for x in h]) for i, j, h, _nb, _tid in records],
            symprec=ctx.symprec,
        )
        _report_orbit_diagnostics(orbits, jump_ids, records, node_lbls)
        if jump_classification == "orbit":
            jump_ids = orbits.ids
            reverse_ids = orbits.reverse_ids
            n_classes = len(set(jump_ids))
    logger.info("Jump classification '%s': %d jump classes", jump_classification, n_classes)

    for (inb, jnb, h, nb, _tid), jump_count in zip(records, jump_ids):
        G.add_edge(
            inb,
            jnb,
            h=encode_jump_vector(h),
            nn_dist=float(nb.nn_distance),
            image=list(map(int, nb.image)),
            glob_jump_cnt=int(jump_count),
        )
    if jump_classification == "orbit":
        # the register only knows triplets, so its class count is replaced
        jr.set_number_of_jumps(n_classes)

    # dump node-link JSON (canonical Stage-1 -> Stage-2 graph exchange format)
    node_link_path = os.path.join(output_dir, "jump_structure.node_link.json")
    node_link_data = _node_link_data(G)
    with open(node_link_path, "w", encoding="utf-8") as fh:
        json.dump(node_link_data, fh, indent=2)

    # retain multiline adjlist as a secondary inspection artifact
    adjlist_path = os.path.join(output_dir, "jump_structure.adjlist")
    nx.write_multiline_adjlist(G, adjlist_path)

    # save refined jump-sites structure used by Stage 2 for symmetry context
    jump_poscar_path = os.path.join(output_dir, "jump_sites.poscar")
    # only the file name goes into the header, so that no local path is written
    pmg.io.vasp.Poscar(
        jump_sites, f"Refined jump sites from {os.path.basename(str(poscar_file))}"
    ).write_file(jump_poscar_path)

    # build and dump energetics template YAML
    energetics_template = prepare_energetics_template(G, reverse_ids=reverse_ids)
    energetics_path = os.path.join(output_dir, "energetics_template.yaml")
    energetics_template["lattice constants"] = {axis: round(val,7) for axis, val in zip("abc", jump_sites.lattice.abc)}
    energetics_template["symbol_alias_map"] = ctx.get_symbol_alias_map()
    remap_pair_keys_in_template(energetics_template)

    dump_data(energetics_path, energetics_template)

    return G, jr, jump_sites



