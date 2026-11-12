import logging
import networkx as nx

from src.common_utils.yaml_with_comments import CommentedData


log = logging.getLogger(__name__)


def prepare_energetics_template(graph: nx.Graph, length_tolerance: float = 0.01) -> dict:
    """Prepare energetics template YAML structure from a jump graph.

    The template has the following (YAML) structure:

    settings:
      length_tolerance: <float>
      units:
        energy: eV
        length: angstrom
        frequency: THz

    site frequencies:
      default_frequency: 10.0
      <SITE_LABEL>: null         # to be filled by user

    sites:
      <SITE_LABEL>: null         # configurational energies to be filled by user

    jumps:
      "A-B":
        - jump_id: [i_forward, i_backward]
          length: <float>
          energy_sp: null
        - ...                    # other distinct lengths for the same pair

    Notes
    -----
    * Jump labels are normalised ("A-B" == "B-A").
    * Each undirected jump (label + representative length) appears exactly once.
    * Forward/backward jumps have different global indices; both are stored as
      a pair [i_forward, i_backward]. Symmetric jumps like "A-A" use the same
      index twice: [i, i].
    * Distinct jumps between the same pair (short vs long A-A, etc.) are
      distinguished by their representative length and stored as separate
      list entries.
    """
    from src.common_utils.label_utils import sort_pair_key

    # HELPER FUNCTIONS
    def fill_per_site_properties(site_properties: list[dict]):
        # Collect unique site labels and initialise placeholders for energies
        seen_labels: list[str] = []
        for node in graph.nodes:
            lbl = graph.nodes[node]["lbl"]
            if lbl not in seen_labels:
                seen_labels.append(lbl)
                for a_property in site_properties:
                    a_property[lbl] = CommentedData(None, commented=True)  # user will provide configurational energy

    def group_jumps(jn_graph: nx.Graph):
        # Intermediate structure for grouping jumps
        # pair_key -> list of {"length": float, "forward": int | None, "backward": int | None}
        grouped: dict[str, list[dict]] = {}

        for start_node, end_node, key, jump_length in jn_graph.edges(data="nn_dist", keys=True):
            start_lbl = jn_graph.nodes[start_node]["lbl"]
            end_lbl = jn_graph.nodes[end_node]["lbl"]

            # Normalised label, order independent ("A-B" == "B-A")
            pair_key = sort_pair_key(start_lbl, end_lbl)

            edge_data = jn_graph.get_edge_data(start_node, end_node, key=key) or {}
            global_idx = edge_data.get("glob_jump_cnt")
            if global_idx is None:
                logging.warning("Edge (%s -> %s, key=%s) is missing 'glob_jump_cnt'", start_lbl, end_lbl, key)
                continue

            # Determine canonical direction w.r.t. sorted labels
            a, b = sorted([start_lbl, end_lbl])
            if start_lbl == end_lbl:
                direction = "symmetric"
            elif (start_lbl, end_lbl) == (a, b):
                direction = "forward"
            else:
                direction = "backward"

            # Find or create a length-group for this pair_key
            groups = grouped.setdefault(pair_key, [])
            group = None
            for g in groups:
                if abs(float(jump_length) - float(g["length"])) < length_tolerance:
                    group = g
                    break
            if group is None:
                group = {"length": float(jump_length), "forward": None, "backward": None}
                groups.append(group)

            # Assign indices according to direction
            if direction == "symmetric":
                # A-A jumps: use the same index for both directions
                if group["forward"] is None:
                    group["forward"] = global_idx
                if group["backward"] is None:
                    group["backward"] = global_idx
            elif direction == "forward":
                if group["forward"] is None:
                    group["forward"] = global_idx
                elif group["forward"] != global_idx:
                    logging.debug("Multiple forward indices for %s (length %.5f): %s, %s",
                                  pair_key, group["length"], group["forward"], global_idx)
            else:  # "backward"
                if group["backward"] is None:
                    group["backward"] = global_idx
                elif group["backward"] != global_idx:
                    logging.debug("Multiple backward indices for %s (length %.5f): %s, %s",
                                  pair_key, group["length"], group["backward"], global_idx)
        return grouped

    def get_yaml_ready(grp_jumps):
        # Convert grp_jumps structure into YAML-ready "jumps" mapping
        jumps: dict[str, list[dict]] = {}
        for pair_key, groups in grp_jumps.items():
            entries: list[dict] = []
            for g in groups:
                f_idx = g["forward"]
                b_idx = g["backward"]

                # Be robust if only one direction was ever seen
                if f_idx is None and b_idx is None:
                    # Should not normally happen, but guard against it
                    logging.warning("No indices recorded for jump %s (length %.5f)", pair_key, g["length"])
                    continue
                if f_idx is None:
                    f_idx = b_idx
                if b_idx is None:
                    b_idx = f_idx

                num_digits = 6
                entries.append({
                    "jump_id": {'forward': int(f_idx), 'backward': int(b_idx)},
                    "length": round(g["length"], num_digits),
                    "energy_sp": CommentedData(None, commented=True),
                })

            # Only add key if we actually have entries
            if entries:
                jumps[pair_key] = entries
        return jumps

    #--------------------------------------#
    # BODY of prepare_energetics_template  #
    #--------------------------------------#

    # Per-site data
    site_energies: dict[str, float | None] = {}
    frequencies: dict[str, float | None] = {"default_frequency": 10.0}
    fill_per_site_properties([site_energies, frequencies])

    # create intermediate grouped jumps structure
    grouped_jumps = group_jumps(graph)

    # convert to YAML-ready structure
    jumps_yaml = get_yaml_ready(grouped_jumps)

    template = {
        "settings": {
            "length_tolerance": float(length_tolerance),
            "units": {
                "energy": "eV",
                "length": "angstrom",
                "frequency": "THz",
            },
        },
        "site frequencies": frequencies,
        "site energies": site_energies,
        "jumps": jumps_yaml,
    }
    return template
