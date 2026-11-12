"""Core graph container for jump networks."""

from __future__ import annotations

import json
from typing import List, Tuple

import networkx as nx


LINKS = "edges"

class JumpNetworkGraph:
    """Own the directed multigraph and edge/node metadata for jumps."""

    def __init__(self, num_sites: int, node_labels: List[str]):
        self.num_sites = num_sites
        self.node_labels = list(node_labels)
        self.graph = nx.MultiDiGraph()
        self.graph.add_nodes_from(range(num_sites))
        for i, label in enumerate(self.node_labels):
            self.graph.nodes[i]["lbl"] = label

    def add_edge(
        self,
        start: int,
        end: int,
        jump_vector: Tuple | List,
        distance: float,
        image: Tuple | List,
        global_jump_index: int,
    ) -> None:
        self.graph.add_edge(
            start,
            end,
            h=list(jump_vector),
            nn_dist=float(distance),
            image=list(image),
            glob_jump_cnt=global_jump_index,
        )

    def get_out_edges(self, node: int):
        for start, end, key, h in self.graph.out_edges(node, data="h", keys=True):
            yield start, end, key, h

    def get_in_edges(self, node: int):
        for start, end, key, h in self.graph.in_edges(node, data="h", keys=True):
            yield start, end, key, h

    def get_edge_data(self, start: int, end: int, key: int) -> dict:
        return dict(self.graph.edges[start, end, key])

    def dump_to_adjlist(self, filename: str) -> None:
        nx.write_multiline_adjlist(self.graph, filename)

    @staticmethod
    def load_from_adjlist(filename: str) -> nx.MultiDiGraph:
        return nx.read_multiline_adjlist(filename, create_using=nx.MultiDiGraph)

    def dump_to_node_link_json(self, filename: str) -> None:
        print(repr(LINKS))
        with open(filename, "w", encoding="utf-8") as fh:
            json.dump(nx.node_link_data(self.graph, edges=LINKS), fh, indent=2)

    @staticmethod
    def load_from_node_link_json(filename: str) -> nx.MultiDiGraph:
        print(repr(LINKS))
        with open(filename, "r", encoding="utf-8") as fh:
            data = json.load(fh)
        return nx.node_link_graph(data, edges=LINKS)
