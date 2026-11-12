#!/usr/bin/env python3
"""Read a NetworkX multiline adjlist produced by jump_network_builder and plot it.

Usage:
    python scripts/plot_jump_adjlist.py path/to/jump_structure.adjlist [--out out.png]

The script extracts `nn_dist` (jump length in angstrom) from edge attributes and labels
edges with the distance. If `--out` is provided the plot is saved to the path instead
of being shown interactively (useful for headless environments / CI).
"""

from __future__ import annotations

import sys
import ast
import re
from pathlib import Path
import argparse

import matplotlib
# use Agg for headless environments so we can save to file without X11
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx

FLOAT_RE = re.compile(r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?\d+)?")


def extract_float_from_repr(s):
    if s is None:
        return None
    if isinstance(s, (float, int)):
        return float(s)
    if isinstance(s, (list, tuple)):
        for el in s:
            try:
                return float(el)
            except Exception:
                continue
        return None
    try:
        val = ast.literal_eval(s)
        if isinstance(val, (float, int)):
            return float(val)
        if isinstance(val, (list, tuple)):
            for el in val:
                if isinstance(el, (float, int)):
                    return float(el)
    except Exception:
        pass
    m = FLOAT_RE.search(str(s))
    if m:
        try:
            return float(m.group(0))
        except Exception:
            return None
    return None


def parse_edge_attrs(raw_attrs):
    attrs = {}
    if isinstance(raw_attrs, str):
        try:
            parsed = ast.literal_eval(raw_attrs)
        except Exception:
            parsed = {}
            nn_m = re.search(r"nn_dist\':?\s*([^,}\]]+)", raw_attrs)
            if nn_m:
                parsed['nn_dist'] = nn_m.group(1)
            h_m = re.search(r"'h'\s*:\s*Matrix\(\[([^\]]+)\]\)", raw_attrs, re.S)
            if h_m:
                nums = FLOAT_RE.findall(h_m.group(1))
                parsed['h'] = [float(n) for n in nums[:3]]
            img_m = re.search(r"image\':?\s*\(([^\)]+)\)", raw_attrs)
            if img_m:
                tup = img_m.group(1)
                nums = FLOAT_RE.findall(tup)
                parsed['image'] = [int(float(n)) for n in nums]
            g_m = re.search(r"'glob_jump_cnt'\s*:\s*([0-9]+)", raw_attrs)
            if g_m:
                parsed['glob_jump_cnt'] = int(g_m.group(1))
        parsed = parsed if isinstance(parsed, dict) else {}
        attrs['nn_dist'] = extract_float_from_repr(parsed.get('nn_dist'))
        if 'h' in parsed:
            attrs['h'] = parsed['h']
        else:
            attrs['h'] = None
        if 'image' in parsed:
            attrs['image'] = parsed['image']
        else:
            attrs['image'] = None
        attrs['glob_jump_cnt'] = int(parsed['glob_jump_cnt']) if 'glob_jump_cnt' in parsed else None

    elif isinstance(raw_attrs, dict):
        attrs['nn_dist'] = float(raw_attrs.get('nn_dist')) if raw_attrs.get('nn_dist') is not None else None
        hraw = raw_attrs.get('h')
        if hraw is None:
            attrs['h'] = None
        else:
            try:
                import sympy as sp
                if isinstance(hraw, sp.Matrix):
                    attrs['h'] = [float(x) for x in list(hraw)]
                else:
                    attrs['h'] = [float(x) for x in hraw]
            except Exception:
                try:
                    attrs['h'] = [float(x) for x in hraw]
                except Exception:
                    attrs['h'] = None
        img = raw_attrs.get('image')
        if img is None:
            attrs['image'] = None
        else:
            try:
                attrs['image'] = [int(x) for x in img]
            except Exception:
                try:
                    attrs['image'] = [int(float(x)) for x in img]
                except Exception:
                    attrs['image'] = None
        attrs['glob_jump_cnt'] = int(raw_attrs.get('glob_jump_cnt')) if raw_attrs.get('glob_jump_cnt') is not None else None
    else:
        attrs = {'nn_dist': None, 'h': None, 'image': None, 'glob_jump_cnt': None}

    return attrs


def build_plot_graph(adjlist_path: Path):
    G = nx.read_multiline_adjlist(str(adjlist_path), create_using=nx.MultiDiGraph)

    H = nx.DiGraph()
    # Use node label if available on the node (networkx stores node attributes if present)
    node_label_map = {}
    try:
        # networkx multiline adjlist stores nodes attributes in a separate section; read nodes attrs if present
        # fallback: use node ids as labels
        for n, data in G.nodes(data=True):
            lbl = data.get('lbl') if isinstance(data, dict) else None
            node_label_map[n] = str(lbl) if lbl else str(n)
    except Exception:
        node_label_map = {n: str(n) for n in G.nodes()}

    for u, v, data in G.edges(data=True):
        parsed = parse_edge_attrs(data)
        nn = parsed.get('nn_dist')
        jl = parsed.get('glob_jump_cnt')
        H.add_node(u)
        H.add_node(v)
        elab = f"{nn:.4f}" if nn is not None else (str(jl) if jl is not None else '')
        if H.has_edge(u, v):
            H[u][v]['labels'].append(elab)
        else:
            H.add_edge(u, v, labels=[elab])

    return H, node_label_map


def plot_graph(H: nx.Graph, node_label_map: dict, out_path: Path | None = None):
    pos = nx.spring_layout(H, seed=42)
    plt.figure(figsize=(10, 8))
    nx.draw_networkx_nodes(H, pos, node_size=600, node_color='lightblue')
    nx.draw_networkx_edges(H, pos, arrowstyle='->', arrowsize=12, connectionstyle='arc3,rad=0.1')
    nx.draw_networkx_labels(H, pos, labels={n: node_label_map.get(n, str(n)) for n in H.nodes()}, font_size=10)

    edge_labels = {}
    for u, v, data in H.edges(data=True):
        lbls = [s for s in data.get('labels', []) if s]
        if not lbls:
            continue
        uniq = sorted(set(lbls))
        edge_labels[(u, v)] = ', '.join(uniq)

    nx.draw_networkx_edge_labels(H, pos, edge_labels=edge_labels, label_pos=0.5, font_size=8, rotate=False)
    plt.title("Jump network")
    plt.axis('off')
    plt.tight_layout()
    if out_path:
        plt.savefig(str(out_path), dpi=200)
        print(f"Saved plot to: {out_path}")
    else:
        plt.show()


def main(argv=None):
    parser = argparse.ArgumentParser(description="Plot a jump_structure.adjlist produced by jump_network_builder")
    parser.add_argument('adjlist', help='Path to jump_structure.adjlist')
    parser.add_argument('--out', '-o', help='Optional output PNG file to save the plot')
    args = parser.parse_args(argv)

    p = Path(args.adjlist)
    if not p.exists():
        print(f"Adjlist not found: {p}")
        return 2

    H, node_label_map = build_plot_graph(p)
    out_path = Path(args.out) if args.out else None
    plot_graph(H, node_label_map, out_path)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

