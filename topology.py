"""Structural summary and hub-removal test for an exported SPOKE subgraph.

Reports node and edge counts, connected components, density, node-label
and relationship-type breakdowns, and the degree distribution. Then
removes the highest-degree nodes (or nodes named with --remove) and
reports how much of the graph stays connected.

Degree here is the number of distinct neighbors inside the exported
subgraph, not the node's degree in the full SPOKE graph.
"""

import argparse
import statistics
from collections import Counter

import networkx as nx

from graph_io import describe, load_graph, resolve_names


def summarize(G, rel_types, top_n):
    n = G.number_of_nodes()
    components = list(nx.connected_components(G))
    largest = max(components, key=len) if components else set()
    degrees = dict(G.degree())
    values = list(degrees.values())
    leaves = sum(1 for d in values if d == 1)

    print("Structure")
    print(f"  nodes:                {n}")
    print(f"  edges:                {G.number_of_edges()} "
          f"(from {sum(rel_types.values())} relationships)")
    print(f"  connected components: {len(components)}")
    print(f"  largest component:    {len(largest)} nodes ({100 * len(largest) / n:.1f}%)")
    print(f"  density:              {nx.density(G):.5f}")

    print("\nNode labels")
    for label, count in Counter(nx.get_node_attributes(G, "label").values()).most_common():
        print(f"  {label}: {count}")

    print("\nRelationship types")
    for rel_type, count in rel_types.most_common():
        print(f"  {rel_type}: {count}")

    print("\nDegree")
    print(f"  mean:   {statistics.mean(values):.2f}")
    print(f"  median: {statistics.median(values)}")
    print(f"  leaves (degree 1): {leaves} ({100 * leaves / n:.1f}%)")
    print(f"  top {top_n}:")
    for key, deg in sorted(degrees.items(), key=lambda kv: -kv[1])[:top_n]:
        print(f"    {describe(G, key)}: {deg}")


def removal_test(G, keys):
    """Remove `keys` and report the connectivity of what remains."""
    H = G.copy()
    H.remove_nodes_from(keys)
    remaining = H.number_of_nodes()

    print(f"\nRemoval test: {', '.join(describe(G, k) for k in keys)}")
    if remaining == 0:
        print("  no nodes remain")
        return

    components = list(nx.connected_components(H))
    largest = max(components, key=len)
    isolated = sum(1 for c in components if len(c) == 1)
    print(f"  remaining nodes:      {remaining}")
    print(f"  connected components: {len(components)}")
    print(f"  largest component:    {len(largest)} nodes ({100 * len(largest) / remaining:.1f}%)")
    print(f"  isolated nodes:       {isolated}")


def plot(G, path, hub_threshold):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    degrees = dict(G.degree())
    hubs = [k for k in G if degrees[k] > hub_threshold]
    rest = [k for k in G if degrees[k] <= hub_threshold]
    pos = nx.spring_layout(G, seed=42, k=0.15, iterations=50)

    fig, ax = plt.subplots(figsize=(10, 10))
    nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.2, width=0.4)
    nx.draw_networkx_nodes(G, pos, nodelist=rest, node_size=15, node_color="lightblue", ax=ax)
    nx.draw_networkx_nodes(G, pos, nodelist=hubs, node_color="crimson", ax=ax,
                           node_size=[degrees[k] * 2 for k in hubs])
    for key in hubs:
        x, y = pos[key]
        ax.text(x, y + 0.05, G.nodes[key]["name"], ha="center", fontsize=10, fontweight="bold")
    ax.axis("off")
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"\nSaved plot to {path} (hubs: degree > {hub_threshold})")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("edges", help="edge-list CSV from export_subgraph.py")
    parser.add_argument("--top", type=int, default=10, help="how many top-degree nodes to list")
    removal = parser.add_mutually_exclusive_group()
    removal.add_argument("--remove-top", type=int, default=2, metavar="K",
                         help="remove the K highest-degree nodes (default 2; 0 to skip)")
    removal.add_argument("--remove", nargs="+", metavar="NAME",
                         help="remove these nodes by name instead")
    parser.add_argument("--plot", metavar="PNG", help="save a network plot to this path")
    parser.add_argument("--hub-threshold", type=int, default=20,
                        help="degree above which nodes are highlighted in the plot")
    args = parser.parse_args()

    G, rel_types = load_graph(args.edges)
    if G.number_of_nodes() == 0:
        raise SystemExit(f"{args.edges} contains no edges.")

    remove_keys = resolve_names(G, args.remove) if args.remove else None
    summarize(G, rel_types, args.top)

    if remove_keys:
        removal_test(G, remove_keys)
    elif args.remove_top > 0:
        top = sorted(G.degree(), key=lambda kv: -kv[1])[:args.remove_top]
        removal_test(G, [key for key, _ in top])

    if args.plot:
        plot(G, args.plot, args.hub_threshold)


if __name__ == "__main__":
    main()
