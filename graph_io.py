"""Shared helpers: Neo4j connection and the edge-list CSV format used by every script.

Every script in this repo reads or writes the same CSV format: one row per
unique SPOKE relationship, with both endpoints identified by their Neo4j
element ID. Keying nodes by element ID rather than by name keeps two
different nodes that happen to share a name (for example a Gene and a
Compound) from being merged into one.
"""

import csv
import getpass
import os
from collections import Counter

import networkx as nx

EDGE_FIELDS = [
    "rel_key", "rel_type",
    "source_key", "source_label", "source_name",
    "target_key", "target_label", "target_name",
]


def connect(uri=None, user=None):
    """Open a Neo4j driver.

    Connection details come from arguments, then the environment variables
    NEO4J_URI and NEO4J_USER. The password is read only from NEO4J_PASSWORD
    or an interactive prompt, so it never appears in code or shell history.
    """
    from neo4j import GraphDatabase

    uri = uri or os.environ.get("NEO4J_URI", "bolt://localhost:7687")
    user = user or os.environ.get("NEO4J_USER", "neo4j")
    password = os.environ.get("NEO4J_PASSWORD") or getpass.getpass(
        f"Neo4j password for {user}@{uri}: "
    )
    driver = GraphDatabase.driver(uri, auth=(user, password))
    driver.verify_connectivity()
    return driver


def write_edges(path, rows):
    """Write edge rows (dicts with EDGE_FIELDS keys) to CSV."""
    with open(path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=EDGE_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def load_graph(path):
    """Load an edge-list CSV into an undirected networkx Graph.

    Returns (G, rel_type_counts). Each node carries `label` and `name`
    attributes. Relationship types are counted once per unique
    relationship. Multiple relationships between the same two nodes
    become a single edge in G, and self-loops are dropped.
    """
    G = nx.Graph()
    rel_types = Counter()
    seen = set()

    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            if row["rel_key"] in seen:
                continue
            seen.add(row["rel_key"])
            rel_types[row["rel_type"]] += 1

            for side in ("source", "target"):
                G.add_node(
                    row[f"{side}_key"],
                    label=row[f"{side}_label"],
                    name=row[f"{side}_name"],
                )
            if row["source_key"] != row["target_key"]:
                G.add_edge(row["source_key"], row["target_key"])

    return G, rel_types


def describe(G, key):
    """Human-readable 'name (Label)' for a node."""
    node = G.nodes[key]
    return f"{node['name']} ({node['label']})"


def resolve_names(G, names):
    """Map node names to keys. Raises if a name is missing or ambiguous."""
    by_name = {}
    for key, data in G.nodes(data=True):
        by_name.setdefault(data["name"], []).append(key)

    keys = []
    for name in names:
        matches = by_name.get(name, [])
        if not matches:
            raise SystemExit(f"No node named {name!r} in this graph.")
        if len(matches) > 1:
            labels = ", ".join(G.nodes[k]["label"] for k in matches)
            raise SystemExit(f"{name!r} matches {len(matches)} nodes ({labels}); ambiguous.")
        keys.append(matches[0])
    return keys
