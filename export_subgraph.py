"""Export the 2-hop neighborhood of a SPOKE Pathway node to an edge-list CSV.

Hop 1 is every node directly connected to the pathway (for a pathway,
mostly its member genes). Hop 2 is every node connected to those.

Three sampling modes:

  full (default)       every hop-2 relationship. Unbiased; use this unless
                       the neighborhood is too large to work with.
  --per-node-limit K   at most K hop-2 relationships per hop-1 node, chosen
                       deterministically (ordered by neighbor name), so
                       every hop-1 node is equally represented.
  --naive-limit N      a single LIMIT N over the whole 2-hop result with no
                       ordering. This reproduces the biased export described
                       in the README and should only be used to demonstrate
                       that bias.

After exporting, the script prints how many hop-2 relationships each hop-1
node received, which makes any sampling imbalance visible immediately.
"""

import argparse
import sys
from collections import Counter

from graph_io import connect, write_edges

HOP1_QUERY = """
MATCH (p:Pathway {identifier: $pathway_id})-[r1]-(n1)
RETURN p AS a, r1 AS r, n1 AS b
"""

HOP2_FULL_QUERY = """
MATCH (p:Pathway {identifier: $pathway_id})--(n1)
WITH DISTINCT p, n1
MATCH (n1)-[r2]-(n2)
WHERE n2 <> p
RETURN n1 AS a, r2 AS r, n2 AS b
"""

HOP2_PER_NODE_QUERY = """
MATCH (p:Pathway {identifier: $pathway_id})--(n1)
WITH DISTINCT p, n1
CALL {
    WITH p, n1
    MATCH (n1)-[r2]-(n2)
    WHERE n2 <> p
    RETURN r2, n2
    ORDER BY n2.name, n2.identifier
    LIMIT $limit
}
RETURN n1 AS a, r2 AS r, n2 AS b
"""

NAIVE_QUERY = """
MATCH (p:Pathway {identifier: $pathway_id})-[r1]-(n1)-[r2]-(n2)
RETURN p, r1, n1, r2, n2
LIMIT $limit
"""


def node_fields(node):
    labels = sorted(node.labels)
    label = labels[0] if labels else "Unknown"
    name = node.get("name") or node.get("identifier") or node.element_id
    return node.element_id, label, str(name)


def edge_row(rel, a, b):
    a_key, a_label, a_name = node_fields(a)
    b_key, b_label, b_name = node_fields(b)
    return {
        "rel_key": rel.element_id, "rel_type": rel.type,
        "source_key": a_key, "source_label": a_label, "source_name": a_name,
        "target_key": b_key, "target_label": b_label, "target_name": b_name,
    }


def export(session, pathway_id, per_node_limit=None, naive_limit=None):
    """Return (rows, hop1_nodes, hop2_counts). rows are unique by relationship."""
    rows = {}
    hop1 = {}
    hop2_counts = Counter()

    def add(rel, a, b):
        rows.setdefault(rel.element_id, edge_row(rel, a, b))

    if naive_limit is not None:
        for rec in session.run(NAIVE_QUERY, pathway_id=pathway_id, limit=naive_limit):
            add(rec["r1"], rec["p"], rec["n1"])
            add(rec["r2"], rec["n1"], rec["n2"])
            hop1[rec["n1"].element_id] = rec["n1"]
            hop2_counts[rec["n1"].element_id] += 1
        # The naive query can only see hop-1 nodes that made it under the
        # limit, so look up the full hop-1 set separately for the report.
        for rec in session.run(HOP1_QUERY, pathway_id=pathway_id):
            hop1.setdefault(rec["b"].element_id, rec["b"])
        return list(rows.values()), hop1, hop2_counts

    for rec in session.run(HOP1_QUERY, pathway_id=pathway_id):
        add(rec["r"], rec["a"], rec["b"])
        hop1[rec["b"].element_id] = rec["b"]

    if per_node_limit is None:
        result = session.run(HOP2_FULL_QUERY, pathway_id=pathway_id)
    else:
        result = session.run(HOP2_PER_NODE_QUERY, pathway_id=pathway_id, limit=per_node_limit)
    for rec in result:
        add(rec["r"], rec["a"], rec["b"])
        hop2_counts[rec["a"].element_id] += 1

    return list(rows.values()), hop1, hop2_counts


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("pathway_id", help="SPOKE Pathway identifier, e.g. WP4629_r127017")
    parser.add_argument("-o", "--out", required=True, help="output edge-list CSV")
    parser.add_argument("--uri", help="Neo4j URI (default: $NEO4J_URI or bolt://localhost:7687)")
    parser.add_argument("--user", help="Neo4j user (default: $NEO4J_USER or neo4j)")
    parser.add_argument("--database", help="Neo4j database name (default: server default)")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--per-node-limit", type=int, metavar="K")
    mode.add_argument("--naive-limit", type=int, metavar="N")
    args = parser.parse_args()

    driver = connect(args.uri, args.user)
    try:
        with driver.session(database=args.database) as session:
            rows, hop1, hop2_counts = export(
                session, args.pathway_id, args.per_node_limit, args.naive_limit
            )
    finally:
        driver.close()

    if not hop1:
        sys.exit(f"No Pathway node with identifier {args.pathway_id!r}, or it has no neighbors.")

    write_edges(args.out, rows)
    print(f"Wrote {len(rows)} relationships to {args.out}")

    print(f"\nHop-2 relationships per hop-1 node ({len(hop1)} hop-1 nodes):")
    for key, node in sorted(hop1.items(), key=lambda kv: -hop2_counts[kv[0]]):
        _, label, name = node_fields(node)
        print(f"  {name} ({label}): {hop2_counts[key]}")
    missing = sum(1 for key in hop1 if hop2_counts[key] == 0)
    if missing:
        print(f"\nWarning: {missing} of {len(hop1)} hop-1 nodes have no hop-2 relationships "
              f"in this export.")


if __name__ == "__main__":
    main()
