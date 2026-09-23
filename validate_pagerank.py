"""Check that networkx PageRank matches Neo4j Graph Data Science (GDS) PageRank.

Loads an exported subgraph into a Neo4j instance that has the GDS plugin,
runs PageRank there and in networkx with the same settings, and compares:

  - top-k overlap
  - Spearman rank correlation across all nodes
  - largest per-node score difference after scaling both to sum to 1
    (GDS does not normalize scores, so raw values are not comparable)

With --seed, both sides run personalized PageRank restarting at the named
nodes, which is the variant PSEVs are built on.

The script only touches nodes labeled ValidationNode, which it deletes and
recreates on every run. Use a local, disposable database, never the shared
SPOKE server.
"""

import argparse

import networkx as nx
from scipy.stats import spearmanr

from graph_io import connect, describe, load_graph, resolve_names

LABEL = "ValidationNode"
GRAPH_NAME = "pagerank_validation"
BATCH = 5000

STANDARD_QUERY = f"""
CALL gds.pageRank.stream($graph, {{
    dampingFactor: $damping, maxIterations: $max_iter, tolerance: $tol
}})
YIELD nodeId, score
RETURN gds.util.asNode(nodeId).key AS key, score
"""

PERSONALIZED_QUERY = f"""
MATCH (s:{LABEL}) WHERE s.key IN $seed_keys
WITH collect(s) AS sources
CALL gds.pageRank.stream($graph, {{
    dampingFactor: $damping, maxIterations: $max_iter, tolerance: $tol,
    sourceNodes: sources
}})
YIELD nodeId, score
RETURN gds.util.asNode(nodeId).key AS key, score
"""


def batches(items, size):
    for i in range(0, len(items), size):
        yield items[i:i + size]


def load_into_neo4j(session, G):
    session.run(f"CREATE CONSTRAINT validation_key IF NOT EXISTS "
                f"FOR (n:{LABEL}) REQUIRE n.key IS UNIQUE").consume()
    session.run(f"MATCH (n:{LABEL}) DETACH DELETE n").consume()

    nodes = [{"key": k, "name": d["name"]} for k, d in G.nodes(data=True)]
    for chunk in batches(nodes, BATCH):
        session.run(f"UNWIND $rows AS row CREATE (:{LABEL} {{key: row.key, name: row.name}})",
                    rows=chunk).consume()

    edges = [{"a": a, "b": b} for a, b in G.edges()]
    for chunk in batches(edges, BATCH):
        session.run(f"""
            UNWIND $rows AS row
            MATCH (a:{LABEL} {{key: row.a}}), (b:{LABEL} {{key: row.b}})
            CREATE (a)-[:LINKED]->(b)
        """, rows=chunk).consume()


def gds_pagerank(session, damping, max_iter, tol, seed_keys):
    session.run("CALL gds.graph.drop($graph, false)", graph=GRAPH_NAME).consume()
    session.run(f"CALL gds.graph.project($graph, '{LABEL}', "
                f"{{LINKED: {{orientation: 'UNDIRECTED'}}}})", graph=GRAPH_NAME).consume()
    try:
        params = {"graph": GRAPH_NAME, "damping": damping, "max_iter": max_iter, "tol": tol}
        if seed_keys:
            result = session.run(PERSONALIZED_QUERY, seed_keys=seed_keys, **params)
        else:
            result = session.run(STANDARD_QUERY, **params)
        return {rec["key"]: rec["score"] for rec in result}
    finally:
        session.run("CALL gds.graph.drop($graph, false)", graph=GRAPH_NAME).consume()


def compare(G, nx_scores, gds_scores, top_k):
    missing = set(nx_scores) ^ set(gds_scores)
    if missing:
        raise SystemExit(f"{len(missing)} nodes present in only one result; graphs differ.")

    keys = list(nx_scores)
    nx_top = sorted(keys, key=lambda k: -nx_scores[k])[:top_k]
    gds_top = sorted(keys, key=lambda k: -gds_scores[k])[:top_k]
    overlap = set(nx_top) & set(gds_top)

    rho, _ = spearmanr([nx_scores[k] for k in keys], [gds_scores[k] for k in keys])
    nx_total, gds_total = sum(nx_scores.values()), sum(gds_scores.values())
    max_diff = max(abs(nx_scores[k] / nx_total - gds_scores[k] / gds_total) for k in keys)

    print(f"\n{'rank':<5}{'networkx':<40}{'GDS':<40}")
    for i, (a, b) in enumerate(zip(nx_top, gds_top), start=1):
        print(f"{i:<5}{describe(G, a):<40}{describe(G, b):<40}")

    print(f"\nTop-{top_k} overlap:           {len(overlap)}/{top_k}")
    print(f"Spearman correlation:     {rho:.4f} (all {len(keys)} nodes)")
    print(f"Max normalized score gap: {max_diff:.2e}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("edges", help="edge-list CSV from export_subgraph.py")
    parser.add_argument("--seed", nargs="+", metavar="NAME",
                        help="run personalized PageRank restarting at these nodes")
    parser.add_argument("--top-k", type=int, default=10)
    parser.add_argument("--damping", type=float, default=0.85)
    parser.add_argument("--max-iter", type=int, default=100)
    parser.add_argument("--tol", type=float, default=1e-8)
    parser.add_argument("--uri", help="Neo4j URI (default: $NEO4J_URI or bolt://localhost:7687)")
    parser.add_argument("--user", help="Neo4j user (default: $NEO4J_USER or neo4j)")
    parser.add_argument("--database", help="Neo4j database with GDS installed")
    args = parser.parse_args()

    G, _ = load_graph(args.edges)
    print(f"Loaded {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    seed_keys = resolve_names(G, args.seed) if args.seed else None
    personalization = {k: 1.0 for k in seed_keys} if seed_keys else None
    mode = f"personalized ({', '.join(args.seed)})" if seed_keys else "standard"
    print(f"Mode: {mode}")

    nx_scores = nx.pagerank(G, alpha=args.damping, personalization=personalization,
                            max_iter=args.max_iter, tol=args.tol)

    driver = connect(args.uri, args.user)
    try:
        with driver.session(database=args.database) as session:
            load_into_neo4j(session, G)
            gds_scores = gds_pagerank(session, args.damping, args.max_iter, args.tol, seed_keys)
    finally:
        driver.close()

    compare(G, nx_scores, gds_scores, args.top_k)


if __name__ == "__main__":
    main()
