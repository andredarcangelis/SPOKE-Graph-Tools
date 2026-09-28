#!/usr/bin/env python3
"""
prepare_inputs.py  (PSEV Artifact Validation Pipeline, input step)

Builds the three input files for psev_artifact_validation.py from real sources.

  spoke-summary  Count nodes per type and edges per type in a SPOKE Neo4j instance, and
                 show the property names on each node type.
  spoke-export   Write nodes.csv and edges.csv, keeping whole node types (default: the 12
                 types used in Nelson et al. 2021) and every edge between kept nodes.
  de-columns     List the fold-change columns and gene ID columns in a GeneLab/OSDR
                 differential expression table.
  build-fc       Convert one or more differential expression tables into fc_long.csv,
                 mapped to SPOKE Gene nodes (with mouse-to-human mapping if needed).

See README.md, "Running on real data", for the full workflow.
"""

import argparse
import csv
import getpass
import os
import re
import sys

import numpy as np
import pandas as pd

PAPER_LABELS = ["Anatomy", "BiologicalProcess", "CellularComponent", "Compound", "Disease",
                "Gene", "MolecularFunction", "Pathway", "PharmacologicalClass", "Protein",
                "SideEffect", "Symptom"]


# ---- SPOKE / Neo4j ----

def connect(args):
    """Return a function that runs a Cypher query and yields records."""
    try:
        from neo4j import GraphDatabase
    except ImportError:
        sys.exit("The Neo4j driver is not installed. Run: pip install neo4j")
    password = os.environ.get("NEO4J_PASSWORD") or getpass.getpass("Neo4j password: ")
    driver = GraphDatabase.driver(args.uri, auth=(args.user, password))
    try:
        driver.verify_connectivity()
    except Exception as e:  # the driver raises several exception types for this
        sys.exit(f"Could not connect to {args.uri} as {args.user}: {e}\n"
                 f"Check the address and port (Bolt, usually 7687), the password, and that "
                 f"this network is allowed to reach the server.")

    def run(query, **params):
        with driver.session(database=args.database) as session:
            yield from session.run(query, **params)

    return run


def spoke_summary(run):
    labels = [r["label"] for r in run("CALL db.labels() YIELD label RETURN label")]
    print("Nodes per type")
    for label in sorted(labels):
        count = next(run(f"MATCH (n:`{label}`) RETURN count(n) AS n"))["n"]
        keys = next(run(f"MATCH (n:`{label}`) RETURN keys(n) AS k LIMIT 1"), {"k": []})["k"]
        mark = "*" if label in PAPER_LABELS else " "
        print(f" {mark} {label:28s} {count:>12,}   properties: {', '.join(sorted(keys))}")
    print("   (* = one of the 12 node types used in Nelson et al. 2021)\n")

    types = [r["t"] for r in run("CALL db.relationshipTypes() YIELD relationshipType AS t RETURN t")]
    print("Edges per type")
    for t in sorted(types):
        count = next(run(f"MATCH ()-[r:`{t}`]->() RETURN count(r) AS n"))["n"]
        print(f"   {t:40s} {count:>12,}")


def parse_filters(items):
    """'Label::cypher condition on n' -> {Label: condition}."""
    out = {}
    for item in items or []:
        if "::" not in item:
            sys.exit(f"--node-filter must look like 'Label::condition', got: {item}")
        label, cond = item.split("::", 1)
        out[label.strip()] = cond.strip()
    return out


SETTINGS_FILE = "export_settings.txt"


def export_settings_text(labels, filters, id_prop, exclude_rels):
    lines = [f"labels: {' '.join(labels)}", f"id_property: {id_prop}",
             f"exclude_edge_types: {' '.join(sorted(exclude_rels))}"]
    lines += [f"node_filter: {k}::{v}" for k, v in sorted(filters.items())]
    return "\n".join(lines) + "\n"


def spoke_export(run, labels, filters, id_prop, name_prop, outdir, exclude_rels):
    """Write nodes.csv and edges.csv. Files are written under temporary names and renamed
    only when the whole export finishes, so an interrupted export never looks complete."""
    os.makedirs(outdir, exist_ok=True)
    unknown = set(filters) - set(labels)
    if unknown:
        sys.exit(f"--node-filter names types that are not being exported: {', '.join(sorted(unknown))}")
    paths = {name: os.path.join(outdir, name) for name in ("nodes.csv", "edges.csv", SETTINGS_FILE)}
    for p in paths.values():
        if os.path.exists(p):
            os.remove(p)

    # A node with several labels is exported once, under the first of its labels that is
    # being kept. Node and edge queries use the same rule, so their IDs always agree.
    canon = "[l IN labels({v}) WHERE l IN $labels][0]"
    keep = set()
    missing_id = 0
    with open(paths["nodes.csv"] + ".tmp", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["node_id", "name", "type"])
        for label in labels:
            where = f"WHERE {filters[label]}" if label in filters else ""
            q = (f"MATCH (n:`{label}`) {where} "
                 f"RETURN n.`{id_prop}` AS id, n.`{name_prop}` AS name, {canon.format(v='n')} AS canon")
            count = 0
            for r in run(q, labels=labels):
                if r["canon"] != label:
                    continue  # exported under another of its labels
                if r["id"] is None:
                    missing_id += 1
                    continue
                # Identifiers repeat across types, so the node ID includes the type.
                nid = f"{label}:{r['id']}"
                if nid in keep:
                    continue
                keep.add(nid)
                w.writerow([nid, "" if r["name"] is None else r["name"], label])
                count += 1
            print(f"  {label:28s} {count:>12,} nodes", flush=True)
    if missing_id:
        print(f"  skipped {missing_id:,} nodes with no '{id_prop}' property")
    if not keep:
        sys.exit(f"No nodes exported. Check the type names and --id-property '{id_prop}' "
                 f"with the spoke-summary command.")

    types = [r["t"] for r in run("CALL db.relationshipTypes() YIELD relationshipType AS t RETURN t")]
    total = 0
    with open(paths["edges.csv"] + ".tmp", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["source", "target", "edge_type"])
        for t in sorted(types):
            if t in exclude_rels:
                continue
            q = (f"MATCH (a)-[r:`{t}`]->(b) "
                 f"WITH a, b, {canon.format(v='a')} AS la, {canon.format(v='b')} AS lb "
                 f"WHERE la IS NOT NULL AND lb IS NOT NULL "
                 f"RETURN la, a.`{id_prop}` AS ia, lb, b.`{id_prop}` AS ib")
            count = 0
            for r in run(q, labels=labels):
                s, d = f"{r['la']}:{r['ia']}", f"{r['lb']}:{r['ib']}"
                # Edges to nodes removed by --node-filter are dropped here.
                if s in keep and d in keep:
                    w.writerow([s, d, t])
                    count += 1
            if count:
                print(f"  {t:40s} {count:>12,} edges", flush=True)
            total += count

    with open(paths[SETTINGS_FILE] + ".tmp", "w") as f:
        f.write(export_settings_text(labels, filters, id_prop, exclude_rels))
    for p in paths.values():
        os.replace(p + ".tmp", p)
    print(f"Wrote {len(keep):,} nodes and {total:,} edges to {outdir}/")


# ---- differential expression ----

ID_COLUMNS = ["ENTREZID", "ENTREZ", "EntrezID", "SYMBOL", "Symbol", "ENSEMBL"]
FLIP_VALUES = {"": False, "false": False, "no": False, "0": False,
               "true": True, "yes": True, "1": True}


def read_table(path, nrows=None):
    sep = "," if path.lower().endswith(".csv") else "\t"
    return pd.read_csv(path, sep=sep, nrows=nrows, low_memory=False)


def de_columns(path):
    df = read_table(path, nrows=5)
    fc_cols = [c for c in df.columns if c.lower().startswith("log2fc")]
    ids = [c for c in df.columns if c in ID_COLUMNS]
    print(f"Gene ID columns: {', '.join(ids) or 'none of ' + ', '.join(ID_COLUMNS)}")
    print("Fold-change columns (copy these into the log2fc_column field of your config):")
    for c in fc_cols:
        print(f"  {c}")
    if not fc_cols:
        print("  none found; columns are: " + ", ".join(df.columns))


def describe_contrast(col):
    """'Log2fc_(Space Flight)v(Ground Control)' -> 'Space Flight vs Ground Control'."""
    m = re.match(r"(?i)log2fc_\((.*)\)v\((.*)\)$", col)
    return f"{m.group(1)} vs {m.group(2)}" if m else col


def clean_entrez(values):
    """GeneLab Entrez fields can be floats ('12345.0') or lists ('12345|67890'); keep the first."""
    s = values.astype(str).str.split(r"[|;,]").str[0].str.strip()
    s = s.str.replace(r"\.0$", "", regex=True)
    return s.where(s.str.fullmatch(r"\d+"), None)


def load_homologs(path):
    """Mouse-to-human pairs as a DataFrame with columns mouse_entrez, mouse_symbol,
    human_entrez, human_symbol. Accepts MGI's HOM_MouseHumanSequence.rpt, or a CSV/TSV
    with columns mouse_entrez, human_entrez (and optionally the symbol columns)."""
    # Read as text: MGI leaves some Entrez IDs blank, and reading numbers would turn the
    # rest into decimals ('12345.0') that no longer match.
    sep = "," if path.lower().endswith(".csv") else "\t"
    df = pd.read_csv(path, sep=sep, dtype=str, keep_default_na=False)
    if "Common Organism Name" in df.columns:
        key = "DB Class Key"
        need = [key, "EntrezGene ID", "Symbol"]
        if not set(need) <= set(df.columns):
            sys.exit(f"{path}: expected MGI columns {need}")
        org = df["Common Organism Name"].str.lower()
        mouse = df[org.str.contains("mouse")][need]
        human = df[org.str.contains("human")][need]
        mouse.columns = [key, "mouse_entrez", "mouse_symbol"]
        human.columns = [key, "human_entrez", "human_symbol"]
        df = mouse.merge(human, on=key).drop(columns=key)
    elif not {"mouse_entrez", "human_entrez"} <= set(df.columns):
        sys.exit("Homolog file not recognized. Use MGI's HOM_MouseHumanSequence.rpt or a table "
                 "with mouse_entrez and human_entrez columns.")
    for col in ("mouse_entrez", "human_entrez"):
        df[col] = clean_entrez(df[col])
    for col in ("mouse_symbol", "human_symbol"):
        if col in df.columns:
            df[col] = df[col].str.strip().replace("", None)
    return df


def build_fc(args):
    config = pd.read_csv(args.config, dtype=str).fillna("")
    required = {"de_file", "log2fc_column", "comparison", "group"}
    if not required <= set(config.columns):
        sys.exit(f"Config needs columns {sorted(required)} (optional: study, flip, padj_column)")
    if config["comparison"].duplicated().any():
        sys.exit("Each comparison name in the config must be unique")

    # SPOKE Gene nodes, keyed the same way as the expression IDs will be.
    nodes = pd.read_csv(args.nodes, dtype=str, keep_default_na=False)
    genes = nodes[nodes["type"] == args.gene_label]
    if args.spoke_gene_key == "entrez":
        gene_lookup = {nid.split(":", 1)[1]: nid for nid in genes["node_id"]}
    else:
        counts = genes["name"].value_counts()
        unique = genes[genes["name"].map(counts) == 1]  # ambiguous symbols are dropped
        gene_lookup = dict(zip(unique["name"], unique["node_id"]))
    print(f"SPOKE {args.gene_label} nodes available for matching: {len(gene_lookup):,}")

    homologs = load_homologs(args.homologs) if args.species == "mouse" else None
    if args.species == "mouse" and homologs is None:
        sys.exit("--homologs is required for mouse data")

    tables = {}
    rows, report = [], []
    for _, c in config.iterrows():
        path = os.path.join(os.path.dirname(os.path.abspath(args.config)), c["de_file"]) \
            if not os.path.isabs(c["de_file"]) else c["de_file"]
        if path not in tables:
            tables[path] = read_table(path)
        de = tables[path]
        if c["log2fc_column"] not in de.columns:
            sys.exit(f"{c['comparison']}: column '{c['log2fc_column']}' not in {c['de_file']}. "
                     f"Run: python prepare_inputs.py de-columns {c['de_file']}")

        id_col = args.id_column or next((x for x in ID_COLUMNS if x in de.columns and x != "ENSEMBL"), None)
        if id_col is None:
            sys.exit(f"{c['de_file']}: no Entrez or symbol column found; set --id-column")
        id_kind = "entrez" if "entrez" in id_col.lower() else "symbol"

        t = pd.DataFrame({"gid": de[id_col], "fc": pd.to_numeric(de[c["log2fc_column"]], errors="coerce")})
        flip_text = c.get("flip", "").strip().lower()
        if flip_text not in FLIP_VALUES:
            sys.exit(f"{c['comparison']}: flip must be blank, true or false (got '{c['flip']}')")
        # Without a study label, a comparison is its own study, so the direction filter
        # never compares it with comparisons from other studies.
        study = c.get("study", "").strip() or c["comparison"]
        r = {"comparison": c["comparison"], "group": c["group"], "study": study,
             "contrast": describe_contrast(c["log2fc_column"]), "flipped": FLIP_VALUES[flip_text],
             "rows": len(t)}

        if args.max_padj is not None:
            padj_col = c.get("padj_column") or re.sub(r"(?i)^log2fc_", "Adj.p.value_", c["log2fc_column"])
            if padj_col not in de.columns:
                sys.exit(f"{c['comparison']}: adjusted p-value column '{padj_col}' not found; "
                         f"set padj_column in the config")
            t = t[pd.to_numeric(de[padj_col], errors="coerce") <= args.max_padj]
        r["after_padj_filter"] = len(t)

        t["gid"] = clean_entrez(t["gid"]) if id_kind == "entrez" else t["gid"].astype(str).str.strip()
        t = t.dropna()
        t = t[np.isfinite(t["fc"])]
        if r["flipped"]:
            t["fc"] = -t["fc"]  # orient every comparison in a group the same way
        t = t.groupby("gid", as_index=False)["fc"].mean()  # e.g. several Ensembl IDs, one Entrez ID
        r["genes_with_id_and_fc"] = len(t)

        # Mouse -> human. One mouse gene mapping to several human genes is ambiguous and
        # dropped; several mouse genes mapping to one human gene are averaged (as in the paper).
        key_out = "human_entrez" if args.spoke_gene_key == "entrez" else "human_symbol"
        if homologs is not None:
            key_in = "mouse_entrez" if id_kind == "entrez" else "mouse_symbol"
            if key_in not in homologs.columns or key_out not in homologs.columns:
                sys.exit(f"Homolog table has no {key_in} / {key_out} columns for this ID type")
            pairs = homologs[[key_in, key_out]].dropna().drop_duplicates()
            n_targets = pairs.groupby(key_in)[key_out].transform("nunique")
            pairs = pairs[n_targets == 1]
            t = t.merge(pairs, left_on="gid", right_on=key_in)[[key_out, "fc"]]
            t.columns = ["gid", "fc"]
            r["mapped_to_human"] = t["gid"].nunique()
        elif id_kind != args.spoke_gene_key:
            sys.exit(f"Human data with {id_kind} IDs cannot match SPOKE genes keyed by "
                     f"{args.spoke_gene_key}; set --spoke-gene-key {id_kind} or --id-column")

        t["node_id"] = t["gid"].map(gene_lookup)
        t = t.dropna(subset=["node_id"]).groupby("node_id", as_index=False)["fc"].mean()
        r["matched_to_spoke"] = len(t)
        report.append(r)
        for nid, fc in zip(t["node_id"], t["fc"]):
            rows.append((c["comparison"], c["group"], r["study"], nid, fc))

    fc = pd.DataFrame(rows, columns=["comparison", "group", "study", "node_id", "log2fc"])

    # Paper's filter: within one study and group, drop a gene whose fold change points in
    # different directions across that study's comparisons (e.g. day 29 up, day 56 down).
    if not args.no_direction_filter:
        keys = [fc["study"], fc["group"], fc["node_id"]]
        up = (fc["log2fc"] > 0).groupby(keys).transform("any")
        down = (fc["log2fc"] < 0).groupby(keys).transform("any")
        mixed = up & down
        dropped = fc[mixed].groupby("comparison").size()
        fc = fc[~mixed]
        for r in report:
            r["dropped_mixed_direction"] = int(dropped.get(r["comparison"], 0))
            r["final_genes"] = int((fc["comparison"] == r["comparison"]).sum())

    fc.to_csv(args.out, index=False)
    rep = pd.DataFrame(report)
    rep.to_csv(args.report, index=False)
    pd.set_option("display.width", 200)
    print("\nComparisons (check that every comparison in a group points the same way):")
    print(rep[["comparison", "group", "study", "contrast", "flipped"]].to_string(index=False))
    print("\nGene counts at each step:")
    print(rep.drop(columns=["group", "study", "contrast", "flipped"]).to_string(index=False))
    print(f"\nWrote {args.out} ({len(fc):,} rows) and {args.report}")


# ---- command line ----

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def neo4j_args(p):
        p.add_argument("--uri", required=True, help="e.g. bolt://HOST:7687")
        p.add_argument("--user", required=True)
        p.add_argument("--database", default=None, help="Neo4j database name, if not the default")

    p = sub.add_parser("spoke-summary", help="count nodes and edges per type")
    neo4j_args(p)

    def export_opts(p):
        p.add_argument("--labels", nargs="+", default=PAPER_LABELS, help="node types to keep")
        p.add_argument("--node-filter", action="append",
                       help="'Label::condition on n', e.g. 'Protein::n.org_ncbi_id = 9606' (repeatable)")
        p.add_argument("--id-property", default="identifier")
        p.add_argument("--exclude-edge-types", nargs="*", default=[])

    p = sub.add_parser("spoke-export", help="write nodes.csv and edges.csv")
    neo4j_args(p)
    export_opts(p)
    p.add_argument("--name-property", default="name")
    p.add_argument("--outdir", default="spoke_export")

    p = sub.add_parser("export-settings",
                       help="print the settings record an export with these options would write")
    export_opts(p)

    p = sub.add_parser("de-columns", help="list fold-change and ID columns in a DE table")
    p.add_argument("de_file")

    p = sub.add_parser("build-fc", help="DE tables -> fc_long.csv")
    p.add_argument("--config", required=True,
                   help="CSV: de_file, log2fc_column, comparison, group[, study, flip, padj_column]")
    p.add_argument("--nodes", required=True, help="nodes.csv from spoke-export")
    p.add_argument("--species", choices=["mouse", "human"], required=True)
    p.add_argument("--homologs", help="MGI HOM_MouseHumanSequence.rpt (mouse data only)")
    p.add_argument("--id-column", default=None, help="gene ID column in the DE tables (default: auto)")
    p.add_argument("--gene-label", default="Gene", help="SPOKE node type for genes")
    p.add_argument("--spoke-gene-key", choices=["entrez", "symbol"], default="entrez",
                   help="what the SPOKE Gene identifier is (entrez) or match on name (symbol)")
    p.add_argument("--max-padj", type=float, default=None,
                   help="keep only genes with adjusted p <= this (default: keep all)")
    p.add_argument("--no-direction-filter", action="store_true")
    p.add_argument("--out", default="fc_long.csv")
    p.add_argument("--report", default="mapping_report.csv")

    args = ap.parse_args()
    if args.cmd == "de-columns":
        de_columns(args.de_file)
    elif args.cmd == "export-settings":
        print(export_settings_text(args.labels, parse_filters(args.node_filter),
                                   args.id_property, set(args.exclude_edge_types)), end="")
    elif args.cmd == "build-fc":
        build_fc(args)
    else:
        run = connect(args)
        if args.cmd == "spoke-summary":
            spoke_summary(run)
        else:
            spoke_export(run, args.labels, parse_filters(args.node_filter), args.id_property,
                         args.name_property, args.outdir, set(args.exclude_edge_types))


if __name__ == "__main__":
    main()
