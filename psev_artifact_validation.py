#!/usr/bin/env python3
"""
psev_artifact_validation.py  (PSEV Artifact Validation Pipeline, scoring step)

Tests which nodes in fold-change-weighted personalized PageRank (PSEV-style) output on a
SPOKE graph reflect the input data, rather than graph structure. Runs on one comparison or
on many comparisons pooled across studies.

For each node it reports:
  z_label  score vs. shuffling fold changes across the input genes (per comparison)
  z_xswap  score vs. degree-preserving rewired graphs (per comparison, optional)
  Stouffer-combined z per group, and Welch's t per test group vs. a control group

Inputs:  nodes.csv (node_id, name, type)
         edges.csv (source, target[, edge_type])
         fc.csv    (comparison, group, node_id, log2fc), one row per gene per comparison
Outputs: node_summary.csv, per_comparison.npz, run_info.json in --outdir

See README.md for the method, its differences from Nelson et al. 2021, and caveats.

Example:
  python psev_artifact_validation.py --nodes nodes.csv --edges edges.csv --fc fc_long.csv
      --control-group ground_vs_baseline --n-perm 200 --n-xswap 20 --outdir results
  (one command; wrapped here for width)
"""

import argparse
import json
import os
import time

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import norm, spearmanr
from scipy.stats import t as tdist

T0 = time.time()


def log(msg):
    print(f"[{time.time() - T0:7.1f}s] {msg}", flush=True)


# ---------------------------------------------------------------- graph + PageRank

def build_transition(src, dst, n, directed=False):
    """P^T (CSR) for the row-stochastic random-walk matrix, plus node degree."""
    # Undirected by default: the walker can traverse any edge in either direction.
    if not directed:
        src, dst = np.concatenate([src, dst]), np.concatenate([dst, src])
    # Parallel edges (same pair, different edge types) are summed, so they weight the walk.
    A = sp.csr_matrix((np.ones(len(src), dtype=np.float32), (src, dst)), shape=(n, n))
    A.sum_duplicates()
    deg = np.asarray(A.sum(axis=1)).ravel()
    inv = np.zeros_like(deg)
    inv[deg > 0] = 1.0 / deg[deg > 0]
    return (sp.diags(inv.astype(np.float32)) @ A).T.tocsr(), deg


def ppr(PT, S, restart=0.1, tol=1e-6, max_iter=300):
    """Personalized PageRank for each column of seed matrix S (n x k), power iteration.
    PPR is linear, so a signed seed vector gives the FC-weighted sum of per-gene PPR vectors.
    Mass reaching nodes with no edges is dropped; every run shares this, so comparisons hold."""
    X = S.copy()
    scale = max(np.abs(S).sum(axis=0).max(), 1e-12)
    for _ in range(max_iter):
        Xn = (1.0 - restart) * (PT @ X) + restart * S
        err = np.abs(Xn - X).sum(axis=0).max() / scale
        X = Xn
        if err < tol:
            break
    return X


def seed_matrix(n, seeds, cols=None):
    cols = range(len(seeds)) if cols is None else cols
    S = np.zeros((n, len(cols)), dtype=np.float32)
    for j, c in enumerate(cols):
        gidx, w = seeds[c]
        S[gidx, j] = w
    return S


# ---------------------------------------------------------------- XSwap

def _keys(s, d, n, directed):
    if directed:
        return s.astype(np.int64) * n + d
    a, b = np.minimum(s, d).astype(np.int64), np.maximum(s, d).astype(np.int64)
    return a * n + b


def _xswap_block(s, d, n, rng, multiplier, directed):
    """Degree-preserving swaps within one edge type, in vectorized rounds.
    Swap (a->b, c->d) into (a->d, c->b). Rejects self-loops, edges that already exist,
    and proposals that collide with each other in the same round."""
    s, d = s.copy(), d.copy()
    E = len(s)
    keys = _keys(s, d, n, directed)
    attempted, accepted = 0, 0
    while attempted < multiplier * E:
        perm = rng.permutation(E)
        h = E // 2
        i, j = perm[:h], perm[h:2 * h]
        k1, k2 = _keys(s[i], d[j], n, directed), _keys(s[j], d[i], n, directed)
        ok = (s[i] != d[j]) & (s[j] != d[i]) & (k1 != k2)
        sk = np.sort(keys)

        def exists(k):
            pos = np.minimum(np.searchsorted(sk, k), len(sk) - 1)
            return sk[pos] == k

        ok &= ~exists(k1) & ~exists(k2)
        cand = np.concatenate([k1[ok], k2[ok]])
        u, c = np.unique(cand, return_counts=True)
        if (c > 1).any():
            oi = np.where(ok)[0]
            dup = u[c > 1]
            ok[oi[np.isin(k1[oi], dup) | np.isin(k2[oi], dup)]] = False
        ii, jj = i[ok], j[ok]
        dj = d[jj].copy()
        d[jj] = d[ii]
        d[ii] = dj
        keys[ii], keys[jj] = k1[ok], k2[ok]
        attempted += h
        accepted += ok.sum()
    return s, d, accepted / max(attempted, 1)


def xswap(src, dst, etype, n, rng, multiplier=10, directed=False):
    """Rewire each edge type separately so node types on each end and every node's
    degree per edge type are preserved."""
    src, dst = src.copy(), dst.copy()
    rates = []
    for t in np.unique(etype):
        m = np.where(etype == t)[0]
        if len(m) < 2:
            continue
        src[m], dst[m], r = _xswap_block(src[m], dst[m], n, rng, multiplier, directed)
        rates.append(r)
    return src, dst, float(np.mean(rates)) if rates else 0.0


# ---------------------------------------------------------------- statistics

def bh(p):
    p = np.asarray(p, dtype=float)
    out = np.full(len(p), np.nan)
    f = np.isfinite(p)
    m = f.sum()
    if m == 0:
        return out
    pf = p[f]
    order = np.argsort(pf)
    ranked = pf[order] * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(ranked[::-1])[::-1]
    tmp = np.empty(m)
    tmp[order] = np.minimum(q, 1.0)
    out[f] = tmp
    return out


def welch(A, B):
    """Row-wise Welch's t-test between A (n x k1) and B (n x k2).
    Treats comparisons as independent samples; comparisons sharing animals or controls are
    not, so p-values are approximate. Low power with few comparisons per group."""
    k1, k2 = A.shape[1], B.shape[1]
    m1, m2 = A.mean(1), B.mean(1)
    v1, v2 = A.var(1, ddof=1), B.var(1, ddof=1)
    se2 = v1 / k1 + v2 / k2
    with np.errstate(divide="ignore", invalid="ignore"):
        t = (m1 - m2) / np.sqrt(se2)
        df = se2 ** 2 / ((v1 / k1) ** 2 / (k1 - 1) + (v2 / k2) ** 2 / (k2 - 1))
    p = 2 * tdist.sf(np.abs(t), df)
    bad = ~(se2 > 0)
    t[bad], p[bad] = 0.0, 1.0
    return t, p


def zfrom(real, s1, s2, N):
    mean = s1 / N
    sd = np.sqrt(np.maximum(s2 / N - mean ** 2, 0))
    return np.where(sd > 0, (real - mean) / np.where(sd > 0, sd, 1), 0.0)


def require(df, cols, path):
    missing = [c for c in cols if c not in df.columns]
    if missing:
        raise SystemExit(f"{path} is missing column(s): {', '.join(missing)}")


def q_within(values_p, types):
    q = np.full(len(values_p), np.nan)
    for t in np.unique(types):
        m = types == t
        q[m] = bh(values_p[m])
    return q


# ---- main ----

def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nodes", required=True)
    ap.add_argument("--edges", required=True)
    ap.add_argument("--fc", required=True)
    ap.add_argument("--outdir", default="results")
    ap.add_argument("--control-group", default=None, help="group used as the Welch reference")
    ap.add_argument("--n-perm", type=int, default=200, help="label shuffles per comparison (0 = skip)")
    ap.add_argument("--n-xswap", type=int, default=0, help="rewired graphs (0 = skip)")
    ap.add_argument("--xswap-multiplier", type=float, default=10, help="swap attempts per edge")
    ap.add_argument("--batch", type=int, default=100)
    ap.add_argument("--restart", type=float, default=0.1)
    ap.add_argument("--tol", type=float, default=1e-6,
                    help="convergence tolerance for PageRank (relative L1 change per iteration)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--directed", action="store_true")
    ap.add_argument("--drop-edge-types", nargs="*", default=[])
    ap.add_argument("--rank-abs", action="store_true",
                    help="rank nodes by |score| instead of signed score (see README)")
    args = ap.parse_args()
    if args.batch < 1:
        raise SystemExit("--batch must be at least 1")
    if args.n_perm < 0 or args.n_xswap < 0:
        raise SystemExit("--n-perm and --n-xswap cannot be negative")
    if not 0 < args.restart < 1:
        raise SystemExit("--restart must be between 0 and 1")
    rng = np.random.default_rng(args.seed)
    os.makedirs(args.outdir, exist_ok=True)

    # ---- graph ----
    nodes = pd.read_csv(args.nodes, dtype=str, keep_default_na=False)
    require(nodes, ["node_id"], args.nodes)
    if "type" not in nodes.columns:
        nodes["type"] = "all"
    if nodes["node_id"].duplicated().any():
        raise SystemExit("nodes.csv has duplicate node_id values")
    index = pd.Index(nodes["node_id"])
    idx = pd.Series(np.arange(len(nodes)), index=index)
    n = len(nodes)
    edges = pd.read_csv(args.edges, dtype=str, keep_default_na=False)
    require(edges, ["source", "target"], args.edges)
    if "edge_type" not in edges.columns:
        edges["edge_type"] = "all"
        if args.n_xswap:
            log("WARNING: no edge_type column; XSwap will mix edge types across node types")
    if args.drop_edge_types:
        before = len(edges)
        edges = edges[~edges["edge_type"].isin(args.drop_edge_types)]
        log(f"Dropped {before - len(edges):,} edges of types {args.drop_edge_types}")
    src = index.get_indexer(edges["source"])
    dst = index.get_indexer(edges["target"])
    ok = (src >= 0) & (dst >= 0)
    if (~ok).any():
        log(f"Skipped {(~ok).sum():,} edges whose endpoints are not in nodes.csv")
    edges, src, dst = edges[ok], src[ok], dst[ok]
    etype = edges["edge_type"].astype(str).values
    PT, deg = build_transition(src, dst, n, args.directed)
    log(f"Graph: {n:,} nodes, {len(edges):,} edges, {len(np.unique(etype))} edge types")

    # ---- comparisons ----
    fc = pd.read_csv(args.fc, dtype={"node_id": str, "comparison": str, "group": str})
    require(fc, ["comparison", "group", "node_id", "log2fc"], args.fc)
    fc["log2fc"] = pd.to_numeric(fc["log2fc"], errors="coerce")
    bad = ~np.isfinite(fc["log2fc"])
    if bad.any():
        log(f"{bad.sum():,} fold-change rows have a missing or non-numeric log2fc and were skipped")
    missing = ~fc["node_id"].isin(index)
    if missing.any():
        log(f"{missing.sum():,} fold-change rows name nodes not in nodes.csv and were skipped")
    fc = fc[~missing & ~bad]
    fc = fc.groupby(["comparison", "group", "node_id"], as_index=False)["log2fc"].mean()
    comps = fc[["comparison", "group"]].drop_duplicates().sort_values(["group", "comparison"])
    if comps["comparison"].duplicated().any():
        raise SystemExit("Each comparison must belong to exactly one group")
    comp_names, comp_groups, seeds = [], [], []
    for c, g in zip(comps["comparison"], comps["group"]):
        sub = fc[fc["comparison"] == c]
        # Paper weights by -log2 FC; the sign flips z and t but not two-sided p or q.
        w = sub["log2fc"].values.astype(np.float64)
        total = np.abs(w).sum()
        if len(w) < 2 or total == 0:
            # One gene cannot be shuffled, and all-zero weights cannot be scaled.
            log(f"WARNING: skipping comparison {c}: {len(w)} usable genes, total |log2fc| = {total:g}")
            continue
        comp_names.append(c)
        comp_groups.append(g)
        seeds.append((index.get_indexer(sub["node_id"]), (w / total).astype(np.float32)))
    C = len(comp_names)
    if C == 0:
        raise SystemExit("No usable comparisons in the fold-change file")
    groups = sorted(set(comp_groups))
    log(f"{C} comparisons in {len(groups)} groups: " +
        ", ".join(f"{g} ({comp_groups.count(g)})" for g in groups))
    # Input genes with no edges keep their seed weight but pass nothing to other nodes.
    isolated = np.unique(np.concatenate([gi[deg[gi] == 0] for gi, _ in seeds]))
    if len(isolated):
        log(f"WARNING: {len(isolated):,} input genes have no edges in this graph and cannot "
            f"affect other nodes' scores")

    # ---- real scores ----
    real = np.zeros((n, C))
    for b in range(0, C, args.batch):
        cols = list(range(b, min(b + args.batch, C)))
        real[:, cols] = ppr(PT, seed_matrix(n, seeds, cols), args.restart, args.tol)
    # Signed ranks follow the paper. --rank-abs ranks by |score|; see README before using.
    rank_src = np.abs(real) if args.rank_abs else real
    ranks = np.column_stack([pd.Series(rank_src[:, j]).rank(method="average").values for j in range(C)])
    log("Real scores done")

    # ---- null 1: label shuffle ----
    z_label = None
    if args.n_perm:
        z_label = np.zeros((n, C))
        for j in range(C):
            gidx, w = seeds[j]
            s1, s2, done = np.zeros(n), np.zeros(n), 0
            while done < args.n_perm:
                k = min(args.batch, args.n_perm - done)
                S = np.zeros((n, k), dtype=np.float32)
                for m in range(k):
                    S[gidx, m] = rng.permutation(w)
                X = ppr(PT, S, args.restart, args.tol).astype(np.float64)
                s1 += X.sum(1)
                s2 += (X ** 2).sum(1)
                done += k
            z_label[:, j] = zfrom(real[:, j], s1, s2, args.n_perm)
            log(f"Label shuffle {j + 1}/{C} ({comp_names[j]})")

    # ---- null 2: XSwap ----
    z_xswap = None
    if args.n_xswap:
        s1, s2 = np.zeros((n, C)), np.zeros((n, C))
        for r in range(args.n_xswap):
            rs, rd, rate = xswap(src, dst, etype, n, rng, args.xswap_multiplier, args.directed)
            PTr, _ = build_transition(rs, rd, n, args.directed)
            for b in range(0, C, args.batch):
                cols = list(range(b, min(b + args.batch, C)))
                X = ppr(PTr, seed_matrix(n, seeds, cols), args.restart, args.tol).astype(np.float64)
                s1[:, cols] += X
                s2[:, cols] += X ** 2
            log(f"XSwap graph {r + 1}/{args.n_xswap} (swap acceptance {rate:.2f})")
        z_xswap = zfrom(real, s1, s2, args.n_xswap)

    # ---- summaries ----
    out = nodes.copy()
    out["degree"] = deg.astype(np.int64)
    # A seed gene's own score is mostly its own fold change, so its z-scores are partly
    # circular. This column shows how many comparisons used the node as an input gene.
    n_input = np.zeros(n, dtype=np.int64)
    for gi, _ in seeds:
        n_input[gi] += 1
    out["input_in_n_comparisons"] = n_input
    types = out["type"].values
    ga = np.array(comp_groups)

    for g in groups:
        cols = ga == g
        k = cols.sum()
        out[f"{g}__mean_rank"] = ranks[:, cols].mean(1)
        for name, Z in (("label", z_label), ("xswap", z_xswap)):
            if Z is None:
                continue
            zc = Z[:, cols].sum(1) / np.sqrt(k)  # Stouffer; assumes independent comparisons
            p = 2 * norm.sf(np.abs(zc))
            out[f"{g}__z_{name}_stouffer"] = zc
            out[f"{g}__q_{name}"] = q_within(p, types)

    welch_cols = []
    if args.control_group:
        if args.control_group not in groups:
            raise SystemExit(f"--control-group {args.control_group} not among groups {groups}")
        ctrl = ga == args.control_group
        for g in groups:
            if g == args.control_group:
                continue
            test = ga == g
            if test.sum() < 2 or ctrl.sum() < 2:
                log(f"Skipping Welch for {g}: needs 2+ comparisons in both groups")
                continue
            for name, M in (("rank", ranks), ("zlabel", z_label)):
                if M is None:
                    continue
                t, p = welch(M[:, test], M[:, ctrl])
                out[f"{g}_vs_ctrl__welch_{name}_t"] = t
                out[f"{g}_vs_ctrl__welch_{name}_q"] = q_within(p, types)
                welch_cols.append(f"{g}_vs_ctrl__welch_{name}_q")

    sort_col = next((c for c in welch_cols if "zlabel" in c), None) or \
        next((c for c in out.columns if c.endswith("__q_label")), None)
    if sort_col:
        out = out.sort_values(sort_col)
    out.to_csv(os.path.join(args.outdir, "node_summary.csv"), index=False)

    # Plain string arrays, so the file opens with np.load() without allow_pickle.
    save = {"comparisons": np.array(comp_names, dtype=str), "groups": ga.astype(str),
            "node_id": nodes["node_id"].to_numpy(dtype=str),
            "raw_score": real.astype(np.float32), "rank": ranks.astype(np.float32)}
    if z_label is not None:
        save["z_label"] = z_label.astype(np.float32)
    if z_xswap is not None:
        save["z_xswap"] = z_xswap.astype(np.float32)
    np.savez_compressed(os.path.join(args.outdir, "per_comparison.npz"), **save)

    # ---- degree-bias diagnostic ----
    seedset = np.unique(np.concatenate([s[0] for s in seeds]))
    mask = (deg > 0)
    mask[seedset] = False
    diag = {"spearman_degree_vs_abs_raw": float(np.mean(
        [spearmanr(deg[mask], np.abs(real[mask, j])).correlation for j in range(C)]))}
    for name, Z in (("z_label", z_label), ("z_xswap", z_xswap)):
        if Z is not None:
            diag[f"spearman_degree_vs_abs_{name}"] = float(np.mean(
                [spearmanr(deg[mask], np.abs(Z[mask, j])).correlation for j in range(C)]))
    info = {"args": vars(args), "comparisons": dict(zip(comp_names, comp_groups)), "diagnostics": diag}
    with open(os.path.join(args.outdir, "run_info.json"), "w") as f:
        json.dump(info, f, indent=2)
    for k, v in diag.items():
        log(f"{k} = {v:.3f}")
    log(f"Wrote results to {args.outdir}/")


if __name__ == "__main__":
    main()
