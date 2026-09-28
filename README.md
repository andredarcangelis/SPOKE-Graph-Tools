# PSEV Artifact Validation Pipeline

Checks whether high-scoring nodes in PSEV-style output on SPOKE reflect the input expression data or graph structure. Personalized PageRank favors well-connected nodes, so a node can rank high because it is a hub or heavily annotated, not because the input genes point to it. This pipeline scores every node against null models that keep that structure and remove the signal.

It runs on a single comparison or on many comparisons grouped by condition (for example, Space vs. Ground across several studies).

## Requirements

Python 3.9+ with numpy, pandas and scipy. Exporting from Neo4j also needs the Neo4j driver.

```
pip install numpy pandas scipy neo4j
```

## Files

| File | Purpose |
|---|---|
| `run_pipeline.sh` | Runs export, fold-change preparation and validation in order, with all settings at the top of the file. |
| `prepare_inputs.py` | Builds the three input files from a SPOKE Neo4j instance and GeneLab/OSDR differential expression tables. |
| `psev_artifact_validation.py` | Scores every node and runs the null models. |
| `tests/` | Four test scripts; see [Tests](#tests). |

## Running on real data

The quickest route is `run_pipeline.sh`. Edit the settings at the top, then run it from the folder that holds your config file:

```
QUICK=1 bash run_pipeline.sh     # short trial: 20 shuffles, no rewiring
bash run_pipeline.sh             # full run
```

Each run gets its own folder, `results/run_<date>_<time>/`, holding the log, a copy of the config and all outputs. Runs never write into each other's folders.

The export is reused only if it finished and was made with the same node types, filters and ID property. If the settings differ, the script stops and shows the difference. Delete `spoke_export/` to export again.

The steps it runs are described below.

### 1. Look at the graph

```
python prepare_inputs.py spoke-summary --uri bolt://HOST:7687 --user USER
```

Prints the number of nodes of each type, the property names on each type, and the number of edges of each type. The password is read from the `NEO4J_PASSWORD` environment variable, or asked for if that is not set. Use this to check three things before exporting: the property holding each node's ID (default `identifier`), whether Gene identifiers are Entrez IDs, and which node types are large enough to need filtering.

### 2. Export a subgraph by node type

```
python prepare_inputs.py spoke-export --uri bolt://HOST:7687 --user USER --outdir spoke_export
```

Keeps every node of the chosen types (default: the 12 types in Nelson et al. 2021) and every edge between kept nodes. Current SPOKE is far larger than the 2021 version, so running on the full graph is not practical. Selecting whole node types keeps the subgraph independent of the input genes. Extracting a neighborhood around the input genes does not, and tends to favor hubs.

Large types can be narrowed with `--node-filter 'Label::condition'`, written in Cypher on `n`, for example to keep only human proteins. Check the property names from step 1 first; the property in the example below is a placeholder.

```
--node-filter 'Protein::n.org_ncbi_id = 9606'
```

A node carrying more than one kept label is exported once, under the first of its labels that is kept, and its edges use the same ID.

The export writes `nodes.csv`, `edges.csv` and `export_settings.txt`. The files appear only when the whole export finishes, so an interrupted export leaves nothing that looks complete. Also record the SPOKE version; results depend on it.

### 3. Build the fold-change file

List the columns in each differential expression table:

```
python prepare_inputs.py de-columns GLDS-244_rna_seq_differential_expression.csv
```

Then write a config CSV with one row per comparison:

```
de_file,log2fc_column,comparison,group,study,flip
GLDS-244_DE.csv,Log2fc_(Space Flight)v(Ground Control),244_sg,space_vs_ground,GLDS-244,
GLDS-244_DE.csv,Log2fc_(Ground Control)v(Basal Control),244_gb,ground_vs_baseline,GLDS-244,
GLDS-288_DE.csv,Log2fc_(Ground Control)v(Space Flight),288_sg,space_vs_ground,GLDS-288,true
```

`flip` flips the sign of a comparison's fold changes so that every comparison in a group points the same way. In the last row the table's contrast is Ground over Space, so it is flipped to Space over Ground. `flip` accepts `true`/`yes`/`1` or `false`/`no`/`0`/blank; anything else stops with an error. `build-fc` prints each comparison's contrast so the orientation can be checked.

`study` groups comparisons for the direction filter below. A comparison with a blank `study` is treated as its own study.

```
python prepare_inputs.py build-fc --config comparisons.csv --nodes spoke_export/nodes.csv \
    --species mouse --homologs HOM_MouseHumanSequence.rpt
```

Mouse data needs a mouse-human homolog table. The paper used HomoloGene, which NCBI no longer updates; this pipeline reads MGI's `HOM_MouseHumanSequence.rpt` (from the MGI website) instead. A mouse gene that maps to several human genes is dropped. Several mouse genes that map to one human gene are averaged, as in the paper. Within each study and group, a gene whose fold change points in different directions across comparisons is dropped, as in the paper; `--no-direction-filter` turns this off. `--max-padj` keeps only genes below an adjusted p-value cutoff (default keeps all, as in the paper).

This writes `fc_long.csv` and `mapping_report.csv`, which counts genes at each step for each comparison. A low match rate to SPOKE usually means the ID type is wrong (see `--id-column` and `--spoke-gene-key`).

### 4. Run the validation

```
python psev_artifact_validation.py --nodes spoke_export/nodes.csv --edges spoke_export/edges.csv \
    --fc fc_long.csv --control-group ground_vs_baseline --n-perm 200 --n-xswap 0 --outdir results
```

Start with `--n-xswap 0` and few comparisons to check timing, then add XSwap. On a 300,000-node, 2.4-million-edge test graph, 40 shuffles took about 8 seconds per comparison and each rewired graph about 25 seconds, with 740 MB peak memory. Time grows roughly with the number of edges, so measure on your own export with `QUICK=1`. `--batch` sets how many walks run at once; lower it if memory runs out.

## Inputs

All CSV.

| File | Columns | Notes |
|---|---|---|
| nodes | `node_id, name, type` | Every node in the graph. |
| edges | `source, target[, edge_type]` | `edge_type` is needed for XSwap to rewire each edge type separately. |
| fold changes | `comparison, group, node_id, log2fc` | Long format, one row per gene per comparison. `node_id` is the SPOKE Gene node after homolog mapping. |

Preprocessing such as homolog mapping, averaging duplicate mappings, and dropping genes that change direction across a study's comparisons should be done before this step. A comparison with fewer than 2 usable genes, or with all fold changes equal to zero, is skipped with a warning.

## psev_artifact_validation.py options

| Option | Default | Meaning |
|---|---|---|
| `--control-group` | none | Group used as the reference for Welch's t. Omit for a single comparison. |
| `--n-perm` | 200 | Label shuffles per comparison. 0 skips. |
| `--n-xswap` | 0 | Rewired graphs. Each one is a full rerun, so use tens, not thousands. |
| `--xswap-multiplier` | 10 | Swap attempts per edge when rewiring. |
| `--restart` | 0.1 | Restart probability, as in Nelson et al. 2021. |
| `--tol` | 1e-6 | PageRank convergence tolerance. |
| `--directed` | off | Walk edges in their stored direction only. Default is undirected. |
| `--drop-edge-types` | none | Remove these edge types before running, for sensitivity analysis (for example, ontology parent-child edges). |
| `--rank-abs` | off | Rank by absolute score instead of signed score. See caveats. |
| `--batch` | 100 | Walks solved together. Lower it if memory is tight. |
| `--seed` | 0 | Random seed. |
| `--outdir` | `results` | Output folder. |

## Outputs

- `node_summary.csv`: one row per node with type, degree, `input_in_n_comparisons` (how many comparisons used the node as an input gene), per-group mean rank, Stouffer-combined z and q for each null, and Welch's t and q for each test group vs. the control.
- `per_comparison.npz`: node × comparison matrices of raw score, rank, `z_label` and `z_xswap`. Opens with `np.load("per_comparison.npz")`.
- `run_info.json`: settings, the comparison-to-group mapping, and degree-bias diagnostics (Spearman correlation between node degree and raw score, `z_label` and `z_xswap`).

## Method

**Scoring.** Personalized PageRank seeded by the input genes, each weighted by its log2 fold change (scaled so the absolute weights sum to 1). PageRank is linear in the seed vector, so one walk with signed weights equals the fold-change-weighted sum of per-gene PageRank vectors.

**Null 1: label shuffle.** For each comparison, the fold-change values are shuffled across the same input genes and the walk is rerun. The gene set and the distribution of fold changes are unchanged. Each node gets `z_label`, its real score compared with its own shuffled scores. This asks whether the specific gene–fold-change pairing matters.

**Null 2: XSwap.** The graph is rewired by degree-preserving edge swaps (Hanhijärvi et al. 2009; used for Hetionet by Himmelstein et al. 2017), separately within each edge type, so every node keeps its degree per edge type and edge endpoints keep their node types. The real seeds are rerun on each rewired graph. Each node gets `z_xswap`. This asks whether the score depends on the actual edges or would appear on any graph with the same degrees.

**Group tests.** Within each group, per-comparison z-scores are combined with Stouffer's method. With a control group named, Welch's t compares each node's values in each test group against the control group, once on raw ranks (as in Nelson et al. 2021) and once on `z_label`. q-values are Benjamini–Hochberg, computed within each node type.

## Differences from Nelson et al. 2021

- The paper z-scores and ranks each gene's PSEV before taking the fold-change-weighted sum. This pipeline sums raw PageRank vectors. Real and null runs use the same operator, so comparisons within this pipeline are consistent, but scores are not directly comparable to the paper's.
- The paper weights by −log2 FC. The sign flips z and t but not two-sided p or q values.
- The paper mapped mouse genes with HomoloGene; this pipeline uses MGI's homology table.
- The paper uses precomputed PSEVs over the full SPOKE graph at the time. Results here depend on the graph export you provide, including its version and any subsetting.

## Caveats

- **Input genes.** An input gene's own score is mostly its own fold change, so its z-scores partly restate the input. Use `input_in_n_comparisons` to separate input genes from the rest when reading Gene results.
- **Independence.** Stouffer and Welch treat comparisons as independent. Comparisons that share animals or control groups are not, so p-values from both are approximate and likely optimistic.
- **Power.** With few comparisons per group, Welch's t has few degrees of freedom and p-values stay modest even when groups separate cleanly. In the synthetic test, 4 vs. 4 comparisons gave q = 0.067 for a node with z ≈ 9 in every test comparison and ≈ 0 in controls.
- **Raw ranks.** A node close to the input genes is reached strongly in every comparison, test or control. Its signed rank swings between top and bottom with the sign of the fold changes reaching it, and its absolute rank stays high throughout. In the synthetic test, Welch on raw ranks did not detect the planted signal with either ranking (q = 0.61 signed, 0.91 absolute), while Welch on `z_label` and Stouffer z did. This is one synthetic setup and has not been checked on real data.
- **Input genes with no edges.** An input gene with no edges in the exported graph cannot affect any other node. The run reports how many there are.
- **XSwap cost.** Each rewired graph requires rewiring every edge type and rerunning every comparison.
- **Dangling mass.** Walk mass reaching a node with no edges is dropped. All runs share this, so it does not bias comparisons, but raw scores do not sum to 1.

## References

- Nelson CA et al. (2021). Knowledge network embedding of transcriptomic data from spaceflown mice uncovers signs and symptoms associated with terrestrial diseases. *Life* 11:42.
- Nelson CA, Butte AJ, Baranzini SE (2019). Integrating biomedical research and electronic health records to create knowledge-based biologically meaningful machine-readable embeddings. *Nat Commun* 10:3045.
- Himmelstein DS et al. (2017). Systematic integration of biomedical knowledge prioritizes drugs for repurposing. *eLife* 6:e26726.
- Hanhijärvi S et al. (2009). Tell me something I don't know: randomization strategies for iterative data mining. *KDD '09*.
