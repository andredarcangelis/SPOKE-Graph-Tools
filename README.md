# SPOKE Topology

Tools for exporting local neighborhoods from the [SPOKE](https://spoke.ucsf.edu/) biomedical knowledge graph and checking their structure: degree distribution, hub nodes, how connectivity depends on those hubs, and whether PageRank computed in Python matches Neo4j's reference implementation.

Work done in Dr. Sylvain Costes' lab at the Trivedi Institute of Space and Global Biomedicine (University of Pittsburgh School of Medicine)

## Background
 
SPOKE links genes, proteins, compounds, diseases, pathways and other biomedical concepts from dozens of public databases into one graph. [Nelson et al. (2021)](https://doi.org/10.3390/life11010042) ranked SPOKE nodes against spaceflight gene expression data using personalized PageRank, reading top-ranked disease and symptom nodes as phenotypes linked to spaceflight.
 
In earlier work applying this approach to two NASA OSDR datasets (OSD-13 and OSD-91), several top-ranked nodes turned out to reflect graph structure rather than the input data. The clearest case was HLA-A, which has 9,273 edges in SPOKE (compared to roughly 17 to 500 for randomly sampled genes) and ranked near the top even when fold-change values were randomly shuffled across genes. Other artifacts came from annotation density and from disease nodes that ranked high through ontology links alone.
 
That motivated looking at SPOKE's structure directly. This repo does that on a small, well-understood test case.
 
## What's here
 
| File | Purpose |
|---|---|
| `export_subgraph.py` | Exports the 2-hop neighborhood of a SPOKE Pathway node to an edge-list CSV |
| `topology.py` | Structural summary of an export, plus a hub-removal test |
| `validate_pagerank.py` | Compares networkx PageRank to Neo4j GDS PageRank on the same graph, standard or personalized |
| `graph_io.py` | Shared Neo4j connection and CSV loading |
 
All three scripts use the same edge-list CSV format: one row per SPOKE relationship, with nodes identified by their Neo4j element ID. Nodes are never keyed by name, since different SPOKE nodes can share a name.
 
## Test case: aerobic glycolysis
 
The test case is the WikiPathways Aerobic Glycolysis pathway (`WP4629_r127017`), which connects to 12 genes in SPOKE. It was chosen because it is small and has nothing to do with spaceflight, so structure can be looked at on its own. Each of the 12 genes has between 283 (PGAM2) and 2,314 (TPI1) connections in the full graph.
 
### Finding 1: an unordered LIMIT silently dropped 10 of 12 genes
 
The first 2-hop export used `LIMIT 1000` with no `ORDER BY`. Only 2 of the 12 genes (GPI and ENO1) appeared in it. This was a query artifact, not a property of the graph:
 
- All 12 genes have hundreds to thousands of connections in SPOKE.
- GPI alone had 626 eligible rows in that query, and 494 appeared in the export, so the limit cut off partway through GPI before any of the next 10 genes were reached.
- Rerunning the same query returned the same two genes, so the bias is deterministic, not random.
Hub-removal results on this export were misleading. Removing GPI and ENO1 left 749 isolated nodes and 0.1% of the graph connected, which looks like extreme hub dependence but only reflects the export having two genes in it.
 
With every gene capped at the same number of neighbors instead (50 per gene), removing the two highest-degree genes (TPI1 and ALDOA) left 91.3% of the remaining graph in one connected component.
 
`export_subgraph.py` reports how many hop-2 relationships each hop-1 node received after every export, so this kind of imbalance shows up immediately. `--naive-limit` reproduces the original biased query for comparison.
 
### Finding 2: networkx PageRank matches GDS
 
On the glycolysis subgraph, standard PageRank in networkx and in Neo4j's Graph Data Science library returned the same top 10 nodes (10/10 overlap). This confirms that PageRank computed locally in Python is mechanically equivalent to the reference implementation.
 
This does not yet cover personalized PageRank, the variant PSEVs are built on. `validate_pagerank.py --seed` runs that comparison.
 
## Setup
 
Requires Python 3.10 or later.
 
```bash
conda create -n spoke_env python=3.11
conda activate spoke_env
pip install -r requirements.txt
```
 
The first line creates an isolated Python environment named `spoke_env`, the second switches your terminal into it, and the third installs the packages listed in `requirements.txt`. Run `conda activate spoke_env` again in each new terminal.
 
### Neo4j access
 
You need two Neo4j instances:
 
- **A SPOKE instance** to export from. If it's on a remote server, forward its port to your machine first:
```bash
  ssh -L 7687:localhost:7687 <user>@<server>
```
 
  This connects to the server and makes its Neo4j port available at `localhost:7687` on your machine. Leave that terminal open.
 
- **A local Neo4j instance with the GDS plugin** for `validate_pagerank.py`, such as a database in Neo4j Desktop. This script deletes and recreates its own nodes on every run, so point it at a disposable database, never the shared SPOKE server.
Connection details are read from environment variables:
 
```bash
export NEO4J_URI=bolt://localhost:7687
export NEO4J_USER=<username>
export NEO4J_PASSWORD=<password>
```
 
`export` sets a variable for the current terminal session only. If you leave out `NEO4J_PASSWORD`, the scripts will prompt for it instead, which keeps it out of your shell history. `--uri`, `--user` and `--database` override these on any script.
 
## Usage
 
### 1. Export a neighborhood
 
```bash
python export_subgraph.py WP4629_r127017 -o glycolysis_full.csv
python export_subgraph.py WP4629_r127017 -o glycolysis_50.csv --per-node-limit 50
python export_subgraph.py WP4629_r127017 -o glycolysis_naive.csv --naive-limit 1000
```
 
The first exports the pathway's full 2-hop neighborhood. The second keeps at most 50 second-hop relationships per gene, chosen at random with a fixed seed (`--seed`, default 0), so the sample is unbiased and the same seed reproduces the same export every time. The third reproduces the biased export from Finding 1.
 
### 2. Analyze structure
 
```bash
python topology.py glycolysis_full.csv --plot glycolysis_full.png
python topology.py glycolysis_full.csv --remove GPI ENO1
```
 
The first prints the structural summary, removes the two highest-degree nodes, and saves a network plot. The second removes specific nodes by name instead. Names that match more than one node are rejected rather than guessed.
 
### 3. Validate PageRank
 
```bash
python validate_pagerank.py glycolysis_full.csv --database <local-db>
python validate_pagerank.py glycolysis_full.csv --database <local-db> --seed TPI1 ALDOA GPI
```
 
The first compares standard PageRank. The second compares personalized PageRank restarting at the named genes. Both print the two top-10 lists side by side, their overlap, the Spearman rank correlation across all nodes, and the largest difference in normalized scores.
 
## Limitations
 
- One pathway has been analyzed so far.
- Degree in `topology.py` is measured inside the exported subgraph, not across all of SPOKE.
- Multiple relationships between the same two nodes are collapsed to one edge for the graph analysis. Relationship-type counts still count each relationship.
- `--per-node-limit` samples uniformly at random per node, using a fixed seed (`--seed`, default 0) so results are reproducible. A different seed will keep a different set of neighbors, so any single export is one random draw, not the full picture.
## References
 
Nelson, C.A., Acuna, A.U., Paul, A.M., Scott, R.T., Butte, A.J., Cekanaviciute, E., Baranzini, S.E., Costes, S.V. (2021). Knowledge Network Embedding of Transcriptomic Data from Spaceflown Mice Uncovers Signs and Symptoms Associated with Terrestrial Diseases. *Life*, 11(1), 42. https://doi.org/10.3390/life11010042
 

