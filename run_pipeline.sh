#!/usr/bin/env bash
# run_pipeline.sh  (PSEV Artifact Validation Pipeline)
#
# Runs the full workflow in order:
#   1. export SPOKE node types to CSV        (reused if a finished export with the same
#                                             settings already exists)
#   2. build fc_long.csv from GeneLab tables (prepare_inputs.py build-fc)
#   3. run the null models                   (psev_artifact_validation.py)
#
# Edit the settings below, then run:   bash run_pipeline.sh
# Quick trial run (few shuffles, no rewiring):   QUICK=1 bash run_pipeline.sh
#
# The Neo4j password is read from NEO4J_PASSWORD if set, otherwise asked for once.

set -euo pipefail   # stop at the first error, including an unset variable

# ------------------------------------------------------------------ settings
NEO4J_URI="${NEO4J_URI:-bolt://localhost:7687}"   # localhost = your SSH tunnel to the server
NEO4J_USER="${NEO4J_USER:-your_username}"         # or: export NEO4J_USER=... before running
NEO4J_DATABASE=""                       # leave empty for the server's default database

# Node types to keep (the 12 from Nelson et al. 2021). Add --node-filter lines below if a
# type is too large; check property names with: python prepare_inputs.py spoke-summary ...
LABELS="Anatomy BiologicalProcess CellularComponent Compound Disease Gene MolecularFunction Pathway PharmacologicClass Protein SideEffect Symptom"
NODE_FILTERS=()                         # e.g. ("Protein::n.org_ncbi_id = 9606")
EXPORT_DIR="spoke_export"

CONFIG="comparisons.csv"                # one row per comparison, see README step 3
SPECIES="mouse"                         # mouse or human
HOMOLOGS="HOM_MouseHumanSequence.rpt"   # only used for mouse
MAX_PADJ=""                             # e.g. 0.05; empty keeps all genes (as in the paper)

CONTROL_GROUP="ground_vs_baseline"      # empty for a single comparison
N_PERM=200
N_XSWAP=20
BATCH=50
SEED=0
RUN_NAME="run_$(date +%Y%m%d_%H%M%S)"   # results go in results/$RUN_NAME
# ------------------------------------------------------------------ end of settings

HERE="$(cd "$(dirname "$0")" && pwd)"
PY="${PYTHON:-python3}"

if [[ "${QUICK:-0}" == "1" ]]; then
    N_PERM=20; N_XSWAP=0; RUN_NAME="${RUN_NAME}_quick"
fi
OUT="results/$RUN_NAME"
i=2
while [[ -e "$OUT" ]]; do               # never write into an earlier run's folder
    OUT="results/${RUN_NAME}_$i"; i=$((i + 1))
done
RUN_NAME="$(basename "$OUT")"
mkdir -p "$OUT"
LOG="$OUT/pipeline.log"
exec > >(tee -a "$LOG") 2>&1           # show output and also save it to the log

echo "== $(date)  $RUN_NAME"
echo "   N_PERM=$N_PERM N_XSWAP=$N_XSWAP SPECIES=$SPECIES CONTROL_GROUP=${CONTROL_GROUP:-none}"

# Check inputs before anything slow starts.
[[ -f "$CONFIG" ]] || { echo "Missing $CONFIG (see README step 3)"; exit 1; }
if [[ "$SPECIES" == "mouse" && ! -f "$HOMOLOGS" ]]; then
    echo "Missing $HOMOLOGS (download HOM_MouseHumanSequence.rpt from MGI)"; exit 1
fi
cp "$CONFIG" "$OUT/"                     # keep the exact config next to the results

# ---- 1. export ----
# The same export options are used for the export and for checking an existing export.
# shellcheck disable=SC2206  # LABELS is split into separate words on purpose
EXPORT_OPTS=(--labels $LABELS)
for f in ${NODE_FILTERS[@]+"${NODE_FILTERS[@]}"}; do EXPORT_OPTS+=(--node-filter "$f"); done

# The export writes export_settings.txt last, so its presence means the export finished.
SETTINGS="$EXPORT_DIR/export_settings.txt"
if [[ -f "$SETTINGS" && -f "$EXPORT_DIR/nodes.csv" && -f "$EXPORT_DIR/edges.csv" ]]; then
    WANT="$("$PY" "$HERE/prepare_inputs.py" export-settings "${EXPORT_OPTS[@]}")"
    if [[ "$WANT" != "$(cat "$SETTINGS")" ]]; then
        echo "The export in $EXPORT_DIR was made with different settings:"
        diff <(echo "$WANT") "$SETTINGS" || true
        echo "Delete $EXPORT_DIR to export again, or change EXPORT_DIR."
        exit 1
    fi
    echo "== 1/3 Using existing export in $EXPORT_DIR (same settings)"
else
    echo "== 1/3 Exporting SPOKE to $EXPORT_DIR"
    if [[ -z "${NEO4J_PASSWORD:-}" ]]; then
        read -rsp "Neo4j password: " NEO4J_PASSWORD; echo
        export NEO4J_PASSWORD
    fi
    ARGS=(--uri "$NEO4J_URI" --user "$NEO4J_USER" --outdir "$EXPORT_DIR")
    [[ -n "$NEO4J_DATABASE" ]] && ARGS+=(--database "$NEO4J_DATABASE")
    "$PY" "$HERE/prepare_inputs.py" spoke-export "${ARGS[@]}" "${EXPORT_OPTS[@]}"
fi

# ---- 2. fold changes ----
echo "== 2/3 Building fold-change file"
ARGS=(--config "$CONFIG" --nodes "$EXPORT_DIR/nodes.csv" --species "$SPECIES"
      --out "$OUT/fc_long.csv" --report "$OUT/mapping_report.csv")
[[ "$SPECIES" == "mouse" ]] && ARGS+=(--homologs "$HOMOLOGS")
[[ -n "$MAX_PADJ" ]] && ARGS+=(--max-padj "$MAX_PADJ")
"$PY" "$HERE/prepare_inputs.py" build-fc "${ARGS[@]}"

# ---- 3. validation ----
echo "== 3/3 Running null models"
ARGS=(--nodes "$EXPORT_DIR/nodes.csv" --edges "$EXPORT_DIR/edges.csv" --fc "$OUT/fc_long.csv"
      --n-perm "$N_PERM" --n-xswap "$N_XSWAP" --batch "$BATCH" --seed "$SEED" --outdir "$OUT")
[[ -n "$CONTROL_GROUP" ]] && ARGS+=(--control-group "$CONTROL_GROUP")
"$PY" "$HERE/psev_artifact_validation.py" "${ARGS[@]}"

echo "== Done $(date). Results in $OUT"
