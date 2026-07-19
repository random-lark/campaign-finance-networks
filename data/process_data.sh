#!/usr/bin/env bash
# Stage 2: build the national static graph CSVs from the cleaned DB.
# Outputs <year>_vertices.csv + <year>_edges.csv into graphs/.
# Usage: bash process_data.sh [YEAR]   (default 2022)   [this script lives in data/]

set -euo pipefail
DATA_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   
GRAPHS_DIR="$(cd "$DATA_DIR/../graphs" && pwd)"
PYTHON="${PYTHON:-python}"
YEAR="${1:-2022}"
( cd "$GRAPHS_DIR" && "$PYTHON" build_graph_data.py "$YEAR" )
echo "Stage 2 complete -> $GRAPHS_DIR/${YEAR}_vertices.csv, ${YEAR}_edges.csv"
