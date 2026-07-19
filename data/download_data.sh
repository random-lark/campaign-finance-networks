#!/usr/bin/env bash
# Stage 1: download FEC bulk data + build the cleaned, committee-only transaction DB.
#   - Downloads cm, cn, ccl, pas2, oth  (and indiv, from which ONLY the 24I/24T earmark
#     legs are kept). The first run per cycle downloads the ~5.2 GB indiv file to extract
#     those legs, then caches them to parquet so re-runs skip it. Needs network + duckdb.
#   - Writes data/duckdb/fec_<year>.duckdb.
# Usage: bash download_data.sh [YEAR]   (default 2022; omit YEAR to build all cycles 2000-2022)

set -euo pipefail
DATA_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"   # this script lives in new_pipeline/data
PYTHON="${PYTHON:-python}"
if [ "$#" -ge 1 ]; then
    ( cd "$DATA_DIR" && "$PYTHON" download.py "$1" )
    echo "Stage 1 complete -> $DATA_DIR/duckdb/fec_${1}.duckdb"
else
    ( cd "$DATA_DIR" && "$PYTHON" download.py )
    echo "Stage 1 complete -> $DATA_DIR/duckdb/ (all cycles 2000-2022)"
fi
