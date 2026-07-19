#!/usr/bin/env python3
"""
Build the static money-flow graph CSVs for a cycle (national).

Reads the cleaned DB (data/duckdb/fec_<year>.duckdb, from download.py) and writes,
into this directory (graphs/):
    <year>_vertices.csv, <year>_edges.csv

Edges:    source_id, target_id, amount, type   (amount summed over the whole cycle
          per (source, target, type) pair — see EDGE_GROUPING note).
Vertices: id, entity_name, zip, state, cmte_tp, org_tp, party, pty_aff,
          total_in, total_out

EDGE_GROUPING note: the graph wants "one edge between a pair" but also a `type`
property. A pair can transact under several FEC codes, so we group by
(source, target, TYPE) — one edge per (pair, code). For the (vast majority of)
single-code pairs this is exactly one edge per pair.

Run via process_data.sh, or directly:  python build_graph_data.py [YEAR]
"""
from __future__ import annotations

import sys
from pathlib import Path
import duckdb

HERE = Path(__file__).resolve().parent
YEAR = sys.argv[1] if len(sys.argv) > 1 else "2022"
DB_PATH = HERE.parent / "data" / "duckdb" / f"fec_{YEAR}.duckdb"
OUT_DIR = HERE                       # write the CSVs directly into graphs/


def main() -> None:
    con = duckdb.connect(str(DB_PATH), read_only=True)

    # One edge per (source, target, type): amount summed over the whole cycle.
    con.execute("""
        CREATE OR REPLACE TEMP TABLE g_edges AS
        SELECT source_id, target_id, type, SUM(amount) AS amount
        FROM transactions
        GROUP BY source_id, target_id, type
    """)

    # Vertex table + in/out money totals derived from the summed edges.
    con.execute("""
        CREATE OR REPLACE TEMP TABLE g_vertices AS
        WITH ids AS (
            SELECT source_id AS id FROM g_edges UNION SELECT target_id FROM g_edges
        ),
        outm AS (SELECT source_id id, SUM(amount) total_out FROM g_edges GROUP BY 1),
        inm  AS (SELECT target_id id, SUM(amount) total_in  FROM g_edges GROUP BY 1)
        SELECT i.id, n.entity_name, n.zip, n.state, n.cmte_tp, n.org_tp, n.party, n.pty_aff,
               COALESCE(inm.total_in, 0)   AS total_in,
               COALESCE(outm.total_out, 0) AS total_out
        FROM ids i
        LEFT JOIN nodes n ON i.id = n.id
        LEFT JOIN inm  ON i.id = inm.id
        LEFT JOIN outm ON i.id = outm.id
    """)

    e_path = OUT_DIR / f"{YEAR}_edges.csv"
    v_path = OUT_DIR / f"{YEAR}_vertices.csv"
    con.execute(f"COPY (SELECT source_id, target_id, amount, type FROM g_edges) "
                f"TO '{e_path}' (HEADER, DELIM ',')")
    con.execute(f"COPY (SELECT id, entity_name, zip, state, cmte_tp, org_tp, party, pty_aff, "
                f"total_in, total_out FROM g_vertices) TO '{v_path}' (HEADER, DELIM ',')")

    n_v = con.execute("SELECT COUNT(*) FROM g_vertices").fetchone()[0]
    n_e = con.execute("SELECT COUNT(*) FROM g_edges").fetchone()[0]
    con.close()
    print(f"[graph {YEAR}] {n_v:,} vertices, {n_e:,} edges  ->  {v_path.name}, {e_path.name}")


if __name__ == "__main__":
    main()
