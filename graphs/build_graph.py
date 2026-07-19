"""
Build an iGraph network from a cycle's CSVs in this directory (graphs/):
    <year>_vertices.csv, <year>_edges.csv    (produced by build_graph_data.py)

Vertex 'name' = the committee/candidate id (the join key used by the edges CSV and
by lpa.py). The human-readable name is attached as 'entity_name' so it does NOT
clobber the id. Edge weights ('amount' etc.) are kept numeric.
"""
from pathlib import Path
import igraph as ig
import pandas as pd

_DIR = Path(__file__).resolve().parent


def _build_graph(edges_csv: Path, vertices_csv: Path) -> ig.Graph:
    """Build a directed iGraph from an edges CSV + a vertices CSV."""

    edges_df = pd.read_csv(edges_csv, dtype={"source_id": str, "target_id": str})
    verts_df = pd.read_csv(vertices_csv, dtype={"id": str})

    # Edge attributes: every edge CSV column except source/target. Must be listed in
    # the tuples and in edge_attrs for TupleList to attach them.
    attr_cols = [c for c in edges_df.columns if c not in ("source_id", "target_id")]
    g = ig.Graph.TupleList(
        edges_df[["source_id", "target_id", *attr_cols]].itertuples(index=False),
        directed=True,
        edge_attrs=attr_cols,
    )

    # Attach vertex attributes in vertex order. The iGraph vertex 'name' stays the id
    # (set by TupleList); never attach a column literally called 'name'.
    ordered = verts_df.set_index("id").reindex(g.vs["name"])
    for col in ordered.columns:
        if col == "name":
            continue
        g.vs[col] = ordered[col].tolist()
    return g


def build_graph(year) -> ig.Graph:
    """Return the static iGraph for cycle `year` from <year>_vertices.csv / <year>_edges.csv."""

    return _build_graph(_DIR / f"{year}_edges.csv", _DIR / f"{year}_vertices.csv")
