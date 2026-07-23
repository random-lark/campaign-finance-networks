#!/usr/bin/env python3
"""
make_appendix_tables.py -- appendix LaTeX tables of empirical p-values and effect sizes (z-scores)
for the three test statistics vs the two weight-only nulls (permuted / dirichlet), per cycle.

For each cycle it rebuilds the absorbing graph from the CURRENT pipeline, computes the real and
unweighted statistic, draws N_PERM=300 null graphs of each kind, and reports:
    p = (#{null >= real} + 1) / (N_PERM + 1)     (one-sided upper tail; hypothesis: real > null)
    z = (real - null_mean) / null_std             (effect size)

Statistics (same definitions as the notebook / run_null_simulations.py):
    polarization  pol       = sum_i rs_i |hR_i - hD_i| / sum_i rs_i
    entropy       ent_diff  = entropy_R - entropy_D
    leakage       leak_diff = mean_{REP rows} hD - mean_{DEM rows} hR

Writes latex_table_polarization.txt, latex_table_entropy.txt, latex_table_leakage.txt.
Usage:  python make_appendix_tables.py      (~45 min)
"""
from __future__ import annotations
import sys, time
from pathlib import Path
import numpy as np
from scipy.sparse import csr_matrix, eye, diags, csc_matrix
from scipy.sparse.linalg import spsolve

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "graphs"))
sys.path.insert(0, str(HERE))
from build_graph import build_graph          # noqa: E402
from lpa import lpa_preprocess                # noqa: E402

YEARS  = [2000, 2002, 2004, 2006, 2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022]
N_PERM = 300
SEED   = 0
# statistic key -> (index in stats() tuple, table filename stem)
TABLES = [("pol", 0, "polarization"), ("ent_diff", 1, "entropy"), ("leak_diff", 2, "leakage")]


def restrict(g):
    cand = [v.index for v in g.vs if v["party"] in ("DEM", "REP")]
    gr = g.copy(); gr.reverse_edges()
    s = gr.add_vertex().index
    gr.add_edges([(s, c) for c in cand])
    keep = sorted(set(gr.subcomponent(s, mode="out")) - {s})
    return g.induced_subgraph(keep)


def precompute(g):
    Y_L, ids_L, ids_U = lpa_preprocess(g)
    id2U = {n: i for i, n in enumerate(ids_U)}
    id2L = {n: i for i, n in enumerate(ids_L)}
    names = g.vs["name"]
    uur = []; uuc = []; uue = []; ulr = []; ulc = []; ule = []
    for e in g.es:
        i = id2U.get(names[e.source])
        if i is None:
            continue
        j = id2U.get(names[e.target])
        if j is not None:
            uur.append(i); uuc.append(j); uue.append(e.index)
        else:
            k = id2L.get(names[e.target])
            if k is not None:
                ulr.append(i); ulc.append(k); ule.append(e.index)
    n2i = {v["name"]: v.index for v in g.vs}
    cp = np.array([g.vs[n2i[c]]["party"] for c in ids_L])
    base = np.asarray(g.es["amount"], float)
    src = np.asarray([e.source for e in g.es])
    tout = np.asarray(g.vs["total_out"], float)[src]
    order = np.argsort(src, kind="stable"); ss = src[order]
    starts = np.concatenate(([0], np.flatnonzero(np.diff(ss)) + 1))
    counts = np.diff(np.append(starts, len(src)))
    return dict(u=len(ids_U), l=len(ids_L), base=base, E=len(base), tout=tout,
                order=order, starts=starts, counts=counts,
                uu=(np.array(uur), np.array(uuc), np.array(uue)),
                ul=(np.array(ulr), np.array(ulc), np.array(ule)),
                dem=np.flatnonzero(cp == "DEM"), rep=np.flatnonzero(cp == "REP"))


def stats(s, amt):
    """(polarization, entropy_diff, leakage_diff) for one weight vector -- one full-H solve."""
    (uur, uuc, uue) = s["uu"]; (ulr, ulc, ule) = s["ul"]; u = s["u"]; l = s["l"]
    TUU = csr_matrix((amt[uue], (uur, uuc)), shape=(u, u))
    TUL = csr_matrix((amt[ule], (ulr, ulc)), shape=(u, l))
    rs = np.asarray(TUU.sum(1)).ravel() + np.asarray(TUL.sum(1)).ravel()
    D = diags(np.where(rs > 0, 1.0 / rs, 0.0))
    H = np.asarray(spsolve((eye(u, format="csc") - (D @ TUU)).tocsc(), csc_matrix(D @ TUL)).todense())
    dem, rep = s["dem"], s["rep"]
    HD = H[:, dem]; HR = H[:, rep]; hD = HD.sum(1); hR = HR.sum(1)
    ld = hD > hR; lr = hR > hD
    W = rs.sum(); pol = float((np.abs(hR - hD) * rs).sum() / W)

    def pe(B, mass, mask, nc):
        if nc < 2:
            return np.nan
        with np.errstate(divide="ignore", invalid="ignore"):
            blogb = np.where(B > 0, B * np.log(B), 0.0).sum(1)
        good = mass > 0; ent = np.zeros(len(mass))
        ent[good] = (np.log(mass[good]) - blogb[good] / mass[good]) / np.log(nc)
        m = mask & good; den = rs[mask].sum()
        return float((ent[m] * rs[m]).sum() / den) if den > 0 else np.nan

    eD = pe(HD, hD, ld, len(dem)); eR = pe(HR, hR, lr, len(rep))
    dl = float(hR[ld].mean()) if ld.any() else np.nan   # DEM -> REP leak
    rl = float(hD[lr].mean()) if lr.any() else np.nan   # REP -> DEM leak
    return (pol, eR - eD, rl - dl)


def null_permuted(s, rng):
    return rng.permutation(s["base"])


def null_dirichlet(s, rng):
    e = rng.standard_exponential(s["E"])[s["order"]]
    gsum = np.add.reduceat(e, s["starts"])
    amt_s = s["tout"][s["order"]] * e / np.repeat(gsum, s["counts"])
    a = np.empty(s["E"]); a[s["order"]] = amt_s
    return a


NULLS = (("perm", null_permuted), ("rand", null_dirichlet))


def write_table(rows, path):
    L = [r"\begin{center}",
         r"\begin{tabular}{|c | c | c | c c | c c|} ",
         r" \hline",
         r" Cycle & Actual & Unweighted & $p_{perm}$ & $z_{perm}$ & $p_{rand}$ & $z_{rand}$ \\ [0.5ex] ",
         r" \hline\hline"]
    for r in rows:
        L.append(f" ${r['year']}$ & ${r['actual']:.3f}$ & ${r['unw']:.3f}$ & "
                 f"${r['p_perm']:.3f}$ & ${r['z_perm']:.3f}$ & "
                 f"${r['p_rand']:.3f}$ & ${r['z_rand']:.3f}$ " + r"\\")
        L.append(r" \hline")
    L += [r"\end{tabular}", r"\end{center}"]
    path.write_text("\n".join(L) + "\n")
    print(f"[wrote] {path.name}")


def main():
    rng = np.random.default_rng(SEED)
    t0 = time.time()
    table_rows = {stem: [] for _, _, stem in TABLES}
    for year in YEARS:
        g = restrict(build_graph(year))
        s = precompute(g)
        real = stats(s, s["base"])
        unw = stats(s, np.ones(s["E"]))
        draws = {}
        for nkey, mk in NULLS:
            arr = np.empty((N_PERM, 3))
            for d in range(N_PERM):
                arr[d] = stats(s, mk(s, rng))
            draws[nkey] = arr
        for key, idx, stem in TABLES:
            row = {"year": year, "actual": real[idx], "unw": unw[idx]}
            for nkey, _ in NULLS:
                col = draws[nkey][:, idx]; col = col[~np.isnan(col)]
                mean = col.mean(); std = col.std()
                row[f"p_{nkey}"] = (np.sum(col >= real[idx]) + 1) / (len(col) + 1)
                row[f"z_{nkey}"] = (real[idx] - mean) / std if std > 0 else np.nan
            table_rows[stem].append(row)
        print(f"[{time.time()-t0:6.0f}s] {year} done", flush=True)
    for _, _, stem in TABLES:
        write_table(table_rows[stem], HERE / f"latex_table_{stem}.txt")
    print(f"[done] 3 tables in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
