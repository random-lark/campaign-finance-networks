#!/usr/bin/env python3
"""
run_null_simulations.py  --  standalone (no notebook) driver for the weight-null
simulations behind the polarization / entropy / leakage figures.

Rebuilds each cycle's absorbing graph from the CURRENT pipeline (so it picks up any
change to the underlying graphs), then runs the SAME two weight-only nulls as the
notebook -- permuted (shuffle the real edge weights) and dirichlet (redistribute
each node's out-strength ~ Dir(1..1)) -- N_PERM=300 draws each, seed 0, full-H
solve, NBINS=100 for the sR-sD density.

Writes three CSVs consumed by the notebook plot cells. All share the uniform
schema  year[, bin_center], real, unweighted, {permuted,dirichlet}_{mean,lo,hi}:

  sims_polarization.csv  (12*100 rows) : per (year, bin_center) sR-sD density
  sims_entropy.csv       (12 rows)     : per year, ent_diff (R-D)
  sims_leakage.csv       (12 rows)     : per year, leak_diff (R-D)

Usage:  python run_null_simulations.py        # ~45 min
"""
from __future__ import annotations
import sys, time
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix, eye, diags, csc_matrix
from scipy.sparse.linalg import spsolve

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "graphs"))
sys.path.insert(0, str(HERE))
from build_graph import build_graph                       # noqa: E402
from lpa import lpa_preprocess, lpa_build_matrix                 # noqa: E402

YEARS  = [2000, 2002, 2004, 2006, 2008, 2010, 2012, 2014, 2016, 2018, 2020, 2022]
N_PERM = 300
NBINS  = 100
SEED   = 0
_BINS  = np.linspace(-1, 1, NBINS + 1)
_CENTERS = 0.5 * (_BINS[:-1] + _BINS[1:])


def restrict(g):
    """Reduce to the absorbing sub-chain: nodes that can reach a DEM/REP candidate."""
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
    return dict(Y_L=Y_L, u=len(ids_U), l=len(ids_L), base=base, E=len(base), tout=tout,
                order=order, starts=starts, counts=counts,
                uu=(np.array(uur), np.array(uuc), np.array(uue)),
                ul=(np.array(ulr), np.array(ulc), np.array(ule)),
                dem=np.flatnonzero(cp == "DEM"), rep=np.flatnonzero(cp == "REP"))


def compute(s, amt):
    """(sR-sD node density over NBINS, ent_diff R-D, leak_diff R-D) for one weight vector."""
    (uur, uuc, uue) = s["uu"]; (ulr, ulc, ule) = s["ul"]; u = s["u"]; l = s["l"]
    TUU = csr_matrix((amt[uue], (uur, uuc)), shape=(u, u))
    TUL = csr_matrix((amt[ule], (ulr, ulc)), shape=(u, l))
    rs = np.asarray(TUU.sum(1)).ravel() + np.asarray(TUL.sum(1)).ravel()
    D = diags(np.where(rs > 0, 1.0 / rs, 0.0))
    H = np.asarray(spsolve((eye(u, format="csc") - (D @ TUU)).tocsc(), csc_matrix(D @ TUL)).todense())
    dem, rep = s["dem"], s["rep"]
    HD = H[:, dem]; HR = H[:, rep]; hD = HD.sum(1); hR = HR.sum(1); ld = hD > hR; lr = hR > hD
    dens = np.histogram(hR - hD, bins=_BINS)[0] / max(u, 1)
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
    dl = float(hR[ld].mean()) if ld.any() else np.nan
    rl = float(hD[lr].mean()) if lr.any() else np.nan
    return dens, (eR - eD), (rl - dl)


def null_permuted(s, rng):
    return rng.permutation(s["base"])


def null_dirichlet(s, rng):
    e = rng.standard_exponential(s["E"])[s["order"]]
    gsum = np.add.reduceat(e, s["starts"])
    amt_s = s["tout"][s["order"]] * e / np.repeat(gsum, s["counts"])
    a = np.empty(s["E"]); a[s["order"]] = amt_s
    return a


NULL_MAKERS = {"permuted": null_permuted, "dirichlet": null_dirichlet}


def main():
    rng = np.random.default_rng(SEED)
    t0 = time.time()
    dens_rows = []              # per (year, bin)
    ent_rows, leak_rows = [], []   # per year
    for year in YEARS:
        g = build_graph(year)
        g = g.connected_components(mode="weak").giant()
        g = restrict(g)
        s = precompute(g)
        d_real, e_real, l_real = compute(s, s["base"])
        d_unw, e_unw, l_unw = compute(s, np.ones(s["E"]))

        # null draws
        agg = {}   # name -> dict of density arrays (per bin)
        ent_row  = {"year": year, "real": e_real, "unweighted": e_unw}
        leak_row = {"year": year, "real": l_real, "unweighted": l_unw}
        for name, mk in NULL_MAKERS.items():
            dens = np.empty((N_PERM, NBINS)); ent = np.empty(N_PERM); leak = np.empty(N_PERM)
            for i in range(N_PERM):
                dd, ee, ll = compute(s, mk(s, rng))
                dens[i] = dd; ent[i] = ee; leak[i] = ll
            agg[name] = dict(mean=dens.mean(0),
                             lo=np.percentile(dens, 2.5, axis=0),
                             hi=np.percentile(dens, 97.5, axis=0))
            ent_row[f"{name}_mean"] = float(np.nanmean(ent))
            ent_row[f"{name}_lo"]   = float(np.nanpercentile(ent, 2.5))
            ent_row[f"{name}_hi"]   = float(np.nanpercentile(ent, 97.5))
            leak_row[f"{name}_mean"] = float(np.nanmean(leak))
            leak_row[f"{name}_lo"]   = float(np.nanpercentile(leak, 2.5))
            leak_row[f"{name}_hi"]   = float(np.nanpercentile(leak, 97.5))
        ent_rows.append(ent_row); leak_rows.append(leak_row)

        for b in range(NBINS):
            dens_rows.append({
                "year": year, "bin_center": _CENTERS[b],
                "real": d_real[b], "unweighted": d_unw[b],
                "permuted_mean": agg["permuted"]["mean"][b],
                "permuted_lo": agg["permuted"]["lo"][b],
                "permuted_hi": agg["permuted"]["hi"][b],
                "dirichlet_mean": agg["dirichlet"]["mean"][b],
                "dirichlet_lo": agg["dirichlet"]["lo"][b],
                "dirichlet_hi": agg["dirichlet"]["hi"][b],
            })
        print(f"[{time.time()-t0:6.0f}s] {year} done", flush=True)

    scalar_cols = ["year", "real", "unweighted",
                   "permuted_mean", "permuted_lo", "permuted_hi",
                   "dirichlet_mean", "dirichlet_lo", "dirichlet_hi"]
    pol_cols = ["year", "bin_center"] + scalar_cols[1:]
    pd.DataFrame(dens_rows)[pol_cols].to_csv(HERE / "sims_polarization.csv", index=False)
    pd.DataFrame(ent_rows)[scalar_cols].to_csv(HERE / "sims_entropy.csv", index=False)
    pd.DataFrame(leak_rows)[scalar_cols].to_csv(HERE / "sims_leakage.csv", index=False)
    print(f"[done] wrote sims_polarization.csv ({len(dens_rows)} rows), "
          f"sims_entropy.csv + sims_leakage.csv ({len(ent_rows)} rows each) in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
