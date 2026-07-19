#!/usr/bin/env python3
"""
Download FEC bulk data and build a cleaned, committee-only transaction database, for
every cycle 2000-2022. One DuckDB per cycle is written to  data/duckdb/fec_<year>.duckdb.

What it does (per cycle)
------------------------
1. Downloads the FEC bulk files. Excludes `oppexp` and most of `indiv`. Each becomes a
   DuckDB table using the FEC's column names.
2. Cleans the transaction tables (pas2, oth) loosely following Helyn
   (https://helyn.com/follow-the-money), plus project-specific exclusions in DROP_CODES. Notes:
       - Dedup on sub_id (FEC's unique key)
       - Amendment handling via file_num, keeping the latest filing per report period
       - Drop memo rows (memo_cd = 'X') from money flow
       - Earmarks: drop 15E and use 24I/24T instead, which come from `indiv.` They are the
         second conduit->recipient legs, so including them does NOT add individuals.
         Always included (the 24I/24T legs are cached to parquet after the one-time indiv
         download, so re-runs are cheap).
       - Drop self-loops. Every committee is its own node.

Output
------
data/duckdb/fec_<year>.duckdb  with:
    raw tables:   cm, cn, ccl, pas2, oth, indiv_earmarks
    transactions: one row per cleaned transaction (source_id, target_id, amount, type, date)
    nodes:        id, entity_name, zip, state, cmte_tp, org_tp, party, pty_aff
                  (committee / candidate)

Run via download_data.sh, or directly:
    python download.py            # all cycles 2000-2022
    python download.py 2016       # a single cycle
"""
from __future__ import annotations

import sys
import shutil
import zipfile
import tempfile
import urllib.request
from pathlib import Path

import duckdb

# ----------------------------------------------------------------------------- Config
HERE = Path(__file__).resolve().parent
DUCKDB_DIR = HERE / "duckdb"                    # One fec_<year>.duckdb per cycle lives here
YEARS = list(range(2000, 2023, 2))              # 2000, 2002, ..., 2022
KEEP_RAW_FILES = False                          # Delete the multi-GB unzipped indiv file after extraction
BASE_URL = "https://www.fec.gov/files/bulk-downloads"

# Restrict the whole pipeline to the 50 states + DC.
# Entities not from these are dropped, along with edges involving them. 
ALLOWED_STATES = {
    "AL","AK","AZ","AR","CA","CO","CT","DE","FL","GA","HI","ID","IL","IN","IA","KS","KY",
    "LA","ME","MD","MA","MI","MN","MS","MO","MT","NE","NV","NH","NJ","NM","NY","NC","ND",
    "OH","OK","OR","PA","RI","SC","SD","TN","TX","UT","VT","VA","WA","WV","WI","WY","DC"}

# FEC bulk file layouts (the files are headerless; these are the official column names).
HEADERS = {
    "cm":  ["CMTE_ID","CMTE_NM","TRES_NM","CMTE_ST1","CMTE_ST2","CMTE_CITY","CMTE_ST",
            "CMTE_ZIP","CMTE_DSGN","CMTE_TP","CMTE_PTY_AFFILIATION","CMTE_FILING_FREQ",
            "ORG_TP","CONNECTED_ORG_NM","CAND_ID"],
    "cn":  ["CAND_ID","CAND_NAME","CAND_PTY_AFFILIATION","CAND_ELECTION_YR","CAND_OFFICE_ST",
            "CAND_OFFICE","CAND_OFFICE_DISTRICT","CAND_ICI","CAND_STATUS","CAND_PCC",
            "CAND_ST1","CAND_ST2","CAND_CITY","CAND_ST","CAND_ZIP"],
    "ccl": ["CAND_ID","CAND_ELECTION_YR","FEC_ELECTION_YR","CMTE_ID","CMTE_TP","CMTE_DSGN","LINKAGE_ID"],
    "pas2":["CMTE_ID","AMNDT_IND","RPT_TP","TRANSACTION_PGI","IMAGE_NUM","TRANSACTION_TP",
            "ENTITY_TP","NAME","CITY","STATE","ZIP_CODE","EMPLOYER","OCCUPATION","TRANSACTION_DT",
            "TRANSACTION_AMT","OTHER_ID","CAND_ID","TRAN_ID","FILE_NUM","MEMO_CD","MEMO_TEXT","SUB_ID"],
    "oth": ["CMTE_ID","AMNDT_IND","RPT_TP","TRANSACTION_PGI","IMAGE_NUM","TRANSACTION_TP",
            "ENTITY_TP","NAME","CITY","STATE","ZIP_CODE","EMPLOYER","OCCUPATION","TRANSACTION_DT",
            "TRANSACTION_AMT","OTHER_ID","TRAN_ID","FILE_NUM","MEMO_CD","MEMO_TEXT","SUB_ID"],
    "indiv":["CMTE_ID","AMNDT_IND","RPT_TP","TRANSACTION_PGI","IMAGE_NUM","TRANSACTION_TP",
            "ENTITY_TP","NAME","CITY","STATE","ZIP_CODE","EMPLOYER","OCCUPATION","TRANSACTION_DT",
            "TRANSACTION_AMT","OTHER_ID","TRAN_ID","FILE_NUM","MEMO_CD","MEMO_TEXT","SUB_ID"],
}

# Zip stem (a 2-digit cycle suffix is appended) and the member filename inside the zip
ZIP_INFO = {"cm":("cm","cm.txt"), "cn":("cn","cn.txt"), "ccl":("ccl","ccl.txt"),
            "pas2":("pas2","itpas2.txt"), "oth":("oth","itoth.txt"), "indiv":("indiv","itcont.txt")}

# This pipeline's excluded edges
DROP_CODES = set()
DROP_CODES |= {"16C","16F","16G","16H","16J","16K","16L","16R","16U",          # Loans (not an explicit donation, and trivial portion of data)
               "20C","20F","20G","20R","22G","22H","22J","22K","22L","22R","22U","22X"}
DROP_CODES |= {"24A","24E"}                                                    # Independent expenditures (again, not an explicit donation)
DROP_CODES |= {"24F","24N","29"}                                               # Communication costs (for/against) + electioneering
                                                                               # Communications money goes to media/vendors, not the candidate;
                                                                               # and 'against' edges point to the enemy. Same category as the 24A/24E IEs above.
                                                                               # But coordinated party expenditures (24C) are kept, as they represent
                                                                               # money from a party to its own candidate.
DROP_CODES |= {"17R","17U","17Y","17Z","20Y","21Y","22Y","22Z","23Y","28L",    # Refunds (reverse edge direction, and not an explicit donation)
               "40Y","40T","40Z","41Y","41T","41Z","42Y","42T","42Z"}
DROP_CODES |= {"15Z","24Z"}                                                    # In-kind (not cash-valued)
DROP_CODES |= {"15E"}                                                          # Earmark recipient-side (I use 24I/24T instead)
DROP_CODES |= {"18G","18H","18J","18K","18L","18U"}                            # Receipt-side codes, filed by recipient. The
                                                                               # pipeline records source=CMTE_ID, so these point
                                                                               # backward and double-count the giver-side 24K/24G;
                                                                               # drop them and keep the correctly-directed giver side.
DROP_CODES |= {"15J","19J","31F","32F","32J"}                                  # JFC allocation MEMO codes ("Memo - recipient committee's
                                                                               # percentage of a contribution given to a joint fundraising
                                                                               # committee"). Properly-flagged ones are dropped by MEMO_CD='X';
                                                                               # a few filers leave MEMO_CD blank, so drop by code too. Verified
                                                                               # to double-count: the real transfer is already captured by a
                                                                               # 24K/24G/24T edge (and these memos point the wrong way).
DROP_LIST = "(" + ",".join(f"'{c}'" for c in sorted(DROP_CODES)) + ")"


# --------------------------------------------------------------------------- Download

def _download_main_txt(family: str, out_path: Path, year: int) -> None:
    """Download <family><yy>.zip and extract only its top-level main .txt (skip indiv/by_date/)."""
    yy = str(year)[-2:]
    stem, member = ZIP_INFO[family]
    url = f"{BASE_URL}/{year}/{stem}{yy}.zip"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    print(f"  downloading {url}")
    with tempfile.NamedTemporaryFile(suffix=".zip", delete=True) as tmp:
        with urllib.request.urlopen(url) as resp:
            shutil.copyfileobj(resp, tmp, length=1 << 20)
        tmp.flush()
        with zipfile.ZipFile(tmp.name) as zf:
            # top-level .txt only -> the main file, never the duplicate by_date/ chunks
            members = [m for m in zf.namelist() if m.lower().endswith(".txt") and "/" not in m]
            pick = member if member in members else (members[0] if members else None)
            if not pick:
                raise RuntimeError(f"No top-level .txt in {url}")
            with zf.open(pick) as src, out_path.open("wb") as dst:
                shutil.copyfileobj(src, dst, length=1 << 20)
    print(f"    -> {out_path} ({out_path.stat().st_size/1e6:.1f} MB)")


def _names_array(family: str) -> str:
    return "[" + ",".join(f"'{c}'" for c in HEADERS[family]) + "]"


def _read_raw(con, family: str, path: Path):
    """Create a raw DuckDB table from a headerless pipe file using FEC column names."""
    names = _names_array(family)
    con.execute(f"""
        CREATE OR REPLACE TABLE {family} AS
        SELECT * FROM read_csv('{path}', delim='|', header=false, names={names},
                               all_varchar=true, ignore_errors=true)
    """)
    n = con.execute(f"SELECT COUNT(*) FROM {family}").fetchone()[0]
    print(f"    table {family}: {n:,} rows")


def _clean_select(tbl: str, target_expr: str) -> str:
    # Amendment handling: within (CMTE_ID, RPT_TP, report-year) keep only the rows from the
    # highest FILE_NUM (the latest filing). report-year is the txn date's year (no RPT_YR
    # column exists). Helyn's "max file_num per period" heuristic.
    return f"""
        WITH base AS (
            SELECT *, TRY_STRPTIME(TRANSACTION_DT,'%m%d%Y') AS dt,
                      TRY_CAST(TRANSACTION_AMT AS DOUBLE)   AS amt,
                      TRY_CAST(FILE_NUM AS BIGINT)          AS fnum
            FROM {tbl}
            WHERE COALESCE(MEMO_CD,'') <> 'X'                       -- drop memo sub-items
              AND TRANSACTION_TP NOT IN {DROP_LIST}                 -- project exclusions
        ),
        latest AS (
            SELECT *, MAX(fnum) OVER (PARTITION BY CMTE_ID, RPT_TP, EXTRACT(year FROM dt)) AS max_fnum
            FROM base
        )
        SELECT CMTE_ID AS source_id,
               {target_expr} AS target_id,
               amt AS amount, TRANSACTION_TP AS type, dt AS date, SUB_ID AS sub_id
        FROM latest
        WHERE (fnum = max_fnum OR fnum IS NULL)
    """


# --------------------------------------------------------------------------- per-cycle build
def build_cycle(year: int) -> None:
    raw_dir = HERE / "raw" / str(year)
    db_path = DUCKDB_DIR / f"fec_{year}.duckdb"
    raw_dir.mkdir(parents=True, exist_ok=True)
    DUCKDB_DIR.mkdir(parents=True, exist_ok=True)
    if db_path.exists():
        db_path.unlink()
    con = duckdb.connect(str(db_path))

    # 1) Download + raw tables -----------------------------------------------------------
    print(f"[1] downloading + loading raw tables for cycle {year}")
    for fam in ["cm", "cn", "ccl", "pas2", "oth"]:
        p = raw_dir / ZIP_INFO[fam][1]
        if not p.exists():
            _download_main_txt(fam, p, year)
        _read_raw(con, fam, p)

    # 2) Master dedup — keep one row per natural key ---------
    con.execute("CREATE OR REPLACE TABLE cm  AS SELECT * FROM cm  QUALIFY row_number() OVER (PARTITION BY CMTE_ID ORDER BY CMTE_ID)=1")
    con.execute("CREATE OR REPLACE TABLE cn  AS SELECT * FROM cn  QUALIFY row_number() OVER (PARTITION BY CAND_ID ORDER BY CAND_ID)=1")

    # 3) Clean transaction tables (Helyn: sub_id dedup, amendments via file_num, memo drop)
    print("[2] cleaning transactions (amendments, memos, exclusions)")
    # oth: committee -> committee (OTHER_ID).  pas2: committee -> recipient committee
    # (OTHER_ID) else the candidate (CAND_ID).
    pas2_target = "COALESCE(NULLIF(OTHER_ID,''), CAND_ID)"
    con.execute(f"CREATE OR REPLACE TABLE clean_oth  AS {_clean_select('oth',  'OTHER_ID')}")
    con.execute(f"CREATE OR REPLACE TABLE clean_pas2 AS {_clean_select('pas2', pas2_target)}")

    # 3b) Earmarks from indiv (24I/24T only)----------------------------------------------
    # The 24I/24T legs are cached to parquet after the (one-time, ~5 GB) indiv download, so
    # re-runs load the cache and skip the download entirely.
    cache = raw_dir / "indiv_earmarks.parquet"
    itcont = raw_dir / "itcont.txt"
    names = _names_array("indiv")

    def _extract_from_itcont():
        con.execute(f"""
            CREATE OR REPLACE TABLE indiv_earmarks AS
            SELECT CMTE_ID AS source_id, OTHER_ID AS target_id,
                   TRY_CAST(TRANSACTION_AMT AS DOUBLE) AS amount,
                   TRANSACTION_TP AS type, TRY_STRPTIME(TRANSACTION_DT,'%m%d%Y') AS date,
                   SUB_ID AS sub_id
            FROM read_csv('{itcont}', delim='|', header=false, names={names},
                          all_varchar=true, ignore_errors=true)
            WHERE TRANSACTION_TP IN ('24I','24T')
              AND COALESCE(OTHER_ID,'') <> ''        -- need a committee recipient
        """)
        con.execute(f"COPY indiv_earmarks TO '{cache}' (FORMAT parquet)")   # Cache for cheap re-runs
        if not KEEP_RAW_FILES:
            itcont.unlink(missing_ok=True)           # Reclaim ~15-20GB

    if itcont.exists():
        print("[2b] extracting 24I/24T earmark legs from indiv (+ caching)")
        _extract_from_itcont()
    elif cache.exists():
        print("[2b] loading cached 24I/24T earmark legs (skipping ~5 GB indiv download)")
        con.execute(f"CREATE OR REPLACE TABLE indiv_earmarks AS SELECT * FROM read_parquet('{cache}')")
    else:
        print("[2b] downloading indiv to extract 24I/24T earmark legs (~5 GB, one-time)")
        _download_main_txt("indiv", itcont, year)
        _extract_from_itcont()
    ne = con.execute("SELECT COUNT(*) FROM indiv_earmarks").fetchone()[0]
    print(f"    indiv_earmarks (24I/24T, committee->committee): {ne:,} rows")

    # 4) Unify -> cleaned transactions, keep only committee/candidate IDs -----------------
    con.execute("""
        CREATE OR REPLACE TABLE tx_all AS
        SELECT * FROM clean_oth
        UNION ALL SELECT * FROM clean_pas2
        UNION ALL SELECT * FROM indiv_earmarks
    """)
    # Helyn's dedup on sub_id requires real FEC ID endpoints; keep only positive amounts
    # (amount > 0 = GROSS INFLOW: FEC negatives -- corrections / reattributions /
    # redesignations, ~4% of rows -- and zeros are dropped, since a flow graph for label
    # propagation needs positive edge weights), and drop self-loops.
    con.execute("""
        CREATE OR REPLACE TABLE transactions AS
        SELECT source_id, target_id, amount, type, date
        FROM (
            SELECT * FROM tx_all
            WHERE source_id IS NOT NULL AND source_id <> ''
              AND target_id IS NOT NULL AND target_id <> ''
              AND amount IS NOT NULL AND amount > 0
            QUALIFY row_number() OVER (PARTITION BY sub_id ORDER BY sub_id) = 1
        )
        WHERE source_id <> target_id          -- drop self-loops
    """)
    print(f"    cleaned transactions: {con.execute('SELECT COUNT(*) FROM transactions').fetchone()[0]:,}")

    # 5) Nodes metadata (committee / candidate) ------------------------------------------
    print("[3] building node metadata table")
    con.execute("""
        CREATE OR REPLACE TABLE node_ids AS
        SELECT DISTINCT source_id AS id FROM transactions
        UNION SELECT DISTINCT target_id FROM transactions
    """)
    con.execute("""
        CREATE OR REPLACE TABLE nodes AS
        SELECT
            n.id,
            COALESCE(cm.CMTE_NM, cn.CAND_NAME, n.id)                                  AS entity_name,
            cm.CMTE_ZIP                                                               AS zip,
            COALESCE(cm.CMTE_ST, cn.CAND_OFFICE_ST,
                     -- fallback for candidate nodes missing from cn: parse the 2-letter
                     -- state out of the candidate id. H/S ids are [office][year][ST][...],
                     -- e.g. H2CA47196 -> CA, S2WI00474 -> WI. (P/president ids have no state.)
                     CASE WHEN substr(n.id, 1, 1) IN ('H', 'S')
                          THEN substr(n.id, 3, 2) END)                               AS state,
            cm.CMTE_TP                                                                AS cmte_tp,
            cm.ORG_TP                                                                 AS org_tp,
            -- 'party' is CANDIDATE-ONLY (from the candidate master); committees get NULL here.
            -- This makes a non-null 'party' a clean marker of candidate (absorbing) nodes for
            -- the absorbing-Markov-chain analysis. Committee affiliation lives in 'pty_aff'.
            cn.CAND_PTY_AFFILIATION                                                   AS party,
            cm.CMTE_PTY_AFFILIATION                                                   AS pty_aff
        FROM node_ids n
        LEFT JOIN cm ON n.id = cm.CMTE_ID
        LEFT JOIN cn ON n.id = cn.CAND_ID
    """)
    print(f"    nodes: {con.execute('SELECT COUNT(*) FROM nodes').fetchone()[0]:,}")

    # 6) RESTRICT to 50 states + DC, and drop unregistered committees ---------------------
    # Exclude (a) nodes whose state is a known non-allowed place (territories / 'US' / foreign),
    # and (b) committee ids (C...) absent from THIS cycle's committee master cm -- i.e. not
    # registered/active this cycle (defunct/ghost committees that only appear as a recipient
    # because another filer named them). Candidate nodes missing from cn are KEPT (their state
    # is parsed from the id above). Drop edges touching excluded ids, then trim to referenced.
    print("[4] restricting to 50 states + DC; dropping committees not in this cycle's cm")
    allowed_sql = "(" + ",".join(f"'{s}'" for s in sorted(ALLOWED_STATES)) + ")"
    con.execute(f"""
        CREATE OR REPLACE TABLE excluded_ids AS
        SELECT id FROM nodes
        WHERE (state IS NOT NULL AND state <> '' AND state NOT IN {allowed_sql})   -- territories / 'US' / foreign
           OR (id LIKE 'C%' AND id NOT IN (SELECT CMTE_ID FROM cm))                -- committee not registered this cycle
    """)
    n_excl = con.execute("SELECT COUNT(*) FROM excluded_ids").fetchone()[0]
    con.execute("""
        CREATE OR REPLACE TABLE transactions AS
        SELECT * FROM transactions
        WHERE source_id NOT IN (SELECT id FROM excluded_ids)
          AND target_id NOT IN (SELECT id FROM excluded_ids)
    """)
    con.execute("""
        CREATE OR REPLACE TABLE nodes AS
        SELECT * FROM nodes
        WHERE id IN (SELECT source_id FROM transactions UNION SELECT target_id FROM transactions)
    """)
    print(f"    dropped {n_excl:,} nodes (non-50+DC or committees absent from cm) -> "
          f"{con.execute('SELECT COUNT(*) FROM transactions').fetchone()[0]:,} transactions, "
          f"{con.execute('SELECT COUNT(*) FROM nodes').fetchone()[0]:,} nodes")

    con.close()
    print(f"[done] wrote {db_path}")


def main() -> None:
    years = [int(sys.argv[1])] if len(sys.argv) > 1 else YEARS
    for year in years:
        print(f"\n========================= cycle {year} =========================")
        build_cycle(year)
    print(f"\n[all done] duckdbs in {DUCKDB_DIR}")

if __name__ == "__main__":
    main()
