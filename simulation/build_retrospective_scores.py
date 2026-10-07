#!/usr/bin/env python3
"""Build retrospective_scores.tsv from fluc Excel + RNAup energy tables.

Prefer ``python -m model.score_rnaup --sites …`` for tile energies. Legacy
defaults still read screening run TSVs that already carry RNAup columns::

    LETM1: runs/20260926T120209Z_3ce6a3cc
    NSD2:  runs/20260926T120220Z_fc8c4409
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SCREENING = ROOT.parent


def block_means_letm1(df: pd.DataFrame, header_row: int, value_rows, dose: str, date: str) -> list[dict]:
    cols = list(df.iloc[header_row].tolist())
    control_idx = next(i for i, c in enumerate(cols) if str(c).strip().lower() == "control")
    records = []
    for r in value_rows:
        row = df.iloc[r]
        ctrl = float(row[control_idx])
        if ctrl <= 0 or pd.isna(ctrl):
            continue
        for i, c in enumerate(cols):
            if pd.isna(c):
                continue
            name = str(c).strip()
            if name.lower() == "control" or name == "5s" or not name.startswith("3l-"):
                continue
            val = row[i]
            if pd.isna(val):
                continue
            records.append({
                "gene": "LETM1",
                "site_short": name,
                "site": f"LETM1-{name}",
                "date": date,
                "dose": dose,
                "fluc": float(val),
                "control": ctrl,
                "fold": float(val) / ctrl,
            })
    return records


def block_means_nsd2(df: pd.DataFrame, header_row: int, value_rows, dose: str, date: str) -> list[dict]:
    cols = list(df.iloc[header_row].tolist())
    control_idx = next(i for i, c in enumerate(cols) if str(c).strip().lower() == "control")
    records = []
    for r in value_rows:
        row = df.iloc[r]
        ctrl = float(row[control_idx])
        if ctrl <= 0 or pd.isna(ctrl):
            continue
        for i, c in enumerate(cols):
            if pd.isna(c) or str(c).strip().lower() == "control":
                continue
            name = str(c).strip()
            val = row[i]
            if pd.isna(val):
                continue
            records.append({
                "gene": "NSD2",
                "site_short": name,
                "site": f"NSD2-{name}",
                "date": date,
                "dose": dose,
                "fluc": float(val),
                "control": ctrl,
                "fold": float(val) / ctrl,
            })
    return records


def summarize(rep: pd.DataFrame) -> pd.DataFrame:
    return (
        rep.groupby(["gene", "site", "site_short"], as_index=False)
        .agg(fold_mean=("fold", "mean"), fold_sd=("fold", "std"), n_replicates=("fold", "count"))
    )


def merge_energy(sum_df: pd.DataFrame, energy_df: pd.DataFrame, rnaup_run: str) -> pd.DataFrame:
    e = energy_df.rename(columns={"name": "site"})
    m = sum_df.merge(
        e[[
            "site", "rnaup_dG_total", "rnaup_dG_duplex", "rnaup_dGu_target",
            "rnaup_dGu_query", "anchor_overlap", "status",
        ]],
        on="site",
        how="left",
    )
    m["dG_rank"] = m["rnaup_dG_total"].rank(method="min", ascending=True).astype(int)
    m["activity_rank"] = m["fold_mean"].rank(method="min", ascending=False).astype(int)
    m["rnaup_run"] = rnaup_run
    m["fluc_source"] = (
        "20ng/well blocks in retrospective xlsx "
        "(fold = fluc/control per replicate, then mean)"
    )
    return m


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--letm1-run", default="runs/20260926T120209Z_3ce6a3cc")
    parser.add_argument("--nsd2-run", default="runs/20260926T120220Z_fc8c4409")
    args = parser.parse_args()

    let_xlsx = next((ROOT / "data/retrospective_study").glob("*LETM1fluc.xlsx"))
    nsd_xlsx = next((ROOT / "data/retrospective_study").glob("*NSD2fluc.xlsx"))
    raw = pd.read_excel(let_xlsx, header=None)
    nsd_raw = pd.read_excel(nsd_xlsx, header=None)

    let_records: list[dict] = []
    let_records += block_means_letm1(raw, 67, range(68, 71), "20ng", "260813")
    let_records += block_means_letm1(raw, 74, range(75, 78), "20ng", "260814")
    let_records += block_means_letm1(raw, 81, range(82, 87), "20ng", "260815")
    let_records += block_means_letm1(raw, 90, range(91, 96), "20ng", "260817")

    nsd_records: list[dict] = []
    nsd_records += block_means_nsd2(nsd_raw, 1, range(2, 7), "20ng", "260814")
    nsd_records += block_means_nsd2(nsd_raw, 20, range(21, 27), "20ng", "260815")

    let_20 = pd.DataFrame(let_records)
    nsd_20 = pd.DataFrame(nsd_records)
    # 3s-e plasmid was not obtained; drop wet-lab rows for that tile.
    nsd_20 = nsd_20.loc[nsd_20["site_short"] != "3s-e"].copy()
    e_let = pd.read_csv(SCREENING / args.letm1_run / "all_binding_sites.tsv", sep="\t")
    e_nsd = pd.read_csv(SCREENING / args.nsd2_run / "all_binding_sites.tsv", sep="\t")

    combined = pd.concat([
        merge_energy(summarize(let_20), e_let, Path(args.letm1_run).name),
        merge_energy(summarize(nsd_20), e_nsd, Path(args.nsd2_run).name),
    ], ignore_index=True)

    out_dir = ROOT / "data" / "retrospective_study"
    combined.to_csv(out_dir / "retrospective_scores.tsv", sep="\t", index=False)
    pd.concat([let_20, nsd_20], ignore_index=True).to_csv(
        out_dir / "retrospective_fluc_replicates_20ng.tsv", sep="\t", index=False
    )
    print(f"wrote {out_dir / 'retrospective_scores.tsv'} ({len(combined)} rows)")


if __name__ == "__main__":
    main()
