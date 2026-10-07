#!/usr/bin/env python3
"""Score BD candidates with ViennaRNA RNAup (moved out of BD_screening).

Typical path: run BLAST/variant filter in BD_screening, then score the pass
table here::

    ../.venv/bin/python -m model.score_rnaup \\
        --utr ../examples/LETM1.fasta \\
        --table ../runs/<id>/candidates.tsv \\
        --out data/retrospective_study/letm1_rnaup.tsv

Or score every wet-lab tile without dropping BLAST failures first::

    ../.venv/bin/python -m model.score_rnaup \\
        --utr ../examples/LETM1.fasta \\
        --sites ../examples/LETM1.sites.fasta \\
        --out letm1_all_tiles_rnaup.tsv
"""

from __future__ import annotations

import argparse
import logging
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq


ENERGY_COLUMNS = [
    "rnaup_dG_total",
    "rnaup_dG_duplex",
    "rnaup_dGu_target",
    "rnaup_dGu_query",
    "interaction_target_start",
    "interaction_target_end",
    "interaction_query_start",
    "interaction_query_end",
    "interaction_utr_start",
    "interaction_utr_end",
    "anchor_overlap",
]


def parse_rnaup_output(text: str) -> dict[str, float | int] | None:
    energy = re.search(
        r"\(\s*(-?\d+(?:\.\d+)?)\s*=\s*(-?\d+(?:\.\d+)?)\s*\+\s*"
        r"(-?\d+(?:\.\d+)?)\s*\+\s*(-?\d+(?:\.\d+)?)\s*\)",
        text,
    )
    coords = re.search(r"(\d+)\s*,\s*(\d+)\s*:\s*(\d+)\s*,\s*(\d+)", text)
    if not energy or not coords:
        return None
    return {
        "rnaup_dG_total": float(energy.group(1)),
        "rnaup_dG_duplex": float(energy.group(2)),
        "rnaup_dGu_target": float(energy.group(3)),
        "rnaup_dGu_query": float(energy.group(4)),
        "interaction_target_start": int(coords.group(1)),
        "interaction_target_end": int(coords.group(2)),
        "interaction_query_start": int(coords.group(3)),
        "interaction_query_end": int(coords.group(4)),
    }


def anchor_interaction(
    parsed: dict[str, float | int],
    slice_utr_start: int,
    candidate_utr_start: int,
    candidate_utr_end: int,
    min_overlap: float = 1.0,
) -> tuple[bool, float, int, int]:
    local_start = min(int(parsed["interaction_target_start"]), int(parsed["interaction_target_end"])) - 1
    local_end = max(int(parsed["interaction_target_start"]), int(parsed["interaction_target_end"]))
    utr_start, utr_end = slice_utr_start + local_start, slice_utr_start + local_end
    overlap = max(0, min(utr_end, candidate_utr_end) - max(utr_start, candidate_utr_start))
    interaction_length = max(1, utr_end - utr_start)
    fraction = overlap / interaction_length
    return fraction >= min_overlap, fraction, utr_start, utr_end


def run_rnaup(
    utr_seq: str,
    bd_seq: str,
    utr_start: int,
    utr_end: int,
    *,
    context: int = 120,
    window: int = 40,
    temperature: float = 37.0,
    include_both: bool = False,
    min_anchor_overlap: float = 1.0,
    rnaup_exe: str | None = None,
) -> dict[str, object]:
    exe = rnaup_exe or shutil.which("RNAup")
    if not exe:
        return {"status": "failed", "failure_reason": "RNAup_not_found"}
    slice_start = max(0, utr_start - context)
    slice_end = min(len(utr_seq), utr_end + context)
    target_rna = utr_seq[slice_start:slice_end].upper().replace("T", "U")
    query_rna = bd_seq.upper().replace("T", "U")
    command = [
        exe, "--interaction_first", "--window", str(window),
        "--temp", str(temperature),
    ]
    if include_both:
        command.append("--include_both")
    try:
        with tempfile.TemporaryDirectory(prefix="bd_design_rnaup_") as work_dir:
            proc = subprocess.run(
                command,
                input=f">target\n{target_rna}\n>query\n{query_rna}\n",
                capture_output=True,
                text=True,
                timeout=120,
                cwd=work_dir,
            )
    except subprocess.TimeoutExpired:
        return {"status": "failed", "failure_reason": "RNAup_timeout"}
    except OSError:
        return {"status": "failed", "failure_reason": "RNAup_execution_error"}
    if proc.returncode:
        return {"status": "failed", "failure_reason": f"RNAup_exit_{proc.returncode}"}
    parsed = parse_rnaup_output(proc.stdout)
    if parsed is None:
        return {"status": "failed", "failure_reason": "RNAup_parse_error"}
    anchored, overlap, hit_start, hit_end = anchor_interaction(
        parsed, slice_start, utr_start, utr_end, min_anchor_overlap
    )
    parsed.update(
        {
            "interaction_utr_start": hit_start,
            "interaction_utr_end": hit_end,
            "anchor_overlap": overlap,
            "status": "eligible" if anchored else "failed",
            "failure_reason": "" if anchored else "off_anchor_interaction",
        }
    )
    return parsed


def _load_utr(path: Path) -> str:
    records = list(SeqIO.parse(str(path), "fasta"))
    if len(records) != 1:
        raise SystemExit(f"--utr must contain exactly one sequence, found {len(records)}")
    return str(records[0].seq).upper().replace("U", "T")


def _sites_from_fasta(utr: str, path: Path) -> pd.DataFrame:
    rows = []
    for rec in SeqIO.parse(str(path), "fasta"):
        sense = str(rec.seq).upper().replace("U", "T")
        hits = [m.start() for m in re.finditer(re.escape(sense), utr)]
        if len(hits) != 1:
            rows.append({
                "name": rec.id,
                "utr_start": pd.NA,
                "utr_end": pd.NA,
                "bd_sequence": str(Seq(sense).reverse_complement()),
                "target_sequence": sense,
                "status": "failed",
                "failure_reason": "not_in_utr" if not hits else "ambiguous",
            })
            continue
        i = hits[0]
        rows.append({
            "name": rec.id,
            "utr_start": i,
            "utr_end": i + len(sense),
            "bd_sequence": str(Seq(sense).reverse_complement()),
            "target_sequence": sense,
            "status": "pending",
            "failure_reason": "",
        })
    return pd.DataFrame(rows)


def _sites_from_table(path: Path, statuses: set[str]) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    needed = {"name", "utr_start", "utr_end", "bd_sequence"}
    missing = needed - set(df.columns)
    if missing:
        raise SystemExit(f"--table missing columns: {sorted(missing)}")
    if "status" in df.columns and statuses:
        df = df[df["status"].isin(statuses)].copy()
    return df.reset_index(drop=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--utr", type=Path, required=True, help="Target 3′UTR FASTA")
    src = parser.add_mutually_exclusive_group(required=True)
    src.add_argument("--table", type=Path, help="Screening TSV (candidates or all_binding_sites)")
    src.add_argument("--sites", type=Path, help="Sense-site FASTA to place on the UTR")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--status", default="pass,pending,eligible",
                        help="With --table, keep these statuses (comma-separated)")
    parser.add_argument("--context", type=int, default=120)
    parser.add_argument("--window", type=int, default=40)
    parser.add_argument("--temp", type=float, default=37.0)
    parser.add_argument("--include-both", action="store_true")
    parser.add_argument("--min-anchor-overlap", type=float, default=1.0)
    parser.add_argument("--rnaup", type=str, default=None, help="Path to RNAup executable")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(message)s")
    utr = _load_utr(args.utr)
    if args.sites:
        table = _sites_from_fasta(utr, args.sites)
    else:
        statuses = {s.strip() for s in args.status.split(",") if s.strip()}
        table = _sites_from_table(args.table, statuses)

    results = []
    for row in table.itertuples(index=False):
        record = {c: getattr(row, c, pd.NA) for c in table.columns}
        if pd.isna(getattr(row, "utr_start", pd.NA)) or pd.isna(getattr(row, "utr_end", pd.NA)):
            record.update({k: pd.NA for k in ENERGY_COLUMNS})
            record["status"] = getattr(row, "status", "failed")
            record["failure_reason"] = getattr(row, "failure_reason", "missing_coordinates")
            results.append(record)
            continue
        scored = run_rnaup(
            utr,
            str(row.bd_sequence),
            int(row.utr_start),
            int(row.utr_end),
            context=args.context,
            window=args.window,
            temperature=args.temp,
            include_both=args.include_both,
            min_anchor_overlap=args.min_anchor_overlap,
            rnaup_exe=args.rnaup,
        )
        record.update(scored)
        results.append(record)
        logging.info("%s  %s  dG=%s", row.name, scored.get("status"), scored.get("rnaup_dG_total"))

    out = pd.DataFrame(results)
    eligible = out["status"].eq("eligible")
    out.loc[eligible, "rank"] = (
        out.loc[eligible, "rnaup_dG_total"].rank(method="min", ascending=True).astype(int)
    )
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, sep="\t", index=False)
    print(f"wrote {args.out}  ({int(eligible.sum())} eligible / {len(out)} sites)")


if __name__ == "__main__":
    main()
