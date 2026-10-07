#!/usr/bin/env python3
"""Fold LETM1 / NSD2 5′ and 3′ UTRs with RNAfold and annotate RNAplot EPS.

Marks: 5′ cap, AUG-adjacent end of 5′UTR, stop-adjacent start of 3′UTR,
polyA-proximal 3′ end, and wet-lab BD tiles (3l / 3s).

Usage (BD_design root)::

    ../.venv/bin/python -m simulation.plot_utr_secondary_structure
"""

from __future__ import annotations

import subprocess
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCREENING = ROOT.parent
OUT = ROOT / "docs" / "figures" / "engineering" / "utr_ss"
CACHE = ROOT / "data" / "retrospective_study"
LETM1_3UTR = SCREENING / "examples" / "LETM1.fasta"
NSD2_3UTR = SCREENING / "examples" / "NSD2.fasta"
LETM1_SITES = SCREENING / "runs/20260926T120209Z_3ce6a3cc/all_binding_sites.tsv"
NSD2_SITES = SCREENING / "runs/20260926T120220Z_fc8c4409/all_binding_sites.tsv"

GENES = {
    "LETM1": {
        "enst": "ENST00000302787",
        "utr5_len": 209,
        "utr3_fa": LETM1_3UTR,
        "sites_tsv": LETM1_SITES,
    },
    "NSD2": {
        "enst": "ENST00000508803",
        "utr5_len": 182,
        "utr3_fa": NSD2_3UTR,
        "sites_tsv": NSD2_SITES,
    },
}


def read_fasta(path: Path) -> str:
    lines = path.read_text().splitlines()
    return "".join(l.strip() for l in lines if l and not l.startswith(">")).upper().replace("T", "U")


def fetch_cdna(enst: str) -> str:
    cache = CACHE / f"{enst}.cdna.fa"
    if cache.is_file():
        return read_fasta(cache)
    url = f"https://rest.ensembl.org/sequence/id/{enst}?type=cdna"
    req = urllib.request.Request(url, headers={"Content-Type": "text/plain", "User-Agent": "bsst"})
    with urllib.request.urlopen(req, timeout=60) as resp:
        seq = resp.read().decode("ascii").strip().upper().replace("T", "U")
    CACHE.mkdir(parents=True, exist_ok=True)
    cache.write_text(f">{enst} cdna\n" + "\n".join(seq[i : i + 60] for i in range(0, len(seq), 60)) + "\n")
    return seq


def load_sites(tsv: Path) -> list[tuple[str, int, int]]:
    rows = []
    lines = tsv.read_text().splitlines()
    hdr = lines[0].split("\t")
    i_name, i_a, i_b = hdr.index("name"), hdr.index("utr_start"), hdr.index("utr_end")
    for line in lines[1:]:
        p = line.split("\t")
        a, b = int(float(p[i_a])), int(float(p[i_b]))
        rows.append((p[i_name].split("-", 1)[-1], a + 1, b))  # 1-based closed
    return rows


def tile_color(short: str) -> tuple[float, float, float]:
    if short.startswith("3s"):
        return (0.50, 0.22, 0.72)
    return (0.78, 0.18, 0.16)


def postscript_marks(
    length: int,
    kind: str,
    sites: list[tuple[str, int, int]],
    start: int = 1,
    full_len: int | None = None,
) -> str:
    """Build RNAplot --post string. `start` is 1-based coord of seq[0] on the parent UTR."""
    cmds: list[str] = []
    parent_len = full_len or (start + length - 1)
    end = start + length - 1

    def clip(i: int, j: int) -> tuple[int, int] | None:
        a, b = max(i, start), min(j, end)
        if a > b:
            return None
        return a - start + 1, b - start + 1

    if kind == "5utr":
        m = clip(1, min(8, length))
        if m:
            cmds.append(f"{m[0]} {m[1]} 9 0.15 0.35 0.85 omark")
        m = clip(max(1, length - 7), length)
        if m:
            cmds.append(f"{m[0]} {m[1]} 9 0.10 0.62 0.28 omark")
        cmds.append("1 cmark")
        cmds.append(f"{length} cmark")
        return " ".join(cmds)

    cmds.append("/Helvetica-Bold 2.4 LabelFont")
    if start == 1:
        m = clip(1, 8)
        if m:
            cmds.append(f"{m[0]} {m[1]} 9 0.10 0.62 0.28 omark")
    if end >= parent_len:
        m = clip(max(1, parent_len - 14), parent_len)
        if m:
            cmds.append(f"{m[0]} {m[1]} 10 0.90 0.50 0.08 omark")
    for short, i, j in sites:
        m = clip(i, j)
        if not m:
            continue
        r, g, b = tile_color(short)
        cmds.append(f"{m[0]} {m[1]} 7 {r:.2f} {g:.2f} {b:.2f} omark")
        mid = (m[0] + m[1]) // 2
        label = short.replace("(", "").replace(")", "")
        k = len([c for c in cmds if c.endswith(" Label")])
        dx = 1.4 if k % 2 == 0 else -6.0
        dy = 2.0 + 1.2 * (k % 3)
        cmds.append(f"0 0 0 setrgbcolor {mid} {dx:.1f} {dy:.1f} ({label}) Label")
    return " ".join(cmds)


def rnafold_plot(seq: str, name: str, post: str, work: Path) -> Path:
    work.mkdir(parents=True, exist_ok=True)
    fold = work / f"{name}.fold"
    folded = subprocess.run(
        ["RNAfold", "--noPS"],
        input=f">{name}\n{seq}\n",
        cwd=work, check=True, capture_output=True, text=True,
    ).stdout
    fold.write_text(folded)
    subprocess.run(
        ["RNAplot", "-f", "eps", "--infile", str(fold), "--post", post],
        cwd=work, check=True, capture_output=True, text=True,
    )
    eps = work / f"{name}_ss.eps"
    if not eps.is_file():
        raise FileNotFoundError(eps)
    png = work / f"{name}_ss.png"
    subprocess.run(
        [
            "gs", "-dSAFER", "-dBATCH", "-dNOPAUSE", "-dEPSCrop",
            "-sDEVICE=png16m", "-r140", f"-sOutputFile={png}", str(eps),
        ],
        check=True, capture_output=True,
    )
    return png


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for gene, meta in GENES.items():
        cdna = fetch_cdna(meta["enst"])
        utr5 = cdna[: meta["utr5_len"]]
        utr3 = read_fasta(meta["utr3_fa"])
        cdna3 = cdna[-len(utr3) :]
        if cdna3 != utr3:
            print(f"WARN {gene}: Ensembl 3′UTR != examples FASTA; plotting the FASTA used for sites")
        sites = load_sites(meta["sites_tsv"])
        print(f"{gene} 5′UTR {len(utr5)} nt, 3′UTR {len(utr3)} nt, {len(sites)} tiles")

        p5 = rnafold_plot(utr5, f"{gene}_5utr", postscript_marks(len(utr5), "5utr", []), OUT)
        print(" wrote", p5.relative_to(ROOT))
        p3 = rnafold_plot(utr3, f"{gene}_3utr", postscript_marks(len(utr3), "3utr", sites), OUT)
        print(" wrote", p3.relative_to(ROOT))

        groups: dict[str, list] = {}
        for short, i, j in sites:
            key = "3s" if short.startswith("3s") else "3l"
            groups.setdefault(key, []).append((short, i, j))
        for key, tiles in groups.items():
            lo = min(t[1] for t in tiles) - 80
            hi = max(t[2] for t in tiles) + 80
            lo = max(1, lo)
            hi = min(len(utr3), hi)
            window = utr3[lo - 1 : hi]
            post = postscript_marks(len(window), "3utr", tiles, start=lo, full_len=len(utr3))
            png = rnafold_plot(window, f"{gene}_3utr_{key}_zoom", post, OUT)
            print(f" wrote {png.relative_to(ROOT)}  window {lo}-{hi}")

    legend = OUT / "README.md"
    legend.write_text(
        """# UTR secondary structures (RNAfold / RNAplot)

MFE structures, 37 °C, ViennaRNA 2.7. Coordinates on 3′UTR match `examples/*.fasta`.

| File | What |
|------|------|
| `LETM1_5utr_ss` / `NSD2_5utr_ss` | Full 5′UTR |
| `LETM1_3utr_ss` / `NSD2_3utr_ss` | Full 3′UTR (crowded; use zooms) |
| `*_3l_zoom_ss` / `*_3s_zoom_ss` | ±80 nt around BD ladders |

Marks (EPS `omark`):

- **blue** — 5′ end / cap (first 8 nt of 5′UTR)
- **green** — AUG-adjacent (last 8 nt of 5′UTR) and stop-adjacent (first 8 nt of 3′UTR)
- **orange** — polyA-proximal (last 15 nt of 3′UTR)
- **red** — 3l BD tiles (name at tile midpoint)
- **purple** — 3s BD tiles (name at tile midpoint)
"""
    )
    print("wrote", legend.relative_to(ROOT))


if __name__ == "__main__":
    main()
