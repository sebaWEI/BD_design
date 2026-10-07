#!/usr/bin/env python3
"""Generate Engineering-cycle figures for the iGEM wiki.

Quantitative panels (C3-B/C/D/E) read
``data/retrospective_study/retrospective_scores.tsv``.
Conceptual panels (C1, C2, C3-A, C3.5) are drawn as SVG with matplotlib.

Usage (from BD_design root)::

    ../.venv/bin/python -m simulation.plot_engineering_figures
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
SCREENING = ROOT.parent  # sibling nested under BD_screening checkout
DATA = ROOT / "data" / "retrospective_study" / "retrospective_scores.tsv"
OUT = ROOT / "docs" / "figures" / "engineering"

# Quiet, print-friendly palette (avoid purple-on-white AI defaults).
INK = "#1c1917"
MUTED = "#57534e"
ACCENT = "#0f766e"
ACCENT2 = "#b45309"
FAIL = "#b91c1c"
OK = "#15803d"
PANEL = "#fafaf9"
LINE = "#d6d3d1"


def _style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 10,
        "figure.facecolor": "white",
        "axes.facecolor": PANEL,
        "axes.edgecolor": INK,
        "axes.labelcolor": INK,
        "xtick.color": INK,
        "ytick.color": INK,
        "text.color": INK,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "savefig.bbox": "tight",
        "savefig.dpi": 200,
    })


def _save(fig: plt.Figure, stem: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("svg", "pdf", "png"):
        path = OUT / f"{stem}.{ext}"
        fig.savefig(path, facecolor="white", bbox_inches="tight", pad_inches=0.2)
        print("wrote", path.relative_to(ROOT))
    plt.close(fig)


def load_scores() -> pd.DataFrame:
    return pd.read_csv(DATA, sep="\t")


def plot_scatter(df: pd.DataFrame) -> None:
    """C3-B / C3-C: ΔG_total vs wet-lab fold (two panels)."""
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2), sharey=False)
    for ax, gene in zip(axes, ("LETM1", "NSD2")):
        g = df[df.gene == gene].copy()
        ax.errorbar(
            g["rnaup_dG_total"], g["fold_mean"],
            yerr=g["fold_sd"].fillna(0),
            fmt="o", color=ACCENT, ecolor=MUTED, elinewidth=0.8,
            capsize=2, markersize=7, zorder=3,
        )
        for _, row in g.iterrows():
            ax.annotate(
                row["site_short"],
                (row["rnaup_dG_total"], row["fold_mean"]),
                textcoords="offset points", xytext=(4, 4),
                fontsize=7, color=MUTED,
            )
        ax.axhline(1.0, color=LINE, lw=1, ls="--", zorder=1)
        ax.set_xlabel(r"RNAup $\Delta G_{\mathrm{total}}$ (kcal/mol)  → more negative")
        ax.set_ylabel("Luciferase fold vs control")
        # Panel label only (wiki caption carries title / stats).
        ax.text(0.98, 0.98, gene, transform=ax.transAxes, ha="right", va="top",
                fontsize=11, fontweight="bold", color=INK)
        # Put more-negative ΔG on the right so “better binding” reads left→right
        # under the (failed) hypothesis; algebraic Spearman still uses raw ΔG.
        lo = float(g["rnaup_dG_total"].min())
        hi = float(g["rnaup_dG_total"].max())
        pad = 0.4
        ax.set_xlim(hi + pad, lo - pad)
    fig.tight_layout()
    _save(fig, "eng_c3_scatter_letm1_nsd2")


def plot_rank_mismatch(df: pd.DataFrame) -> None:
    """C3-D: order by ΔG vs order by activity."""
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 6.5), gridspec_kw={"hspace": 0.45})
    cmap = plt.cm.tab10

    for col, gene in enumerate(("LETM1", "NSD2")):
        g = df[df.gene == gene].copy()
        colors = {s: cmap(i % 10) for i, s in enumerate(sorted(g["site_short"]))}

        by_dg = g.sort_values("rnaup_dG_total")  # most negative first
        by_act = g.sort_values("fold_mean", ascending=False)

        ax_top = axes[0, col]
        ax_bot = axes[1, col]
        y = np.arange(len(g))

        ax_top.barh(
            y, by_dg["rnaup_dG_total"].abs(),
            color=[colors[s] for s in by_dg["site_short"]],
            edgecolor="none", height=0.7,
        )
        ax_top.set_yticks(y)
        ax_top.set_yticklabels(by_dg["site_short"])
        ax_top.invert_yaxis()
        ax_top.set_xlabel(r"|$\Delta G_{\mathrm{total}}$| (kcal/mol)")
        ax_top.set_title(f"{gene}: ranked by ΔG (best → worst)")

        ax_bot.barh(
            y, by_act["fold_mean"],
            color=[colors[s] for s in by_act["site_short"]],
            edgecolor="none", height=0.7,
        )
        ax_bot.set_yticks(y)
        ax_bot.set_yticklabels(by_act["site_short"])
        ax_bot.invert_yaxis()
        ax_bot.set_xlabel("Luciferase fold vs control")
        ax_bot.set_title(f"{gene}: ranked by wet-lab activity (high → low)")
        ax_bot.axvline(1.0, color=LINE, lw=1, ls="--")

    fig.suptitle(
        "Cycle 3 — Same tiles, mismatched order (ΔG ranking ≠ activity ranking)",
        fontsize=12, y=0.98,
    )
    _save(fig, "eng_c3_rank_mismatch")


def plot_energy_terms(df: pd.DataFrame) -> None:
    """C3-E: alternative energy axes still fail to track activity."""
    terms = [
        ("rnaup_dG_total", r"$\Delta G_{\mathrm{total}}$"),
        ("rnaup_dG_duplex", r"$\Delta G_{\mathrm{duplex}}$"),
        ("rnaup_dGu_target", r"$\Delta G_{\mathrm{open,target}}$"),
        ("rnaup_dGu_query", r"$\Delta G_{\mathrm{open,query}}$"),
    ]
    fig, axes = plt.subplots(2, 4, figsize=(11, 5.5), sharey="row")
    for row, gene in enumerate(("LETM1", "NSD2")):
        g = df[df.gene == gene]
        for col, (key, label) in enumerate(terms):
            ax = axes[row, col]
            ax.scatter(g["fold_mean"], g[key], color=ACCENT, s=36, zorder=3)
            rho, p = stats.spearmanr(g["fold_mean"], g[key])
            ax.set_title(f"{gene}\n{label}\nρ={rho:.2f}", fontsize=9)
            if row == 1:
                ax.set_xlabel("Fold vs control")
            if col == 0:
                ax.set_ylabel("kcal/mol")
            ax.axvline(1.0, color=LINE, lw=0.8, ls="--")
    fig.suptitle(
        "Cycle 3 — Reweighting individual energy terms does not rescue the correlation",
        fontsize=12, y=1.02,
    )
    fig.tight_layout()
    _save(fig, "eng_c3_energy_terms")


def plot_c1_overconservative() -> None:
    """C1: deleting every BLAST hit is too broad."""
    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5)
    ax.axis("off")
    ax.set_title("Cycle 1 — Sequence similarity alone is over-conservative", loc="left", pad=12)

    def mrna(y: float, label: str, regions: list[tuple[float, float, str, str]]) -> None:
        ax.text(0.15, y + 0.55, label, fontsize=9, color=MUTED, va="bottom")
        for x0, x1, name, color in regions:
            ax.add_patch(FancyBboxPatch(
                (x0, y), x1 - x0, 0.45, boxstyle="round,pad=0.02,rounding_size=0.05",
                facecolor=color, edgecolor=INK, linewidth=0.8,
            ))
            ax.text((x0 + x1) / 2, y + 0.22, name, ha="center", va="center", fontsize=8, color=INK)

    mrna(3.6, "Target mRNA", [
        (1.0, 2.2, "5′UTR", "#e7e5e4"),
        (2.2, 5.5, "CDS", "#d6d3d1"),
        (5.5, 8.8, "3′UTR (BD acts here)", "#99f6e4"),
    ])
    mrna(2.2, "Off-target A", [
        (1.0, 2.2, "5′UTR", "#e7e5e4"),
        (2.2, 5.5, "CDS", "#fecaca"),
        (5.5, 8.8, "3′UTR", "#e7e5e4"),
    ])
    mrna(0.8, "Off-target B", [
        (1.0, 2.2, "5′UTR", "#e7e5e4"),
        (2.2, 5.5, "CDS", "#e7e5e4"),
        (5.5, 8.8, "3′UTR", "#fecaca"),
    ])

    # BLAST hit brackets
    ax.annotate("", xy=(3.5, 2.65), xytext=(6.5, 3.6),
                arrowprops=dict(arrowstyle="-", color=FAIL, lw=1.6, connectionstyle="arc3,rad=-0.2"))
    ax.text(4.6, 3.15, "BLAST hit\nin CDS", color=FAIL, fontsize=8, ha="center")
    ax.annotate("", xy=(6.8, 1.25), xytext=(6.8, 3.6),
                arrowprops=dict(arrowstyle="-", color=ACCENT2, lw=1.6))
    ax.text(7.5, 2.4, "BLAST hit\nin 3′UTR", color=ACCENT2, fontsize=8)

    # Rule boxes
    ax.add_patch(FancyBboxPatch(
        (0.2, 4.35), 4.4, 0.5, boxstyle="round,pad=0.02,rounding_size=0.05",
        facecolor="#fee2e2", edgecolor=FAIL, linewidth=1.0,
    ))
    ax.text(2.4, 4.6, "Old rule: any BLAST hit → drop site  ✗", ha="center", va="center",
            fontsize=9, color=FAIL, fontweight="bold")
    ax.add_patch(FancyBboxPatch(
        (5.0, 4.35), 4.7, 0.5, boxstyle="round,pad=0.02,rounding_size=0.05",
        facecolor="#ecfdf5", edgecolor=OK, linewidth=1.0,
    ))
    ax.text(7.35, 4.6, "Learn: only UTR-relevant hits matter", ha="center", va="center",
            fontsize=9, color=OK, fontweight="bold")

    ax.text(5.0, 0.15,
            "A CDS match is sequence similarity, not the interaction a 3′UTR-acting element set out to avoid.",
            ha="center", fontsize=8, color=MUTED)
    _save(fig, "eng_c1_blast_overconservative")


def plot_c2_utr_strand_variants() -> None:
    """C2: region gate, strand, variants, readable table."""
    fig = plt.figure(figsize=(10, 5.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.2, 0.85], hspace=0.55, wspace=0.35)
    fig.suptitle("Cycle 2 — UTR context, same-strand only, optional variants, readable summary",
                 fontsize=12, y=0.98)

    # Panel 1: region
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.set_xlim(0, 10)
    ax1.set_ylim(0, 3)
    ax1.axis("off")
    ax1.set_title("Region gate", fontsize=10, loc="left")
    for x0, x1, lab, fc in [
        (0.5, 3.0, "5′UTR", "#99f6e4"),
        (3.0, 7.0, "CDS — keep", "#e7e5e4"),
        (7.0, 9.5, "3′UTR", "#99f6e4"),
    ]:
        ax1.add_patch(FancyBboxPatch((x0, 1.2), x1 - x0, 0.7, boxstyle="round,pad=0.02,rounding_size=0.04",
                                     facecolor=fc, edgecolor=INK, lw=0.8))
        ax1.text((x0 + x1) / 2, 1.55, lab, ha="center", va="center", fontsize=8)
    ax1.text(5, 0.4, "Drop only if hit overlaps 5′/3′UTR\nand length ≥ N nt",
             ha="center", fontsize=8, color=MUTED)

    # Panel 2: strand
    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_xlim(0, 10)
    ax2.set_ylim(0, 3)
    ax2.axis("off")
    ax2.set_title("Same-strand only", fontsize=10, loc="left")
    ax2.annotate("", xy=(8, 2.1), xytext=(2, 2.1),
                 arrowprops=dict(arrowstyle="->", color=OK, lw=2))
    ax2.text(5, 2.35, "sstart < send  → consider", ha="center", color=OK, fontsize=8)
    ax2.annotate("", xy=(2, 0.9), xytext=(8, 0.9),
                 arrowprops=dict(arrowstyle="->", color=FAIL, lw=2))
    ax2.text(5, 1.15, "sstart > send  → ignore", ha="center", color=FAIL, fontsize=8)
    ax2.text(5, 0.25, "Reverse hit is the BD sequence elsewhere,\nnot a bindable sense site",
             ha="center", fontsize=7, color=MUTED)

    # Panel 3: variants
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.set_xlim(0, 10)
    ax3.set_ylim(0, 3)
    ax3.axis("off")
    ax3.set_title("Optional --variants", fontsize=10, loc="left")
    ax3.add_patch(FancyBboxPatch((0.8, 1.3), 8.4, 0.55, boxstyle="round,pad=0.02,rounding_size=0.04",
                                 facecolor="#fef3c7", edgecolor=ACCENT2, lw=0.8))
    ax3.text(5, 1.55, "candidate genomic window", ha="center", va="center", fontsize=8)
    for x, lab in [(2.5, "SNP"), (5.0, "indel"), (7.5, "SNP")]:
        ax3.plot(x, 2.2, marker="D", color=ACCENT2, markersize=9)
        ax3.text(x, 2.55, lab, ha="center", fontsize=7, color=ACCENT2)
    ax3.text(5, 0.45, "Off by default\ncommon_all overlap → drop",
             ha="center", fontsize=8, color=MUTED)

    # Bottom: fake readable table
    ax4 = fig.add_subplot(gs[1, :])
    ax4.axis("off")
    ax4.set_title("Readable blast_matches summary (Learn)", fontsize=10, loc="left")
    table = ax4.table(
        cellText=[
            ["site", "matched_gene", "region", "aligned_nt", "drops_site"],
            ["tile-A", "GENE_X", "3′UTR", "21", "yes"],
            ["tile-A", "GENE_Y", "— (CDS)", "28", "no"],
            ["tile-B", "GENE_Z", "5′UTR", "14", "no"],
        ],
        cellLoc="center",
        loc="center",
        colWidths=[0.16, 0.22, 0.18, 0.16, 0.16],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(9)
    table.scale(1, 1.6)
    for j in range(5):
        table[0, j].set_facecolor("#d6d3d1")
        table[0, j].set_text_props(fontweight="bold")
    table[1, 4].set_text_props(color=FAIL, fontweight="bold")
    table[2, 4].set_text_props(color=OK)
    table[3, 4].set_text_props(color=OK)

    _save(fig, "eng_c2_utr_strand_variants")


def plot_c3a_rnaup_scheme() -> None:
    """C3-A: opening + duplex → ΔG_total."""
    fig, ax = plt.subplots(figsize=(9.2, 4.2))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 5)
    ax.axis("off")
    ax.set_title("Cycle 3 — RNAup thermodynamic picture (exploratory ranking)", loc="left")

    # native fold
    ax.add_patch(FancyBboxPatch((0.4, 2.4), 3.2, 1.6, boxstyle="round,pad=0.03,rounding_size=0.08",
                                facecolor="#e7e5e4", edgecolor=INK, lw=1))
    ax.text(2.0, 3.5, "Target mRNA", ha="center", fontsize=9, fontweight="bold")
    ax.text(2.0, 3.0, "local secondary\nstructure", ha="center", fontsize=8, color=MUTED)

    ax.annotate("", xy=(4.4, 3.2), xytext=(3.7, 3.2),
                arrowprops=dict(arrowstyle="->", color=INK, lw=1.5))
    ax.text(4.05, 3.55, "open", fontsize=8, color=ACCENT2, ha="center")

    # open + duplex
    ax.add_patch(FancyBboxPatch((4.5, 2.2), 3.6, 2.0, boxstyle="round,pad=0.03,rounding_size=0.08",
                                facecolor="#ccfbf1", edgecolor=ACCENT, lw=1.2))
    ax.text(6.3, 3.7, "Duplex region", ha="center", fontsize=9, fontweight="bold", color=ACCENT)
    ax.text(6.3, 3.15, "40 nt sense site\n+ 120 nt context / side", ha="center", fontsize=8, color=MUTED)
    ax.plot([5.0, 8.0], [2.55, 2.55], color=ACCENT, lw=3)
    ax.plot([5.0, 8.0], [2.75, 2.75], color=ACCENT2, lw=3)
    ax.text(6.3, 2.35, "BD  (antisense)", ha="center", fontsize=7, color=ACCENT2)

    ax.annotate("", xy=(9.0, 3.2), xytext=(8.2, 3.2),
                arrowprops=dict(arrowstyle="->", color=INK, lw=1.5))

    ax.add_patch(FancyBboxPatch((9.1, 2.3), 2.6, 1.8, boxstyle="round,pad=0.03,rounding_size=0.08",
                                facecolor="#fff7ed", edgecolor=ACCENT2, lw=1.2))
    ax.text(10.4, 3.55, r"$\Delta G_{\mathrm{total}}$", ha="center", fontsize=11, fontweight="bold")
    ax.text(10.4, 2.85, "rank sites\n(more negative\nfirst)", ha="center", fontsize=8, color=MUTED)

    ax.text(
        6.0, 1.2,
        r"$\Delta G_{\mathrm{total}}=\Delta G_{\mathrm{duplex}}+\Delta G_{\mathrm{open,target}}+\Delta G_{\mathrm{open,query}}$",
        ha="center", fontsize=12,
    )
    ax.text(
        6.0, 0.45,
        "Retrospective Cycle 3 scored every wet-lab tile (correlation analysis),\n"
        "not only sites that pass today’s BLAST / --variants shortlist.",
        ha="center", fontsize=8, color=MUTED,
    )
    _save(fig, "eng_c3_rnaup_scheme")


def plot_c35_filter_not_predictor() -> None:
    """C3.5: reject predictor path; keep filter path."""
    fig, ax = plt.subplots(figsize=(10, 4.6))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 5)
    ax.axis("off")
    ax.set_title("Cycle 3.5 — BSST claims a filter, not an upregulation predictor", loc="left")

    # Left: predictor (rejected)
    ax.add_patch(FancyBboxPatch((0.3, 0.6), 5.2, 3.8, boxstyle="round,pad=0.04,rounding_size=0.1",
                                facecolor="#fef2f2", edgecolor=FAIL, lw=1.5, linestyle="--"))
    ax.text(2.9, 4.05, "Predictor path — not claimed", ha="center", fontsize=10,
            color=FAIL, fontweight="bold")
    steps_l = ["sites", "RNAup ΔG", "rank", "“predicted\nupregulation”"]
    for i, lab in enumerate(steps_l):
        x = 0.7 + i * 1.15
        ax.add_patch(FancyBboxPatch((x, 2.2), 1.0, 1.1, boxstyle="round,pad=0.02,rounding_size=0.06",
                                    facecolor="white", edgecolor=FAIL, lw=0.9))
        ax.text(x + 0.5, 2.75, lab, ha="center", va="center", fontsize=7, color=FAIL)
        if i < len(steps_l) - 1:
            ax.annotate("", xy=(x + 1.1, 2.75), xytext=(x + 1.0, 2.75),
                        arrowprops=dict(arrowstyle="->", color=FAIL, lw=1))
    # big X
    ax.plot([1.2, 4.6], [1.2, 3.6], color=FAIL, lw=3, alpha=0.5)
    ax.plot([1.2, 4.6], [3.6, 1.2], color=FAIL, lw=3, alpha=0.5)
    ax.text(2.9, 1.0, "failed wet-lab correlation", ha="center", fontsize=8, color=FAIL)

    # Right: filter
    ax.add_patch(FancyBboxPatch((5.9, 0.6), 5.8, 3.8, boxstyle="round,pad=0.04,rounding_size=0.1",
                                facecolor="#ecfdf5", edgecolor=OK, lw=1.5))
    ax.text(8.8, 4.05, "Filter path — shipped product", ha="center", fontsize=10,
            color=OK, fontweight="bold")
    steps_r = [
        (6.2, "sites"),
        (7.35, "same-strand\nUTR BLAST"),
        (8.7, "optional\n--variants"),
        (10.0, "shortlist"),
        (11.15, "wet\nlab"),
    ]
    for i, (x, lab) in enumerate(steps_r):
        ax.add_patch(FancyBboxPatch((x, 2.15), 1.05, 1.2, boxstyle="round,pad=0.02,rounding_size=0.06",
                                    facecolor="white", edgecolor=OK, lw=0.9))
        ax.text(x + 0.52, 2.75, lab, ha="center", va="center", fontsize=6.5, color=OK)
        if i < len(steps_r) - 1:
            ax.annotate("", xy=(steps_r[i + 1][0], 2.75), xytext=(x + 1.05, 2.75),
                        arrowprops=dict(arrowstyle="->", color=OK, lw=1))
    # optional RNAup bypass
    ax.annotate(
        "RNAup\noptional",
        xy=(10.5, 2.15), xytext=(8.8, 1.15),
        fontsize=7, color=MUTED, ha="center",
        arrowprops=dict(arrowstyle="->", color=MUTED, lw=1, ls="--",
                        connectionstyle="arc3,rad=0.3"),
    )
    ax.text(8.8, 0.85, "pre-screen · decision support · not a black-box activity model",
            ha="center", fontsize=8, color=MUTED)

    _save(fig, "eng_c35_filter_not_predictor")


RNA_COFOLD = {
    ("LETM1", "3l-1"): -56.62, ("LETM1", "3l-2"): -62.30, ("LETM1", "3l-3"): -51.25,
    ("LETM1", "3l-4"): -62.78, ("LETM1", "3l-5"): -63.85, ("LETM1", "3l-6"): -57.27,
    ("LETM1", "3l-7"): -53.17, ("LETM1", "3l-8"): -62.56, ("LETM1", "3l-9"): -60.45,
    ("NSD2", "3l-3"): -47.58, ("NSD2", "3l-4"): -51.95, ("NSD2", "3l-8"): -53.65,
    ("NSD2", "3l-9"): -54.64, ("NSD2", "3s-f"): -59.91, ("NSD2", "3s-g"): -57.03,
    ("NSD2", "3s-h"): -63.34, ("NSD2", "3s-i"): -60.44,
}

MM_PBSA = {
    ("LETM1", "3l-1"): -14860.75, ("LETM1", "3l-2"): -15623.12, ("LETM1", "3l-3"): -15717.80,
    ("LETM1", "3l-4"): -15404.06, ("LETM1", "3l-5"): -15284.33, ("LETM1", "3l-6"): -15183.77,
    ("LETM1", "3l-7"): -14730.90, ("LETM1", "3l-8"): -14673.87, ("LETM1", "3l-9"): -14923.07,
    ("NSD2", "3l-3"): -16529.23, ("NSD2", "3l-4"): -16502.95, ("NSD2", "3l-8"): -16674.64,
    ("NSD2", "3l-9"): -16700.23, ("NSD2", "3s-f"): -15002.44, ("NSD2", "3s-g"): -15397.68,
    ("NSD2", "3s-h"): -15169.72, ("NSD2", "3s-i"): -15686.02,
}


def _site_distance_table(scores: pd.DataFrame) -> pd.DataFrame:
    let = pd.read_csv(SCREENING / "runs/20260926T120209Z_3ce6a3cc/all_binding_sites.tsv", sep="\t")
    nsd = pd.read_csv(SCREENING / "runs/20260926T120220Z_fc8c4409/all_binding_sites.tsv", sep="\t")
    pos = pd.concat([let, nsd], ignore_index=True)[["name", "utr_end"]].rename(columns={"name": "site"})
    utr_len = {"LETM1": 2945, "NSD2": 3283}
    df = scores.merge(pos, on="site")
    df["dist_3p"] = df["gene"].map(utr_len) - df["utr_end"]
    df["rna_cofold"] = [
        RNA_COFOLD.get((row.gene, row.site_short), np.nan) for row in df.itertuples()
    ]
    df["mmpbsa"] = [
        MM_PBSA.get((row.gene, row.site_short), np.nan) for row in df.itertuples()
    ]
    df["cluster"] = np.where(
        df.gene.eq("LETM1"),
        "LETM1-3l",
        np.where(df.site_short.str.startswith("3l"), "NSD2-3l", "NSD2-3s"),
    )
    return df


def _ols_pred(X: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, float, float]:
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ beta
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot else np.nan
    rho = float(stats.spearmanr(pred, y).statistic) if len(y) >= 3 else np.nan
    return beta, pred, r2, rho


def _pred_vs_obs(ax, pred, y, sd, sites, xlabel: str, note: str, cluster: str) -> None:
    ax.errorbar(pred, y, yerr=sd, fmt="o", color=ACCENT, ecolor=MUTED,
                elinewidth=0.8, capsize=2, markersize=7, zorder=3)
    lo = float(min(pred.min(), y.min()) - 0.2)
    hi = float(max(pred.max(), y.max()) + 0.2)
    ax.plot([lo, hi], [lo, hi], color=LINE, lw=1, ls="--")
    for lab, x0, y0 in zip(sites, pred, y):
        ax.annotate(lab, (x0, y0), textcoords="offset points", xytext=(4, 4),
                    fontsize=7, color=MUTED)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Luciferase fold vs control")
    ax.text(0.98, 0.98, cluster, transform=ax.transAxes, ha="right", va="top",
            fontsize=10, fontweight="bold")
    if note:
        ax.text(0.02, 0.02, note, transform=ax.transAxes, ha="left", va="bottom",
                fontsize=8, color=MUTED)


def _fit_linear_ed(energy: np.ndarray, dist: np.ndarray, y: np.ndarray):
    X = np.column_stack([energy, dist, np.ones(len(y))])
    return _ols_pred(X, y)


def _loo_pred(energy: np.ndarray, dist: np.ndarray, y: np.ndarray) -> np.ndarray:
    n = len(y)
    loo = np.empty(n)
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        Xtr = np.column_stack([energy[mask], dist[mask], np.ones(n - 1)])
        beta, *_ = np.linalg.lstsq(Xtr, y[mask], rcond=None)
        loo[i] = beta[0] * energy[i] + beta[1] * dist[i] + beta[2]
    return loo


def _r2(y: np.ndarray, pred: np.ndarray) -> float:
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    return 1 - ss_res / ss_tot if ss_tot else np.nan


def plot_energy_plus_dist_clusters(df: pd.DataFrame) -> None:
    """Per cluster: RNAup / MM-PBSA / RNAcofold as linear aE+wd+b (predicted vs observed)."""
    table = _site_distance_table(df)
    clusters = ["LETM1-3l", "NSD2-3l", "NSD2-3s"]
    energies = [
        ("RNAup", "rnaup_dG_total"),
        ("MM/PBSA", "mmpbsa"),
        ("RNAcofold", "rna_cofold"),
    ]
    fig = plt.figure(figsize=(12.2, 13.2))
    gs = fig.add_gridspec(4, 3, height_ratios=[1.0, 1.0, 1.0, 1.15], hspace=0.38, wspace=0.28)
    axes = [[fig.add_subplot(gs[r, c]) for c in range(3)] for r in range(3)]
    rows_txt = []
    print("\n=== linear aE+wd+b (predicted vs observed) ===")
    for row, cluster in enumerate(clusters):
        g = table[table.cluster.eq(cluster)].reset_index(drop=True)
        y = g["fold_mean"].to_numpy(float)
        d = g["dist_3p"].to_numpy(float)
        sd = g["fold_sd"].fillna(0).to_numpy(float)
        sites = g["site_short"].tolist()
        n = len(y)
        for col, (label, ecol) in enumerate(energies):
            E = g[ecol].to_numpy(float)
            _, pred, r2, rho = _fit_linear_ed(E, d, y)
            loo = _loo_pred(E, d, y)
            loo_r2 = _r2(y, loo)
            loo_mae = float(np.mean(np.abs(y - loo)))
            _pred_vs_obs(
                axes[row][col], pred, y, sd, sites,
                f"{label} predicted fold",
                "",
                cluster,
            )
            rows_txt.append([
                cluster,
                label,
                str(n),
                f"{r2:.3f}",
                f"{loo_r2:.3f}",
                f"{loo_mae:.3f}",
            ])
            print(
                f"{cluster:10s} {label:10s}  r²={r2:.3f}  LOO r²={loo_r2:.3f}  "
                f"LOO MAE={loo_mae:.3f}  ρ={rho:.3f}"
            )
    ax_tab = fig.add_subplot(gs[3, :])
    ax_tab.axis("off")
    col_labels = ["Cluster", "Energy", "n", "R²", "LOO R²", "LOO MAE"]
    tbl = ax_tab.table(
        cellText=rows_txt,
        colLabels=col_labels,
        cellLoc="center",
        bbox=[0.0, 0.22, 1.0, 0.78],
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)
    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor(LINE)
        if r == 0:
            cell.set_facecolor("#e7e5e4")
            cell.set_text_props(fontweight="bold", color=INK)
        else:
            cell.set_facecolor("white")
            cell.set_text_props(color=INK)
    ax_tab.text(
        0.0, 0.08,
        "Linear model  fold = aE + w d + b.  "
        "LOO: refit on n−1 sites, predict the held-out site.  "
        "LOO R² uses the full-sample mean in the denominator.",
        transform=ax_tab.transAxes, ha="left", va="bottom", fontsize=8, color=MUTED,
    )
    _save(fig, "eng_c3_energy_plus_dist_clusters")


def write_readme(df: pd.DataFrame) -> None:
    rho_l, p_l = stats.spearmanr(
        df.loc[df.gene == "LETM1", "rnaup_dG_total"],
        df.loc[df.gene == "LETM1", "fold_mean"],
    )
    rho_n, p_n = stats.spearmanr(
        df.loc[df.gene == "NSD2", "rnaup_dG_total"],
        df.loc[df.gene == "NSD2", "fold_mean"],
    )
    text = f"""# Engineering figures

Generated by `simulation/plot_engineering_figures.py` from
`data/retrospective_study/retrospective_scores.tsv`.

## Files

| File | Cycle | Role |
|------|-------|------|
| `eng_c1_blast_overconservative` | 1 | Conceptual: any-BLAST-hit rule is too broad |
| `eng_c2_utr_strand_variants` | 2 | Conceptual: UTR gate, strand, optional variants, readable table |
| `eng_c3_rnaup_scheme` | 3 | Conceptual: RNAup ΔG decomposition |
| `eng_c3_scatter_letm1_nsd2` | 3 | Quantitative: ΔG vs luciferase fold |
| `eng_c3_rank_mismatch` | 3 | Quantitative: ΔG order ≠ activity order |
| `eng_c3_energy_terms` | 3 | Quantitative: individual energy terms still fail |
| `eng_c3_energy_plus_dist_clusters` | 3 | RNAup / MM-PBSA / RNAcofold linear E+d; in-sample and LOO table |
| `eng_c35_filter_not_predictor` | 3.5 | Conceptual: filter path vs rejected predictor path |

Each stem is written as `.svg`, `.pdf`, and `.png`.

## Quantitative data notes

- Luciferase: 20 ng/well blocks from the retrospective Excel sheets; fold = sample/control per replicate, then mean (± SD across replicates).
- LETM1: all nine `3l-1`…`3l-9` tiles (n=16 replicates pooled across 20 ng dates).
- NSD2: tiles present in the fluc sheet (`3l-3/4/8/9`, `3s-f`…`3s-i`; 11 wells/tile). `3s-e` omitted (plasmid not obtained).
- RNAup: `bsst filter --sites … --offtarget-min-length 100 --rnaup` so BLAST does not drop tiles before scoring (retrospective correlation needs every measured tile).
- Spearman (ΔG_total vs fold; hypothesis “more negative → higher fold” needs ρ < 0):
  LETM1 ρ={rho_l:.2f} (p={p_l:.3f}); NSD2 ρ={rho_n:.2f} (p={p_n:.3f}).

Regenerate::

```bash
../.venv/bin/python -m simulation.plot_engineering_figures
```
"""
    path = OUT / "README.md"
    path.write_text(text)
    print("wrote", path.relative_to(ROOT))


def main() -> None:
    _style()
    df = load_scores()
    plot_scatter(df)
    plot_rank_mismatch(df)
    plot_energy_terms(df)
    plot_c3a_rnaup_scheme()
    plot_c35_filter_not_predictor()
    plot_c1_overconservative()
    plot_c2_utr_strand_variants()
    plot_energy_plus_dist_clusters(df)
    write_readme(df)


if __name__ == "__main__":
    main()
