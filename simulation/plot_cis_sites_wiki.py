#!/usr/bin/env python3
"""Wiki-style combined cis ConstrainedRing figure (Moderna + TMSB10).

Layout::

    [Moderna profile] [TMSB10 profile]   ← same plot type on one row
    [Moderna fit]     [TMSB10 fit]

Y-axis labels only on the leftmost panel of each row.

Usage (BD_design root)::

    ../.venv/bin/python -m simulation.plot_cis_sites_wiki
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from simulation.fit_cis_sites import (  # noqa: E402
    OUT_FIG,
    UTR3_NT,
    UTR5_NT,
    CDS_NT,
    scores_constrained,
)

DATA = ROOT / "data"

WIKI = {
    "blue": "#169bb7",
    "green": "#059669",
    "ink": "#1a3348",
    "muted": "#7a8896",
    "line": "#dfe7ed",
    "white": "#ffffff",
}

ORDER = ("Moderna", "TMSB10")


def _style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Inter", "Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 10,
        "axes.titlesize": 13,
        "axes.titleweight": "normal",
        "axes.titlelocation": "left",
        "axes.labelsize": 11,
        "axes.labelcolor": WIKI["ink"],
        "axes.edgecolor": WIKI["line"],
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.facecolor": WIKI["white"],
        "figure.facecolor": WIKI["white"],
        "xtick.color": WIKI["ink"],
        "ytick.color": WIKI["ink"],
        "xtick.labelsize": 9.5,
        "ytick.labelsize": 9.5,
        "text.color": WIKI["ink"],
        "grid.color": WIKI["line"],
        "grid.linewidth": 0.7,
        "legend.frameon": False,
        "legend.fontsize": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.facecolor": WIKI["white"],
    })


def _save(fig: plt.Figure, stem: str) -> None:
    OUT_FIG.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        path = OUT_FIG / f"{stem}.{ext}"
        fig.savefig(path, facecolor=WIKI["white"], bbox_inches="tight", pad_inches=0.2)
        print("wrote", path.relative_to(ROOT))
    plt.close(fig)


def _load_con(backbone: str) -> tuple[pd.Series, pd.DataFrame]:
    params = pd.read_csv(DATA / "cis_fit_params.csv")
    row = params[(params.backbone == backbone) & (params.model == "ConstrainedRing")].iloc[0]
    tab = pd.read_csv(DATA / "cis_fit_predictions.csv")
    tab = tab[tab.backbone == backbone].copy()
    return row, tab


def _draw_profile(
    ax: plt.Axes,
    backbone: str,
    con: pd.Series,
    tab: pd.DataFrame,
    *,
    ylabel: bool,
    show_legend: bool,
) -> None:
    utr3 = float(UTR3_NT[backbone])
    mrna = UTR5_NT + CDS_NT + utr3
    x = tab["utr_pos_1based"].to_numpy(float)
    y = tab["fold_mean"].to_numpy(float)
    yerr = tab["fold_sd"].fillna(0).to_numpy(float)

    x_max = float(x.max()) + 8.0
    grid_pos = np.linspace(0.0, min(x_max, utr3 - 1.0), 300)
    grid_dist = utr3 - grid_pos
    y_fit = (
        float(con.a)
        * scores_constrained(grid_dist, mrna, float(con.kuhn_nt), float(con.n_5a), float(con.kappa))
        + float(con.b)
    )

    ax.yaxis.grid(True, which="major")
    ax.xaxis.grid(False)
    ax.axhline(1.0, color=WIKI["muted"], lw=0.9, ls=(0, (3, 2.5)), zorder=1)
    ax.fill_between(
        grid_pos, 1.0, y_fit, where=(y_fit >= 1.0),
        color=WIKI["blue"], alpha=0.10, interpolate=True, zorder=2, linewidth=0,
    )
    ax.plot(grid_pos, y_fit, color=WIKI["blue"], lw=2.0, zorder=3, label="Constrained ring")
    ax.errorbar(
        x, y, yerr=yerr, fmt="none",
        ecolor=WIKI["green"], elinewidth=1.0, capsize=2.5, capthick=1.0, zorder=4,
    )
    ax.scatter(
        x, y, s=36, color=WIKI["green"], edgecolors=WIKI["ink"], linewidths=0.6,
        zorder=5, label="Observed (mean ± SD)",
    )
    ymax = max(float((y + yerr).max()), float(y_fit.max())) + 0.25
    ax.set_xlim(-2.0, x_max)
    ax.set_ylim(0.0, ymax)
    ax.set_xlabel("3'UTR position (nt)")
    ax.set_ylabel("Luciferase fold" if ylabel else "")
    ax.set_title(backbone)
    if show_legend:
        ax.legend(loc="upper right", handlelength=1.6)


def _draw_fit(
    ax: plt.Axes,
    backbone: str,
    tab: pd.DataFrame,
    *,
    ylabel: bool,
) -> None:
    y = tab["fold_mean"].to_numpy(float)
    yerr = tab["fold_sd"].fillna(0).to_numpy(float)
    pred = tab["pred_constrained"].to_numpy(float)

    ax.yaxis.grid(True, which="major")
    ax.xaxis.grid(False)
    lo = 0.0
    hi = float(max(pred.max(), (y + yerr).max()) + 0.35)
    ax.plot([lo, hi], [lo, hi], color=WIKI["muted"], lw=1.0, ls=(0, (4, 2.5)), zorder=1)
    ax.errorbar(
        pred, y, yerr=yerr, fmt="none",
        ecolor=WIKI["green"], elinewidth=1.0, capsize=2.5, capthick=1.0, zorder=3,
    )
    ax.scatter(
        pred, y, s=36, color=WIKI["green"], edgecolors=WIKI["ink"], linewidths=0.6, zorder=4,
    )
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Predicted fold")
    ax.set_ylabel("Observed fold" if ylabel else "")
    ax.text(
        0.04, 0.96, r"$y = x$", transform=ax.transAxes,
        ha="left", va="top", fontsize=9, color=WIKI["muted"],
    )


def main() -> None:
    _style()
    loaded = {b: _load_con(b) for b in ORDER}

    fig, axes = plt.subplots(
        2, 2, figsize=(7.2, 6.6),
        gridspec_kw={"wspace": 0.22, "hspace": 0.40},
    )

    # Row 0: positional profiles
    for j, backbone in enumerate(ORDER):
        con, tab = loaded[backbone]
        _draw_profile(
            axes[0, j], backbone, con, tab,
            ylabel=(j == 0),
            show_legend=(j == 0),
        )

    # Row 1: parity fits
    for j, backbone in enumerate(ORDER):
        _, tab = loaded[backbone]
        _draw_fit(axes[1, j], backbone, tab, ylabel=(j == 0))

    fig.suptitle(
        r"Constrained-ring model on cis UTRs · $TE \propto C_{\mathrm{eff}}$",
        fontsize=13, fontweight="normal", color=WIKI["ink"],
        x=0.01, ha="left", y=1.01,
    )
    _save(fig, "cis_fit_wiki")


if __name__ == "__main__":
    main()
