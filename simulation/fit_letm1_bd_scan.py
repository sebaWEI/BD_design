#!/usr/bin/env python3
"""Fit LETM1 1018 BD-scan folds with GaussianRing / ConstrainedRing.

Uses fold_vs_blank from data/1018_bd_scan.csv. Excludes plate controls and
dead wells: blank, 5s, -4, -5, -7, -9.

Site coordinate: distance from the 40 nt window 3′ end to the major polyA
site on LETM1-3l.gbk (PAS = 5371). Transcript length = 5371 nt.
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
    KAPPA_GRID,
    KUHN_NT_GRID,
    N5A_SCALE,
    FitResult,
    loo_free,
    loo_pos,
    ols_free,
    ols_pos,
    r2_score,
    scores_constrained,
    scores_ring,
)

DATA = ROOT / "data"
OUT_FIG = ROOT / "docs" / "figures" / "model"
OUT_TAB = DATA

EXCLUDE = {"blank", "5s", "-4", "-5", "-7", "-9"}
MRNA_NT = 5371.0  # LETM1-3l.gbk
PAS = 5371

# pekinghsc/docs/wiki-plot-spec.md
WIKI = {
    "blue": "#169bb7",
    "green": "#059669",
    "green_deep": "#047857",
    "ink": "#1a3348",
    "muted": "#7a8896",
    "line": "#dfe7ed",
    "paper": "#f6f9fa",
    "white": "#ffffff",
    "note": "#0e7490",
}


def _style() -> None:
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Inter", "Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 10,
        "axes.titlesize": 13,
        "axes.titleweight": "medium",
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


def load_table() -> pd.DataFrame:
    df = pd.read_csv(DATA / "1018_bd_scan.csv")
    df["sample"] = df["sample"].astype(str)
    keep = df[~df["sample"].isin(EXCLUDE)].copy()
    keep = keep[keep["n_replicates"].fillna(0).astype(int) > 0]
    keep = keep[keep["mrna_start"].notna()]
    keep["dist_3p_nt"] = keep["dist_to_pas_nt"].astype(float)
    keep["fold_mean"] = keep["fold_vs_blank"].astype(float)
    keep["fold_sd"] = keep["fold_vs_blank_sd"].astype(float)
    keep["mrna_nt"] = MRNA_NT
    keep["site_mid"] = (keep["mrna_start"].astype(float) + keep["mrna_end"].astype(float)) / 2.0
    return keep.sort_values("dist_3p_nt", ascending=False).reset_index(drop=True)


def fit(tab: pd.DataFrame) -> tuple[FitResult, FitResult, pd.DataFrame]:
    dist = tab["dist_3p_nt"].to_numpy(float)
    y = tab["fold_mean"].to_numpy(float)
    mrna_nt = float(tab["mrna_nt"].iloc[0])
    n5a_grid = np.unique(np.round(N5A_SCALE, 4))

    best_ring: FitResult | None = None
    ring_pack = None
    con_cands: list[tuple[FitResult, np.ndarray, np.ndarray, np.ndarray]] = []

    for kuhn in KUHN_NT_GRID:
        for n_5a in n5a_grid:
            if n_5a <= 0:
                continue
            s = scores_ring(dist, mrna_nt, kuhn, n_5a)
            if not np.all(np.isfinite(s)) or s.max() - s.min() < 1e-18:
                continue
            beta, pred = ols_free(s, y)
            r2 = r2_score(y, pred)
            loo = loo_free(s, y)
            cand = FitResult(
                backbone="LETM1", model="GaussianRing",
                kuhn_nt=kuhn, n_5a=float(n_5a), kappa=0.0, k_match=0.0,
                a=float(beta[0]), b=float(beta[1]),
                r2=r2, r2_loo=r2_score(y, loo),
                rho=float(np.corrcoef(s, y)[0, 1]),
                n_sites=len(y), a_nonneg=bool(beta[0] >= 0),
            )
            if best_ring is None or r2 > best_ring.r2:
                best_ring = cand
                ring_pack = (s, pred, loo)

            for kappa in KAPPA_GRID:
                s = scores_constrained(dist, mrna_nt, kuhn, n_5a, kappa)
                if not np.all(np.isfinite(s)) or s.max() - s.min() < 1e-18:
                    continue
                beta, pred, ok = ols_pos(s, y)
                if not ok:
                    continue
                r2 = r2_score(y, pred)
                loo = loo_pos(s, y)
                cand = FitResult(
                    backbone="LETM1", model="ConstrainedRing",
                    kuhn_nt=kuhn, n_5a=float(n_5a), kappa=float(kappa),
                    k_match=float(kappa * 2.0 * 0.616),
                    a=float(beta[0]), b=float(beta[1]),
                    r2=r2, r2_loo=r2_score(y, loo),
                    rho=float(np.corrcoef(s, y)[0, 1]),
                    n_sites=len(y), a_nonneg=True,
                )
                con_cands.append((cand, s, pred, loo))

    if best_ring is None or ring_pack is None or not con_cands:
        raise RuntimeError("no valid LETM1 fit")

    r2_max = max(c.r2 for c, *_ in con_cands)
    near = [
        (c, s, p, l) for c, s, p, l in con_cands
        if c.r2 >= r2_max - 0.03 and c.r2_loo >= 0.0
    ]
    if not near:
        near = [(c, s, p, l) for c, s, p, l in con_cands if c.r2 >= r2_max - 0.03]
    if not near:
        near = con_cands

    def _key(item: tuple[FitResult, np.ndarray, np.ndarray, np.ndarray]):
        c = item[0]
        sigma_n = 1.0 / np.sqrt(max(2.0 * c.kappa, 1e-12))
        return (c.r2_loo, sigma_n, -c.kappa)

    best_con, s_c, pred_c, loo_c = max(near, key=_key)
    out = tab[["sample", "mrna_start", "mrna_end", "dist_3p_nt", "site_mid",
               "fold_mean", "fold_sd"]].copy()
    out["score_ring"] = ring_pack[0]
    out["pred_ring"] = ring_pack[1]
    out["score_constrained"] = s_c
    out["pred_constrained"] = pred_c
    out["loo_constrained"] = loo_c
    return best_ring, best_con, out


def plot(con: FitResult, tab: pd.DataFrame) -> None:
    """Wiki figure per pekinghsc/docs/wiki-plot-spec.md."""
    d_lo = float(tab["dist_3p_nt"].min()) - 40.0
    d_hi = float(tab["dist_3p_nt"].max()) + 40.0
    grid = np.linspace(max(d_lo, 40.0), d_hi, 400)
    y_fit = con.a * scores_constrained(grid, MRNA_NT, con.kuhn_nt, con.n_5a, con.kappa) + con.b

    # full-width twin panel
    fig, (ax, axp) = plt.subplots(
        1, 2, figsize=(7.2, 4.2),
        gridspec_kw={"width_ratios": [1.45, 1.0], "wspace": 0.32},
    )

    # --- left: positional profile ---
    ax.yaxis.grid(True, which="major")
    ax.xaxis.grid(False)
    ax.axhline(1.0, color=WIKI["muted"], lw=0.9, ls=(0, (3, 2.5)), zorder=1)
    ax.fill_between(
        grid, 1.0, y_fit, where=(y_fit >= 1.0),
        color=WIKI["blue"], alpha=0.10, interpolate=True, zorder=2, linewidth=0,
    )
    ax.plot(
        grid, y_fit, color=WIKI["blue"], lw=2.0, zorder=3,
        label="Constrained ring",
    )
    ax.errorbar(
        tab["dist_3p_nt"], tab["fold_mean"], yerr=tab["fold_sd"].fillna(0),
        fmt="none", ecolor=WIKI["green"], elinewidth=1.0, capsize=2.5,
        capthick=1.0, zorder=4,
    )
    ax.scatter(
        tab["dist_3p_nt"], tab["fold_mean"],
        s=36, color=WIKI["green"], edgecolors=WIKI["ink"], linewidths=0.6,
        zorder=5, label="Observed (mean ± SD)",
    )
    ax.set_xlim(d_hi, max(d_lo, 40.0))
    ymax = max(float(tab["fold_mean"].max() + tab["fold_sd"].max()), float(y_fit.max())) + 0.5
    ax.set_ylim(0.0, ymax)
    ax.set_xlabel("Distance to polyA (nt)")
    ax.set_ylabel("Luciferase fold vs blank")
    ax.set_title("Positional profile")
    ax.legend(loc="upper left", handlelength=1.6)

    # --- right: parity ---
    axp.yaxis.grid(True, which="major")
    axp.xaxis.grid(False)
    lo = 0.0
    hi = float(max(tab["pred_constrained"].max(), tab["fold_mean"].max()) + 0.7)
    axp.plot([lo, hi], [lo, hi], color=WIKI["muted"], lw=1.0, ls=(0, (4, 2.5)), zorder=1)
    axp.errorbar(
        tab["pred_constrained"], tab["fold_mean"], yerr=tab["fold_sd"].fillna(0),
        fmt="none", ecolor=WIKI["green"], elinewidth=1.0, capsize=2.5,
        capthick=1.0, zorder=3,
    )
    axp.scatter(
        tab["pred_constrained"], tab["fold_mean"],
        s=36, color=WIKI["green"], edgecolors=WIKI["ink"], linewidths=0.6, zorder=4,
    )
    axp.set_xlim(lo, hi)
    axp.set_ylim(lo, hi)
    axp.set_aspect("equal", adjustable="box")
    axp.set_xlabel("Predicted fold")
    axp.set_ylabel("Observed fold")
    axp.set_title("Fit")
    axp.text(
        0.04, 0.96, r"$y = x$", transform=axp.transAxes,
        ha="left", va="top", fontsize=9, color=WIKI["muted"],
    )

    fig.suptitle(
        r"Constrained-ring model on LETM1 · $TE \propto C_{\mathrm{eff}}$",
        fontsize=13, fontweight="medium", color=WIKI["ink"],
        x=0.01, ha="left", y=1.02,
    )
    _save(fig, "letm1_bd_scan_fit")


def main() -> None:
    _style()
    tab = load_table()
    print(f"sites used ({len(tab)}): {', '.join(tab['sample'].tolist())}")
    ring, con, pred = fit(tab)
    print(
        f"ConstrainedRing:  R2={con.r2:.3f} LOO={con.r2_loo:.3f} "
        f"rho={con.rho:.3f} a={con.a:.3g} "
        f"kuhn={con.kuhn_nt} n5a={con.n_5a:.2f} kappa={con.kappa:.4g} "
        f"k_match={con.k_match:.4g}"
    )
    OUT_TAB.mkdir(parents=True, exist_ok=True)
    pred_path = OUT_TAB / "fit_letm1_bd_scan_predictions.csv"
    params_path = OUT_TAB / "fit_letm1_bd_scan_params.csv"
    pred_out = pred.drop(columns=["score_ring", "pred_ring"], errors="ignore")
    pred_out.to_csv(pred_path, index=False)
    pd.DataFrame([con.__dict__]).to_csv(params_path, index=False)
    print("wrote", pred_path.relative_to(ROOT))
    print("wrote", params_path.relative_to(ROOT))
    plot(con, pred)
    _ = ring  # fitted but not plotted


if __name__ == "__main__":
    main()
