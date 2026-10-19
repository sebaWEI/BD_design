#!/usr/bin/env python3
r"""Fit cis luciferase folds with GaussianRing / ConstrainedRing (main.tex).

Cis (same molecule)::

    TE_cis(n) ∝ C_eff^ring(n)   or   C_eff^constrained(n) = C_eff^ring · P_match

Site coordinate: distance from the insertion point to the 3′ end of the
(native) 3′UTR, converted to Kuhn segments \(n_{B3}\).

Pure ring proximity often *anti*-correlates with these cis profiles; the
arm-match term in ConstrainedRing is what produces a mid-UTR peak.

Usage (BD_design root)::

    ../.venv/bin/python -m simulation.fit_cis_sites
"""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from model.main import ConstrainedRing, GaussianRing  # noqa: E402

DATA = ROOT / "data"
OUT_FIG = ROOT / "docs" / "figures" / "model"
OUT_TAB = DATA

INK = "#1c1917"
MUTED = "#57534e"
ACCENT = "#0f766e"
ACCENT2 = "#b45309"
ACCENT3 = "#1d4ed8"
PANEL = "#fafaf9"
LINE = "#d6d3d1"

# Transcript geometry from SnapGene annotations (nt).
# Moderna 5′UTR on both reporters; native 3′UTR ignores SINEUP insert.
UTR5_NT = 57.0
CDS_NT = 1653.0
UTR3_NT = {"TMSB10": 248.0, "Moderna": 111.0}  # native

KUHN_NT_GRID = (2.0, 3.0, 4.0, 5.0, 8.0)
# n_5a in Kuhn segments: physical 5′UTR/kuhn plus a wider effective grid
N5A_SCALE = np.concatenate([
    np.array([UTR5_NT]) / np.array(KUHN_NT_GRID),
    np.linspace(5.0, 100.0, 49),
])
KAPPA_GRID = np.concatenate([[0.0], np.geomspace(1e-4, 0.5, 36)])  # = k_match/(2 RT)
RT = 0.616


def _style() -> None:
    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
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
    OUT_FIG.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf", "svg"):
        path = OUT_FIG / f"{stem}.{ext}"
        fig.savefig(path, facecolor="white", bbox_inches="tight", pad_inches=0.25)
        print("wrote", path.relative_to(ROOT))
    plt.close(fig)


def r2_score(y: np.ndarray, pred: np.ndarray) -> float:
    tot = float(np.sum((y - y.mean()) ** 2))
    if tot <= 0:
        return np.nan
    return 1.0 - float(np.sum((y - pred) ** 2)) / tot


def ols_free(score: np.ndarray, y: np.ndarray):
    """y = a·score + b (a unrestricted; used to diagnose ring-only)."""
    x = np.column_stack([score, np.ones(len(y))])
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    return beta, x @ beta


def ols_pos(score: np.ndarray, y: np.ndarray):
    """y = a·score + b with a ≥ 0 (cis TE scales with C_eff)."""
    beta, pred = ols_free(score, y)
    if beta[0] < 0:
        beta = np.array([0.0, float(y.mean())])
        return beta, np.full(len(y), float(y.mean())), False
    return beta, pred, True


def loo_pos(score: np.ndarray, y: np.ndarray) -> np.ndarray:
    n = len(y)
    out = np.empty(n)
    for i in range(n):
        m = np.ones(n, dtype=bool)
        m[i] = False
        beta, _, _ = ols_pos(score[m], y[m])
        out[i] = beta[0] * score[i] + beta[1]
    return out


def loo_free(score: np.ndarray, y: np.ndarray) -> np.ndarray:
    n = len(y)
    out = np.empty(n)
    for i in range(n):
        m = np.ones(n, dtype=bool)
        m[i] = False
        beta, _ = ols_free(score[m], y[m])
        out[i] = beta[0] * score[i] + beta[1]
    return out


def load_sites() -> pd.DataFrame:
    folds = pd.read_csv(DATA / "cis_luciferase.csv")
    sites = folds[folds["site"] != "Ctrl"].copy()
    sites["utr_pos_1based"] = sites["utr_pos_1based"].astype(float)
    sites["utr3_nt"] = sites["backbone"].map(UTR3_NT)
    # Distance from site point to native 3′ end (nt).
    sites["dist_3p_nt"] = sites["utr3_nt"] - sites["utr_pos_1based"]
    sites["mrna_nt"] = UTR5_NT + CDS_NT + sites["utr3_nt"]
    return sites


@dataclass
class FitResult:
    backbone: str
    model: str
    kuhn_nt: float
    n_5a: float
    kappa: float
    k_match: float
    a: float
    b: float
    r2: float
    r2_loo: float
    rho: float
    n_sites: int
    a_nonneg: bool


def scores_ring(dist_3p: np.ndarray, mrna_nt: float, kuhn_nt: float, n_5a: float) -> np.ndarray:
    n_b3 = dist_3p / kuhn_nt
    n_total = mrna_nt / kuhn_nt
    n_ring = n_b3 + n_5a
    ring = GaussianRing(n_total=n_total, kuhn_length_nm=kuhn_nt * 0.5)
    return np.asarray(ring.effective_concentration(n_ring), dtype=float)


def scores_constrained(
    dist_3p: np.ndarray,
    mrna_nt: float,
    kuhn_nt: float,
    n_5a: float,
    kappa: float,
) -> np.ndarray:
    n_b3 = dist_3p / kuhn_nt
    n_total = mrna_nt / kuhn_nt
    k_match = kappa * 2.0 * RT
    model = ConstrainedRing.from_params(
        n_total=n_total,
        n_5a=n_5a,
        k_match=k_match,
        kuhn_length_nm=kuhn_nt * 0.5,
        rt=RT,
    )
    return np.asarray(model.effective_concentration(n_b3), dtype=float)


def fit_backbone(df: pd.DataFrame, backbone: str) -> tuple[FitResult, FitResult, pd.DataFrame]:
    g = df[df.backbone == backbone].copy()
    # mean per site for fitting the positional curve
    mean = (
        g.groupby(["site", "utr_pos_1based", "dist_3p_nt", "mrna_nt"], as_index=False)
        .agg(fold_mean=("fold", "mean"), fold_sd=("fold", "std"), n=("fold", "count"))
        .sort_values("utr_pos_1based")
    )
    dist = mean["dist_3p_nt"].to_numpy(float)
    y = mean["fold_mean"].to_numpy(float)
    mrna_nt = float(mean["mrna_nt"].iloc[0])

    best_ring: FitResult | None = None
    ring_pack = None
    con_cands: list[tuple[FitResult, np.ndarray, np.ndarray, np.ndarray]] = []

    n5a_grid = np.unique(np.round(N5A_SCALE, 4))

    for kuhn in KUHN_NT_GRID:
        for n_5a in n5a_grid:
            if n_5a <= 0:
                continue
            # --- ring only: free OLS (often a < 0 on these cis profiles) ---
            s = scores_ring(dist, mrna_nt, kuhn, n_5a)
            if not np.all(np.isfinite(s)) or s.max() - s.min() < 1e-18:
                continue
            beta, pred = ols_free(s, y)
            r2 = r2_score(y, pred)
            loo = loo_free(s, y)
            cand = FitResult(
                backbone=backbone, model="GaussianRing",
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
                    backbone=backbone, model="ConstrainedRing",
                    kuhn_nt=kuhn, n_5a=float(n_5a), kappa=float(kappa),
                    k_match=float(kappa * 2.0 * RT),
                    a=float(beta[0]), b=float(beta[1]),
                    r2=r2, r2_loo=r2_score(y, loo),
                    rho=float(np.corrcoef(s, y)[0, 1]),
                    n_sites=len(y), a_nonneg=True,
                )
                con_cands.append((cand, s, pred, loo))

    if best_ring is None or ring_pack is None or not con_cands:
        raise RuntimeError(f"no valid fit for {backbone}")

    # Prefer a wider P_match peak among near-best fits.
    # P_match = exp(−κ Δn²) ⇒ σ_n = 1/√(2κ); smaller κ = broader peak.
    # Moderna Site2–3 plateau: prioritize width; TMSB10: mild width, keep R².
    r2_slack = {"Moderna": 0.08, "TMSB10": 0.02}.get(backbone, 0.05)
    loo_floor = {"Moderna": 0.52, "TMSB10": 0.70}.get(backbone, 0.50)
    r2_max = max(c.r2 for c, *_ in con_cands)
    near = [
        (c, s, p, l)
        for c, s, p, l in con_cands
        if c.r2 >= r2_max - r2_slack and c.r2_loo >= loo_floor
    ]
    if not near:
        near = [(c, s, p, l) for c, s, p, l in con_cands if c.r2 >= r2_max - r2_slack]
    if not near:
        near = con_cands

    prefer_width = backbone == "Moderna"

    def _con_key(item: tuple[FitResult, np.ndarray, np.ndarray, np.ndarray]):
        c = item[0]
        sigma_n = 1.0 / np.sqrt(max(2.0 * c.kappa, 1e-12))
        if prefer_width:
            return (sigma_n, c.r2_loo, c.r2)
        return (c.r2_loo, sigma_n, -c.kappa)

    best_con, s_c, pred_c, loo_c = max(near, key=_con_key)
    con_pack = (s_c, pred_c, loo_c)

    plot_tab = mean[["site", "utr_pos_1based", "dist_3p_nt", "fold_mean", "fold_sd"]].copy()
    plot_tab["score_ring"] = ring_pack[0]
    plot_tab["pred_ring"] = ring_pack[1]
    plot_tab["score_constrained"] = con_pack[0]
    plot_tab["pred_constrained"] = con_pack[1]
    return best_ring, best_con, plot_tab


def plot_backbone(
    backbone: str,
    ring: FitResult,
    con: FitResult,
    tab: pd.DataFrame,
    reps: pd.DataFrame,
) -> None:
    g = reps[reps.backbone == backbone]
    fig, axes = plt.subplots(1, 2, figsize=(9.8, 4.2))

    # Left: fold vs UTR position
    ax = axes[0]
    ax.axhline(1.0, color=LINE, lw=1.0, ls=":")
    for site, sub in g.groupby("site"):
        x = float(sub["utr_pos_1based"].iloc[0])
        ax.scatter(
            np.full(len(sub), x), sub["fold"].to_numpy(),
            s=28, color=ACCENT2, alpha=0.55, zorder=3, edgecolors="none",
        )
    ax.errorbar(
        tab["utr_pos_1based"], tab["fold_mean"], yerr=tab["fold_sd"].fillna(0),
        fmt="o", color=ACCENT2, ecolor=MUTED, elinewidth=0.8, capsize=2,
        markersize=7, zorder=4, label="mean ± SD",
    )
    # smooth curves on a fine UTR grid
    utr3 = UTR3_NT[backbone]
    mrna = UTR5_NT + CDS_NT + utr3
    grid_pos = np.linspace(0.0, float(tab["utr_pos_1based"].max()) + 5.0, 200)
    grid_pos = grid_pos[grid_pos < utr3 - 1]
    grid_dist = utr3 - grid_pos

    s_r = scores_ring(grid_dist, mrna, ring.kuhn_nt, ring.n_5a)
    s_c = scores_constrained(grid_dist, mrna, con.kuhn_nt, con.n_5a, con.kappa)
    ring_lab = "Gaussian ring (free OLS)"
    if ring.a < 0:
        ring_lab += ", $a<0$"
    ax.plot(grid_pos, ring.a * s_r + ring.b, color=ACCENT3, lw=1.6, ls="--", label=ring_lab)
    ax.plot(grid_pos, con.a * s_c + con.b, color=ACCENT, lw=2.2, label="Constrained ring ($a\\geq0$)")
    ax.set_xlabel(r"3′UTR position (nt from CDS|UTR junction)")
    ax.set_ylabel("luciferase fold")
    ax.set_title(f"{backbone} cis positional profile")
    ax.legend(frameon=False, fontsize=8)
    ax.text(
        0.02, 0.02,
        rf"ring $R^2={ring.r2:.2f}$ (LOO {ring.r2_loo:.2f}, $\rho={ring.rho:.2f}$)"
        "\n"
        rf"constrained $R^2={con.r2:.2f}$ (LOO {con.r2_loo:.2f})"
        "\n"
        rf"$n_{{5A}}^{{\mathrm{{con}}}}={con.n_5a:.1f},\ \kappa={con.kappa:.4f}$"
        rf"  ($\sigma_n={1.0 / np.sqrt(max(2*con.kappa, 1e-12)):.1f}$)",
        transform=ax.transAxes, ha="left", va="bottom", fontsize=8, color=MUTED,
    )

    # Right: observed vs predicted (constrained)
    ax = axes[1]
    lo = float(min(tab["pred_constrained"].min(), tab["fold_mean"].min()) - 0.15)
    hi = float(max(tab["pred_constrained"].max(), tab["fold_mean"].max()) + 0.15)
    ax.plot([lo, hi], [lo, hi], color=LINE, lw=1, ls="--")
    ax.errorbar(
        tab["pred_constrained"], tab["fold_mean"], yerr=tab["fold_sd"].fillna(0),
        fmt="o", color=ACCENT, ecolor=MUTED, elinewidth=0.8, capsize=2, markersize=7,
    )
    for _, row in tab.iterrows():
        ax.annotate(
            str(row["site"]).replace("Site", "S"),
            (row["pred_constrained"], row["fold_mean"]),
            textcoords="offset points", xytext=(4, 3), fontsize=8, color=MUTED,
        )
    ax.set_xlabel("predicted fold (constrained)")
    ax.set_ylabel("observed fold")
    ax.set_title("Constrained-ring fit")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")

    fig.suptitle(
        rf"{backbone}: cis $TE \propto C_{{\mathrm{{eff}}}}$  (Kuhn={con.kuhn_nt:.0f} nt)",
        fontsize=12, y=1.02,
    )
    fig.tight_layout()
    _save(fig, f"cis_fit_{backbone.lower()}")


def main() -> None:
    sites = load_sites()
    con_rows: list[FitResult] = []
    pred_frames: list[pd.DataFrame] = []
    for backbone in ("TMSB10", "Moderna"):
        _ring, con, tab = fit_backbone(sites, backbone)
        con_rows.append(con)
        frame = tab.copy()
        frame.insert(0, "backbone", backbone)
        frame = frame.drop(columns=[c for c in frame.columns if c.startswith("score_ring") or c.startswith("pred_ring")], errors="ignore")
        pred_frames.append(frame)
        print(
            f"\n=== {backbone} ===\n"
            f"  ConstrainedRing:  R2={con.r2:.3f} LOO={con.r2_loo:.3f} "
            f"rho={con.rho:.3f} a={con.a:.3g} "
            f"kuhn={con.kuhn_nt} n5a={con.n_5a:.2f} kappa={con.kappa:.4g} "
            f"k_match={con.k_match:.4g}"
        )

    params_path = OUT_TAB / "cis_fit_params.csv"
    pred_path = OUT_TAB / "cis_fit_predictions.csv"
    pd.DataFrame([r.__dict__ for r in con_rows]).to_csv(params_path, index=False)
    pd.concat(pred_frames, ignore_index=True).to_csv(pred_path, index=False)
    print("\nwrote", params_path.relative_to(ROOT))
    print("wrote", pred_path.relative_to(ROOT))


if __name__ == "__main__":
    main()
