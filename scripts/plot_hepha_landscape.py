#!/usr/bin/env python3
"""Illustrative 3D landscape of the joint HEPHA score in mRNA_model-2.pdf §3.5.

    S_HEPHA(n, ipTM) ∝ C_ring_eff(n) P_match(n) exp(-ΔG_BD / RT) Φ_ED(ipTM)

ΔG_BD is held constant (no sequence data). Φ_ED is the capture–release
form P(1−P) ∝ k_on k_off / (k_on + k_off)², with k_off = k0 exp(−β ipTM).
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap, LightSource
from scipy.stats import chi2, spearmanr

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
OUT = ROOT / "docs" / "figures" / "engineering"

LETM1_MRNA_NT = 5371.0
LETM1_UTR5_NT = 209.0
LETM1_KUHN_NT = 4.0
N5A_CONTOUR = LETM1_UTR5_NT / LETM1_KUHN_NT

# Visualisation-only polymer / kinetics (not fitted).
N_RING = 120.0         # Kuhn segments on the closed mRNA
N5A = 20.0             # 5′ → AUG; match peak sits mid-axis
KUHN_NM = 2.0
RMIN_NM, RMAX_NM = 1.0, 8.0
KAPPA = 0.012          # k_match / (2 RT)
DG_OVER_RT = 0.0       # constant affinity; overall scale only
KON_R = 1.0            # k_on [R]
BETA = 5.0             # k_off = K0_OFF exp(-BETA * ipTM)
K0_OFF = KON_R * np.exp(BETA * 0.75)  # peak at ipTM ≈ 0.75


def c_ring_eff(n_b3: np.ndarray) -> np.ndarray:
    n_seg = n_b3 + N5A
    n_tot = N_RING
    n_seg = np.clip(n_seg, 1e-6, n_tot - 1e-6)
    s2 = (KUHN_NM ** 2) * n_seg * (n_tot - n_seg) / n_tot
    penc = np.clip(
        chi2.cdf(3.0 * RMAX_NM ** 2 / s2, df=3)
        - chi2.cdf(3.0 * RMIN_NM ** 2 / s2, df=3),
        0.0,
        1.0,
    )
    v_m3 = (4.0 / 3.0 * np.pi * (RMAX_NM ** 3 - RMIN_NM ** 3)) * 1e-27
    na = 6.02214076e23
    return penc / (na * v_m3)


def p_match(n_b3: np.ndarray) -> np.ndarray:
    return np.exp(-KAPPA * (N5A - n_b3) ** 2)


def phi_ed(iptm: np.ndarray) -> np.ndarray:
    koff = K0_OFF * np.exp(-BETA * iptm)
    return KON_R * koff / (KON_R + koff) ** 2


def plot_illustrative_landscape() -> None:
    n = np.linspace(1.0, 40.0, 161)
    iptm = np.linspace(0.20, 1.00, 161)
    N, I = np.meshgrid(n, iptm)

    s_bd = c_ring_eff(n) * p_match(n) * np.exp(-DG_OVER_RT)
    phi = phi_ed(iptm)
    Z = np.outer(phi, s_bd)
    Z = Z / Z.max()

    j = int(np.argmax(s_bd))
    i = int(np.argmax(phi))
    n_star, iptm_star = float(n[j]), float(iptm[i])
    z_star = float(Z[i, j])

    cmap = LinearSegmentedColormap.from_list(
        "hepha",
        ["#0c4a6e", "#0f766e", "#65a30d", "#fbbf24", "#fff7ed"],
    )
    ls = LightSource(azdeg=310, altdeg=42)
    rgb = ls.shade(Z, cmap=cmap, vert_exag=0.55, blend_mode="soft")

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.labelsize": 11,
        "figure.facecolor": "white",
    })
    fig = plt.figure(figsize=(9.4, 5.6))
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(
        N, I, Z,
        facecolors=rgb,
        rstride=2, cstride=2,
        linewidth=0, antialiased=True, shade=False,
        zorder=1,
    )
    ax.contour(
        N, I, Z,
        zdir="z", offset=0.0,
        levels=np.linspace(0.15, 1.0, 8),
        cmap=cmap, linewidths=0.7, alpha=0.85,
        zorder=0,
    )
    ax.scatter(
        [n_star], [iptm_star], [z_star],
        c="#fff7ed", edgecolors="#1c1917", linewidths=0.9,
        s=70, depthshade=False, zorder=10,
    )
    ax.plot(
        [n_star, n_star], [iptm_star, iptm_star], [0, z_star],
        color="#1c1917", lw=0.9, ls=":", zorder=4,
    )

    ax.set_xlabel(r"$n_{B3}$ (Kuhn)", labelpad=12)
    ax.set_ylabel("ipTM", labelpad=10)
    ax.set_zlabel(r"$S_{\mathrm{HEPHA}}$", labelpad=8)
    ax.set_xlim(n.min(), n.max())
    ax.set_ylim(iptm.min(), iptm.max())
    ax.set_zlim(0.0, 1.08)
    ax.view_init(elev=24, azim=-48)
    ax.xaxis.pane.set_facecolor("#fafaf9")
    ax.yaxis.pane.set_facecolor("#fafaf9")
    ax.zaxis.pane.set_facecolor("#f8fafc")
    ax.xaxis.pane.set_edgecolor("#e7e5e4")
    ax.yaxis.pane.set_edgecolor("#e7e5e4")
    ax.zaxis.pane.set_edgecolor("#e7e5e4")
    ax.set_title("Joint BD–ED design landscape (illustrative)", pad=10)
    ax.text2D(
        0.02, 0.98,
        rf"$n^\ast={n_star:.1f},\ \mathrm{{ipTM}}^\ast={iptm_star:.2f}$"
        "\n"
        rf"$n_{{5A}}={N5A:.0f},\ N={N_RING:.0f}$",
        transform=ax.transAxes,
        color="#44403c",
        fontsize=9,
        va="top",
    )
    mappable = cm.ScalarMappable(cmap=cmap)
    mappable.set_clim(0, 1)
    cb = fig.colorbar(mappable, ax=ax, shrink=0.55, pad=0.10, aspect=18)
    cb.set_label(r"$S_{\mathrm{HEPHA}}$ (normalized)")
    fig.subplots_adjust(left=0.04, right=0.90, bottom=0.06, top=0.90)

    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        path = OUT / f"eng_c3_hepha_landscape.{ext}"
        fig.savefig(path, dpi=220 if ext == "png" else None, facecolor="white")
        print("wrote", path.relative_to(ROOT))
    plt.close(fig)


def _r2(y: np.ndarray, pred: np.ndarray) -> float:
    tot = float(np.sum((y - y.mean()) ** 2))
    return 1.0 - float(np.sum((y - pred) ** 2)) / tot if tot else np.nan


def _ols_pos(score: np.ndarray, y: np.ndarray):
    """y = a score + b with a ≥ 0 (formula says S_BD scales fold, not anti-scales)."""
    x = np.column_stack([score, np.ones(len(y))])
    beta, *_ = np.linalg.lstsq(x, y, rcond=None)
    if beta[0] < 0:
        beta = np.array([0.0, float(y.mean())])
        pred = np.full(len(y), float(y.mean()))
        return beta, pred, False
    return beta, x @ beta, True


def _loo_ab(score: np.ndarray, y: np.ndarray) -> np.ndarray:
    n = len(y)
    out = np.empty(n)
    for i in range(n):
        m = np.ones(n, dtype=bool)
        m[i] = False
        beta, _, _ = _ols_pos(score[m], y[m])
        out[i] = beta[0] * score[i] + beta[1]
    return out


def fit_letm1_sbd():
    from scripts.plot_engineering_figures import _site_distance_table, load_scores
    from scripts.try_closedloop_occupancy import c_ring_eff_molar

    g = _site_distance_table(load_scores())
    g = g[g.cluster == "LETM1-3l"].copy()
    g["n_b3"] = g["dist_3p"].to_numpy(float) / LETM1_KUHN_NT
    g = g.sort_values("n_b3")
    n = g["n_b3"].to_numpy(float)
    dist = g["dist_3p"].to_numpy(float)
    y = g["fold_mean"].to_numpy(float)
    sd = g["fold_sd"].fillna(0).to_numpy(float)
    sites = g["site_short"].to_numpy()
    c = c_ring_eff_molar(dist, LETM1_MRNA_NT, utr5_nt=LETM1_UTR5_NT)

    best = None
    for n5a in np.linspace(40.0, 200.0, 161):
        for kap in np.concatenate([[0.0], np.geomspace(1e-5, 0.08, 36)]):
            p = np.exp(-kap * (n5a - n) ** 2)
            score = c * p
            beta, pred, ok = _ols_pos(score, y)
            if not ok:
                continue
            r2 = _r2(y, pred)
            cand = dict(n5a=float(n5a), kappa=float(kap), beta=beta, pred=pred,
                        score=score, r2=r2, c=c, p=p)
            if best is None or r2 > best["r2"]:
                best = cand
    if best is None:
        raise RuntimeError("no positive-slope LETM1 fit")
    best["loo"] = _loo_ab(best["score"], y)
    best["r2_loo"] = _r2(y, best["loo"])
    best["n"] = n
    best["y"] = y
    best["sd"] = sd
    best["sites"] = sites
    best["dist"] = dist
    best["rho"] = float(spearmanr(best["score"], y).statistic)
    return best


def plot_letm1_fit() -> None:
    from scripts.try_closedloop_occupancy import c_ring_eff_molar

    fit = fit_letm1_sbd()
    n, y, sd, sites = fit["n"], fit["y"], fit["sd"], fit["sites"]
    n5a, kap = fit["n5a"], fit["kappa"]
    a, b0 = fit["beta"]

    n_grid = np.linspace(28.0, 78.0, 400)
    dist_grid = n_grid * LETM1_KUHN_NT
    c_grid = c_ring_eff_molar(dist_grid, LETM1_MRNA_NT, utr5_nt=LETM1_UTR5_NT)
    p_grid = np.exp(-kap * (n5a - n_grid) ** 2)
    y_grid = a * c_grid * p_grid + b0

    p_contour = np.exp(-0.012 * (N5A_CONTOUR - n_grid) ** 2)
    s_contour = c_grid * p_contour
    # show shape only (may anti-correlate)
    s_c = (s_contour - s_contour.min()) / (s_contour.max() - s_contour.min() + 1e-12)
    s_c = y.min() + s_c * (y.max() - y.min())

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.labelsize": 11,
        "figure.facecolor": "white",
        "axes.facecolor": "#fafaf9",
        "axes.edgecolor": "#1c1917",
        "axes.spines.top": False,
        "axes.spines.right": False,
    })
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2))

    ax = axes[0]
    ax.plot(n_grid, y_grid, color="#0f766e", lw=2.0, zorder=2,
            label=r"fit $a\,C_{\mathrm{eff}}P_{\mathrm{match}}+b$")
    ax.plot(n_grid, s_c, color="#a8a29e", lw=1.2, ls="--", zorder=1,
            label=fr"contour $n_{{5A}}={N5A_CONTOUR:.0f}$ (shape)")
    ax.axvline(N5A_CONTOUR, color="#d6d3d1", lw=1, ls=":")
    ax.errorbar(n, y, yerr=sd, fmt="o", color="#b45309", ecolor="#78716c",
                elinewidth=0.8, capsize=2, markersize=7, zorder=3)
    for lab, x0, y0 in zip(sites, n, y):
        ax.annotate(lab, (x0, y0), textcoords="offset points",
                    xytext=(4, 4), fontsize=8, color="#44403c")
    ax.set_xlabel(r"$n_{B3}$ (Kuhn)  $= d_{3'}/4$")
    ax.set_ylabel("luciferase fold (20 ng)")
    ax.set_title("LETM1-3l")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.set_xlim(28, 78)
    ax.text(
        0.98, 0.03,
        rf"$n_{{5A}}^{{\mathrm{{eff}}}}={n5a:.0f}$ (contour {N5A_CONTOUR:.0f})"
        "\n"
        rf"$\kappa={kap:.4f},\ R^2={fit['r2']:.2f},\ R^2_{{\mathrm{{LOO}}}}={fit['r2_loo']:.2f}$",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=8, color="#44403c",
    )

    ax = axes[1]
    lo = float(min(fit["pred"].min(), y.min()) - 0.3)
    hi = float(max(fit["pred"].max(), y.max()) + 0.3)
    ax.plot([lo, hi], [lo, hi], color="#d6d3d1", lw=1, ls="--")
    ax.errorbar(fit["pred"], y, yerr=sd, fmt="o", color="#b45309",
                ecolor="#78716c", elinewidth=0.8, capsize=2, markersize=7)
    for lab, x0, y0 in zip(sites, fit["pred"], y):
        ax.annotate(lab, (x0, y0), textcoords="offset points",
                    xytext=(4, 3), fontsize=8, color="#44403c")
    ax.set_xlabel("predicted fold")
    ax.set_ylabel("observed fold")
    ax.set_title("in-sample")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal", adjustable="box")

    fig.suptitle(
        r"LETM1 fit of $S_{BD}\propto C^{\mathrm{ring}}_{\mathrm{eff}}P_{\mathrm{match}}$  ($a\geq 0$)",
        fontsize=12, y=1.02,
    )
    fig.tight_layout()
    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        path = OUT / f"eng_c3_hepha_letm1_fit.{ext}"
        fig.savefig(path, dpi=200 if ext == "png" else None, facecolor="white",
                    bbox_inches="tight", pad_inches=0.2)
        print("wrote", path.relative_to(ROOT))
    print(
        f"LETM1: n5A_eff={n5a:.1f} kappa={kap:.5f} a={a:.3g} b={b0:.3f} "
        f"R2={fit['r2']:.3f} LOO={fit['r2_loo']:.3f} rho={fit['rho']:.2f}"
    )
    plt.close(fig)


if __name__ == "__main__":
    plot_letm1_fit()
