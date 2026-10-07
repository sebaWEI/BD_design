#!/usr/bin/env python3
"""Try occupancy × 3D closed-loop contact on the retrospective tiles.

Also the trans site score from model/mRNA_model.pdf §3.2::

    Site_score(n) ∝ C_ring_eff(n) * exp(-ΔG_BD(n) / RT)

Usage (repo root)::

    ../.venv/bin/python -m simulation.try_closedloop_occupancy
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2, spearmanr

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from simulation.plot_engineering_figures import (  # noqa: E402
    ACCENT,
    LINE,
    MUTED,
    _save,
    _site_distance_table,
    _style,
    load_scores,
)

RT = 0.616  # kcal/mol at 37 °C
NA = 6.02214076e23
CLUSTERS = ["LETM1-3l", "NSD2-3l", "NSD2-3s"]
ENERGIES = [
    ("RNAup", "rnaup_dG_total", True),
    ("MM/PBSA", "mmpbsa", False),
    ("RNAcofold", "rna_cofold", True),
]

# Closed-loop Gaussian ring (mRNA_model.pdf §3).
# Kuhn length ~2 nm, ~0.5 nm/nt → 4 nt/segment.
KUHN_NT = 4.0
KUHN_NM = 2.0
RMIN_NM = 1.0
RMAX_NM = 8.0
# spliced mRNA: Ensembl lookup (LETM1-201 5371 nt; NSD2-218 7560 nt).
# n = (dist_3p + 5′UTR) / Kuhn  — arc BD→polyA→cap→AUG on the closed ring.
MRNA_NT = {"LETM1-3l": 5371, "NSD2-3l": 7560, "NSD2-3s": 7560}
UTR5_NT = {"LETM1-3l": 209, "NSD2-3l": 182, "NSD2-3s": 182}


def r2_score(y: np.ndarray, pred: np.ndarray) -> float:
    ss_res = float(np.sum((y - pred) ** 2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    return 1 - ss_res / ss_tot if ss_tot else np.nan


def ols(X: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return beta, X @ beta


def contact(L: np.ndarray) -> np.ndarray:
    return np.power(np.clip(L, 1.0, None), -1.5)


def theta_langmuir(energy: np.ndarray, g_half: float) -> np.ndarray:
    x = np.clip((energy - g_half) / RT, -60, 60)
    return 1.0 / (1.0 + np.exp(x))


def fit_linear_el(energy: np.ndarray, L: np.ndarray, y: np.ndarray):
    X = np.column_stack([energy, L, np.ones(len(y))])
    beta, pred = ols(X, y)
    return beta, pred, r2_score(y, pred)


def fit_contact_only(L: np.ndarray, y: np.ndarray):
    X = np.column_stack([contact(L), np.ones(len(y))])
    beta, pred = ols(X, y)
    return beta, pred, r2_score(y, pred)


def fit_theta_P(energy: np.ndarray, L: np.ndarray, y: np.ndarray, g_grid: np.ndarray):
    P = contact(L)
    best = None
    ones = np.ones(len(y))
    for g_half in g_grid:
        th = theta_langmuir(energy, g_half)
        score = th * P
        if score.max() - score.min() < 1e-15:
            continue
        X = np.column_stack([score, ones])
        beta, pred = ols(X, y)
        r2 = r2_score(y, pred)
        cand = dict(g_half=g_half, beta=beta, pred=pred, r2=r2, theta=th, P=P)
        if best is None or r2 > best["r2"]:
            best = cand
    return best


def loo_linear(energy: np.ndarray, L: np.ndarray, y: np.ndarray) -> np.ndarray:
    n = len(y)
    X = np.column_stack([energy, L, np.ones(n)])
    out = np.empty(n)
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        b, *_ = np.linalg.lstsq(X[mask], y[mask], rcond=None)
        out[i] = X[i] @ b
    return out


def loo_contact_only(L: np.ndarray, y: np.ndarray) -> np.ndarray:
    n = len(y)
    X = np.column_stack([contact(L), np.ones(n)])
    out = np.empty(n)
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        b, *_ = np.linalg.lstsq(X[mask], y[mask], rcond=None)
        out[i] = X[i] @ b
    return out


def loo_theta_P(
    energy: np.ndarray, L: np.ndarray, y: np.ndarray, g_grid: np.ndarray
) -> np.ndarray:
    n = len(y)
    out = np.empty(n)
    P = contact(L)
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        fit = fit_theta_P(energy[mask], L[mask], y[mask], g_grid)
        if fit is None:
            out[i] = y[mask].mean()
            continue
        th_i = theta_langmuir(np.array([energy[i]]), fit["g_half"])[0]
        out[i] = fit["beta"][0] * th_i * P[i] + fit["beta"][1]
    return out


def g_half_grid(energy: np.ndarray) -> np.ndarray:
    lo, hi = float(energy.min()), float(energy.max())
    pad = max(3.0, hi - lo)
    return np.linspace(lo - pad, hi + pad, 81)


def ring_s2_nm2(n_seg: np.ndarray, n_total: float, kuhn_nm: float) -> np.ndarray:
    """<R²> = b² n (N-n) / N  for a Gaussian ring."""
    n = np.clip(n_seg, 1e-6, n_total - 1e-6)
    return (kuhn_nm ** 2) * n * (n_total - n) / n_total


def p_ring_enc(s2: np.ndarray, rmin: float, rmax: float) -> np.ndarray:
    """P(Rmin ≤ R ≤ Rmax) from χ²_3 of the Cartesian Gaussians."""
    hi = chi2.cdf(3.0 * rmax ** 2 / s2, df=3)
    lo = chi2.cdf(3.0 * rmin ** 2 / s2, df=3)
    return np.clip(hi - lo, 0.0, 1.0)


def c_ring_eff_molar(
    dist_3p_nt: np.ndarray,
    mrna_nt: float,
    kuhn_nt: float = KUHN_NT,
    kuhn_nm: float | None = None,
    rmin: float = RMIN_NM,
    rmax: float = RMAX_NM,
    n_bridge: float = 0.0,
    utr5_nt: float = 0.0,
) -> np.ndarray:
    if kuhn_nm is None:
        kuhn_nm = kuhn_nt * 0.5
    n_seg = (dist_3p_nt + utr5_nt) / kuhn_nt + n_bridge
    n_total = mrna_nt / kuhn_nt + n_bridge
    s2 = ring_s2_nm2(n_seg, n_total, kuhn_nm)
    penc = p_ring_enc(s2, rmin, rmax)
    v_nm3 = 4.0 / 3.0 * np.pi * (rmax ** 3 - rmin ** 3)
    v_m3 = v_nm3 * 1e-27
    return penc / (NA * v_m3)


E_REF = -58.0  # kcal/mol; only a numerical shift, cancels into OLS slope


def trans_site_score(energy, dist_3p_nt, mrna_nt, utr5_nt) -> np.ndarray:
    """Site_score ∝ C_ring_eff(n) exp(-ΔG/RT), n = L + u5."""
    c_eff = c_ring_eff_molar(dist_3p_nt, mrna_nt, utr5_nt=utr5_nt)
    return c_eff * np.exp(-(energy - E_REF) / RT)


def fit_trans(energy, L, y, mrna_nt, utr5_nt):
    s = trans_site_score(energy, L, mrna_nt, utr5_nt)
    X = np.column_stack([s, np.ones(len(y))])
    beta, pred = ols(X, y)
    return beta, pred, r2_score(y, pred), s


def loo_trans(energy, L, y, mrna_nt, utr5_nt) -> np.ndarray:
    n = len(y)
    out = np.empty(n)
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        s_tr = trans_site_score(energy[mask], L[mask], mrna_nt, utr5_nt)
        X = np.column_stack([s_tr, np.ones(n - 1)])
        b, *_ = np.linalg.lstsq(X, y[mask], rcond=None)
        s_i = trans_site_score(energy[i : i + 1], L[i : i + 1], mrna_nt, utr5_nt)[0]
        out[i] = b[0] * s_i + b[1]
    return out


KUHN_GRID = (3.0, 4.0, 8.0, 16.0)
GAMMA_GRID = (0.75, 1.0, 1.5, 2.0)
TEFF_GRID = (np.inf, 32.0, 12.0, 6.0)  # kcal/mol; inf = drop energy
BRIDGE_GRID = (0.0, 10.0)


def damped_score(energy, L, mrna_nt, utr5_nt, kuhn_nt, gamma, t_eff, n_bridge) -> np.ndarray:
    c_eff = c_ring_eff_molar(
        L, mrna_nt, kuhn_nt=kuhn_nt, n_bridge=n_bridge, utr5_nt=utr5_nt,
    )
    pos = np.power(np.clip(c_eff, 1e-30, None), gamma)
    if np.isinf(t_eff):
        return pos
    return pos * np.exp(-(energy - E_REF) / t_eff)


def fit_damped_grid(energy, L, y, mrna_nt, utr5_nt, with_energy: bool):
    teffs = TEFF_GRID if with_energy else (np.inf,)
    best = None
    for kuhn in KUHN_GRID:
        for gamma in GAMMA_GRID:
            for t_eff in teffs:
                for n_bridge in BRIDGE_GRID:
                    s = damped_score(
                        energy, L, mrna_nt, utr5_nt, kuhn, gamma, t_eff, n_bridge,
                    )
                    if not np.all(np.isfinite(s)) or s.max() - s.min() < 1e-18:
                        continue
                    X = np.column_stack([s, np.ones(len(y))])
                    beta, pred = ols(X, y)
                    r2 = r2_score(y, pred)
                    cand = dict(
                        kuhn=kuhn, gamma=gamma, t_eff=t_eff, n_bridge=n_bridge,
                        beta=beta, pred=pred, r2=r2, score=s,
                    )
                    if best is None or r2 > best["r2"]:
                        best = cand
    return best


def loo_damped(energy, L, y, mrna_nt, utr5_nt, with_energy: bool) -> np.ndarray:
    n = len(y)
    out = np.empty(n)
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        fit = fit_damped_grid(
            energy[mask], L[mask], y[mask], mrna_nt, utr5_nt, with_energy,
        )
        if fit is None:
            out[i] = y[mask].mean()
            continue
        s_i = damped_score(
            energy[i : i + 1], L[i : i + 1], mrna_nt, utr5_nt,
            fit["kuhn"], fit["gamma"], fit["t_eff"], fit["n_bridge"],
        )[0]
        out[i] = fit["beta"][0] * s_i + fit["beta"][1]
    return out


XSTAR_GRID = np.linspace(80, 450, 38)
SIGMA_GRID = (40.0, 80.0, 120.0, 180.0, 250.0)


def parallel_score(energy, L, x_star, sigma, t_eff) -> np.ndarray:
    pos = np.exp(-0.5 * ((L - x_star) / sigma) ** 2)
    if np.isinf(t_eff):
        return pos
    return pos * np.exp(-(energy - E_REF) / t_eff)


def fit_parallel_grid(energy, L, y, with_energy: bool):
    teffs = (np.inf, 32.0, 12.0) if with_energy else (np.inf,)
    best = None
    for x_star in XSTAR_GRID:
        for sigma in SIGMA_GRID:
            for t_eff in teffs:
                s = parallel_score(energy, L, x_star, sigma, t_eff)
                if s.max() - s.min() < 1e-15:
                    continue
                X = np.column_stack([s, np.ones(len(y))])
                beta, pred = ols(X, y)
                r2 = r2_score(y, pred)
                cand = dict(
                    x_star=float(x_star), sigma=sigma, t_eff=t_eff,
                    beta=beta, pred=pred, r2=r2, score=s,
                )
                if best is None or r2 > best["r2"]:
                    best = cand
    return best


def loo_parallel(energy, L, y, with_energy: bool) -> np.ndarray:
    n = len(y)
    out = np.empty(n)
    for i in range(n):
        mask = np.ones(n, dtype=bool)
        mask[i] = False
        fit = fit_parallel_grid(energy[mask], L[mask], y[mask], with_energy)
        s_i = parallel_score(
            energy[i : i + 1], L[i : i + 1],
            fit["x_star"], fit["sigma"], fit["t_eff"],
        )[0]
        out[i] = fit["beta"][0] * s_i + fit["beta"][1]
    return out


def main() -> None:
    table = _site_distance_table(load_scores())
    print("RT = {:.3f} kcal/mol (37 C)".format(RT))
    print("P_contact = L^(-3/2);  fold = a * theta(ΔG) * P + b")
    print("MM/PBSA: theta fixed at 1 (not a binding ΔG)\n")
    hdr = (
        f"{'cluster':10s} {'energy':10s} n  "
        f"{'lin_r2':>7s} {'lin_loo':>7s}  "
        f"{'P_r2':>7s} {'P_loo':>7s}  "
        f"{'thP_r2':>7s} {'thP_loo':>7s}  ΔG½  <θ>"
    )
    print(hdr)
    print("-" * len(hdr))
    for cluster in CLUSTERS:
        g = table[table.cluster.eq(cluster)].reset_index(drop=True)
        y = g["fold_mean"].to_numpy(float)
        L = g["dist_3p"].to_numpy(float)
        n = len(y)
        _, _, p_r2 = fit_contact_only(L, y)
        p_loo = r2_score(y, loo_contact_only(L, y))
        for label, col, use_theta in ENERGIES:
            E = g[col].to_numpy(float)
            _, _, lin_r2 = fit_linear_el(E, L, y)
            lin_loo = r2_score(y, loo_linear(E, L, y))
            if use_theta:
                grid = g_half_grid(E)
                fit = fit_theta_P(E, L, y, grid)
                thp_r2 = fit["r2"] if fit else np.nan
                thp_loo = r2_score(y, loo_theta_P(E, L, y, grid)) if fit else np.nan
                g_half = fit["g_half"] if fit else np.nan
                th_mean = float(fit["theta"].mean()) if fit else np.nan
                extra = f"{g_half:6.1f}  {th_mean:.2f}"
            else:
                thp_r2 = p_r2
                thp_loo = p_loo
                extra = "   n/a  1.00"
            print(
                f"{cluster:10s} {label:10s} {n}  "
                f"{lin_r2:7.3f} {lin_loo:7.3f}  "
                f"{p_r2:7.3f} {p_loo:7.3f}  "
                f"{thp_r2:7.3f} {thp_loo:7.3f}  {extra}"
            )
        print()
    print("lin = aE + wL + b;  P = a L^(-3/2) + b;  thP = a θ P + b")
    print("If all ΔG << ΔG½ then θ≈1 and thP collapses to P.")

    print("\n=== trans Site_score = C_ring_eff(n) exp(-ΔG/RT)  (mRNA_model.pdf) ===")
    print(
        f"Kuhn={KUHN_NT:.0f} nt ({KUHN_NM:.1f} nm),  "
        f"capture {RMIN_NM:.0f}–{RMAX_NM:.0f} nm,  n = (dist_3p + 5'UTR) / Kuhn,  N = mRNA/Kuhn"
    )
    hdr2 = (
        f"{'cluster':10s} {'energy':10s}  "
        f"{'trans_r2':>8s} {'trans_loo':>9s}  rho(score,fold)  "
        f"Ceff span  bind span"
    )
    print(hdr2)
    print("-" * 88)
    for cluster in CLUSTERS:
        g = table[table.cluster.eq(cluster)].reset_index(drop=True)
        y = g["fold_mean"].to_numpy(float)
        L = g["dist_3p"].to_numpy(float)
        mrna_nt = MRNA_NT[cluster]
        utr5 = UTR5_NT[cluster]
        n_arc = L + utr5
        c_eff = c_ring_eff_molar(L, mrna_nt, utr5_nt=utr5)
        print(
            f"  [{cluster}] N={mrna_nt} nt, u5={utr5} nt, "
            f"n_nt={n_arc.min():.0f}–{n_arc.max():.0f} (BD→polyA→AUG)"
        )
        for label, col, use_theta in ENERGIES:
            if not use_theta:
                continue
            E = g[col].to_numpy(float)
            _, _, r2, score = fit_trans(E, L, y, mrna_nt, utr5)
            loo_r2 = r2_score(y, loo_trans(E, L, y, mrna_nt, utr5))
            rho = float(spearmanr(score, y).statistic)
            bind = np.exp(-(E - E.mean()) / RT)
            print(
                f"{cluster:10s} {label:10s}  "
                f"{r2:8.3f} {loo_r2:9.3f}  {rho:14.3f}  "
                f"{c_eff.max()/c_eff.min():6.2f}x  "
                f"{bind.max()/bind.min():8.1f}x"
            )
        print()
    print("rho: Spearman of Site_score vs fold (model wants rho > 0).")

    print("\n=== damped energy × learnable polymer  score = C_eff^γ exp(-ΔG/T_eff) ===")
    print("grids: Kuhn nt", KUHN_GRID, " γ", GAMMA_GRID, " T_eff", TEFF_GRID, " n_bridge", BRIDGE_GRID)
    print(f"{'cluster':10s} {'energy':10s}  r2    loo    Kuhn  γ    T_eff  bridge  bind_span  rho")
    print("-" * 92)
    for cluster in CLUSTERS:
        g = table[table.cluster.eq(cluster)].reset_index(drop=True)
        y = g["fold_mean"].to_numpy(float)
        L = g["dist_3p"].to_numpy(float)
        mrna_nt = MRNA_NT[cluster]
        utr5 = UTR5_NT[cluster]
        dist_fit = fit_damped_grid(
            np.zeros_like(y), L, y, mrna_nt, utr5, with_energy=False,
        )
        dist_loo = r2_score(y, loo_damped(np.zeros_like(y), L, y, mrna_nt, utr5, False))
        te = "inf"
        print(
            f"{cluster:10s} {'(dist only)':10s}  {dist_fit['r2']:5.3f} {dist_loo:6.3f}  "
            f"{dist_fit['kuhn']:4.0f}  {dist_fit['gamma']:.2f}  {te:>5s}  "
            f"{dist_fit['n_bridge']:4.0f}     n/a     "
            f"{float(spearmanr(dist_fit['score'], y).statistic):+.3f}"
        )
        for label, col, use_theta in ENERGIES:
            if not use_theta:
                continue
            E = g[col].to_numpy(float)
            fit = fit_damped_grid(E, L, y, mrna_nt, utr5, with_energy=True)
            loo_r2 = r2_score(y, loo_damped(E, L, y, mrna_nt, utr5, True))
            te = "inf" if np.isinf(fit["t_eff"]) else f"{fit['t_eff']:.0f}"
            if np.isinf(fit["t_eff"]):
                span = 1.0
            else:
                bind = np.exp(-(E - E.mean()) / fit["t_eff"])
                span = float(bind.max() / bind.min())
            rho = float(spearmanr(fit["score"], y).statistic)
            print(
                f"{cluster:10s} {label:10s}  {fit['r2']:5.3f} {loo_r2:6.3f}  "
                f"{fit['kuhn']:4.0f}  {fit['gamma']:.2f}  {te:>5s}  "
                f"{fit['n_bridge']:4.0f}  {span:8.1f}x  {rho:+.3f}"
            )
        print()
    print("T_eff=inf drops energy. Larger T_eff damps exp(-ΔG/T) toward 1.")
    print("n_bridge = extra Kuhn segments for the protein 5′–3′ bridge.")

    print("\n=== parallel juxtaposition: P(L) = exp(-(L-x*)^2 / 2σ^2) ===")
    print("L = dist_3p; x* = 3'UTR register aligned with 5' initiation along the protein.")
    print(f"{'cluster':10s} {'energy':10s}  r2    loo     x*    σ    T_eff  rho")
    print("-" * 72)
    for cluster in CLUSTERS:
        g = table[table.cluster.eq(cluster)].reset_index(drop=True)
        y = g["fold_mean"].to_numpy(float)
        L = g["dist_3p"].to_numpy(float)
        dist_fit = fit_parallel_grid(np.zeros_like(y), L, y, with_energy=False)
        dist_loo = r2_score(y, loo_parallel(np.zeros_like(y), L, y, False))
        print(
            f"{cluster:10s} {'(dist only)':10s}  {dist_fit['r2']:5.3f} "
            f"{dist_loo:6.3f}  {dist_fit['x_star']:5.0f}  {dist_fit['sigma']:4.0f}    inf  "
            f"{float(spearmanr(dist_fit['score'], y).statistic):+.3f}"
        )
        E = g["rnaup_dG_total"].to_numpy(float)
        fit = fit_parallel_grid(E, L, y, with_energy=True)
        loo_r2 = r2_score(y, loo_parallel(E, L, y, True))
        te = "inf" if np.isinf(fit["t_eff"]) else f"{fit['t_eff']:.0f}"
        print(
            f"{cluster:10s} {'RNAup':10s}  {fit['r2']:5.3f} {loo_r2:6.3f}  "
            f"{fit['x_star']:5.0f}  {fit['sigma']:4.0f}  {te:>5s}  "
            f"{float(spearmanr(fit['score'], y).statistic):+.3f}"
        )
        print()
    print("If x* is outside the ladder, points sit on one flank: monotonic in dist_3p.")
    plot_y_vs_n(table)


def plot_y_vs_n(table) -> None:
    """Wet-lab fold and ring C_eff vs n = L+u5 and vs L/u5."""
    _style()
    fig, axes = plt.subplots(3, 4, figsize=(16.4, 10.4))
    for row, cluster in enumerate(CLUSTERS):
        g = table[table.cluster.eq(cluster)].reset_index(drop=True)
        y = g["fold_mean"].to_numpy(float)
        sd = g["fold_sd"].fillna(0).to_numpy(float)
        L = g["dist_3p"].to_numpy(float)
        u5 = float(UTR5_NT[cluster])
        n_nt = L + u5
        ratio = L / u5
        c_eff = c_ring_eff_molar(L, MRNA_NT[cluster], utr5_nt=u5)
        sites = g["site_short"].tolist()
        dG = g["rnaup_dG_total"].to_numpy(float)
        xs = [
            (n_nt, "n (nt)  BD → polyA → AUG", "n"),
            (ratio, "(BD–polyA) / (5′–AUG)", "L/u5"),
        ]
        for k, (x, xlabel, tag) in enumerate(xs):
            rho_y = float(spearmanr(x, y).statistic)
            rho_c = float(spearmanr(x, c_eff).statistic)
            ax = axes[row, 2 * k]
            ax.scatter(x, y, c=dG, cmap="coolwarm", s=42, zorder=3, edgecolors="white",
                       vmin=-63, vmax=-49)
            ax.errorbar(x, y, yerr=sd, fmt="none", ecolor=MUTED, elinewidth=0.8,
                        capsize=2, zorder=2)
            for lab, x0, y0 in zip(sites, x, y):
                ax.annotate(lab, (x0, y0), textcoords="offset points", xytext=(3, 3),
                            fontsize=7, color=MUTED)
            ax.set_xlabel(xlabel)
            ax.set_ylabel("Luciferase fold")
            ax.text(0.98, 0.98, cluster, transform=ax.transAxes, ha="right", va="top",
                    fontsize=10, fontweight="bold")
            ax.text(0.98, 0.02, f"ρ({tag}, fold)={rho_y:+.2f}", transform=ax.transAxes,
                    ha="right", va="bottom", fontsize=8, color=MUTED)

            ax2 = axes[row, 2 * k + 1]
            ax2.scatter(x, c_eff, c=dG, cmap="coolwarm", s=42, zorder=3, edgecolors="white",
                        vmin=-63, vmax=-49)
            order = np.argsort(x)
            ax2.plot(x[order], c_eff[order], color=LINE, lw=1.2, zorder=1)
            for lab, x0, y0 in zip(sites, x, c_eff):
                ax2.annotate(lab, (x0, y0), textcoords="offset points", xytext=(3, 3),
                             fontsize=7, color=MUTED)
            ax2.set_xlabel(xlabel)
            ax2.set_ylabel(r"$C^{\mathrm{ring}}_{\mathrm{eff}}$ (M)")
            ax2.text(0.98, 0.98, cluster, transform=ax.transAxes, ha="right", va="top",
                     fontsize=10, fontweight="bold")
            ax2.text(0.02, 0.02, f"ρ({tag}, C_eff)={rho_c:+.2f}  (model: −)",
                     transform=ax.transAxes, ha="left", va="bottom", fontsize=8, color=MUTED)
    fig.tight_layout()
    _save(fig, "eng_c3_y_vs_n")


if __name__ == "__main__":
    main()
