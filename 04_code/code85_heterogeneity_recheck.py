# -*- coding: utf-8 -*-
"""
code85_heterogeneity_recheck.py  (seed-deterministic)
================================================
Task: close constant-registry gap G3 (ATM activity distribution across cells) and test whether "the routes agree".

Route ledger:
  code 81: ideal ceil-mapping generative model, I(D;fate)=1.325 bits (upper-bound reference)
  code 82: homogeneous ensemble, dispersion ratio 0.01 (under-dispersed; shows pure counting noise is insufficient)
  code 84: mechanistic chain + repair stochasticity, sigma=0, dispersion ratio 8.6(10Gy)-36.0(2.5Gy) (over-dispersed at low dose)
  code 85 (this file): inject ATM gain heterogeneity g~lognormal(sigma) onto the code-84 skeleton,
    with sigma not fitted to the dispersion ratio but anchored to independent literature:
      Lopez-Pujol 2023 (Int J Radiat Biol): per-cell gamma-H2AX foci counts,
      variance/mean (Fano) slope 4.07-4.75 (four lymphocyte subtypes, 1-2 Gy).
      In this model the foci equivalent is F = g*N_DSB, Poisson compounded with lognormal:
      Fano = 1 + CV_g^2*mu, mu(2Gy)=70 -> Fano~4 requires CV_g~0.21.
    The dispersion-ratio landing point at this sigma is a zero-fit prediction of literature value 5.8 (Loewer 2010).

Usage:
  python code85_heterogeneity_recheck.py scan          # sigma scan, N=2000
  python code85_heterogeneity_recheck.py final 0.25    # final check at given sigma*, N=6000
"""
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["svg.fonttype"] = "none"
plt.rcParams["font.size"] = 9.5

SEED = 20260925
ROOT = Path(__file__).resolve().parent.parent
OUT_JSON = ROOT / "结果" / "代码85_异质性复核_结果.json"
OUT_PNG = ROOT / "结果" / "代码85_异质性复核_图.png"
OUT_SVG = ROOT / "结果" / "代码85_异质性复核_图.svg"

# ---- constants kept byte-identical to code 84 ----
DOSES = np.array([0.1, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0])
K0_DSB = 35.0
R_FAST = 0.35
FRAC_SLOW = 0.10
T_HALF_SLOW = 20.0
R_SLOW = np.log(2) / T_HALF_SLOW
T_PERIOD = 5.5
JIT_IPI = 0.085
W_PULSE = 3.5
W = 48.0
N_MAX = 12
Q21, TAU21 = 100.0, 10.0
QPU, TAUPU, PUMA_NTHR = 100.0, 4.0, 3
TH_PUMA, R0_FATE, TH_P21 = 800.0, 0.15, 1500.0
LAM0 = 0.003
A_MIN = 10.0
LIT_DISPERSION = 5.8
LIT_FANO_BAND = (4.07, 4.75)   # Lopez-Pujol 2023 four-subtype slopes
CODE81_END2END = 1.325
CODE84_END2END = 0.817


def mi_discrete(x, y):
    x = np.asarray(x); y = np.asarray(y)
    Kx = int(x.max()) + 1; Ky = int(y.max()) + 1
    N = len(x)
    cxy = np.bincount(x * Ky + y, minlength=Kx * Ky).reshape(Kx, Ky) / N
    cx = cxy.sum(1, keepdims=True); cy = cxy.sum(0, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        term = cxy * np.log2(cxy / (cx @ cy))
    return float(np.nansum(term))


def qbin(v, nbin):
    v = np.asarray(v, dtype=float)
    edges = np.unique(np.quantile(v, np.linspace(0, 1, nbin + 1)))
    if len(edges) <= 2:
        return np.zeros(len(v), dtype=int)
    return np.clip(np.digitize(v, edges[1:-1]), 0, len(edges) - 2)


def first_pulse_time(taus, lam0, rng):
    """Byte-identical to code 84; lam0 is the per-cell effective hazard rate."""
    if len(taus) == 0:
        return None
    ts = np.sort(taus)
    u = rng.exponential(1.0)
    lam_total = lam0 * ts.sum()
    if lam_total <= u:
        return None
    cs = np.cumsum(ts)
    n = len(ts)
    lam_knots = lam0 * (cs + ts * (n - np.arange(n) - 1))
    lam_knots = np.concatenate([[0.0], lam_knots])
    t_knots = np.concatenate([[0.0], ts])
    k = min(int(np.searchsorted(lam_knots, u, side="right")) - 1, n - 1)
    slope = lam0 * (n - k)
    if slope <= 0:
        return None
    t1 = t_knots[k] + (u - lam_knots[k]) / slope
    return t1 if t1 <= W else None


def simulate_dose(D, n, rng, sigma_ln):
    """Code-84 skeleton + per-cell ATM gain g~lognormal(sigma) (multiplies only the trigger hazard rate)."""
    n_dsb = rng.poisson(K0_DSB * D, n)
    if sigma_ln > 0:
        g = np.exp(sigma_ln * rng.standard_normal(n) - 0.5 * sigma_ln ** 2)
    else:
        g = np.ones(n)
    lam_cell = LAM0 * g
    taus_list = []
    for i in range(n):
        ni = n_dsb[i]
        if ni == 0:
            taus_list.append(np.empty(0)); continue
        slow = rng.random(ni) < FRAC_SLOW
        rates = np.where(slow, R_SLOW, R_FAST)
        taus_list.append(rng.exponential(1.0 / rates))
    t1 = np.full(n, np.nan)
    times = np.full((n, N_MAX), np.nan)
    A_at = lambda taus, t: int(np.sum(taus > t))
    for i in range(n):
        tp = first_pulse_time(taus_list[i], lam_cell[i], rng)
        if tp is None:
            continue
        t1[i] = tp
        times[i, 0] = tp
        for j in range(1, N_MAX):
            interval = T_PERIOD * (1.0 + JIT_IPI * rng.standard_normal())
            interval = max(interval, 0.3 * T_PERIOD)
            tp = times[i, j - 1] + interval
            if tp > W:
                break
            if A_at(taus_list[i], tp) < A_MIN:
                break
            times[i, j] = tp
    active = ~np.isnan(times)
    n_pulse = active.sum(1)
    responded = n_pulse > 0

    def integrate_tau(q, tau, mask):
        with np.errstate(invalid="ignore"):
            contrib = np.where(mask, q * tau * (1.0 - np.exp(-(W - times) / tau)), 0.0)
        return np.nansum(contrib, axis=1)
    pidx = np.arange(N_MAX)[None, :]
    m_p21 = integrate_tau(Q21, TAU21, active)
    puma_mask = active & (pidx >= PUMA_NTHR - 1)
    m_puma = integrate_tau(QPU, TAUPU, puma_mask)
    ratio = m_puma / (m_p21 + 1.0)
    fate = np.where((m_puma >= TH_PUMA) & (ratio >= R0_FATE), 2,
                    np.where(m_p21 >= TH_P21, 1, 0))
    return dict(n_dsb=n_dsb, g=g, t1=t1, times=times, n_pulse=n_pulse,
                responded=responded, m_p21=m_p21, m_puma=m_puma, fate=fate)


def run_ensemble(sigma_ln, n_cells, seed):
    rng = np.random.default_rng(seed)
    parts = []
    for D in DOSES:
        p = simulate_dose(D, n_cells, rng, sigma_ln)
        parts.append(p)
    out = {}
    # anchor quantities
    p10 = parts[-1]
    out["count_10Gy"] = float(p10["n_pulse"].mean())
    t1v = p10["t1"]; t1v = t1v[~np.isnan(t1v)]
    out["t1_10Gy_mean"] = float(t1v.mean())
    out["resp_frac"] = [float(p["responded"].mean()) for p in parts]
    # Fano anchor (2 Gy equivalent foci count F = g*N_DSB)
    p2 = parts[4]
    F = p2["g"] * p2["n_dsb"]
    out["fano_2Gy"] = float(F.var() / F.mean())
    out["cv_g_implied_by_fano4"] = float(np.sqrt((4.0 - 1.0) / (K0_DSB * 2.0)))
    # dispersion ratio (responder-conditional convention, same as code 84)
    disp = {}
    for di, D in enumerate(DOSES):
        p = parts[di]
        t1d = p["t1"]; t1d = t1d[~np.isnan(t1d)]
        ipi_sd_cell = []
        tt = p["times"]
        for i in range(tt.shape[0]):
            tv = tt[i][~np.isnan(tt[i])]
            if len(tv) >= 3:
                ipi_sd_cell.append(np.std(np.diff(tv)))
        if len(t1d) > 10 and ipi_sd_cell:
            disp[f"{D}Gy"] = dict(t1_mean=float(t1d.mean()), t1_sd=float(t1d.std()),
                                  ipi_sd=float(np.mean(ipi_sd_cell)),
                                  dispersion_ratio=float(np.std(t1d) / np.mean(ipi_sd_cell)))
    out["dispersion"] = disp
    # end-to-end information (same convention as code 81/84)
    d_idx = np.concatenate([np.full(n_cells, i) for i in range(len(DOSES))])
    E = {k: np.concatenate([p[k] for p in parts]) for k in
         ["n_dsb", "t1", "times", "n_pulse", "responded", "m_p21", "m_puma", "fate"]}
    y5 = E["fate"]
    out["I_D_fate"] = mi_discrete(d_idx, y5)
    return out


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "scan"
    log = lambda *a: print(*a, flush=True)

    if mode == "scan":
        sigmas = [0.0, 0.1, 0.2, 0.3, 0.4]
        N = 2000
        log(f"[code85] scan mode: sigma in {sigmas}, N={N}, seed {SEED}")
        results = {}
        for s in sigmas:
            log(f"  σ={s} ...")
            results[f"{s}"] = run_ensemble(s, N, SEED + int(round(s * 1000)))
            r = results[f"{s}"]
            dr = {k: round(v["dispersion_ratio"], 2) for k, v in r["dispersion"].items()}
            log(f"    <N>10={r['count_10Gy']:.2f}, t1@10={r['t1_10Gy_mean']:.2f}h, "
                f"Fano(2Gy)={r['fano_2Gy']:.2f}, I(D;fate)={r['I_D_fate']:.3f}, dispersion_ratio={dr}")
        payload = dict(meta=dict(script="代码85_异质性复核.py", mode="scan", N=N, seed=SEED,
                                 lit_dispersion=LIT_DISPERSION,
                                 lit_fano_band=LIT_FANO_BAND),
                       scan=results)
        # keep existing final results if present
        if OUT_JSON.exists():
            try:
                old = json.loads(OUT_JSON.read_text(encoding="utf-8"))
                if "final" in old:
                    payload["final"] = old["final"]
            except Exception:
                pass
        OUT_JSON.write_text(json.dumps(payload, ensure_ascii=False, indent=2),
                            encoding="utf-8")
        log("JSON written:", OUT_JSON)
        make_figure(results, None)
    else:
        sigma_star = float(sys.argv[2]) if len(sys.argv) > 2 else 0.25
        N = 6000
        log(f"[code85] final mode: sigma*={sigma_star}, N={N}, seed {SEED + 777}")
        fin = run_ensemble(sigma_star, N, SEED + 777)
        dr = {k: round(v["dispersion_ratio"], 2) for k, v in fin["dispersion"].items()}
        log(f"  <N>10={fin['count_10Gy']:.2f}, t1@10={fin['t1_10Gy_mean']:.2f}h, "
            f"Fano(2Gy)={fin['fano_2Gy']:.2f}, I(D;fate)={fin['I_D_fate']:.3f}, dispersion_ratio={dr}")
        payload = {}
        if OUT_JSON.exists():
            payload = json.loads(OUT_JSON.read_text(encoding="utf-8"))
        payload["final"] = dict(sigma_star=sigma_star, N=N, **fin)
        payload["meta"]["mode"] = "scan+final"
        OUT_JSON.write_text(json.dumps(payload, ensure_ascii=False, indent=2),
                            encoding="utf-8")
        log("JSON updated:", OUT_JSON)
        make_figure(payload.get("scan"), payload["final"])
    log("[code85] done.")


def make_figure(scan, final):
    fig, axes = plt.subplots(2, 2, figsize=(12.5, 8.6))
    if scan:
        sigmas = sorted(float(s) for s in scan.keys())
        # (a) Fano anchor
        ax = axes[0, 0]
        fanos = [scan[f"{s}"]["fano_2Gy"] for s in sigmas]
        ax.plot(sigmas, fanos, "o-", color="#33527a")
        ax.axhspan(LIT_FANO_BAND[0], LIT_FANO_BAND[1], color="#0b6b3a", alpha=0.15)
        ax.axhline(4.0, color="#0b6b3a", ls="--", lw=1.0)
        ax.text(0.02, 4.05, "literature Fano 4.07--4.75 (Lopez-Pujol 2023)",
                fontsize=8, color="#0b6b3a")
        cv_imp = scan[f"{sigmas[0]}"]["cv_g_implied_by_fano4"]
        ax.axvline(cv_imp, color="#8c1d18", ls=":", lw=1.0)
        ax.text(cv_imp + 0.01, ax.get_ylim()[0] + 0.3,
                f"implied CV_g={cv_imp:.2f}", fontsize=8, color="#8c1d18")
        ax.set_xlabel(r"$\sigma_{ln}$ of ATM gain $g$"); ax.set_ylabel("Fano factor at 2 Gy")
        ax.set_title("(a) Independent anchor: foci count Fano factor", fontsize=10.5)
        # (b) dispersion ratio
        ax = axes[0, 1]
        for D, c in [(2.5, "#d9863d"), (5.0, "#33527a"), (10.0, "#0b6b3a")]:
            ys = [scan[f"{s}"]["dispersion"].get(f"{D}Gy", {}).get("dispersion_ratio", np.nan)
                  for s in sigmas]
            ax.plot(sigmas, ys, "o-", color=c, label=f"{D} Gy")
        ax.axhline(LIT_DISPERSION, color="#8c1d18", ls="--", lw=1.0)
        ax.text(0.02, LIT_DISPERSION + 0.6, "literature 5.8", fontsize=8, color="#8c1d18")
        if final:
            ax.plot([final["sigma_star"]],
                    [final["dispersion"].get("10.0Gy", final["dispersion"].get("10Gy", {})).get("dispersion_ratio", np.nan)],
                    "*", ms=16, color="#8c1d18")
        ax.set_yscale("log")
        ax.set_xlabel(r"$\sigma_{ln}$"); ax.set_ylabel("dispersion ratio SD($t_1$)/SD(IPI)")
        ax.legend(fontsize=8)
        ax.set_title("(b) Dispersion ratio vs gain heterogeneity (zero-fit prediction)", fontsize=10.5)
        # (c) anchor stability
        ax = axes[1, 0]
        cnt = [scan[f"{s}"]["count_10Gy"] for s in sigmas]
        t1m = [scan[f"{s}"]["t1_10Gy_mean"] for s in sigmas]
        ax.plot(sigmas, cnt, "o-", color="#33527a", label="<N> at 10 Gy")
        ax.axhspan(0, 7, color="#33527a", alpha=0.08)
        ax.set_xlabel(r"$\sigma_{ln}$"); ax.set_ylabel("pulse count", color="#33527a")
        ax2 = ax.twinx()
        ax2.plot(sigmas, t1m, "s--", color="#c05640", label="$t_1$ at 10 Gy")
        ax2.axhspan(1.5, 3.0, color="#c05640", alpha=0.10)
        ax2.set_ylabel("first-pulse delay (h)", color="#c05640")
        h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="center right")
        ax.set_title("(b2) Anchor stability: counting law and $t_1$ vs $\\sigma_{ln}$", fontsize=10.5)
        # (d) information
        ax = axes[1, 1]
        info = [scan[f"{s}"]["I_D_fate"] for s in sigmas]
        ax.plot(sigmas, info, "o-", color="#33527a")
        ax.axhline(CODE81_END2END, color="#0b6b3a", ls=":", lw=1.0)
        ax.text(0.02, CODE81_END2END + 0.02, "code 81: 1.325", fontsize=8, color="#0b6b3a")
        ax.axhline(CODE84_END2END, color="#d9863d", ls=":", lw=1.0)
        ax.text(0.02, CODE84_END2END - 0.06, "code 84: 0.817", fontsize=8, color="#d9863d")
        if final:
            ax.plot([final["sigma_star"]], [final["I_D_fate"]], "*", ms=16, color="#8c1d18")
        ax.set_xlabel(r"$\sigma_{ln}$"); ax.set_ylabel("I(D; fate) (bits)")
        ax.set_title("(c) End-to-end information vs gain heterogeneity", fontsize=10.5)
    fig.suptitle("Code 85: ATM gain heterogeneity, route reconciliation (G3 closure test)",
                 fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(OUT_PNG, dpi=200, bbox_inches="tight")
    fig.savefig(OUT_SVG, bbox_inches="tight")
    print("figure written:", OUT_PNG, flush=True)


if __name__ == "__main__":
    main()
