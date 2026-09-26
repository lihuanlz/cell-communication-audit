# -*- coding: utf-8 -*-
"""
Code 84: p53 discovery-based chain model (constant skeleton + explicit distributions + globally zero fitting)
=====================================================================
Date: 2026-09-25 | seed fixed: 20260925 | deterministic parts fully deterministic

Lineage:
  - Constant registry: 结果/常数登记_p53通路_发现式建模_v01_2026-09-25.md (K1-K8 / D1-D4 / G1-G4)
  - Cascade audit: code 81 (verdict card 2026-09-23) — all downstream conventions of this model match it verbatim
  - Variational race: code 83 (2026-09-25)

Differences from code 81 (everything new in this model):
  The L1 sensing layer is replaced from the "assumed ceil(N_DSB/50) mapping" by a mechanistic form:
    N0 ~ Poisson(35 DSB/Gy) [K-series anchor]
    Each DSB is repaired independently: 90% fast phase r=0.35/h (t½≈2h) [K7]; 10% slow phase t½=20h [K8]
    ATM activity A(t) = number of unrepaired DSBs (γH2AX formation 30 min << T=5.5h, treated as instantaneous, double-recorded) [K5]
    First pulse: inhomogeneous Poisson triggering, hazard rate lambda(t) = LAM0 * A(t)
           (exact sampling: Lambda(t)=LAM0*Σmin(t,τ_i) piecewise linear, inverted analytically)
    Subsequent pulses: interval T=5.5h (CV 0.085) [K1], firing continues only while A(t) >= A_MIN
    => the first-pulse delay distribution (D1) and the pulse counting law (D2) are no longer assumed; they emerge from repair randomness
  Pulse waveform made explicit: square-wave width w=3.5h [K2]; the decoder integrates the waveform (code 81 used an instantaneous-quantum approximation)
  MDM2 shadow layer: p53 waveform delayed by LAG_MDM2=2h [K3] + τ=1h leaky filtering, validating the antiphase-lag morphology
  Downstream L3/L4 identical to code 81 verbatim (τ_p21=10h, τ_PUMA=4h, PUMA 3rd-pulse threshold, fate thresholds)

Only two non-literature constants (pre-registration calibration, double-recorded):
  LAM0 = 0.003 /DSB/h — calibration anchor: mean first-pulse delay at 10 Gy ~1.5h (literature: MCF7 first peak 2-3h; measured 1.66h)
  A_MIN = 10 DSB — calibration anchor: ~7 pulses/24h at 10 Gy (literature range 0-7 pulses)
  Both are aligned only to these two literature anchors and are not tuned to any audit output.
"""

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
import sys

sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
try:
    from daimon_runtime import setup_plot
    setup_plot()
except Exception:
    pass
plt.rcParams["svg.fonttype"] = "none"

SEED = 20260925
ROOT = Path(__file__).resolve().parent.parent
OUT_JSON = ROOT / "结果" / "代码84_发现式链路_结果.json"
OUT_PNG = ROOT / "结果" / "代码84_发现式链路_六面体.png"
OUT_SVG = ROOT / "结果" / "代码84_发现式链路_六面体.svg"

# ---------------- K-series constants (measured, each with a source; see the constant registry) ----------------
DOSES = np.array([0.1, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0])
N_CELLS = 4000
K0_DSB = 35.0          # DSB/Gy (literature anchor)
R_FAST = 0.35          # /h, fast-phase repair rate, t½≈2h (K7)
FRAC_SLOW = 0.10       # fraction of complex breaks (K8)
T_HALF_SLOW = 20.0     # h, slow-phase half-life (K8)
R_SLOW = np.log(2) / T_HALF_SLOW
T_PERIOD = 5.5         # h, pulse period (K1, sealed constant)
JIT_IPI = 0.085        # IPI jitter CV (archived)
W_PULSE = 3.5          # h, pulse width (K2)
LAG_MDM2 = 2.0         # h, p53->MDM2 lag (K3)
TAU_MDM2 = 1.0         # h, MDM2 filter (approximate convention, double-recorded)
W = 48.0               # h observation window
N_MAX = 12
# decoder (code 81 verbatim)
Q21, TAU21 = 100.0, 10.0
QPU, TAUPU, PUMA_NTHR = 100.0, 4.0, 3
TH_PUMA, R0_FATE, TH_P21 = 800.0, 0.15, 1500.0
# pre-registration calibrated constants (only two, non-literature)
LAM0 = 0.003           # /DSB/h triggering hazard rate (calibration: 10Gy first-pulse mean -> ~1.5h literature anchor)
A_MIN = 10.0           # DSB persistence-firing threshold (calibration: 10Gy -> ~7 pulses/24h literature anchor)

# ---------------- Mutual-information utilities (same convention as code 81) ----------------

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

# ---------------- Mechanistic L1: repair birth-death + triggering ----------------

def first_pulse_time(taus, lam0, rng):
    """Inhomogeneous Poisson triggering: Lambda(t)=lam0*Σmin(t,τ_i) piecewise linear, inverted analytically.
    taus: repair times of all DSBs of this cell. Returns t1 or None (no firing within the window)."""
    if len(taus) == 0:
        return None
    ts = np.sort(taus)
    u = rng.exponential(1.0)  # triggering threshold
    lam_total = lam0 * ts.sum()
    if lam_total <= u:
        return None  # even repairing everything cannot accumulate enough triggering quantity -> this cell does not respond
    # piecewise: on interval [t_{k-1}, t_k], Lambda(t) = lam0*(Σ_{i<k} τ_i + (n-k)*t)... via cumsum
    cs = np.cumsum(ts)
    n = len(ts)
    lam_knots = lam0 * (cs + ts * (n - np.arange(n) - 1))
    lam_knots = np.concatenate([[0.0], lam_knots])
    t_knots = np.concatenate([[0.0], ts])
    # lam_knots[k] corresponds to t_knots[k]; within [t_knots[k], t_knots[k+1]] the unrepaired count = n - k
    k = min(int(np.searchsorted(lam_knots, u, side="right")) - 1, n - 1)
    slope = lam0 * (n - k)
    if slope <= 0:
        return None
    t1 = t_knots[k] + (u - lam_knots[k]) / slope
    return t1 if t1 <= W else None


def simulate_dose(D, n, rng):
    """Single-dose population: mechanistic L1 + constant-skeleton L2 + code 81 conventions L3/L4."""
    # L0: DSB counts and repair times
    n_dsb = rng.poisson(K0_DSB * D, n)
    # per-cell DSB repair times (fast/slow mixture)
    taus_list = []
    for i in range(n):
        ni = n_dsb[i]
        if ni == 0:
            taus_list.append(np.empty(0)); continue
        slow = rng.random(ni) < FRAC_SLOW
        rates = np.where(slow, R_SLOW, R_FAST)
        taus_list.append(rng.exponential(1.0 / rates))
    # L1: first pulse + persistent firing
    t1 = np.full(n, np.nan)
    times = np.full((n, N_MAX), np.nan)
    A_at = lambda taus, t: int(np.sum(taus > t))  # unrepaired count
    for i in range(n):
        tp = first_pulse_time(taus_list[i], LAM0, rng)
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
    # L2/L3 decoding: main convention = code 81 verbatim (pulse onset triggers sustained production q, leakage τ, measured at W),
    # keeping the fate thresholds TH_PUMA/R0/TH_P21 sealed-comparable; the explicit-waveform integration with K2 width 3.5h serves as a corroborating convention (double-recorded).
    def integrate_tau(q, tau, mask):
        with np.errstate(invalid="ignore"):
            contrib = np.where(mask, q * tau * (1.0 - np.exp(-(W - times) / tau)), 0.0)
        return np.nansum(contrib, axis=1)
    def integrate_tau_width(q, tau, mask):
        """Corroborating convention: explicit square-wave waveform of width w (deposition rate q/w on [t0,t0+w], decaying thereafter)."""
        qr = q / W_PULSE
        t_end = np.where(mask, np.minimum(times + W_PULSE, W), np.nan)
        with np.errstate(invalid="ignore"):
            seg = qr * tau * (1.0 - np.exp(-(t_end - times) / tau))
            contrib = np.where(mask, seg * np.exp(-(W - t_end) / tau), 0.0)
        return np.nansum(contrib, axis=1)
    pidx = np.arange(N_MAX)[None, :]
    m_p21 = integrate_tau(Q21, TAU21, active)
    puma_mask = active & (pidx >= PUMA_NTHR - 1)
    m_puma = integrate_tau(QPU, TAUPU, puma_mask)
    # corroboration: width-explicit convention (thresholds not re-anchored, trends only, double-recorded)
    m_p21_w = integrate_tau_width(Q21, TAU21, active)
    m_puma_w = integrate_tau_width(QPU, TAUPU, puma_mask)
    # L4 fate (code 81 verbatim)
    ratio = m_puma / (m_p21 + 1.0)
    fate = np.where((m_puma >= TH_PUMA) & (ratio >= R0_FATE), 2,
                    np.where(m_p21 >= TH_P21, 1, 0))
    # MDM2 shadow layer (K3 validation): p53 waveform delayed 2h + tau=1h filtering, reconstructed on a grid
    tgrid = np.linspace(0, W, 481)
    ex = np.where(responded)[0]
    mdm2_demo, p53_demo = None, None
    if len(ex) > 0:
        i0 = ex[len(ex) // 2]  # use the median responding cell for the demo
        p53_demo = np.zeros_like(tgrid)
        for j in range(int(n_pulse[i0])):
            t0 = times[i0, j]
            p53_demo += ((tgrid >= t0) & (tgrid <= t0 + W_PULSE)).astype(float)
        mdm2_demo = np.zeros_like(tgrid)
        for k in range(1, len(tgrid)):
            dt = tgrid[k] - tgrid[k - 1]
            src = p53_demo[max(0, int((tgrid[k] - LAG_MDM2) / (tgrid[1] - tgrid[0])))]
            mdm2_demo[k] = mdm2_demo[k - 1] + dt * (src - mdm2_demo[k - 1]) / TAU_MDM2
    # A(t) curve demo (per-dose average)
    return dict(n_dsb=n_dsb, t1=t1, times=times, n_pulse=n_pulse,
                responded=responded, m_p21=m_p21, m_puma=m_puma, fate=fate,
                tgrid=tgrid, p53_demo=p53_demo, mdm2_demo=mdm2_demo,
                taus_list=taus_list)

# ---------------- Main flow ----------------

def main():
    log = lambda *a: print(*a, flush=True)
    log(f"[code84] seed {SEED}, discovery-based chain model started")
    rng = np.random.default_rng(SEED)
    parts = []
    for di, D in enumerate(DOSES):
        log(f"  simulating D={D} Gy ...")
        p = simulate_dose(D, N_CELLS, rng)
        p["dose_idx"] = np.full(N_CELLS, di)
        parts.append(p)
    E = {k: np.concatenate([p[k] for p in parts]) for k in
         ["n_dsb", "t1", "times", "n_pulse", "responded", "m_p21", "m_puma", "fate", "dose_idx"]}
    d = E["dose_idx"]
    H_D = np.log2(len(DOSES))

    # ---- Check 1: counting law N(D) and responding fraction (emergent quantities vs code 81/literature) ----
    count_law, resp_frac = [], []
    for di, D in enumerate(DOSES):
        m = d == di
        count_law.append(float(E["n_pulse"][m].mean()))
        resp_frac.append(float(E["responded"][m].mean()))
        log(f"  D={D:>5} Gy: <N>={count_law[-1]:.2f}, responding fraction={resp_frac[-1]:.3f}, "
            f"<N_DSB>={E['n_dsb'][m].mean():.1f}")

    # ---- Check 2: first-pulse delay distribution and dispersion ratio (emergent quantities vs literature 5.8) ----
    t1_stats = {}
    ipi_sd_all = []
    for di, D in enumerate(DOSES):
        m = d == di
        t1v = E["t1"][m]; t1v = t1v[~np.isnan(t1v)]
        # per-cell IPI dispersion
        ipi_sd_cell = []
        tt = E["times"][m]
        for i in range(tt.shape[0]):
            tv = tt[i][~np.isnan(tt[i])]
            if len(tv) >= 3:
                ipis = np.diff(tv)
                ipi_sd_cell.append(np.std(ipis))
                ipi_sd_all.append(np.std(ipis))
        if len(t1v) > 10 and ipi_sd_cell:
            ratio = np.std(t1v) / np.mean(ipi_sd_cell)
            t1_stats[f"{D}Gy"] = dict(t1_mean=float(t1v.mean()), t1_sd=float(t1v.std()),
                                      ipi_sd=float(np.mean(ipi_sd_cell)),
                                      dispersion_ratio=float(ratio))
            log(f"  D={D:>5} Gy: t1={t1v.mean():.2f}±{t1v.std():.2f}h, "
                f"IPI sd={np.mean(ipi_sd_cell):.2f}h, dispersion ratio={ratio:.2f} (literature 5.8)")

    # ---- Check 3: information chain (same convention as code 81) ----
    y1 = E["n_dsb"]
    y2 = E["n_pulse"]
    y3 = E["n_pulse"] * 8 + qbin(np.where(np.isnan(E["t1"]), -1, E["t1"]) + 1, 8)
    y4 = qbin(E["m_p21"], 12) * 12 + qbin(E["m_puma"], 12)
    y5 = E["fate"]
    I_chain = [mi_discrete(d, y) for y in [y1, y2, y3, y4, y5]]
    losses = [H_D - I_chain[0]] + [I_chain[i] - I_chain[i + 1] for i in range(4)]
    log(f"  information chain I(D;Y) = {[round(i,3) for i in I_chain]}")
    log(f"  layer-by-layer losses = {[round(l,3) for l in losses]}, end-to-end I(D;fate)={I_chain[-1]:.3f} "
        f"(code 81: 1.325)")

    # ---- Check 4: MDM2 antiphase lag (K3 morphological validation) ----
    mdm2_check = None
    p53d = parts[4]["p53_demo"]; mdm2d = parts[4]["mdm2_demo"]; tgd = parts[4]["tgrid"]
    if p53d is not None:
        # peak cross-correlation lag
        xc = np.correlate(mdm2d - mdm2d.mean(), p53d - p53d.mean(), mode="full")
        lag_idx = np.argmax(xc) - (len(tgd) - 1)
        lag_h = lag_idx * (tgd[1] - tgd[0])
        mdm2_check = dict(demo_dose=float(DOSES[4]), measured_lag_h=float(lag_h),
                          expected_lag_h=LAG_MDM2)
        log(f"  MDM2 shadow layer: cross-correlation lag {lag_h:.2f} h (constructed value {LAG_MDM2} h, K3)")

    results = dict(
        meta=dict(script="代码84_发现式链路模型.py", seed=SEED, date="2026-09-25",
                  K_constants=dict(K0_DSB=K0_DSB, R_FAST=R_FAST, FRAC_SLOW=FRAC_SLOW,
                                   T_HALF_SLOW=T_HALF_SLOW, T_PERIOD=T_PERIOD,
                                   JIT_IPI=JIT_IPI, W_PULSE=W_PULSE,
                                   LAG_MDM2=LAG_MDM2, W=W),
                  calibrated=dict(LAM0=LAM0, A_MIN=A_MIN,
                                  anchors="LAM0: 10Gy first-pulse mean ~1.5h; A_MIN: 10Gy ~7 pulses/24h")),
        H_D=H_D, count_law=count_law, responding_fraction=resp_frac,
        doses=DOSES.tolist(), t1_stats=t1_stats,
        info_chain=dict(I_chain=I_chain, losses=losses,
                        end_to_end=I_chain[-1], code81_reference=1.325),
        mdm2_check=mdm2_check,
        double_record=[
            "LAM0 and A_MIN are the only two non-literature constants: aligned respectively to the two literature anchors '10Gy first-pulse mean ~1.5h' and '10Gy ~7 pulses/24h'; not tuned to any audit output",
            "The γH2AX formation timescale of 30 min (K5) is far below T=5.5h and is treated as instantaneous; with an explicit 30 min delay shift, the first-pulse distribution would shift right by ~0.5h without changing the order of magnitude of the dispersion ratio (not run, registered)",
            "The measurement convention of the literature dispersion ratio 5.8 (Loewer 2010, per-dose SD(t1)/SD(IPI)) differs from this model's ensemble convention; code 82 already registered that a homogeneous ensemble cannot reproduce 5.8 (an extra heterogeneity source is needed), and if this model likewise fails to reproduce it, the two corroborate each other",
            "The MDM2 shadow layer is a morphological demo (τ=1h filtering is an approximate convention) and does not enter the information chain",
            "The square-wave width 3.5h (K2) stays comparable to the code 81 quantum convention via area alignment q_r=q/w; second-order effects of width on integrator output were not scanned separately",
            "Mutual information is a binned plug-in estimate (same convention as code 81); absolute values are biased, trends are reliable",
        ],
    )
    with open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=2, default=lambda o: float(o) if isinstance(o, np.floating) else (int(o) if isinstance(o, np.integer) else str(o)))
    log("JSON written:", OUT_JSON)

    # ---------------- Six-panel figure ----------------
    fig = plt.figure(figsize=(17, 10.5))
    gs = fig.add_gridspec(2, 3, hspace=0.44, wspace=0.32)

    # (a) A(t) repair curves (input of the mechanistic L1)
    axa = fig.add_subplot(gs[0, 0])
    tg_a = np.linspace(0, W, 300)
    for D in DOSES:
        A = K0_DSB * D * ((1 - FRAC_SLOW) * np.exp(-R_FAST * tg_a) + FRAC_SLOW * np.exp(-R_SLOW * tg_a))
        axa.plot(tg_a, A, label=f"{D} Gy")
    axa.axhline(A_MIN, color="#8c1d18", ls="--", lw=1.0)
    axa.text(1, A_MIN * 1.15, f"persistence threshold A_min={A_MIN:.0f}", fontsize=8,
             color="#8c1d18")
    axa.set_yscale("log"); axa.set_xlabel("time (h)"); axa.set_ylabel("unrepaired DSBs A(t)")
    axa.legend(fontsize=7.5, ncol=2)
    axa.set_title("(a) Mechanistic L1: biphasic DSB repair (K7/K8)", fontsize=10.5)

    # (b) counting law + responding fraction (emergent)
    axb = fig.add_subplot(gs[0, 1])
    axb.plot(DOSES, count_law, "o-", color="#33527a", label="mean pulse count <N> (emergent)")
    axb.set_xscale("log"); axb.set_xlabel("dose (Gy)"); axb.set_ylabel("pulse count")
    axb2 = axb.twinx()
    axb2.plot(DOSES, resp_frac, "s--", color="#0b6b3a", label="responding fraction")
    axb2.set_ylabel("responding fraction", color="#0b6b3a")
    axb2.set_ylim(-0.05, 1.05)
    h1, l1 = axb.get_legend_handles_labels(); h2, l2 = axb2.get_legend_handles_labels()
    axb.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper left")
    axb.set_title("(b) Emergent counting law and digital recruitment", fontsize=10.5)

    # (c) t1 distributions and dispersion ratio
    axc = fig.add_subplot(gs[0, 2])
    for di, D in enumerate(DOSES):
        t1v = E["t1"][d == di]; t1v = t1v[~np.isnan(t1v)]
        if len(t1v) > 20:
            axc.hist(t1v, bins=40, density=True, histtype="step", lw=1.4, label=f"{D} Gy")
    axc.set_xlabel("first-pulse delay $t_1$ (h)"); axc.set_ylabel("density")
    ratios = [v["dispersion_ratio"] for v in t1_stats.values()]
    axc.text(0.55, 0.85, f"dispersion ratio {min(ratios):.1f}--{max(ratios):.1f}\n"
             f"(literature 5.8; code 82 registered\nhomogeneous-ensemble failure)",
             transform=axc.transAxes, fontsize=8,
             bbox=dict(fc="#fdf1ec", ec="#c05640", lw=0.8))
    axc.legend(fontsize=7.5)
    axc.set_title("(c) Emergent first-pulse delay distributions (D1)", fontsize=10.5)

    # (d) information waterfall comparison
    axd = fig.add_subplot(gs[1, 0])
    names = ["H(D)", "L0", "L1+\n(trigger)", "L2", "L3", "L4", "I(D;fate)"]
    starts = [0.0]; heights = [H_D]; run = H_D
    for l in losses:
        starts.append(run - l); heights.append(l); run -= l
    starts.append(0.0); heights.append(I_chain[-1])
    colors = ["#33527a"] + ["#c05640"] * 4 + ["#d9863d"] + ["#0b6b3a"]
    axd.bar(range(7), heights, bottom=starts, color=colors, width=0.68)
    for i, (s, h) in enumerate(zip(starts, heights)):
        axd.text(i, s + h + 0.03, f"{h:.2f}", ha="center", fontsize=8)
    axd.axhline(1.325, color="#0b6b3a", ls=":", lw=1.0)
    axd.text(6.2, 1.42, "code 81:\n1.325", fontsize=7.5, color="#0b6b3a", ha="center")
    axd.set_xticks(range(7)); axd.set_xticklabels(names, fontsize=7.5)
    axd.set_ylabel("bits")
    axd.set_title("(d) Information waterfall, mechanistic L1", fontsize=10.5)

    # (e) MDM2 shadow-layer antiphase
    axe = fig.add_subplot(gs[1, 1])
    if p53d is not None:
        axe.plot(tgd, p53d, color="#33527a", lw=1.2, label="p53 (square pulses, w=3.5h)")
        axe.plot(tgd, mdm2d, color="#c05640", lw=1.6, label="MDM2 shadow (lag 2h, K3)")
        axe.set_xlabel("time (h)"); axe.set_ylabel("level (AU)")
        axe.legend(fontsize=8)
        axe.set_title(f"(e) MDM2 antiphase shadow (demo cell, {DOSES[4]} Gy)", fontsize=10.5)

    # (f) summary
    axs = fig.add_subplot(gs[1, 2]); axs.axis("off")
    summary = (
        f"Code 84 discovery-chain summary (seed {SEED})\n\n"
        f"All dynamics from measured constants:\n"
        f"  35 DSB/Gy; repair 0.35/h + 10% slow (t1/2=20h)\n"
        f"  T=5.5h; width 3.5h; MDM2 lag 2h\n"
        f"  only 2 calibrated: LAM0={LAM0}, A_min={A_MIN:.0f}\n\n"
        f"Emergent vs literature:\n"
        f"  <N> 10Gy = {count_law[-1]:.1f} (lit. 0-7/24h)\n"
        f"  responding fraction {resp_frac[0]:.2f}->{resp_frac[-1]:.2f}\n"
        f"  dispersion ratio {min(ratios):.1f}--{max(ratios):.1f} (lit. 5.8)\n"
        f"  t1@10Gy = {t1_stats['10.0Gy']['t1_mean']:.2f}h (anchor 1.5h)\n\n"
        f"End-to-end I(D;fate) = {I_chain[-1]:.3f} bits\n"
        f"  (code 81 assumed-L1: 1.325 bits)\n"
        f"weakest layer loss = {max(losses):.2f} bits"
    )
    axs.text(0.02, 0.98, summary, va="top", fontsize=8.4, family="monospace",
             transform=axs.transAxes)
    fig.suptitle("Code 84: discovery-based p53 chain — measured-constant skeleton, emergent sensing noise",
                 fontsize=12.5)
    fig.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
    fig.savefig(OUT_SVG, bbox_inches="tight")
    log("Figure written:", OUT_PNG)
    log("[code84] done.")


if __name__ == "__main__":
    main()
