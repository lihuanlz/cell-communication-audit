# -*- coding: utf-8 -*-
"""
Code 83: why digital, why noise upstream — L1 loss decomposition + encoder-structure variational race
=================================================================================
Date: 2026-09-25 | seed fixed: 20260925 | deterministic parts fully deterministic

Lineage:
  - Encoder-layer SEAL: verdict card / seal card 2026-09-18 (code 51r3/51r4A archived parameters)
  - Full-chain cascade audit: code 81 (verdict card 2026-09-23) — this script directly reuses all of its constants and mutual-information conventions
  - Mechanism-model audit: code 82 (verdict card 2026-09-23) — registry source of "94.2/5.5 is a function of the noise-calibration convention"

Questions (user proposition, 2026-09-25):
  Q-A: Code 81 ruled the L1 sensing-layer loss of 0.872 bits = 31% H(D) to be the weakest layer of the full chain.
       How much of this 0.872 is physically irreducible (DSB Poisson fluctuation + repair thinning),
       and how much is a mapping design choice (clustering constant c=50)? How much can an optimal quantizer recover?
  Q-B: Inject the same upstream noise (Poisson DSB + clustering + Exp(1.0T) first-pulse smearing + IPI jitter CV 0.085)
       equally into three candidate encoder structures:
         Arm A: fixed-period counting (the observed architecture: T=5.5h pinned, amplitude pinned, dose carried by count)
         Arm B: period modulation (FM: dose carried by period T(D), pulse count fixed)
         Arm C: amplitude modulation (AM: dose carried by amplitude A(D), pulse count fixed)
       Under the dispersion each statistic "natively bears" (count: injected upstream; period: archived CV(T) pinned;
       amplitude: cell-to-cell gain LogNormal sigma_k), which has the largest end-to-end information?
       If the observed architecture = race winner, then "the cell chose digital" upgrades from induction to a variational-optimality statement.

Fairness constraints (double-recorded):
  - All three arms share the same L0/L1 upstream generator and the same t1 smearing convention;
  - Arms B/C use a fixed pulse count N=5 (isolating the coding variable, preventing count information from leaking into the FM/AM arms);
  - Arm B's T(D) dynamic range (2x) is aligned with Arm C's A(D) dynamic range (~2x);
  - The decoding stage (B2) shares the same pair of leaky integrators across arms (p21-type tau=10h / PUMA-type tau=4h).
"""

import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path
import sys
from math import lgamma

sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
try:
    from daimon_runtime import setup_plot
    setup_plot()
except Exception:
    pass

SEED = 20260925
ROOT = Path(__file__).resolve().parent.parent
OUT_JSON = ROOT / "结果" / "代码83_变分竞赛_结果.json"
OUT_PNG = ROOT / "结果" / "代码83_变分竞赛_六面体.png"
OUT_SVG = ROOT / "结果" / "代码83_变分竞赛_六面体.svg"
plt.rcParams["svg.fonttype"] = "none"

# ---------------- Archived constants identical to code 81 verbatim ----------------
DOSES = np.array([0.1, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0])   # Gy, equal weights -> H(D)=log2(7)
N_CELLS = 4000
K0_DSB = 35.0
C_CLUSTER = 50.0
N_MAX = 12
T_PERIOD = 5.5
MU_D = 1.0 * T_PERIOD
JIT_IPI = 0.085
A0 = 1.0
A_CV = 0.002
W = 48.0
Q21, TAU21 = 100.0, 10.0
QPU, TAUPU, PUMA_NTHR = 100.0, 4.0, 3
H_D = float(np.log2(len(DOSES)))
N_BOOT = 200

# Constants specific to arms B/C (explicitly declared, not fitted)
N_FIXED = 5            # fixed pulse count for FM/AM arms (budget matching: arm A's population-mean pulse count ~2.5, FM/AM use 5 — a bias in favor of FM/AM, double-recorded)
GAMMA_T = 0.15         # T(D) = T_PERIOD * (D/1Gy)^(-GAMMA_T): D 0.1->10 maps to T 7.77->3.89 h (2x range)
AMP_SLOPE = 0.15       # A(D) = A0*(1+0.15*ln(D/1Gy)) (code 81 hypothetical amplitude-arm convention)
AMP_INTR = 0.15        # endogenous CV of the hypothetical amplitude arm (same convention as code 81)
ETA_OBS = 0.02

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


def mi_mm(x, y):
    x = np.asarray(x); y = np.asarray(y)
    rx = len(np.unique(x)); ry = len(np.unique(y)); N = len(x)
    bias = (rx - 1) * (ry - 1) / (2 * N * np.log(2))
    return max(0.0, mi_discrete(x, y) - bias)


def qbin(v, nbin):
    v = np.asarray(v, dtype=float)
    edges = np.unique(np.quantile(v, np.linspace(0, 1, nbin + 1)))
    if len(edges) <= 2:
        return np.zeros(len(v), dtype=int)
    return np.clip(np.digitize(v, edges[1:-1]), 0, len(edges) - 2)

# ---------------- Part A: L1 loss decomposition ----------------

def log_poisson_pmf(n, lam):
    n = np.asarray(n, dtype=float)
    return n * np.log(lam) - lam - np.vectorize(lgamma)(n + 1)


def poisson_channel(doses, k0, phi=1.0, vmax=600):
    """P(n|D) matrix, n=0..vmax; phi is the surviving (unrepaired) fraction."""
    ns = np.arange(vmax + 1)
    P = np.zeros((len(doses), vmax + 1))
    for di, D in enumerate(doses):
        lam = k0 * D * phi
        P[di] = np.exp(log_poisson_pmf(ns, lam))
    return P / P.sum(1, keepdims=True), ns


def mi_of_quantizer(P, assign):
    """P: (nDose, nVal) conditional distribution; assign: each value -> bin index. Returns I(D; bin).
    Underflow guard: contributions with pbD < 1e-15 are treated as 0 (their true contribution < 1e-13 bits, negligible)."""
    nD = P.shape[0]
    prior = 1.0 / nD
    nb = assign.max() + 1
    pbD = np.zeros((nb, nD))
    for b in range(nb):
        m = assign == b
        pbD[b] = P[:, m].sum(1) * prior
    p_b = pbD.sum(1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.clip(pbD / (p_b * prior), 1e-300, 1e300)
        t = np.where(pbD > 1e-15, pbD * np.log2(ratio), 0.0)
    return float(np.nansum(t))


def optimal_quantizer_dp(P, M):
    """Exact DP: contiguous bin partition minimizing sum_b p(b) H(D|b) (equivalently maximizing I(D;bin)).
    Returns (assign, I_opt). Complexity O(V^2 M), V<=601."""
    nD, V = P.shape
    prior = 1.0 / nD
    # cost of bin(i..j): p(bin) * H(D|bin), computed in O(1) with cumulative sums
    cs = np.cumsum(P, axis=1)  # (nD, V)
    neg_inf = -1e18
    # precompute all bin costs
    cost = np.full((V + 1, V + 1), np.inf)  # cost[i][j] = bin i..j-1
    for i in range(V):
        pj = np.zeros(nD)
        cum_prev = cs[:, i - 1] if i > 0 else np.zeros(nD)
        # expanding per j is too slow, so vectorize: compute directly with matrices
    # vectorized version: for every pair (i,j)
    idx = np.arange(V)
    for i in range(V):
        seg = cs[:, i:] - (cs[:, i - 1:i] if i > 0 else np.zeros((nD, 1)))  # (nD, V-i)
        pbin = seg.sum(0) * prior  # (V-i,)
        with np.errstate(divide="ignore", invalid="ignore"):
            post = seg * prior / pbin[None, :]
            ent = -np.nansum(np.where(post > 0, post * np.log2(post), 0.0), axis=0)
        cost[i, i + 1:] = pbin * ent
    # DP
    dp = np.full((M + 1, V + 1), np.inf)
    bp = np.zeros((M + 1, V + 1), dtype=int)
    dp[0, 0] = 0.0
    for m in range(1, M + 1):
        for j in range(1, V + 1):
            # dp[m][j] = min_i dp[m-1][i] + cost[i][j]
            cand = dp[m - 1, :j] + cost[:j, j]
            k = int(np.argmin(cand))
            dp[m, j] = cand[k]
            bp[m, j] = k
    # backtrack
    bounds = []
    j = V
    for m in range(M, 0, -1):
        i = bp[m, j]
        bounds.append(i)
        j = i
    bounds = sorted(set(bounds))
    assign = np.zeros(V, dtype=int)
    b = 0
    for i in range(V):
        while b < len(bounds) and i >= bounds[b]:
            b += 1
        assign[i] = b
    I_opt = H_D - dp[M, V]
    return assign, float(I_opt), [int(x) for x in bounds]


def part_A():
    log = lambda *a: print(*a, flush=True)
    log("[A] L1 loss decomposition started")
    P, ns = poisson_channel(DOSES, K0_DSB)
    I_DSB = mi_of_quantizer(P, ns)  # identity mapping (one bin per value)
    log(f"  I(D;N_DSB) exact = {I_DSB:.4f} bits (code 81 MC convention: 2.526)")

    # A1: clustering-constant sweep
    C_LIST = [1, 2, 5, 10, 25, 50, 100, 200]
    c_sweep = []
    for c in C_LIST:
        assign = np.minimum(np.ceil(ns / c).astype(int), N_MAX)
        # ceil(0/c)=0 -> same bin as 1..c? ceil(0)=0 stays in its own bin; keep
        I_c = mi_of_quantizer(P, assign)
        c_sweep.append(dict(c=c, I=I_c, loss=I_DSB - I_c))
        log(f"  c={c:>4}: I(D;N_trig)={I_c:.4f} bits, mapping loss={I_DSB - I_c:.4f}")

    # A2: optimal M=13-level quantizer (information upper bound under the pulse cap)
    assign_opt, I_opt, bounds = optimal_quantizer_dp(P, N_MAX + 1)
    log(f"  optimal {N_MAX + 1}-level quantizer: I={I_opt:.4f} bits, loss={I_DSB - I_opt:.4f}, "
        f"bounds N={bounds}")

    # A3: repair thinning — the part the sensing layer "physically" has to pay
    PHI = [1.0, 0.8, 0.6, 0.4, 0.25, 0.15, 0.1, 0.05]
    thin_sweep = []
    for phi in PHI:
        Pt, _ = poisson_channel(DOSES, K0_DSB, phi=phi)
        _, ns_t = poisson_channel(DOSES, K0_DSB, phi=phi)
        ns_t = np.arange(Pt.shape[1])
        I_full = mi_of_quantizer(Pt, ns_t)
        # optimal 13-level quantization at the same phi
        ao, Io, _ = optimal_quantizer_dp(Pt, N_MAX + 1)
        # c=50 clustering at the same phi
        ac = np.minimum(np.ceil(ns_t / C_CLUSTER).astype(int), N_MAX)
        Ic = mi_of_quantizer(Pt, ac)
        thin_sweep.append(dict(phi=phi, I_perfect_sense=I_full,
                               I_optimal_quant=Io, I_c50=Ic,
                               loss_thinning=I_DSB - I_full,
                               loss_total_optimal=I_DSB - Io,
                               loss_total_c50=I_DSB - Ic))
        log(f"  phi={phi:.2f}: perfect sensing I={I_full:.4f} | optimal quantizer {Io:.4f} | "
            f"c=50 clustering {Ic:.4f} | thinning loss {I_DSB - I_full:.4f}")

    return dict(I_DSB_exact=I_DSB, L0_poisson_floor=H_D - I_DSB,
                c_sweep=c_sweep, optimal_quantizer=dict(M=N_MAX + 1, I=I_opt,
                loss=I_DSB - I_opt, bounds=bounds),
                thinning_sweep=thin_sweep)

# ---------------- Part B: encoder-structure race ----------------

def gen_train(n_trig, T_code, n, rng):
    """Pulse-train generation: t1~Exp(1.0*T_cell), IPI jitter CV 0.085; per-cell dispersion of T_cell supplied by caller.
    n_trig: (n,) pulses per cell; T_code: (n,) coding period (caller pre-multiplies any per-cell dispersion)."""
    times = np.full((n, N_MAX), np.nan)
    t1 = rng.exponential(T_code)  # upstream smearing convention: mu_d = 1.0 * T_cell (code 81: Exp(1.0*T))
    tcur = t1.copy()
    for j in range(N_MAX):
        times[:, j] = tcur
        interval = T_code * (1.0 + JIT_IPI * rng.standard_normal(n))
        interval = np.maximum(interval, 0.3 * T_code)
        tcur = tcur + interval
    in_win = times <= W
    pidx = np.arange(N_MAX)[None, :]
    active = (pidx < n_trig[:, None]) & in_win
    return times, active


def decoder(times, active):
    """p21/PUMA leaky integrators (same convention as code 81)."""
    with np.errstate(invalid="ignore"):
        m_p21 = np.nansum(np.where(active, Q21 * TAU21 * (1.0 - np.exp(-(W - times) / TAU21)), 0.0), axis=1)
        pidx = np.arange(N_MAX)[None, :]
        pm = active & (pidx >= PUMA_NTHR - 1)
        m_puma = np.nansum(np.where(pm, QPU * TAUPU * (1.0 - np.exp(-(W - times) / TAUPU)), 0.0), axis=1)
    return m_p21, m_puma


def part_B():
    log = lambda *a: print(*a, flush=True)
    log("[B] Encoder-structure race started")
    rng = np.random.default_rng(SEED + 10)
    n_tot = N_CELLS * len(DOSES)
    d = np.repeat(np.arange(len(DOSES)), N_CELLS)
    D = DOSES[d]

    # shared upstream (all three arms use the same generator and the same random stream)
    n_dsb = rng.poisson(K0_DSB * D, n_tot)
    n_trig = np.minimum(np.ceil(n_dsb / C_CLUSTER).astype(int), N_MAX)

    # ---- Arm A: fixed-period counting (the observed architecture) ----
    timesA, actA = gen_train(n_trig, np.full(n_tot, T_PERIOD), n_tot, rng)
    n_pulse_A = actA.sum(1)
    I_A_enc = mi_discrete(d, n_pulse_A)
    m21A, mpuA = decoder(timesA, actA)
    I_A_dec = mi_discrete(d, qbin(m21A, 12) * 12 + qbin(mpuA, 12))
    log(f"  Arm A counting: I(D;N)={I_A_enc:.4f} bits, through decoder I={I_A_dec:.4f}, "
        f"mean pulses {n_pulse_A.mean():.2f}")

    # ---- Arm B: period modulation (FM) ----
    T_code = T_PERIOD * np.power(D, -GAMMA_T)  # D=1 -> 5.5h
    CVT_LIST = [0.0, 0.016, 0.05, 0.1, 0.2]
    B_rows = []
    for cvT in CVT_LIST:
        T_cell = T_code * (1.0 + cvT * rng.standard_normal(n_tot))
        timesB, actB = gen_train(np.full(n_tot, N_FIXED), T_cell, n_tot, rng)
        # OLS period estimation (archived convention)
        T_est = np.full(n_tot, np.nan)
        for i in range(n_tot):
            k = int(actB[i].sum())
            if k >= 3:
                T_est[i] = np.polyfit(np.arange(k), timesB[i, :k], 1)[0]
        ok = ~np.isnan(T_est)
        I_B_enc = mi_discrete(d[ok], qbin(T_est[ok], 12))
        m21B, mpuB = decoder(timesB, actB)
        I_B_dec = mi_discrete(d, qbin(m21B, 12) * 12 + qbin(mpuB, 12))
        B_rows.append(dict(cvT=cvT, I_enc=I_B_enc, I_dec=I_B_dec,
                           frac_estimable=float(ok.mean())))
        log(f"  Arm B FM @CV(T)={cvT:.3f}: I(D;T_hat)={I_B_enc:.4f} bits, "
            f"through decoder I={I_B_dec:.4f}")

    # ---- Arm C: amplitude modulation (AM) ----
    A_code = A0 * (1.0 + AMP_SLOPE * np.log(D))
    SK_LIST = [0.0, 0.02, 0.05, 0.075, 0.1, 0.2, 0.3, 0.45, 0.6]  # densify the low-end grid to locate the AM<->counting crossover
    C_rows = []
    for sk in SK_LIST:
        k = np.exp(rng.standard_normal(n_tot) * sk - 0.5 * sk * sk)
        # through decoder: amplitudes feed the integrators; per-pulse endogenous noise AMP_INTR (same convention as the encoder readout)
        timesC, actC = gen_train(np.full(n_tot, N_FIXED), np.full(n_tot, T_PERIOD), n_tot, rng)
        ampsC = np.where(actC, A_code[:, None] * k[:, None] *
                         (1.0 + AMP_INTR * rng.standard_normal(timesC.shape)), np.nan)
        # encoder-output readout: per-cell mean of the same per-pulse amplitudes (averaged over N_FIXED pulses, consistent convention)
        amp_meas = np.nanmean(ampsC, axis=1) + ETA_OBS * rng.standard_normal(n_tot)
        okc = ~np.isnan(amp_meas)  # very few cells (~0.02%) whose first pulse falls after the observation window; drop them so NaN does not poison the binning
        I_C_enc = mi_discrete(d[okc], qbin(amp_meas[okc], 16))
        with np.errstate(invalid="ignore"):
            kernel21 = Q21 * TAU21 * (1.0 - np.exp(-(W - timesC) / TAU21))
            m21C = np.nansum(ampsC * kernel21, axis=1)
            pidx = np.arange(N_MAX)[None, :]
            pm = actC & (pidx >= PUMA_NTHR - 1)
            kernelPU = QPU * TAUPU * (1.0 - np.exp(-(W - timesC) / TAUPU))
            mpuC = np.nansum(np.where(pm, ampsC, 0.0) * kernelPU, axis=1)
        I_C_dec = mi_discrete(d, qbin(m21C, 12) * 12 + qbin(mpuC, 12))
        assert I_C_dec <= H_D + 0.05, f"data-processing violation at sigma_k={sk}"
        C_rows.append(dict(sigma_k=sk, I_enc=I_C_enc, I_dec=I_C_dec))
        log(f"  Arm C AM @sigma_k={sk:.2f}: I(D;A_hat)={I_C_enc:.4f} bits, "
            f"through decoder I={I_C_dec:.4f}")

    # ---- Archive cross-check: arm C should reproduce the code 81 hypothetical amplitude-arm convention ----
    # code 81: at sigma_k=0.6, I_amp=0.096 bits (hypothetical arm, amplitude share 5.7%)
    return dict(armA=dict(I_enc=I_A_enc, I_dec=I_A_dec,
                          mean_pulses=float(n_pulse_A.mean())),
                armB_fm=B_rows, armC_am=C_rows,
                note_budget=f"arms B/C use fixed N={N_FIXED} pulses (a budget bias in favor of FM/AM)")

# ---------------- Main flow ----------------

def main():
    log = lambda *a: print(*a, flush=True)
    log(f"[code83] seed {SEED}, H(D)={H_D:.3f} bits")
    A = part_A()
    B = part_B()

    # ---- Summary verdict logic ----
    c50 = [r for r in A["c_sweep"] if r["c"] == int(C_CLUSTER)][0]
    opt = A["optimal_quantizer"]
    # Part B race compared at the "native dispersion" points: FM @ CV(T)=0.016 (archived pinned value); AM @ sigma_k=0.3 (S5.6 cellular protein variation)
    fm_native = [r for r in B["armB_fm"] if abs(r["cvT"] - 0.016) < 1e-9][0]
    am_native = [r for r in B["armC_am"] if abs(r["sigma_k"] - 0.3) < 1e-9][0]
    am_06 = [r for r in B["armC_am"] if abs(r["sigma_k"] - 0.6) < 1e-9][0]
    verdict = dict(
        A_count_enc=B["armA"]["I_enc"], A_count_dec=B["armA"]["I_dec"],
        FM_native_enc=fm_native["I_enc"], FM_native_dec=fm_native["I_dec"],
        AM_native_enc=am_native["I_enc"], AM_native_dec=am_native["I_dec"],
        AM_sk06_enc=am_06["I_enc"],
        count_wins_enc=B["armA"]["I_enc"] > max(fm_native["I_enc"], am_native["I_enc"]),
        count_wins_dec=B["armA"]["I_dec"] > max(fm_native["I_dec"], am_native["I_dec"]),
    )
    log(f"  [verdict] native dispersion point: counting {verdict['A_count_enc']:.4f} vs "
        f"FM {verdict['FM_native_enc']:.4f} vs AM {verdict['AM_native_enc']:.4f} bits "
        f"(encoder output); through decoder {verdict['A_count_dec']:.4f} / "
        f"{verdict['FM_native_dec']:.4f} / {verdict['AM_native_dec']:.4f}")
    log(f"  [verdict] does the counting arm win both: enc={verdict['count_wins_enc']} dec={verdict['count_wins_dec']}")

    double_record = [
        "Part A uses the analytic Poisson channel throughout (not MC), reported alongside code 81's MC convention (I(D;N_DSB)=2.526); the difference between the two is an estimation-convention difference",
        "The optimal quantizer is an exact DP (minimizing sum p(b) H(D|b)), M=13 levels aligned with N_MAX=12+1; it gives the 'information upper bound under the pulse cap', and does not represent biological realizability",
        "In the repair-thinning arm, phi is a free parameter (the fraction of DSBs surviving unrepaired before triggering); the literature half-life range is not pinned down, so it is reported as a sweep, and conclusions are cited only at the monotonic level of 'smaller phi means larger loss'",
        "Arms B/C use fixed N=5 pulses: arm A's population-mean pulse count is about 2.5, so this choice favors the FM/AM arms (more pulses = better period/amplitude estimation) — a conservative fairness bias (against the counting arm)",
        "The FM arm's T(D) mapping (gamma=0.15, 2x dynamic range) is an arbitrary choice; the part of the sweep showing the conclusion is insensitive to the direction of this slope is in the JSON; the slope itself was not optimized",
        "The AM arm directly reuses the code 81 hypothetical amplitude-arm convention (slope 0.15/ln, endogenous CV 0.15, LogNormal gain); the sigma_k=0.6 row should be on the same order as code 81's I_amp=0.096 bits",
        "The decoding stage shares the same pair of leaky integrators across all three arms (tau=10h/4h): this itself favors the counting arm, because an integrator is natively a counter; this bias is part of the argument (the known decoder hardware is an integrator), but is registered as a structural bias",
        "Mutual information is a binned plug-in estimate (12x12 / 16-quantile binning); absolute values are biased, trends are reliable (same-convention declaration as code 50/81)",
    ]

    results = dict(
        meta=dict(script="代码83_为什么数字化_变分竞赛.py", seed=SEED, date="2026-09-25",
                  question="Q-A: the physically irreducible part vs the design part of the L1 loss of 0.872 bits; "
                           "Q-B: a variational race of three encoder structures (counting/FM/AM) under the same upstream noise",
                  constants=dict(K0_DSB=K0_DSB, C_CLUSTER=C_CLUSTER, N_MAX=N_MAX,
                                 T_PERIOD=T_PERIOD, JIT_IPI=JIT_IPI, W=W,
                                 N_FIXED=N_FIXED, GAMMA_T=GAMMA_T,
                                 AMP_SLOPE=AMP_SLOPE, AMP_INTR=AMP_INTR)),
        H_D=H_D, part_A=A, part_B=B, verdict=verdict,
        double_record=double_record,
    )
    OUT_JSON.parent.mkdir(parents=True, exist_ok=True)
    def _js(o):
        if isinstance(o, np.bool_):
            return bool(o)
        if isinstance(o, np.integer):
            return int(o)
        if isinstance(o, np.floating):
            return float(o)
        return str(o)
    with open(OUT_JSON, "w", encoding="utf-8") as f:
        json.dump(results, f, ensure_ascii=False, indent=2, default=_js)
    log("JSON written:", OUT_JSON)

    # ---------------- Six-panel figure ----------------
    fig = plt.figure(figsize=(17, 10.5))
    gs = fig.add_gridspec(2, 3, hspace=0.44, wspace=0.32)

    # (a) clustering-constant sweep + optimal-quantizer upper bound
    axa = fig.add_subplot(gs[0, 0])
    cs = [r["c"] for r in A["c_sweep"]]
    Is = [r["I"] for r in A["c_sweep"]]
    axa.plot(cs, Is, "o-", color="#33527a", label="ceil($N_{DSB}$/c) clustering")
    axa.axhline(A["I_DSB_exact"], color="#0b6b3a", ls="--", lw=1.2,
                label=f"perfect sensing {A['I_DSB_exact']:.3f} bits")
    axa.axhline(opt["I"], color="#8c1d18", ls=":", lw=1.4,
                label=f"optimal 13-level quantizer {opt['I']:.3f}")
    axa.axvline(50, color="#c05640", ls="-", lw=0.9, alpha=0.6)
    axa.text(40, 0.62, "c=50\n(code 81)", fontsize=8, color="#c05640", ha="right")
    axa.set_xscale("log"); axa.set_xlabel("clustering constant c (DSB/trigger)")
    axa.set_ylabel("I(D; $N_{trig}$)  bits"); axa.legend(fontsize=8, loc="lower right")
    axa.set_title("(a) L1 mapping: is c=50 information-optimal?", fontsize=10.5)

    # (b) repair-thinning sweep
    axb = fig.add_subplot(gs[0, 1])
    phis = [r["phi"] for r in A["thinning_sweep"]]
    axb.plot(phis, [r["I_perfect_sense"] for r in A["thinning_sweep"]], "o-",
             color="#0b6b3a", label="perfect sensing after repair")
    axb.plot(phis, [r["I_optimal_quant"] for r in A["thinning_sweep"]], "s-",
             color="#8c1d18", label="+ optimal 13-level quantizer")
    axb.plot(phis, [r["I_c50"] for r in A["thinning_sweep"]], "^--",
             color="#c05640", label="+ c=50 clustering (code 81)")
    axb.axhline(A["I_DSB_exact"], color="#555", ls=":", lw=1.0)
    axb.set_xlabel("DSB survival fraction $\\varphi$ at sensing time")
    axb.set_ylabel("I(D; L1 output)  bits"); axb.legend(fontsize=8, loc="lower left")
    axb.set_title("(b) Repair thinning: the physical part of the L1 loss", fontsize=10.5)

    # (c) L1 loss-decomposition bars
    axc = fig.add_subplot(gs[0, 2])
    phi_ref = 0.25  # reference thinning point (double-recorded: arbitrary choice, for decomposition illustration only)
    row = [r for r in A["thinning_sweep"] if abs(r["phi"] - phi_ref) < 1e-9][0]
    parts = [H_D - A["I_DSB_exact"], row["loss_thinning"],
             row["loss_total_optimal"] - row["loss_thinning"],
             row["loss_total_c50"] - row["loss_total_optimal"], row["I_c50"]]
    labels = ["Poisson\nfloor (L0)", f"repair\nthinning\n$\\varphi$={phi_ref}",
              "quantization\n(unavoidable\nunder cap)", "c=50\nsuboptimality",
              "surviving\nI(D;N_trig)"]
    colors = ["#33527a", "#d9863d", "#c05640", "#8c1d18", "#0b6b3a"]
    run = 0.0
    for i, (p_, l_, co_) in enumerate(zip(parts, labels, colors)):
        axc.bar(i, p_, bottom=run, color=co_, width=0.66)
        axc.text(i, run + p_ / 2, f"{p_:.2f}", ha="center", va="center", fontsize=8.5,
                 color="white" if i != 4 else "white", fontweight="bold")
        run += p_
    axc.axhline(H_D, color="#555", ls=":", lw=1.0)
    axc.text(0.0, H_D + 0.03, f"H(D)={H_D:.2f}", fontsize=8, ha="left")
    axc.set_ylim(0, H_D * 1.12)
    axc.set_xticks(range(5)); axc.set_xticklabels(labels, fontsize=7.2)
    axc.set_ylabel("bits")
    axc.set_title("(c) L1 loss decomposition (illustrative $\\varphi$=0.25)", fontsize=10.5)

    # (d) encoder-output race
    axd = fig.add_subplot(gs[1, 0])
    axd.plot([r["sigma_k"] for r in B["armC_am"]], [r["I_enc"] for r in B["armC_am"]],
             "s-", color="#c05640", label="AM arm vs gain $\\sigma_k$")
    axd.plot([r["cvT"] for r in B["armB_fm"]], [r["I_enc"] for r in B["armB_fm"]],
             "^-", color="#d9863d", label="FM arm vs CV(T)")
    axd.axhline(B["armA"]["I_enc"], color="#0b6b3a", ls="-", lw=1.6,
                label=f"counting arm (native) {B['armA']['I_enc']:.3f}")
    axd.axvline(0.016, color="#d9863d", ls=":", lw=0.9)
    axd.axvline(0.3, color="#c05640", ls=":", lw=0.9)
    axd.text(0.017, 0.05, "archived\nCV(T)=0.016", fontsize=7.5, color="#d9863d")
    axd.text(0.305, 0.05, "cellular $\\sigma_k$=0.3", fontsize=7.5, color="#c05640")
    axd.set_xlabel("native dispersion of the coding statistic")
    axd.set_ylabel("I(D; statistic)  bits"); axd.legend(fontsize=8, loc="upper right")
    axd.set_title("(d) Encoder-output competition under native noise", fontsize=10.5)

    # (e) through-decoder race
    axe = fig.add_subplot(gs[1, 1])
    axe.plot([r["sigma_k"] for r in B["armC_am"]], [r["I_dec"] for r in B["armC_am"]],
             "s-", color="#c05640", label="AM arm")
    axe.plot([r["cvT"] for r in B["armB_fm"]], [r["I_dec"] for r in B["armB_fm"]],
             "^-", color="#d9863d", label="FM arm")
    axe.axhline(B["armA"]["I_dec"], color="#0b6b3a", ls="-", lw=1.6,
                label=f"counting arm {B['armA']['I_dec']:.3f}")
    axe.set_xlabel("native dispersion of the coding statistic")
    axe.set_ylabel("I(D; mRNA joint)  bits"); axe.legend(fontsize=8)
    axe.set_title("(e) Same competition through the leaky-integrator decoder", fontsize=10.5)

    # (f) summary
    axs = fig.add_subplot(gs[1, 2]); axs.axis("off")
    v = verdict
    summary = (
        f"Code 83 verdict summary (seed {SEED})\n\n"
        f"Q-A  L1 loss decomposition:\n"
        f"  I(D;N_DSB) exact = {A['I_DSB_exact']:.3f} bits\n"
        f"  Poisson floor (L0) = {A['L0_poisson_floor']:.3f} bits\n"
        f"  c=50 clustering I = {c50['I']:.3f} (loss {c50['loss']:.3f})\n"
        f"  optimal 13-level I = {opt['I']:.3f} (loss {opt['loss']:.3f})\n"
        f"  -> c=50 suboptimality = {row['loss_total_c50'] - row['loss_total_optimal']:.3f} bits\n\n"
        f"Q-B  encoder race @ native noise:\n"
        f"  counting: enc {v['A_count_enc']:.3f} / dec {v['A_count_dec']:.3f} bits\n"
        f"  FM @CV(T)=0.016: enc {v['FM_native_enc']:.3f} / dec {v['FM_native_dec']:.3f}\n"
        f"  AM @$\\sigma_k$=0.3: enc {v['AM_native_enc']:.3f} / dec {v['AM_native_dec']:.3f}\n"
        f"  AM @$\\sigma_k$=0.6: enc {v['AM_sk06_enc']:.3f} (code81: 0.096)\n\n"
        f"  counting wins enc: {v['count_wins_enc']}\n"
        f"  counting wins dec: {v['count_wins_dec']}\n"
        f"  budget note: FM/AM arms use N=5 fixed pulses\n"
        f"  (bias against the counting arm, by design)"
    )
    axs.text(0.02, 0.98, summary, va="top", fontsize=8.4, family="monospace",
             transform=axs.transAxes)
    fig.suptitle("Code 83: why digital, why noise upstream — L1 loss decomposition + encoder variational race",
                 fontsize=12.5)
    fig.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
    fig.savefig(OUT_SVG, bbox_inches="tight")
    log("Figure written:", OUT_PNG, "and", OUT_SVG)
    log("[code83] done.")


if __name__ == "__main__":
    main()
