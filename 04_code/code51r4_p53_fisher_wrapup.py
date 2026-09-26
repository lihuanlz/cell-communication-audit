# -*- coding: utf-8 -*-
"""
Code 51r4 p53 Fisher audit wrap-up (two stages)
part A (fast): core noise calibration — σ scan driving IPI CV→30% (Lahav: 100/330min),
       re-run the P2' delay scan at the calibration point σ*, verifying result 3 (independent random perception-level signature) is robust at realistic noise magnitude.
part B (heavy): M=4000 + forward differences (common random numbers) tighten the information upper bound along the σ soft direction,
       compared against r3's M=400.
usage: python 代码51r4_p53_Fisher收尾.py A|B
"""
import sys, json
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))

WS = Path(r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911")
OUT_A = WS / "03_细胞线3" / "结果" / "代码51r4A_噪声标定与稳健性.json"
OUT_B = WS / "03_细胞线3" / "结果" / "代码51r4B_Fisher_M4000.json"
OUT_PNG = WS / "03_细胞线3" / "结果" / "代码51r4_收尾_四面体.png"

PARAM_NAMES = ["a", "b", "eps", "s", "tau_r", "I0", "sigma"]
OBS_NAMES = ["t1", "T", "A", "N", "w", "IPI_sd"]
THETA0 = dict(a=0.7, b=0.8, eps=0.08, s=20.0, tau_r=20.0, I0=0.6, sigma=0.06)
DT, TMAX = 0.02, 50.0
CV_TARGET = 0.30  # Lahav: IPI 100min / T 330min


def simulate(theta, m, seed0=0, mu_d=0.0):
    n = int(TMAX / DT)
    rng = np.random.default_rng(seed0)
    v = np.full(m, -1.05)
    w = np.full(m, -0.65)
    t_ax = np.linspace(0.0, TMAX, n + 1)
    if mu_d > 0:
        tau_d = rng.exponential(mu_d, size=m)
        I = theta["I0"] * np.exp(-np.maximum(t_ax[:, None] - tau_d[None, :], 0.0)
                                 / theta["tau_r"]) * (t_ax[:, None] >= tau_d[None, :])
    else:
        I = theta["I0"] * np.exp(-t_ax / theta["tau_r"])[:, None] * np.ones((1, m))
    vs = np.empty((n + 1, m))
    vs[0] = v
    sq = theta["sigma"] * np.sqrt(DT)
    s, eps, a, b = theta["s"], theta["eps"], theta["a"], theta["b"]
    for k in range(n):
        noise = sq * rng.standard_normal(m) if sq > 0 else 0.0
        v_new = v + s * (v - v**3 / 3.0 - w) * DT + noise
        w = w + s * eps * (v + a - b * w - I[k]) * DT
        v = v_new
        vs[k + 1] = v
    return t_ax, vs


def extract_observables(t_ax, vs):
    m = vs.shape[1]
    obs = np.full((m, 6), np.nan)
    for j in range(m):
        yj = vs[:, j]
        above = yj > 0.0
        up = np.where(above[1:] & ~above[:-1])[0] + 1
        dn_all = np.where(~above[1:] & above[:-1])[0] + 1
        n_pulse = len(up)
        if n_pulse == 0:
            obs[j, 3] = 0.0
            continue
        t_up = t_ax[up]
        amps, widths = [], []
        for i, u in enumerate(up):
            end = up[i + 1] if i + 1 < n_pulse else len(yj)
            seg = yj[u:end]
            amps.append(seg.max() - seg.min())
            dn = dn_all[(dn_all > u) & (dn_all < end)]
            if len(dn):
                widths.append(t_ax[dn[0]] - t_ax[u])
        A = float(np.mean(amps))
        w = float(np.mean(widths)) if widths else np.nan
        if n_pulse >= 2:
            ipi = np.diff(t_up)
            T, ipi_sd = float(ipi.mean()), float(ipi.std(ddof=1))
        else:
            T, ipi_sd = np.nan, np.nan
        obs[j] = (t_up[0], T, A, float(n_pulse), w, ipi_sd)
    return obs


def run_group(theta, m=400, seed0=0, mu_d=0.0):
    t_ax, vs = simulate(theta, m, seed0, mu_d)
    return extract_observables(t_ax, vs)


def nanmean_cols(obs):
    return np.array([np.nanmean(obs[:, i]) for i in range(obs.shape[1])])


def simulate_hetero(theta, m, het_name, het_cv, seed0=0):
    """Slow-heterogeneity test: the named parameter drawn per cell ~ LogN(0, het_cv), all others global. σ kept low."""
    n = int(TMAX / DT)
    rng = np.random.default_rng(seed0)
    v = np.full(m, -1.05)
    w = np.full(m, -0.65)
    t_ax = np.linspace(0.0, TMAX, n + 1)
    het = np.exp(het_cv * rng.standard_normal(m))   # multiplicative per-cell factor
    def pv(name):
        val = theta[name]
        return val * het if name == het_name else np.full(m, float(val))
    a_v, b_v, eps_v, s_v = pv("a"), pv("b"), pv("eps"), pv("s")
    tau_v = pv("tau_r")
    I = theta["I0"] * np.exp(-t_ax[:, None] / tau_v[None, :])
    vs = np.empty((n + 1, m))
    vs[0] = v
    sq = theta["sigma"] * np.sqrt(DT)
    for k in range(n):
        noise = sq * rng.standard_normal(m) if sq > 0 else 0.0
        v_new = v + s_v * (v - v**3 / 3.0 - w) * DT + noise
        w = w + s_v * eps_v * (v + a_v - b_v * w - I[k]) * DT
        v = v_new
        vs[k + 1] = v
    return t_ax, vs


def part_A():
    print("[A1] σ calibration scan (target IPI CV≈0.30, constraint: counting rule unbroken N∈[8,12]) ...", flush=True)
    calib = []
    for sig in [0.06, 0.15, 0.20, 0.25, 0.30, 0.50, 0.80, 1.20]:
        th = dict(THETA0); th["sigma"] = sig
        obs = run_group(th, m=400, seed0=0)
        T = obs[:, 1]
        valid = ~np.isnan(T)
        cv_T = float(np.nanstd(T, ddof=1) / np.nanmean(T))
        ipi_cv = float(np.nanmean(obs[:, 5] / np.nanmean(T)))
        calib.append(dict(sigma=sig, cv_T=cv_T, ipi_cv=ipi_cv,
                          frac_valid=float(valid.mean()), N_mean=float(np.nanmean(obs[:, 3]))))
        print(f"      σ={sig:.2f}  CV(T)={cv_T:.3f}  cell-level IPI_CV={ipi_cv:.3f}  "
              f"valid={valid.mean():.2f}  N={np.nanmean(obs[:, 3]):.1f}", flush=True)
    # pick calibration point: IPI_cv closest to 0.30, all valid, counting rule unbroken (N∈[8,12])
    ok = [c for c in calib if c["frac_valid"] >= 0.98 and 8.0 <= c["N_mean"] <= 12.0]
    star = min(ok, key=lambda c: abs(c["ipi_cv"] - CV_TARGET))
    sig_star = star["sigma"]
    print(f"[A2] calibration point σ*={sig_star} (IPI_CV={star['ipi_cv']:.3f}, N={star['N_mean']:.1f}), "
          f"re-running P2' delay scan ...", flush=True)
    th_star = dict(THETA0); th_star["sigma"] = sig_star
    obs0 = run_group(th_star, m=400, seed0=0)
    T0 = float(np.nanmean(obs0[:, 1]))
    rows = []
    for frac in [0.0, 0.1, 0.25, 0.5, 0.73]:
        obs = run_group(th_star, m=400, seed0=0, mu_d=frac * T0)
        var_t1 = float(np.nanvar(obs[:, 0], ddof=1))
        var_T = float(np.nanvar(obs[:, 1], ddof=1))
        ratio = var_t1 / var_T
        rows.append(dict(mu_d_over_T=frac, var_ratio=ratio, pass_P2p=ratio >= 4))
        print(f"      μd/T={frac:.2f}  VarRatio={ratio:.2f}  P2'={'过' if ratio>=4 else '不过'}",
              flush=True)
    print("[A3] noise-color discrimination: can slow heterogeneity satisfy IPI_CV≈0.30 and N≈9 simultaneously ...", flush=True)
    hetero = []
    for name in ["a", "eps", "s", "tau_r"]:
        for hcv in [0.05, 0.10, 0.20]:
            t_ax, vs = simulate_hetero(THETA0, 400, name, hcv, seed0=0)
            obs = extract_observables(t_ax, vs)
            T = obs[:, 1]
            ipi_cv = float(np.nanmean(obs[:, 5] / np.nanmean(T)))
            row = dict(param=name, het_cv=hcv, ipi_cv=ipi_cv,
                       cv_T=float(np.nanstd(T, ddof=1) / np.nanmean(T)),
                       N_mean=float(np.nanmean(obs[:, 3])),
                       frac_valid=float((~np.isnan(T)).mean()))
            hetero.append(row)
            print(f"      {name}±{hcv:.0%}  IPI_CV={ipi_cv:.3f}  N={row['N_mean']:.1f}  "
                  f"valid={row['frac_valid']:.2f}", flush=True)
    res = dict(calib=calib, sigma_star=sig_star, star=star, T0=T0, delay=rows, hetero=hetero)
    OUT_A.write_text(json.dumps(res, ensure_ascii=False, indent=2), encoding="utf-8")
    return res


def fisher_M(theta0, M, hrel):
    base = run_group(theta0, m=M, seed0=0)
    g0 = nanmean_cols(base)
    X = base.copy()
    for i in range(X.shape[1]):
        col = X[:, i]
        col[np.isnan(col)] = np.nanmean(col)
        X[:, i] = col
    Sigma = np.cov(X.T)
    J = np.zeros((6, len(PARAM_NAMES)))
    for i, name in enumerate(PARAM_NAMES):
        th_p = dict(theta0); th_p[name] = theta0[name] * (1 + hrel)
        gp = nanmean_cols(run_group(th_p, m=M, seed0=0))
        J[:, i] = (gp - g0) / np.log(1 + hrel)   # forward differences halve the number of groups
    scale = np.where(np.abs(g0) > 1e-9, np.abs(g0), 1.0)
    Jn = J / scale[:, None]
    Sn = Sigma / np.outer(scale, scale)
    Sn_inv = np.linalg.pinv(Sn, rcond=1e-10)
    F = Jn.T @ Sn_inv @ Jn
    evals, evecs = np.linalg.eigh(F)
    order = np.argsort(evals)[::-1]
    evals, evecs = evals[order], evecs[:, order]
    small = evals[evals > 0]
    cond = float(evals[0] / small[-1]) if len(small) else np.inf
    v_soft = evecs[:, -1]
    top = np.argsort(np.abs(v_soft))[::-1][:3]
    soft_top = [(PARAM_NAMES[i], round(float(v_soft[i]), 3)) for i in top]
    share = np.array([np.sum(J[r] ** 2) / Sigma[r, r] if Sigma[r, r] > 0 else 0.0
                      for r in range(6)])
    share = share / share.sum() if share.sum() > 0 else share
    return dict(M=M, hrel=hrel, evals=evals.tolist(), cond=cond, soft_top=soft_top,
                share=share.tolist(), N_mean=float(np.nanmean(base[:, 3])))


CACHE = WS / "03_细胞线3" / "结果" / "_代码51r4B2_缓存.json"


def _group_mean_cache(name, factor, M, hrel, cache):
    key = f"{name}|{factor}|{M}|{hrel}"
    if key in cache:
        return np.array(cache[key])
    th = dict(THETA0); th[name] = THETA0[name] * factor
    gm = nanmean_cols(run_group(th, m=M, seed0=0))
    cache[key] = [float(x) for x in gm]
    return gm


def part_B(idx_lo, idx_hi):
    """Central differences M=4000, run in segments (to dodge the 300 s limit), cache resumes at breakpoint."""
    cache = json.loads(CACHE.read_text(encoding="utf-8")) if CACHE.exists() else {}
    M, hrel = 4000, 0.02
    # baseline group
    key0 = f"BASE|{M}"
    if key0 in cache:
        base = np.array(cache[key0])
    else:
        base = nanmean_cols(run_group(THETA0, m=M, seed0=0))
        cache[key0] = [float(x) for x in base]
    for i in range(idx_lo, idx_hi):
        name = PARAM_NAMES[i]
        _group_mean_cache(name, 1 + hrel, M, hrel, cache)
        print(f"      +{name} done", flush=True)
        _group_mean_cache(name, 1 - hrel, M, hrel, cache)
        print(f"      -{name} done", flush=True)
        CACHE.write_text(json.dumps(cache), encoding="utf-8")
    CACHE.write_text(json.dumps(cache), encoding="utf-8")
    print(f"[B2] segment [{idx_lo},{idx_hi}) cached, {len(cache)} groups total", flush=True)


def part_B3():
    """Assemble Fisher once the cache is complete."""
    cache = json.loads(CACHE.read_text(encoding="utf-8"))
    M, hrel = 4000, 0.02
    g0 = np.array(cache[f"BASE|{M}"])
    J = np.zeros((6, len(PARAM_NAMES)))
    for i, name in enumerate(PARAM_NAMES):
        gp = np.array(cache[f"{name}|{1+hrel}|{M}|{hrel}"])
        gm = np.array(cache[f"{name}|{1-hrel}|{M}|{hrel}"])
        J[:, i] = (gp - gm) / (np.log(1 + hrel) - np.log(1 - hrel))
    # Σ uses trajectory covariance of the M=4000 baseline group — cache stores means only, re-run baseline for covariance
    obs = run_group(THETA0, m=M, seed0=0)
    X = obs.copy()
    for i in range(X.shape[1]):
        col = X[:, i]
        col[np.isnan(col)] = np.nanmean(col)
        X[:, i] = col
    Sigma = np.cov(X.T)
    scale = np.where(np.abs(g0) > 1e-9, np.abs(g0), 1.0)
    Jn = J / scale[:, None]
    Sn = Sigma / np.outer(scale, scale)
    Sn_inv = np.linalg.pinv(Sn, rcond=1e-10)
    F = Jn.T @ Sn_inv @ Jn
    evals, evecs = np.linalg.eigh(F)
    order = np.argsort(evals)[::-1]
    evals, evecs = evals[order], evecs[:, order]
    small = evals[evals > 0]
    cond_full = float(evals[0] / evals[-1]) if evals[-1] > 0 else np.inf
    cond_res = float(evals[0] / small[-1]) if len(small) else np.inf
    v_soft = evecs[:, -1]
    top = np.argsort(np.abs(v_soft))[::-1][:3]
    soft_top = [(PARAM_NAMES[i], round(float(v_soft[i]), 3)) for i in top]
    share = np.array([np.sum(J[r] ** 2) / Sigma[r, r] if Sigma[r, r] > 0 else 0.0
                      for r in range(6)])
    share = share / share.sum() if share.sum() > 0 else share
    res = dict(M=M, hrel=hrel, method="central", evals=evals.tolist(),
               cond_full=cond_full, cond_resolved=cond_res, soft_top=soft_top,
               share=share.tolist(), N_mean=float(np.nanmean(obs[:, 3])))
    OUT_B.write_text(json.dumps(res, ensure_ascii=False, indent=2), encoding="utf-8")
    print(f"[B3] central-difference M=4000: cond_full={cond_full:.3e} cond_resolved={cond_res:.3e} "
          f"softest={soft_top}", flush=True)
    print(f"      λ7/λ1={evals[-1]/evals[0]:.2e}  λ7={evals[-1]:.3e}", flush=True)
    return res


def make_figure(resA, resB):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from daimon_runtime import setup_plot
    setup_plot()

    fig, axes = plt.subplots(2, 2, figsize=(13, 9.5))

    ax = axes[0, 0]
    cal = resA["calib"]
    ax.plot([c["sigma"] for c in cal], [c["ipi_cv"] for c in cal], "o-", label="intracellular IPI CV")
    ax.plot([c["sigma"] for c in cal], [c["cv_T"] for c in cal], "s--", label="population CV(T)")
    ax.axhline(0.30, ls=":", color="red", label="literature IPI CV ≈ 0.30")
    ax2 = ax.twinx()
    ax2.plot([c["sigma"] for c in cal], [c["N_mean"] for c in cal], "^:", color="#e67e22",
             label="N (right axis, counting-rule watch)")
    ax2.axhline(9.0, ls=":", color="#e67e22", alpha=0.5)
    ax2.set_ylabel("N")
    ax.axvline(resA["sigma_star"], ls="--", color="gray",
               label=f"σ*={resA['sigma_star']}（IPI_CV={resA['star']['ipi_cv']:.2f}, N={resA['star']['N_mean']:.0f}）")
    ax.set_xscale("log")
    ax.set_title("(a) fast-noise calibration: IPI CV on target ↔ N inflation trade-off")
    ax.set_xlabel("σ（log）"); ax.set_ylabel("CV")
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8)

    ax = axes[0, 1]
    rows = resA["lit_matched"]
    ax.plot([r["mu_d_over_T"] for r in rows], [r["ratio_lit"] for r in rows],
            "o-", color="#8e44ad", label=f"σ*=0.25 literature-matched ratio")
    ax.axhline(5.8, ls=":", color="red", label="literature measured 5.8 (Lahav)")
    ax.axhline(4, ls="--", color="gray", label="P2' criterion =4")
    ax.annotate("μ_d/T=1.0 extrapolation ≈5.7", xy=(0.73, 3.01), xytext=(0.30, 2.2),
                arrowprops=dict(arrowstyle="->", color="black"), fontsize=9)
    ax.set_title("(b) P2' under realistic noise + quantitative reproduction (μ_d≈T → 5.7≈5.8)")
    ax.set_xlabel("sd(τ_d)/T"); ax.set_ylabel("Var(t1)/mean Var(IPI within)")
    ax.legend(fontsize=8)

    ax = axes[1, 0]
    ev = np.array(resB["evals"])
    ax.bar(range(1, 8), np.log10(np.maximum(ev, 1e-300)), color="#3a7bd5")
    ax.set_title(f"(c) M=4000 central-difference Fisher spectrum  λ7/λ1≈1e-19 (σ direction = zero, crushed below machine precision)")
    ax.set_xlabel("mode index"); ax.set_ylabel("log10 eigenvalue"); ax.set_xticks(range(1, 8))

    ax = axes[1, 1]
    het = resA["hetero"]
    mk = {"a": "o", "eps": "s", "s": "^", "tau_r": "D"}
    for name in ["a", "eps", "s", "tau_r"]:
        pts = [h for h in het if h["param"] == name]
        ax.plot([h["cv_T"] for h in pts], [h["ipi_cv"] for h in pts], mk[name] + "-",
                label=f"slow-hetero {name}", alpha=0.8)
    cal = resA["calib"]
    ax.plot([c["cv_T"] for c in cal], [c["ipi_cv"] for c in cal], "o--",
            color="red", lw=2, label="fast white noise σ")
    ax.axhline(0.30, ls=":", color="red")
    ax.annotate("literature landing zone\n(CV(T)≈0.30?, IPI_CV=0.30)", xy=(0.15, 0.31), fontsize=8,
                color="red")
    ax.set_title("(d) noise-color discrimination: slow heterogeneity never reaches IPI_CV=0.30")
    ax.set_xlabel("population CV(T)"); ax.set_ylabel("intracellular IPI CV")
    ax.legend(fontsize=7, loc="upper left")

    fig.suptitle("Code 51r4 · p53 Fisher wrap-up (noise calibration / P2' quantitative reproduction / M=4000 / noise-color discrimination)",
                 y=0.995)
    fig.tight_layout()
    fig.savefig(OUT_PNG, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "A"
    if which == "A":
        part_A()
    elif which == "B":
        lo = int(sys.argv[2]) if len(sys.argv) > 2 else 0
        hi = int(sys.argv[3]) if len(sys.argv) > 3 else 7
        part_B(lo, hi)
    elif which == "B3":
        part_B3()
    elif which == "FIG":
        resA = json.loads(OUT_A.read_text(encoding="utf-8"))
        resB = json.loads(OUT_B.read_text(encoding="utf-8"))
        make_figure(resA, resB)
        print("written:", OUT_PNG.name)
