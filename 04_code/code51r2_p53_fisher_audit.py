# -*- coding: utf-8 -*-
"""
Code 51r2 p53 pulse-generator Fisher identifiability audit (pathfinding round 2)
Fixes: (1) pulse-width w extraction off-by-one (in r1 the w row was dead -> the three zero directions were an artifact)
       (2) k_obs moved out of the dynamics parameter vector (absolute-scale degeneracy belongs to the population level;
           tested separately with a sigma_k corruption test, same convention as code 50 leg 1)
       (3) P2 (t1 dispersed / T precise) changed to a threshold scan: the proposition is refined to
           "holds only in the near-threshold noise-triggered regime"; literature facts locate the p53 operating point in reverse
Outputs: 7-dynamics-parameter Fisher spectrum / heatmap / channel information shares, I0-scan Var(t1)/Var(T),
         counting law, I(D;N), sigma_k corruption control
"""
import sys, json
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))

WS = Path(r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911")
OUT_JSON = WS / "03_细胞线3" / "结果" / "代码51r2_p53_Fisher审计_结果.json"
OUT_PNG = WS / "03_细胞线3" / "结果" / "代码51r2_p53_Fisher审计_六面体.png"

PARAM_NAMES = ["a", "b", "eps", "s", "tau_r", "I0", "sigma"]
OBS_NAMES = ["t1", "T", "A", "N", "w", "IPI_sd"]
THETA0 = dict(a=0.7, b=0.8, eps=0.08, s=20.0, tau_r=20.0, I0=0.6, sigma=0.06)
DT, TMAX = 0.02, 50.0
M = 400
IC = 0.33  # FHN excitability threshold (memo: Ic≈0.33)


def simulate(theta, m=M, seed0=0):
    n = int(TMAX / DT)
    rng = np.random.default_rng(seed0)
    v = np.full(m, -1.05)
    w = np.full(m, -0.65)
    t_ax = np.linspace(0.0, TMAX, n + 1)
    I = theta["I0"] * np.exp(-t_ax / theta["tau_r"])
    vs = np.empty((n + 1, m))
    vs[0] = v
    sq = theta["sigma"] * np.sqrt(DT)
    s, eps, a, b = theta["s"], theta["eps"], theta["a"], theta["b"]
    for k in range(n):
        noise = sq * rng.standard_normal(m) if sq > 0 else 0.0
        v = v + s * (v - v**3 / 3.0 - w) * DT + noise
        w = w + s * eps * (v + a - b * w - I[k]) * DT
        vs[k + 1] = v
    return t_ax, vs


def extract_observables(t_ax, vs):
    """(t1, T, A, N, w, IPI_sd). Crossing-index fix: up/down aligned with +1."""
    m = vs.shape[1]
    obs = np.full((m, 6), np.nan)
    for j in range(m):
        yj = vs[:, j]
        above = yj > 0.0
        up = np.where(above[1:] & ~above[:-1])[0] + 1   # index of the first True
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


def run_group(theta, m=M, seed0=0):
    t_ax, vs = simulate(theta, m, seed0)
    return extract_observables(t_ax, vs)


def nanmean_cols(obs):
    return np.array([np.nanmean(obs[:, i]) for i in range(obs.shape[1])])


def fisher_audit(hrel):
    base = run_group(THETA0, seed0=0)
    g0 = nanmean_cols(base)
    X = base.copy()
    for i in range(X.shape[1]):
        col = X[:, i]
        col[np.isnan(col)] = np.nanmean(col)
        X[:, i] = col
    Sigma = np.cov(X.T)
    p = len(PARAM_NAMES)
    J = np.zeros((6, p))
    for i, name in enumerate(PARAM_NAMES):
        th_p = dict(THETA0); th_m = dict(THETA0)
        th_p[name] = THETA0[name] * (1 + hrel)
        th_m[name] = THETA0[name] * (1 - hrel)
        gp = nanmean_cols(run_group(th_p, seed0=0))
        gm = nanmean_cols(run_group(th_m, seed0=0))
        J[:, i] = (gp - gm) / (np.log(1 + hrel) - np.log(1 - hrel))
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
    # Channel information share (diagonal approximation): share_r ∝ Σ_i J[r,i]^2 / Sigma[r,r]
    share = np.array([np.sum(J[r] ** 2) / Sigma[r, r] if Sigma[r, r] > 0 else 0.0
                      for r in range(6)])
    share = share / share.sum() if share.sum() > 0 else share
    return dict(hrel=hrel, g0=g0.tolist(), evals=evals.tolist(), evecs=evecs.tolist(),
                cond=cond, share=share.tolist(),
                var_t1=float(np.nanvar(base[:, 0], ddof=1)),
                var_T=float(np.nanvar(base[:, 1], ddof=1)),
                N_mean=float(np.nanmean(base[:, 3])),
                Sigma_diag=np.diag(Sigma).tolist())


def threshold_scan():
    """I0 from near-threshold to far above threshold: behavior of Var(t1)/Var(T) -> regime localization of P2."""
    rows = []
    for I0 in [0.34, 0.36, 0.40, 0.45, 0.50, 0.60, 0.72]:
        th = dict(THETA0); th["I0"] = I0
        obs = run_group(th, m=400, seed0=0)
        n_valid = int(np.sum(~np.isnan(obs[:, 1])))
        var_t1 = float(np.nanvar(obs[:, 0], ddof=1))
        var_T = float(np.nanvar(obs[:, 1], ddof=1)) if n_valid > 10 else np.nan
        ratio = var_t1 / var_T if var_T and var_T > 0 else np.nan
        rows.append(dict(I0=I0, supercrit=I0 / IC, N_mean=float(np.nanmean(obs[:, 3])),
                         frac_multi=n_valid / 400, var_t1=var_t1, var_T=var_T,
                         ratio=ratio, sd_t1=float(np.sqrt(var_t1)),
                         sd_T=float(np.sqrt(var_T)) if var_T == var_T else np.nan))
    return rows


def channel_capacity():
    doses = np.array([0.30, 0.36, 0.42, 0.50, 0.60, 0.72, 0.85, 1.00, 1.20])
    counts = []
    for d in doses:
        th = dict(THETA0); th["I0"] = float(d)
        obs = run_group(th, m=300, seed0=1000)
        counts.append(obs[:, 3])
    n_max = int(max(c.max() for c in counts)) + 1
    P = np.zeros((len(doses), n_max))
    for i, c in enumerate(counts):
        for k in range(n_max):
            P[i, k] = np.mean(c == k)
    Pn = P.mean(axis=0)
    I_DN = 0.0
    for i in range(len(doses)):
        for k in range(n_max):
            if P[i, k] > 1e-12 and Pn[k] > 1e-12:
                I_DN += (1 / len(doses)) * P[i, k] * np.log2(P[i, k] / Pn[k])
    return dict(doses=doses.tolist(), N_means=[float(c.mean()) for c in counts],
                N_sds=[float(c.std(ddof=1)) for c in counts],
                I_DN_bits=float(I_DN), H_D_bits=float(np.log2(len(doses))))


def gain_corruption():
    """sigma_k corruption: per-cell k_i~LogN(0,sigma_k) multiplied onto the observation y=k_i*v.
    Watch the cross-cell CV of A inflate with sigma_k while the CV of T stays immune — minimal control with the same convention as code 50 leg 1."""
    rows = []
    for sk in [0.0, 0.2, 0.4, 0.6]:
        rng = np.random.default_rng(5)
        t_ax, vs = simulate(THETA0, m=400, seed0=0)
        ks = np.exp(sk * rng.standard_normal(400))
        obs = extract_observables(t_ax, vs * ks[None, :])
        A = obs[:, 2]; T = obs[:, 1]
        rows.append(dict(sigma_k=sk,
                         cv_A=float(np.nanstd(A, ddof=1) / np.nanmean(A)),
                         cv_T=float(np.nanstd(T, ddof=1) / np.nanmean(T))))
    return rows


def make_figure(res, scan, cap, gain):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from daimon_runtime import setup_plot
    setup_plot()

    fig, axes = plt.subplots(2, 3, figsize=(17, 9.5))

    ax = axes[0, 0]
    t_ax, vs = simulate(THETA0, m=6, seed0=7)
    for j in range(3):
        ax.plot(t_ax, vs[:, j] + 3.2 * j, lw=0.8)
    ax.set_title("(a) FHN pulse-train example (3 traces, vertically offset)")
    ax.set_xlabel("time (au)"); ax.set_ylabel("v (offset)")

    ax = axes[0, 1]
    ev = np.array(res["evals"])
    ax.bar(range(1, 8), np.log10(np.maximum(ev, 1e-300)), color="#3a7bd5")
    ax.set_title(f"(b) Fisher eigenspectrum (7 dynamics parameters) cond={res['cond']:.1e}")
    ax.set_xlabel("mode index"); ax.set_ylabel("log10 eigenvalue"); ax.set_xticks(range(1, 8))

    ax = axes[0, 2]
    E = np.abs(np.array(res["evecs"]))
    im = ax.imshow(E, aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_yticks(range(7)); ax.set_yticklabels(PARAM_NAMES)
    ax.set_xticks(range(7)); ax.set_xticklabels([f"v{i+1}" for i in range(7)])
    ax.set_title("(c) eigenvector composition (v1 stiffest -> v7 softest)")
    fig.colorbar(im, ax=ax, fraction=0.046)

    ax = axes[1, 0]
    sh = np.array(res["share"])
    ax.bar(OBS_NAMES, sh, color=["#c0392b" if o == "A" else "#27ae60" for o in OBS_NAMES])
    ax.set_title("(d) channel information share (diagonal approx.): how dead is amplitude A")
    ax.set_ylabel("information share")

    ax = axes[1, 1]
    sc = [r for r in scan if r["ratio"] == r["ratio"]]
    ax.plot([r["supercrit"] for r in sc], [r["ratio"] for r in sc], "o-", color="#8e44ad",
            label="Var(t1)/Var(T)")
    ax.axhline(4, ls="--", color="gray", label="P2 criterion =4 (literature (240/100)²≈5.8)")
    ax.axvline(1.0, ls=":", color="red", label="threshold Ic")
    ax.set_title("(e) dispersion only near threshold: first-pulse/interval variance ratio vs supercriticality")
    ax.set_xlabel("I0 / Ic (supercritical multiple)"); ax.set_ylabel("Var(t1)/Var(T)")
    ax.legend()

    ax = axes[1, 2]
    d = np.array(cap["doses"]); nm = np.array(cap["N_means"]); nsd = np.array(cap["N_sds"])
    ax.errorbar(np.log(d), nm, yerr=nsd, fmt="o-", capsize=3, label="N(I0)")
    mask = nm > 0
    cf = np.polyfit(np.log(d[mask]), nm[mask], 1)
    ax.plot(np.log(d), np.polyval(cf, np.log(d)), "--", color="gray",
            label=f"log-law slope={cf[0]:.2f}")
    ax2 = ax.twinx()
    g = gain
    ax2.plot([r["sigma_k"] for r in g], [r["cv_A"] for r in g], "s--", color="#c0392b",
             label="CV(A) vs σ_k (right axis)")
    ax2.plot([r["sigma_k"] for r in g], [r["cv_T"] for r in g], "^--", color="#27ae60",
             label="CV(T) vs σ_k (right axis)")
    ax2.set_ylabel("cross-cell CV")
    ax.set_xlabel("ln I0 (sigma_k panel: read bottom categories)"); ax.set_ylabel("N")
    ax.set_title(f"(f) counting law + gain corruption control  I(D;N)={cap['I_DN_bits']:.2f}/{cap['H_D_bits']:.2f} bits")
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, fontsize=8, loc="upper left")

    fig.suptitle("Code 51r2 · p53 pulse-generator Fisher audit round 2 (w fixed, 7 dynamics parameters)", y=0.995)
    fig.tight_layout()
    fig.savefig(OUT_PNG, bbox_inches="tight")
    plt.close(fig)


def main():
    print("[1/4] Fisher h=2% ...", flush=True)
    res2 = fisher_audit(0.02)
    print(f"      cond={res2['cond']:.3e} N={res2['N_mean']:.1f} "
          f"share={dict(zip(OBS_NAMES, [round(s,3) for s in res2['share']]))}", flush=True)
    print("[2/4] Fisher h=5% ...", flush=True)
    res5 = fisher_audit(0.05)
    print(f"      cond={res5['cond']:.3e}", flush=True)
    print("[3/4] threshold scan ...", flush=True)
    scan = threshold_scan()
    for r in scan:
        rs = f"{r['ratio']:.2f}" if r["ratio"] == r["ratio"] else "NA"
        print(f"      I0={r['I0']:.2f} (x{r['supercrit']:.2f}Ic) N={r['N_mean']:.1f} "
              f"multi={r['frac_multi']:.2f} sd(t1)={r['sd_t1']:.3f} ratio={rs}", flush=True)
    print("[4/4] channel capacity + gain corruption ...", flush=True)
    cap = channel_capacity()
    gain = gain_corruption()
    print(f"      I(D;N)={cap['I_DN_bits']:.3f}/{cap['H_D_bits']:.3f} bits", flush=True)
    for r in gain:
        print(f"      σ_k={r['sigma_k']:.1f}  CV(A)={r['cv_A']:.3f}  CV(T)={r['cv_T']:.3f}", flush=True)

    out = dict(theta0=THETA0, param_names=PARAM_NAMES, obs_names=OBS_NAMES,
               fisher_h2=res2, fisher_h5=res5, scan=scan, capacity=cap, gain=gain)
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    make_figure(res2, scan, cap, gain)
    print("Written:", OUT_JSON.name, OUT_PNG.name, flush=True)


if __name__ == "__main__":
    main()
