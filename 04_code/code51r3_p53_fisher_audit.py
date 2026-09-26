# -*- coding: utf-8 -*-
"""
代码51r3 p53脉冲发生器 Fisher 审计（探路三轮）
修正：①积分器恢复同步欧拉（r2 误为半隐式，N 漂移 9.0→8.3）
      ②新增上游随机感知延迟级 τ_d~Exp(μ_d)：I(t)=I0·exp(-(t-τ_d)/τ_r)·1[t≥τ_d]
P2'（跑前钉死）：Var(t1)/Var(T) ≥ 4 当且仅当存在上游随机延迟且其离散度
      sd(τ_d)/T ≳ O(0.1)；振荡器核心自身的触发离散度远不够（r2 实测 sd(t1)≈0.5%·T）。
      文献标定：t1 sd ≈ 0.73·T（240/330min）。
同时保留：7参数Fisher（h=2%/5%双跑+最软方向鉴定）、通道信息份额、
      计数律、I(D;N)、σ_k 腐蚀对照。
"""
import sys, json
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))

WS = Path(r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911")
OUT_JSON = WS / "03_细胞线3" / "结果" / "代码51r3_p53_Fisher审计_结果.json"
OUT_PNG = WS / "03_细胞线3" / "结果" / "代码51r3_p53_Fisher审计_六面体.png"

PARAM_NAMES = ["a", "b", "eps", "s", "tau_r", "I0", "sigma"]
OBS_NAMES = ["t1", "T", "A", "N", "w", "IPI_sd"]
THETA0 = dict(a=0.7, b=0.8, eps=0.08, s=20.0, tau_r=20.0, I0=0.6, sigma=0.06)
DT, TMAX = 0.02, 50.0
M = 400
IC = 0.33


def simulate(theta, m=M, seed0=0, mu_d=0.0):
    """同步 Euler-Maruyama；mu_d>0 时每细胞独立延迟 τ_d~Exp(mu_d)。"""
    n = int(TMAX / DT)
    rng = np.random.default_rng(seed0)
    v = np.full(m, -1.05)
    w = np.full(m, -0.65)
    t_ax = np.linspace(0.0, TMAX, n + 1)
    if mu_d > 0:
        tau_d = rng.exponential(mu_d, size=m)          # 每细胞感知延迟
        I = theta["I0"] * np.exp(-np.maximum(t_ax[:, None] - tau_d[None, :], 0.0)
                                 / theta["tau_r"]) * (t_ax[:, None] >= tau_d[None, :])
    else:
        tau_d = np.zeros(m)
        I = theta["I0"] * np.exp(-t_ax / theta["tau_r"])[:, None] * np.ones((1, m))
    vs = np.empty((n + 1, m))
    vs[0] = v
    sq = theta["sigma"] * np.sqrt(DT)
    s, eps, a, b = theta["s"], theta["eps"], theta["a"], theta["b"]
    for k in range(n):
        noise = sq * rng.standard_normal(m) if sq > 0 else 0.0
        v_new = v + s * (v - v**3 / 3.0 - w) * DT + noise   # 同步更新：w 用旧 v
        w = w + s * eps * (v + a - b * w - I[k]) * DT
        v = v_new
        vs[k + 1] = v
    return t_ax, vs, tau_d


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


def run_group(theta, m=M, seed0=0, mu_d=0.0):
    t_ax, vs, _ = simulate(theta, m, seed0, mu_d)
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
    J = np.zeros((6, len(PARAM_NAMES)))
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
    v_soft = evecs[:, -1]
    top = np.argsort(np.abs(v_soft))[::-1][:3]
    soft_top = [(PARAM_NAMES[i], round(float(v_soft[i]), 3)) for i in top]
    share = np.array([np.sum(J[r] ** 2) / Sigma[r, r] if Sigma[r, r] > 0 else 0.0
                      for r in range(6)])
    share = share / share.sum() if share.sum() > 0 else share
    return dict(hrel=hrel, g0=g0.tolist(), evals=evals.tolist(), evecs=evecs.tolist(),
                cond=cond, soft_top=soft_top, share=share.tolist(),
                N_mean=float(np.nanmean(base[:, 3])))


def delay_scan():
    """P2'：μ_d/T ∈ {0, 0.1, 0.25, 0.5, 0.73, 1.0}，判据 Var(t1)/Var(T)≥4。"""
    # 先用无延迟组估计 T0 用于归一
    obs0 = run_group(THETA0, seed0=0)
    T0 = float(np.nanmean(obs0[:, 1]))
    var_T0 = float(np.nanvar(obs0[:, 1], ddof=1))
    rows = []
    for frac in [0.0, 0.1, 0.25, 0.5, 0.73, 1.0]:
        mu_d = frac * T0
        obs = run_group(THETA0, seed0=0, mu_d=mu_d)
        var_t1 = float(np.nanvar(obs[:, 0], ddof=1))
        var_T = float(np.nanvar(obs[:, 1], ddof=1))
        cv_t1 = float(np.nanstd(obs[:, 0], ddof=1) / np.nanmean(obs[:, 0]))
        cv_T = float(np.nanstd(obs[:, 1], ddof=1) / np.nanmean(obs[:, 1]))
        rows.append(dict(mu_d_over_T=frac, var_ratio=var_t1 / var_T,
                         cv_t1=cv_t1, cv_T=cv_T, cv_ratio=cv_t1 / cv_T,
                         sd_t1_over_T=float(np.sqrt(var_t1)) / T0,
                         pass_P2p=(var_t1 / var_T) >= 4))
    return dict(T0=T0, var_T0=var_T0, rows=rows)


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
    rows = []
    for sk in [0.0, 0.2, 0.4, 0.6]:
        rng = np.random.default_rng(5)
        t_ax, vs, _ = simulate(THETA0, m=400, seed0=0)
        ks = np.exp(sk * rng.standard_normal(400))
        obs = extract_observables(t_ax, vs * ks[None, :])
        A = obs[:, 2]; T = obs[:, 1]
        rows.append(dict(sigma_k=sk,
                         cv_A=float(np.nanstd(A, ddof=1) / np.nanmean(A)),
                         cv_T=float(np.nanstd(T, ddof=1) / np.nanmean(T))))
    return rows


def make_figure(res, delay, cap, gain):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from daimon_runtime import setup_plot
    setup_plot()

    fig, axes = plt.subplots(2, 3, figsize=(17, 9.5))

    ax = axes[0, 0]
    t_ax, vs, _ = simulate(THETA0, m=6, seed0=7)
    for j in range(3):
        ax.plot(t_ax, vs[:, j] + 3.2 * j, lw=0.8)
    ax.set_title("(a) FHN 脉冲串（同步积分，N=%.1f）" % res["N_mean"])
    ax.set_xlabel("时间 (au)"); ax.set_ylabel("v（错位）")

    ax = axes[0, 1]
    ev = np.array(res["evals"])
    ax.bar(range(1, 8), np.log10(np.maximum(ev, 1e-300)), color="#3a7bd5")
    ax.set_title("(b) Fisher 特征谱  cond(h2)=%.1e" % res["cond"])
    ax.set_xlabel("模式序号"); ax.set_ylabel("log10 特征值"); ax.set_xticks(range(1, 8))

    ax = axes[0, 2]
    E = np.abs(np.array(res["evecs"]))
    im = ax.imshow(E, aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_yticks(range(7)); ax.set_yticklabels(PARAM_NAMES)
    ax.set_xticks(range(7)); ax.set_xticklabels([f"v{i+1}" for i in range(7)])
    ax.set_title("(c) 特征向量（最软 v7：%s）" % ", ".join(f"{n} {v}" for n, v in res["soft_top"][:2]))
    fig.colorbar(im, ax=ax, fraction=0.046)

    ax = axes[1, 0]
    sh = np.array(res["share"])
    ax.bar(OBS_NAMES, sh, color=["#c0392b" if o == "A" else "#27ae60" for o in OBS_NAMES])
    ax.set_title("(d) 通道信息份额：A=%.1f%%，时间类合计=%.1f%%"
                 % (100 * sh[2], 100 * (sh[0] + sh[1] + sh[3] + sh[4])))
    ax.set_ylabel("信息份额")

    ax = axes[1, 1]
    rows = delay["rows"]
    ax.plot([r["mu_d_over_T"] for r in rows], [r["var_ratio"] for r in rows],
            "o-", color="#8e44ad", label="Var(t1)/Var(T)")
    ax.axhline(4, ls="--", color="gray", label="P2' 判据 =4")
    ax.axvline(0.73, ls=":", color="red", label="文献 t1 sd ≈ 0.73·T")
    ax.set_yscale("log")
    ax.set_title("(e) P2'：上游随机感知延迟 → t1弥散/IPI精确 签名")
    ax.set_xlabel("sd(τ_d) / T"); ax.set_ylabel("Var(t1)/Var(T)（log）")
    ax.legend()

    ax = axes[1, 2]
    d = np.array(cap["doses"]); nm = np.array(cap["N_means"]); nsd = np.array(cap["N_sds"])
    ax.errorbar(np.log(d), nm, yerr=nsd, fmt="o-", capsize=3, label="N(I0)")
    mask = nm > 0
    cf = np.polyfit(np.log(d[mask]), nm[mask], 1)
    ax.plot(np.log(d), np.polyval(cf, np.log(d)), "--", color="gray",
            label="对数律 斜率=%.2f" % cf[0])
    ax.set_xlabel("ln I0"); ax.set_ylabel("N")
    ax.set_title("(f) 计数律 + I(D;N)=%.2f/%.2f bits；σ_k: CV(A)%.2f→%.2f, CV(T)钉死%.3f"
                 % (cap["I_DN_bits"], cap["H_D_bits"],
                    gain[0]["cv_A"], gain[-1]["cv_A"], gain[-1]["cv_T"]))
    ax.legend(loc="upper left", fontsize=8)

    fig.suptitle("代码51r3 · p53 Fisher 审计三轮（同步积分 + 上游随机感知级）", y=0.995)
    fig.tight_layout()
    fig.savefig(OUT_PNG, bbox_inches="tight")
    plt.close(fig)


def main():
    print("[1/4] Fisher h=2%% / h=5%% ...", flush=True)
    res2 = fisher_audit(0.02)
    res5 = fisher_audit(0.05)
    print(f"      N={res2['N_mean']:.1f} cond(h2)={res2['cond']:.2e} cond(h5)={res5['cond']:.2e}",
          flush=True)
    print(f"      最软方向(h2)={res2['soft_top']}  (h5)={res5['soft_top']}", flush=True)
    print(f"      share={dict(zip(OBS_NAMES, [round(s,3) for s in res2['share']]))}", flush=True)
    print("[2/4] P2' 延迟扫描 ...", flush=True)
    delay = delay_scan()
    print(f"      T0={delay['T0']:.3f}", flush=True)
    for r in delay["rows"]:
        print(f"      μd/T={r['mu_d_over_T']:.2f}  VarRatio={r['var_ratio']:.2f}  "
              f"CV(t1)/CV(T)={r['cv_ratio']:.2f}  P2'={'过' if r['pass_P2p'] else '不过'}",
              flush=True)
    print("[3/4] 信道容量 ...", flush=True)
    cap = channel_capacity()
    print(f"      I(D;N)={cap['I_DN_bits']:.3f}/{cap['H_D_bits']:.3f} bits", flush=True)
    print("[4/4] 增益腐蚀 ...", flush=True)
    gain = gain_corruption()
    for r in gain:
        print(f"      σ_k={r['sigma_k']:.1f}  CV(A)={r['cv_A']:.3f}  CV(T)={r['cv_T']:.3f}",
              flush=True)

    out = dict(theta0=THETA0, param_names=PARAM_NAMES, obs_names=OBS_NAMES,
               fisher_h2=res2, fisher_h5=res5, delay=delay, capacity=cap, gain=gain)
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    make_figure(res2, delay, cap, gain)
    print("已写出：", OUT_JSON.name, OUT_PNG.name, flush=True)


if __name__ == "__main__":
    main()
