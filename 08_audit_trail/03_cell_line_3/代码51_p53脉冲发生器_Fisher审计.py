# -*- coding: utf-8 -*-
"""
代码51 p53脉冲发生器 Fisher 可辨识性审计（探路）
谱系：p53可行性备忘录 §7 步骤2-3（已规划未执行）→ 本脚本接跑
模型：FHN约化 + 损伤指数修复（手稿方法节参数套：s=20, tau_r=20）
输出：Fisher特征谱 / 特征向量-参数热图 / 零方向对仿射群轨道cos / t1 vs T 方差比 / I(D;N)信道容量
纪律：预言跑前钉死在探路卡（P1: |cos|>=0.99 且 cond>=1e6；P2: Var(t1)/Var(T)>=4）
"""
import sys, json
from pathlib import Path
import numpy as np

sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))

WS = Path(r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911")
OUT_JSON = WS / "03_细胞线3" / "结果" / "代码51_p53_Fisher审计_结果.json"
OUT_PNG = WS / "03_细胞线3" / "结果" / "代码51_p53_Fisher审计_四面体.png"

# ---------- 模型 ----------
# v' = s*(v - v^3/3 - w) + sigma*xi
# w' = s*eps*(v + a - b*w - I(t)),  I(t)=I0*exp(-t/tau_r)
PARAM_NAMES = ["a", "b", "eps", "s", "tau_r", "I0", "sigma", "k_obs"]
THETA0 = dict(a=0.7, b=0.8, eps=0.08, s=20.0, tau_r=20.0, I0=0.6, sigma=0.06, k_obs=1.0)
DT, TMAX = 0.02, 50.0
M = 400  # 轨迹数/参数组


def simulate(theta, m=M, seed0=0):
    """向量化 Euler-Maruyama。返回 v 轨迹 (n_steps+1, m) 与时间轴。"""
    n = int(TMAX / DT)
    rng = np.random.default_rng(seed0)
    v = np.full(m, -1.05)  # 静息附近
    w = np.full(m, -0.65)
    t_ax = np.linspace(0.0, TMAX, n + 1)
    I = theta["I0"] * np.exp(-t_ax / theta["tau_r"])
    vs = np.empty((n + 1, m))
    vs[0] = v
    sq = theta["sigma"] * np.sqrt(DT)
    s, eps, a, b = theta["s"], theta["eps"], theta["a"], theta["b"]
    for k in range(n):
        noise = sq * rng.standard_normal(m) if sq > 0 else 0.0
        dv = s * (v - v**3 / 3.0 - w) * DT + noise
        dw = s * eps * (v + a - b * w - I[k]) * DT
        v = v + dv
        w = w + dw
        vs[k + 1] = v
    return t_ax, vs


def extract_observables(t_ax, vs, k_obs):
    """从 y=k_obs*v 提取 (t1, T, A, N, w, IPI_sd)。返回 (m,6)，不足2脉冲的T/IPI记NaN。"""
    y = k_obs * vs
    m = y.shape[1]
    obs = np.full((m, 6), np.nan)
    thr = 0.0
    for j in range(m):
        yj = y[:, j]
        above = yj > thr
        up = np.where((~above[:-1]) & (above[1:]))[0]
        if len(up) == 0:
            obs[j, 3] = 0.0
            continue
        t_up = t_ax[up]
        t1 = t_up[0]
        n_pulse = len(up)
        # 每脉冲：峰、谷（上穿点到下次上穿/终点窗口）、宽度
        amps, widths = [], []
        for i, u in enumerate(up):
            end = up[i + 1] if i + 1 < len(up) else len(yj) - 1
            seg = yj[u:end]
            amps.append(seg.max() - seg.min())
            dn = np.where((above[u:end])[1:] & ~(above[u:end])[:-1])[0]
            if len(dn):
                widths.append(t_ax[u + dn[0]] - t_ax[u])
        A = float(np.mean(amps)) if amps else np.nan
        w = float(np.mean(widths)) if widths else np.nan
        if n_pulse >= 2:
            ipi = np.diff(t_up)
            T = float(ipi.mean())
            ipi_sd = float(ipi.std(ddof=1)) if len(ipi) > 1 else 0.0
        else:
            T, ipi_sd = np.nan, np.nan
        obs[j] = (t1, T, A, float(n_pulse), w, ipi_sd)
    return obs


def run_group(theta, m=M, seed0=0):
    t_ax, vs = simulate(theta, m, seed0)
    return extract_observables(t_ax, vs, theta["k_obs"])


def nanmean_cols(obs):
    return np.array([np.nanmean(obs[:, i]) for i in range(obs.shape[1])])


def fisher_audit(hrel):
    """对数扰动 hrel，共同随机数（seed0=0..M-1 同源）中心差分。"""
    base = run_group(THETA0, seed0=0)
    g0 = nanmean_cols(base)
    # 协方差（NaN 成对删除近似：用 mask 填充列均值后估协方差——400条中NaN极少）
    X = base.copy()
    for i in range(X.shape[1]):
        col = X[:, i]
        col[np.isnan(col)] = np.nanmean(col)
        X[:, i] = col
    Sigma = np.cov(X.T)
    J = np.zeros((6, 8))
    for i, name in enumerate(PARAM_NAMES):
        th_p = dict(THETA0); th_m = dict(THETA0)
        th_p[name] = THETA0[name] * (1 + hrel)
        th_m[name] = THETA0[name] * (1 - hrel)
        gp = nanmean_cols(run_group(th_p, seed0=0))
        gm = nanmean_cols(run_group(th_m, seed0=0))
        dlog = np.log(1 + hrel) - np.log(1 - hrel)
        J[:, i] = (gp - gm) / dlog
    # Fisher（相对信息：可观测量除以其均值量级，避免单位主导）
    scale = np.where(np.abs(g0) > 1e-9, np.abs(g0), 1.0)
    Jn = J / scale[:, None]
    Sn = Sigma / np.outer(scale, scale)
    evals_S = np.linalg.eigvalsh(Sn)
    tol = evals_S.max() * 1e-10
    Sn_inv = np.linalg.pinv(Sn, rcond=tol if tol > 0 else 1e-12)
    F = Jn.T @ Sn_inv @ Jn
    evals, evecs = np.linalg.eigh(F)
    order = np.argsort(evals)[::-1]
    evals, evecs = evals[order], evecs[:, order]
    cond = evals[0] / evals[-1] if evals[-1] > 0 else np.inf
    # 零方向 vs k_obs 轴
    e_kobs = np.zeros(8); e_kobs[PARAM_NAMES.index("k_obs")] = 1.0
    v_min = evecs[:, -1]
    cos_kobs = abs(float(v_min @ e_kobs))
    # P2：t1 与 T 的跨细胞方差
    var_t1 = float(np.nanvar(base[:, 0], ddof=1))
    var_T = float(np.nanvar(base[:, 1], ddof=1))
    return dict(hrel=hrel, g0=g0.tolist(), evals=evals.tolist(), evecs=evecs.tolist(),
                cond=float(cond), cos_kobs=cos_kobs, var_t1=var_t1, var_T=var_T,
                var_ratio=var_t1 / var_T if var_T > 0 else np.inf,
                N_mean=float(np.nanmean(base[:, 3])), Sigma_diag=np.diag(Sigma).tolist())


def channel_capacity():
    """I0 扫描 → N 分布 → I(D;N)（均匀剂量先验，N天然离散）。"""
    doses = np.array([0.30, 0.36, 0.42, 0.50, 0.60, 0.72, 0.85, 1.00, 1.20])
    counts = []
    for d in doses:
        th = dict(THETA0); th["I0"] = float(d)
        obs = run_group(th, m=300, seed0=1000)
        counts.append(obs[:, 3])
    n_max = int(max(c.max() for c in counts)) + 1
    # P(N|D)
    P = np.zeros((len(doses), n_max))
    for i, c in enumerate(counts):
        for k in range(n_max):
            P[i, k] = np.mean(c == k)
    Pn = P.mean(axis=0)  # 均匀先验
    I_DN = 0.0
    for i in range(len(doses)):
        for k in range(n_max):
            if P[i, k] > 1e-12 and Pn[k] > 1e-12:
                I_DN += (1 / len(doses)) * P[i, k] * np.log2(P[i, k] / Pn[k])
    H_D = np.log2(len(doses))
    return dict(doses=doses.tolist(), N_means=[float(c.mean()) for c in counts],
                N_sds=[float(c.std(ddof=1)) for c in counts], I_DN_bits=float(I_DN), H_D_bits=float(H_D))


def make_figure(res2, res5, cap):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from daimon_runtime import setup_plot
    setup_plot()

    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    # (a) 示例轨迹 + 脉冲统计
    t_ax, vs = simulate(THETA0, m=6, seed0=7)
    ax = axes[0, 0]
    for j in range(3):
        ax.plot(t_ax, vs[:, j] + 3.2 * j, lw=0.8)
    ax.set_title("(a) FHN 脉冲串（3 条示例，纵向错位）")
    ax.set_xlabel("时间 (au)"); ax.set_ylabel("y（错位）")

    # (b) Fisher 特征谱
    ax = axes[0, 1]
    ev = np.array(res2["evals"])
    ax.bar(range(1, 9), np.log10(np.maximum(ev, 1e-300)), color="#3a7bd5")
    ax.set_title(f"(b) Fisher 特征谱  cond={res2['cond']:.2e}")
    ax.set_xlabel("模式序号"); ax.set_ylabel("log10 特征值")
    ax.set_xticks(range(1, 9))

    # (c) 特征向量-参数热图
    ax = axes[1, 0]
    E = np.abs(np.array(res2["evecs"]))
    im = ax.imshow(E, aspect="auto", cmap="viridis", vmin=0, vmax=1)
    ax.set_yticks(range(8)); ax.set_yticklabels(PARAM_NAMES)
    ax.set_xticks(range(8)); ax.set_xticklabels([f"v{i+1}" for i in range(8)])
    ax.set_title(f"(c) 特征向量组成（v8=最软）  cos(v8, k_obs)={res2['cos_kobs']:.4f}")
    ax.set_xlabel("模式（v1最硬 → v8最软）"); ax.set_ylabel("参数")
    fig.colorbar(im, ax=ax, fraction=0.046)

    # (d) 计数律 + 信道容量
    ax = axes[1, 1]
    d = np.array(cap["doses"]); nm = np.array(cap["N_means"]); nsd = np.array(cap["N_sds"])
    ax.errorbar(np.log(d), nm, yerr=nsd, fmt="o-", capsize=3, label="N(I₀) 均值±sd")
    mask = nm > 0
    if mask.sum() >= 2:
        cf = np.polyfit(np.log(d[mask]), nm[mask], 1)
        ax.plot(np.log(d), np.polyval(cf, np.log(d)), "--", color="gray",
                label=f"对数律拟合 斜率={cf[0]:.2f}")
    ax.set_title(f"(d) 计数律与信道容量  I(D;N)={cap['I_DN_bits']:.2f}/{cap['H_D_bits']:.2f} bits")
    ax.set_xlabel("ln I₀"); ax.set_ylabel("脉冲数 N")
    ax.legend()

    fig.suptitle("代码51 · p53 脉冲发生器 Fisher 审计（FHN代表元，手稿参数套）", y=0.995)
    fig.tight_layout()
    fig.savefig(OUT_PNG, bbox_inches="tight")
    plt.close(fig)


def main():
    print("[1/3] Fisher 审计 h=±2% ...", flush=True)
    res2 = fisher_audit(0.02)
    print(f"      cond={res2['cond']:.3e}  cos(v_min,k_obs)={res2['cos_kobs']:.4f}  "
          f"Var(t1)/Var(T)={res2['var_ratio']:.2f}  N_mean={res2['N_mean']:.1f}", flush=True)
    print("[2/3] Fisher 审计 h=±5%（稳健性）...", flush=True)
    res5 = fisher_audit(0.05)
    print(f"      cond={res5['cond']:.3e}  cos(v_min,k_obs)={res5['cos_kobs']:.4f}", flush=True)
    print("[3/3] 信道容量 I(D;N) ...", flush=True)
    cap = channel_capacity()
    print(f"      I(D;N)={cap['I_DN_bits']:.3f} bits / H(D)={cap['H_D_bits']:.3f} bits", flush=True)

    out = dict(theta0=THETA0, param_names=PARAM_NAMES, dt=DT, tmax=TMAX, M=M,
               fisher_h2=res2, fisher_h5=res5, capacity=cap,
               verdict=dict(P1_cos=res2["cos_kobs"] >= 0.99, P1_cond=res2["cond"] >= 1e6,
                            P2_var_ratio=res2["var_ratio"] >= 4))
    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    make_figure(res2, res5, cap)
    print("已写出：", OUT_JSON.name, OUT_PNG.name, flush=True)


if __name__ == "__main__":
    main()
