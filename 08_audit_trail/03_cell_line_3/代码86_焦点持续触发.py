# -*- coding: utf-8 -*-
"""
代码86_焦点持续触发.py  (种子确定性)
================================================
任务：路线对账的最终章。代码85 已证伪"ATM 增益异质性补弥散比"路线；
本件换机制口径——触发信号不是活 DSB 数 A(t)（快相 t½≈2h，代码84 口径），
而是 γH2AX 焦点数 F(t)（实测单相衰减 t½≈11h：dSTORM 2Gy 50 foci@30min ->
12 foci@24h，2^(-24/11)=0.22 ≈ 12/50=0.24）。

生物学依据：ATM 在焦点处被激活，焦点在 DSB 重接合后仍长期存在
（K5 形成 30 min；Foray 学派核穿梭模型同向）。

由此产生一个零拟合恒等式预言：
  若触发是（近似）无记忆 Poisson 过程，则 CV(t1)≈1，
  弥散比 = SD(t1)/SD(IPI) ≈ mean(t1)/SD(IPI) ≈ 2.3h/0.4h ≈ 5.8，
  即文献"弥散比 5.8"不是独立常数，而是无记忆触发 + 精确振荡器的折叠签名。
  同时预言 弥散比(D) ∝ 1/D（低剂量更散），供实验证伪。

常数表：全部同代码84/85；唯一校准常数 LAM0'（/focus/h）只对齐
t1@10Gy ≈ 2.2h 文献锚（MCF7 首峰 2-3h），弥散比不参与校准。

用法：
  python 代码86_焦点持续触发.py calib    # LAM0' 三点校准，N=1000
  python 代码86_焦点持续触发.py final    # 定稿，N=4000，σ∈{0, 0.2}
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
OUT_JSON = ROOT / "结果" / "代码86_焦点持续触发_结果.json"
OUT_PNG = ROOT / "结果" / "代码86_焦点持续触发_图.png"
OUT_SVG = ROOT / "结果" / "代码86_焦点持续触发_图.svg"

DOSES = np.array([0.1, 0.25, 0.5, 1.0, 2.5, 5.0, 10.0])
K0_DSB = 35.0
R_FAST = 0.35
FRAC_SLOW = 0.10
T_HALF_SLOW = 20.0
R_SLOW = np.log(2) / T_HALF_SLOW
T_HALF_FOCI = 11.0            # h，γH2AX 焦点衰减（dSTORM 文献锚，见头注）
R_FOCI = np.log(2) / T_HALF_FOCI
T_PERIOD = 5.5
JIT_IPI = 0.085
W = 48.0
N_MAX = 12
Q21, TAU21 = 100.0, 10.0
QPU, TAUPU, PUMA_NTHR = 100.0, 4.0, 3
TH_PUMA, R0_FATE, TH_P21 = 800.0, 0.15, 1500.0
A_MIN = 10.0                  # 终止阈值仍用活 DSB 数（代码84 封卷口径）
T1_ANCHOR = 2.2               # h，t1@10Gy 校准锚（MCF7 首峰 2-3h 中带）
LIT_DISPERSION = 5.8
IPI_SD_REF = 0.4              # h，文献 IPI 精度量级


def mi_discrete(x, y):
    x = np.asarray(x); y = np.asarray(y)
    Kx = int(x.max()) + 1; Ky = int(y.max()) + 1
    N = len(x)
    cxy = np.bincount(x * Ky + y, minlength=Kx * Ky).reshape(Kx, Ky) / N
    cx = cxy.sum(1, keepdims=True); cy = cxy.sum(0, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        term = cxy * np.log2(cxy / (cx @ cy))
    return float(np.nansum(term))


def first_pulse_time(taus, lam0, rng):
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


def simulate_dose(D, n, rng, lam0_foci, sigma_ln):
    """触发走焦点 F(t)（单相 t½=11h），终止走活 DSB A(t)（双相，代码84 口径）。"""
    n_dsb = rng.poisson(K0_DSB * D, n)
    if sigma_ln > 0:
        g = np.exp(sigma_ln * rng.standard_normal(n) - 0.5 * sigma_ln ** 2)
    else:
        g = np.ones(n)
    lam_cell = lam0_foci * g
    taus_dsb, taus_foci = [], []
    for i in range(n):
        ni = n_dsb[i]
        if ni == 0:
            taus_dsb.append(np.empty(0)); taus_foci.append(np.empty(0)); continue
        slow = rng.random(ni) < FRAC_SLOW
        rates = np.where(slow, R_SLOW, R_FAST)
        taus_dsb.append(rng.exponential(1.0 / rates))
        taus_foci.append(rng.exponential(1.0 / R_FOCI, ni))
    t1 = np.full(n, np.nan)
    times = np.full((n, N_MAX), np.nan)
    A_at = lambda taus, t: int(np.sum(taus > t))
    for i in range(n):
        tp = first_pulse_time(taus_foci[i], lam_cell[i], rng)  # 触发：焦点
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
            if A_at(taus_dsb[i], tp) < A_MIN:                  # 终止：活 DSB
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
    return dict(n_dsb=n_dsb, t1=t1, times=times, n_pulse=n_pulse,
                responded=responded, m_p21=m_p21, m_puma=m_puma, fate=fate)


def run_ensemble(lam0_foci, sigma_ln, n_cells, seed):
    rng = np.random.default_rng(seed)
    parts = [simulate_dose(D, n_cells, rng, lam0_foci, sigma_ln) for D in DOSES]
    out = dict(lam0_foci=lam0_foci, sigma_ln=sigma_ln, N=n_cells)
    p10 = parts[-1]
    out["count_10Gy"] = float(p10["n_pulse"].mean())
    t1v = p10["t1"]; t1v = t1v[~np.isnan(t1v)]
    out["t1_10Gy_mean"] = float(t1v.mean())
    out["t1_10Gy_sd"] = float(t1v.std())
    out["cv_t1_10Gy"] = float(t1v.std() / t1v.mean())
    out["resp_frac"] = [float(p["responded"].mean()) for p in parts]
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
            ipi_sd = float(np.mean(ipi_sd_cell))
            disp[f"{D}Gy"] = dict(t1_mean=float(t1d.mean()), t1_sd=float(t1d.std()),
                                  cv_t1=float(t1d.std() / t1d.mean()), ipi_sd=ipi_sd,
                                  dispersion_ratio=float(t1d.std() / ipi_sd),
                                  identity_check=float(t1d.mean() / ipi_sd))
    out["dispersion"] = disp
    d_idx = np.concatenate([np.full(n_cells, i) for i in range(len(DOSES))])
    fate = np.concatenate([p["fate"] for p in parts])
    out["I_D_fate"] = mi_discrete(d_idx, fate)
    return out


def main():
    mode = sys.argv[1] if len(sys.argv) > 1 else "calib"
    log = lambda *a: print(*a, flush=True)

    if mode == "calib":
        log(f"[代码86] calib：LAM0' 三点扫描，N=1000，种子 {SEED}")
        grid = {}
        for lam in [0.0008, 0.0013, 0.0020]:
            r = run_ensemble(lam, 0.0, 1000, SEED + int(lam * 1e7))
            grid[f"{lam}"] = dict(t1_10Gy_mean=r["t1_10Gy_mean"],
                                  count_10Gy=r["count_10Gy"],
                                  dispersion_10Gy=r["dispersion"].get("10.0Gy", {}))
            log(f"  LAM0'={lam}: t1@10={r['t1_10Gy_mean']:.2f}h, "
                f"<N>10={r['count_10Gy']:.2f}, "
                f"弥散比@10={r['dispersion'].get('10.0Gy', {}).get('dispersion_ratio', float('nan')):.2f}")
        payload = {}
        if OUT_JSON.exists():
            payload = json.loads(OUT_JSON.read_text(encoding="utf-8"))
        payload["calib"] = grid
        payload.setdefault("meta", dict(script="代码86_焦点持续触发.py", seed=SEED,
                                        t_half_foci_h=T_HALF_FOCI,
                                        anchor="LAM0' 只对齐 t1@10Gy~2.2h"))
        OUT_JSON.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
        log("JSON 已写出:", OUT_JSON)
        return

    # final
    lam_star = float(sys.argv[2]) if len(sys.argv) > 2 else 0.0013
    N = 4000
    log(f"[代码86] final：LAM0'={lam_star}，N={N}，σ ∈ {{0, 0.2}}")
    finals = {}
    for s in [0.0, 0.2]:
        r = run_ensemble(lam_star, s, N, SEED + 555 + int(s * 1000))
        finals[f"{s}"] = r
        dr = {k: round(v["dispersion_ratio"], 2) for k, v in r["dispersion"].items()}
        cv = {k: round(v["cv_t1"], 2) for k, v in r["dispersion"].items()}
        log(f"  σ={s}: <N>10={r['count_10Gy']:.2f}, t1@10={r['t1_10Gy_mean']:.2f}±"
            f"{r['t1_10Gy_sd']:.2f}h (CV={r['cv_t1_10Gy']:.2f}), "
            f"I(D;fate)={r['I_D_fate']:.3f}")
        log(f"    弥散比={dr}")
        log(f"    CV(t1)={cv}")
    payload = {}
    if OUT_JSON.exists():
        payload = json.loads(OUT_JSON.read_text(encoding="utf-8"))
    payload["final"] = finals
    payload["meta"] = dict(script="代码86_焦点持续触发.py", seed=SEED, mode="final",
                           lam0_foci=lam_star, t_half_foci_h=T_HALF_FOCI,
                           anchors="LAM0' 只对齐 t1@10Gy~2.2h；t½(foci)=11h 为文献锚",
                           identity_prediction="CV(t1)~1 时 弥散比~mean(t1)/SD(IPI)")
    OUT_JSON.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    log("JSON 已更新:", OUT_JSON)
    make_figure(finals, lam_star)
    log("[代码86] 完成。")


def make_figure(finals, lam_star):
    fig, axes = plt.subplots(2, 2, figsize=(13.5, 9.2))
    r0 = finals["0.0"]
    # (a) 弥散比 vs 剂量（两条路线对比）
    ax = axes[0, 0]
    common = [k for k in r0["dispersion"].keys() if k in finals["0.2"]["dispersion"]]
    doses_d = [float(k.replace("Gy", "")) for k in common]
    for s, c, lab in [("0.0", "#33527a", "code 86 (foci trigger, g fixed)"),
                      ("0.2", "#0b6b3a", "code 86 + ATM gain (Fano-anchored)")]:
        r = finals[s]
        ys = [r["dispersion"][k]["dispersion_ratio"] for k in common]
        ax.plot(doses_d, ys, "o-", color=c, label=lab)
    ax.axhline(LIT_DISPERSION, color="#8c1d18", ls="--", lw=1.0)
    ax.text(0.03, 0.06, "literature 5.8", fontsize=8, color="#8c1d18",
            transform=ax.transAxes)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("dose (Gy)"); ax.set_ylabel("dispersion ratio")
    ax.legend(fontsize=8, loc="upper right")
    ax.set_title("(a) Dispersion ratio vs dose: zero-fit landing", fontsize=10.5)
    # (b) CV(t1) 恒等式验证
    ax = axes[0, 1]
    cvs = [r0["dispersion"][k]["cv_t1"] for k in common]
    idc = [r0["dispersion"][k]["identity_check"] for k in common]
    rat = [r0["dispersion"][k]["dispersion_ratio"] for k in common]
    ax.plot(doses_d, cvs, "o-", color="#33527a", label="CV($t_1$)")
    ax.axhline(1.0, color="#8c1d18", ls="--", lw=1.0)
    ax.text(0.03, 0.06, "memoryless trigger: CV=1", fontsize=8, color="#8c1d18",
            transform=ax.transAxes)
    ax.set_xscale("log")
    ax.set_xlabel("dose (Gy)"); ax.set_ylabel("CV($t_1$)")
    ax.set_ylim(0.8, 1.45)
    ax2 = ax.twinx()
    ax2.plot(doses_d, np.array(rat) / np.array(idc), "s--", color="#d9863d",
             label="ratio / identity")
    ax2.axhline(1.0, color="#d9863d", ls=":", lw=0.8)
    ax2.set_ylabel("ratio / [mean($t_1$)/SD(IPI)]", color="#d9863d")
    ax2.set_ylim(0.8, 1.45)
    h1, l1 = ax.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=8, loc="upper right")
    ax.set_title("(b) Identity test: dispersion = mean($t_1$)/SD(IPI) when CV=1", fontsize=10.5)
    # (c) 锚定量摘要
    ax = axes[1, 0]
    ax.axis("off")
    ax.text(0.5, 0.55, f"t1 at 10 Gy: {r0['t1_10Gy_mean']:.2f} ± {r0['t1_10Gy_sd']:.2f} h\n"
            f"CV = {r0['cv_t1_10Gy']:.2f} (exponential: 1.00)\n"
            f"count law <N> = {r0['count_10Gy']:.2f} (anchor ~7)\n"
            f"responding fraction 0.1->10 Gy:\n  "
            + " -> ".join(f"{v:.2f}" for v in r0["resp_frac"]),
            transform=ax.transAxes, ha="center", va="center", fontsize=10,
            bbox=dict(fc="#f4f6fa", ec="#33527a", lw=0.8))
    ax.set_title("(c) Anchors at LAM0' = " + f"{lam_star}", fontsize=10.5)
    # (d) 路线对账总表
    ax = axes[1, 1]
    ax.axis("off")
    lines = [
        "Route reconciliation ledger (dispersion ratio, 10 Gy)",
        "",
        "code 81 ideal map ......... n/a (assumed ceil mapping)",
        "code 82 homogeneous ....... 0.01  (under-dispersed)",
        "code 84 live-DSB hazard ... 8.6   (fat tail, CV~2)",
        "code 85 + ATM gain ........ 10.2  (falsified direction)",
        f"code 86 foci trigger ...... {r0['dispersion']['10.0Gy']['dispersion_ratio']:.1f}   (zero-fit)",
        f"code 86 + gain sigma 0.2 .. {finals['0.2']['dispersion']['10.0Gy']['dispersion_ratio']:.1f}",
        "literature ................ 5.8",
        "",
        f"identity: ratio ~ mean(t1)/SD(IPI) = "
        f"{r0['dispersion']['10.0Gy']['identity_check']:.1f}",
        f"I(D;fate): sigma0 = {finals['0.0']['I_D_fate']:.3f}, "
        f"sigma0.2 = {finals['0.2']['I_D_fate']:.3f} bits",
        "(code 81: 1.325; code 84: 0.817)",
    ]
    ax.text(0.02, 0.98, "\n".join(lines), transform=ax.transAxes, va="top",
            fontsize=9, family="monospace")
    ax.set_title("(d) Four-route ledger", fontsize=10.5)
    fig.suptitle("Code 86: persistent-foci trigger closes the dispersion account", fontsize=12)
    fig.subplots_adjust(left=0.07, right=0.97, top=0.90, bottom=0.08,
                        hspace=0.34, wspace=0.42)
    fig.savefig(OUT_PNG, dpi=200, bbox_inches="tight")
    fig.savefig(OUT_SVG, bbox_inches="tight")
    print("图已写出:", OUT_PNG, flush=True)


if __name__ == "__main__":
    main()
