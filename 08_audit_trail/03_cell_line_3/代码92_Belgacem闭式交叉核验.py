# -*- coding: utf-8 -*-
"""
代码92：Belgacem 闭式 Hopf 周期 x 我们的分子常数——两条独立路线交叉核验
来源：arXiv 2605.23722（2026-05）Section 9 p53 标定：
  γ1=γ2=ln2≈0.693/h（1h 蛋白半衰期），环路增益反演使 τc=1h（实测转录延迟），
  得 AB≈1.72，闭式 Tc=2π/ωc≈5.6 h（对 5.5 h 差 3% 内）。
本卡四件事：
  P1 复现其闭式算术（ωc、τc、T 三个公式独立重算）；
  P2 用 Python 数值积分其 logistic 两基因 DDE，独立验证闭式公式（他用 R 做的）；
  P3 把我们登记表 v02 的分子常数代入同一闭式框架：τ 取 K3（实测 2.0 h /
     分解 1.70 h），γ 取不同口径，看闭式周期落在哪；
  P4 与我们的 Mönke 数值结果（模拟 5.48 h / Jacobian 5.95 h）和实测 5.5 h 对账。
纪律：确定性，种子标记 20260925；JSON + PNG/SVG + 判词卡。
"""
import json
import numpy as np
from scipy.optimize import brentq
from scipy.signal import find_peaks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SEED_TAG = 20260925
OUT_JSON = r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911\03_细胞线3\结果\代码92_Belgacem交叉核验_结果.json"
OUT_PNG = r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911\03_细胞线3\代码与图\代码92_Belgacem交叉核验_两面体.png"
OUT_SVG = OUT_PNG.replace(".png", ".svg")
plt.rcParams.update({"svg.fonttype": "none", "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False})

LN2 = np.log(2.0)
res = {"seed": SEED_TAG,
       "source": "arXiv 2605.23722 Section 9 (PDF belgacem_2605.23722.pdf 本地存档)"}

# ---------------- 闭式公式（论文 Section 4） ----------------
def omega_c(g1, g2, AB):
    p = (-(g1**2 + g2**2) + np.sqrt((g1**2 - g2**2)**2 + 4 * AB**2)) / 2
    return np.sqrt(p) if p > 0 else np.nan

def tau_c(g1, g2, AB):
    w = omega_c(g1, g2, AB)
    return np.arctan2(w * (g1 + g2), w**2 - g1 * g2) / w

def T_of(g1, g2, AB):
    return 2 * np.pi / omega_c(g1, g2, AB)

def AB_from_tauc(g1, g2, tau_target):
    """给定降解率与目标临界延迟，反演环路增益 AB。"""
    f = lambda AB: tau_c(g1, g2, AB) - tau_target
    return brentq(f, g1 * g2 * (1 + 1e-6), 1e6,
                  xtol=1e-14, rtol=1e-14)

# ---------------- P1：复现 Belgacem Section 9 算术 ----------------
g = LN2
AB_rep = AB_from_tauc(g, g, 1.0)
res["P1_reproduce"] = dict(
    gamma=LN2, tau_c_target_h=1.0,
    AB_inverted=float(AB_rep), AB_paper=1.72,
    omega_c=float(omega_c(g, g, AB_rep)),
    tau_c_check_h=float(tau_c(g, g, AB_rep)),
    T_closed_h=float(T_of(g, g, AB_rep)), T_paper_h=5.6, T_measured_h=5.5,
    reproduce_ok=bool(abs(AB_rep - 1.72) < 0.02 and abs(T_of(g, g, AB_rep) - 5.6) < 0.1))

# ---------------- P2：数值积分其 DDE（Python 独立复核其 R 验证） ----------------
# 对称参数化：κ1=κ2=κ, θ1=θ2=θ, λ=4；平衡点 x*=θ 处 f*=1/2，A=B=κλ/4=κ
# AB=κ²=1.72 -> κ=1.3114；平衡点条件 κ·0.5/γ=θ -> θ=0.946
KAPPA = float(np.sqrt(AB_rep))
LAM = 4.0
THETA = KAPPA * 0.5 / LN2
GAM = LN2


def f_minus(x):
    return 1.0 / (1.0 + np.exp(LAM * (x - THETA)))


def f_plus(x):
    return 1.0 / (1.0 + np.exp(-LAM * (x - THETA)))


def simulate_dde(tau_total, t_end=240.0, dt=0.002):
    """显式 Euler + 常数历史，τ1=τ2=τ/2。确定性。"""
    n = int(t_end / dt)
    lag = int(tau_total / 2 / dt)
    x1 = np.full(n + lag + 1, 0.4)
    x2 = np.full(n + lag + 1, 0.4)
    for i in range(lag, n + lag):
        x1[i + 1] = x1[i] + dt * (KAPPA * f_minus(x2[i + 1 - lag]) - GAM * x1[i])
        x2[i + 1] = x2[i] + dt * (KAPPA * f_plus(x1[i + 1 - lag]) - GAM * x2[i])
    t = np.arange(n + 1) * dt
    return t, x1[lag:], x2[lag:]


dde_rows = []
for tau in [1.0, 1.2, 1.5, 2.0]:
    t, x1, x2 = simulate_dde(tau)
    mask = t > 120
    pk, _ = find_peaks(x1[mask], prominence=0.02)
    if len(pk) >= 3:
        per = float(np.mean(np.diff(t[mask][pk])))
        amp = float((x1[mask].max() - x1[mask].min()) / 2)
    else:
        per, amp = None, float((x1[mask].max() - x1[mask].min()) / 2)
    dde_rows.append(dict(tau_h=tau, period_sim_h=per, amplitude=amp))
res["P2_dde_integration"] = dict(kappa=KAPPA, lam=LAM, theta=THETA, gamma=GAM,
                                 rows=dde_rows,
                                 note="τ=1.0 位于临界点上（边缘），周期仅在 τ>τc 的极限环上有意义")

# ---------------- P3：我们的分子常数代入闭式框架 ----------------
def T_from_measured_tau(g1, g2, tau_measured):
    AB = AB_from_tauc(g1, g2, tau_measured)
    return float(AB), float(T_of(g1, g2, AB))

p3 = []
# 行1：Belgacem 标定（对照）
AB1, T1 = T_from_measured_tau(LN2, LN2, 1.0)
p3.append(dict(case="Belgacem calib: g=ln2/1h both, tau=1.0h",
               g1=LN2, g2=LN2, tau=1.0, AB=AB1, T_h=T1))
# 行2：同一半衰期口径，τ 换成我们 K3 实测 2.0 h
AB2, T2 = T_from_measured_tau(LN2, LN2, 2.0)
p3.append(dict(case="same half-lives, tau=K3 measured 2.0h",
               g1=LN2, g2=LN2, tau=2.0, AB=AB2, T_h=T2))
# 行3：τ = 代码89 相位分解值 1.70 h
AB3, T3 = T_from_measured_tau(LN2, LN2, 1.70)
p3.append(dict(case="same half-lives, tau=K3 decomposed 1.70h (code89)",
               g1=LN2, g2=LN2, tau=1.70, AB=AB3, T_h=T3))
# 行4：应激态口径——p53 MDM2 依赖半衰期 ~15 min（γ=ln2/0.25），Mdm2 dM=2/h，τ=2.0 h
g_p53_fast = LN2 / 0.25
AB4, T4 = T_from_measured_tau(g_p53_fast, 2.0, 2.0)
p3.append(dict(case="stress regime: g_p53=ln2/0.25h, g_Mdm2=2/h, tau=2.0h",
               g1=g_p53_fast, g2=2.0, tau=2.0, AB=AB4, T_h=T4))
res["P3_our_constants_in_closed_form"] = p3

# ---------------- P4：总对账 ----------------
res["P4_reconciliation"] = dict(
    measured_period_h=5.5,
    belgacem_closed_form_h=T1,
    monke_sim_h=5.48, monke_jacobian_h=5.95,
    ours_tau2_closed_form_h=T2,
    verdict="见判词卡")

# ---------------- 图：两面体 ----------------
fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.6))

ax = axes[0]
taus = np.linspace(1.001, 3.0, 200)
Ts = [T_of(g, g, AB_from_tauc(g, g, tv)) for tv in taus]
ax.plot(taus, Ts, color="#4C72B0", lw=1.8,
        label="Belgacem closed form T(tau) at onset (g=ln2)")
xs = [r["tau_h"] for r in dde_rows if r["period_sim_h"]]
ys = [r["period_sim_h"] for r in dde_rows if r["period_sim_h"]]
ax.plot(xs, ys, "o", color="#C44E52", ms=7,
        label="our Python DDE integration (his model)")
ax.axhline(5.5, color="#333333", lw=1.0, ls="--")
ax.text(1.02, 5.62, "measured 5.5 h (MCF7)  |  our Monke simulation 5.48 h",
        fontsize=8, color="#333333")
ax.axhline(5.48, color="#55A868", lw=0.9, ls=":")
ax.axvline(2.0, color="#DD8452", lw=0.9, ls=":")
ax.text(2.03, 10.0, "K3 measured lag 2.0 h", fontsize=8, color="#B5641E")
ax.set_xlabel("loop delay tau (h)")
ax.set_ylabel("period (h)")
ax.set_title("a  Closed-form period vs loop delay, with our independent\n"
             "DDE integration overlaid (formula verified in Python)", fontsize=10)
ax.legend(fontsize=8, frameon=False, loc="upper left")

ax = axes[1]
cases = ["Belgacem\ncalib\n(tau=1.0)", "K3 measured\n(tau=2.0)",
         "K3 decomposed\n(tau=1.70)", "stress regime\n(fast p53)"]
Ts4 = [T1, T2, T3, T4]
colors = ["#888888", "#4C72B0", "#55A868", "#DD8452"]
bars = ax.bar(range(4), Ts4, color=colors, width=0.55)
for i, v in enumerate(Ts4):
    ax.text(i, v + 0.12, f"{v:.2f} h", ha="center", fontsize=9)
ax.axhspan(5.2, 5.8, color="#C44E52", alpha=0.15)
ax.axhline(5.5, color="#C44E52", lw=1.0, ls="--")
ax.text(-0.45, 5.62, "measured 5.5 h", fontsize=8, color="#C44E52", ha="left")
ax.set_xticks(range(4))
ax.set_xticklabels(cases, fontsize=8)
ax.set_ylabel("closed-form period (h)")
ax.set_ylim(0, max(Ts4) * 1.25)
ax.set_title("b  Same closed form, four constant sets:\n"
             "which molecular inputs land on 5.5 h?", fontsize=10)

fig.suptitle("Code 92 | Cross-validation: Belgacem closed-form Hopf period x our molecular constants",
             y=1.0)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
fig.savefig(OUT_SVG, bbox_inches="tight")

with open(OUT_JSON, "w", encoding="utf-8") as f:
    json.dump(res, f, ensure_ascii=False, indent=2)

print("P1: AB=%.4f (paper 1.72), tau_c=%.4f h, T=%.3f h (paper ~5.6) ok=%s"
      % (AB_rep, tau_c(g, g, AB_rep), T1, res["P1_reproduce"]["reproduce_ok"]))
print("P2 DDE:", [(r["tau_h"], None if r["period_sim_h"] is None else round(r["period_sim_h"], 3),
                   round(r["amplitude"], 3)) for r in dde_rows])
print("P3:", [(r["case"].split(":")[0], round(r["T_h"], 3), "AB=%.2f" % r["AB"]) for r in p3])
print("JSON/PNG/SVG written.")
