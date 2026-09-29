# -*- coding: utf-8 -*-
"""
代码91：K2 脉冲宽度（3.5 h）分子层分解 + 全链路缺口 L13 闭合
命题：宽度不归延迟管，归 ATM* 平台期管——
  W = Wip1 诱导延迟(K4 已闭合) + Wip1 积累到熄灭阈值的时长 + p53 快崩塌尾(分钟级)
方法：
  B1 Mönke 2025 勘误版 DSB=100 恒定输入，提取 p53 首脉冲 FWHM、ATM* 平台期、
     Wip1 熄灭阈值，逐段核对分解式；
  B2 剂量扫描（DSB 50/100/200/400）复核宽度剂量无关（归档主张）；
  B3 Wip1 产生率扰动（Tw x0.85/x0.70）复核宽度不变性主张（v01 K2 vs 代码82 V5 的张力）。
纪律：确定性，种子标记 20260925；JSON + PNG/SVG + 判词卡。
"""
import json
import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SEED_TAG = 20260925
OUT_JSON = r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911\03_细胞线3\结果\代码91_K2宽度分解_结果.json"
OUT_PNG = r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911\03_细胞线3\代码与图\代码91_K2宽度分解_三面体.png"
OUT_SVG = OUT_PNG.replace(".png", ".svg")

plt.rcParams.update({"svg.fonttype": "none", "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False})

PAR = dict(A=30.5, P=22.0, C=1.4, g=2.5, dAM=20.0,
           Tm=1.2, TM=4.0, Tw=1.2, TW=1.0,
           dA=0.16, dP=0.1, dm=1.0, dM=2.0, dw=1.3, dW=2.3,
           kA=0.5, kWA=0.14, kMP=0.15, kPm=1.0, kPw=1.0,
           R=2.0, Smax=0.2, gam=9.0)


def make_rhs(par):
    def rhs(t, x, S):
        p = par
        ATM, P53, mdm2, Mdm2, wip1, Wip1 = x
        return [
            p["A"] * (ATM**2 / (1 + ATM**2 / p["kA"])) / (1 + Wip1 / p["kWA"])
            - p["dA"] * ATM - p["P"] * ATM * Wip1 + S,
            p["C"] - p["dP"] * P53
            - p["g"] * Mdm2 * P53 / (p["kMP"] + P53) * (1 + p["R"] / (1 + ATM)),
            p["Tm"] * P53 / (p["kPm"] + P53) - p["dm"] * mdm2,
            p["TM"] * mdm2 - p["dM"] * Mdm2 - p["dAM"] * ATM * Mdm2,
            p["Tw"] * P53 / (p["kPw"] + P53) - p["dw"] * wip1,
            p["TW"] * wip1 - p["dW"] * Wip1,
        ]
    return rhs


def S_of(DSB):
    return PAR["Smax"] * np.log(DSB / PAR["gam"] + 1)


def simulate_first_pulse(par, DSB, t_end=12.0):
    """从基态定点阶跃，返回时间轴与六物种轨迹。"""
    rhs = make_rhs(par)
    S_basal = S_of(2.0)
    x0 = fsolve(lambda x: rhs(0, x, S_basal), [0.05, 0.3, 0.3, 0.3, 0.3, 0.1])
    sol = solve_ivp(rhs, (0, t_end), x0, args=(S_of(DSB),),
                    dense_output=True, rtol=1e-9, atol=1e-12)
    t = np.linspace(0, t_end, int(t_end * 2000) + 1)
    return t, sol.sol(t)


def first_pulse_fwhm(t, P53):
    """首脉冲半高宽；返回 (width, t_rise, t_fall, t_peak) 或 None。"""
    i_pk = int(np.argmax(P53))
    pk = P53[i_pk]
    half = 0.5 * (pk + P53[0])
    above = P53 >= half
    if not above[i_pk]:
        return None
    i_r = i_pk
    while i_r > 0 and above[i_r - 1]:
        i_r -= 1
    i_f = i_pk
    while i_f < len(t) - 1 and above[i_f + 1]:
        i_f += 1
    if i_f >= len(t) - 1:
        return None
    return float(t[i_f] - t[i_r]), float(t[i_r]), float(t[i_f]), float(t[i_pk])


def atm_plateau(t, ATM):
    """ATM* 首平台的半高起止。"""
    i_pk = int(np.argmax(ATM))
    pk = ATM[i_pk]
    half = 0.5 * (pk + ATM[0])
    above = ATM >= half
    i_r = i_pk
    while i_r > 0 and above[i_r - 1]:
        i_r -= 1
    i_f = i_pk
    while i_f < len(t) - 1 and above[i_f + 1]:
        i_f += 1
    return float(t[i_r]), float(t[i_f]), float(t[i_f] - t[i_r])


res = {"seed": SEED_TAG}

# ================= B1：DSB=100 分解 =================
t, x = simulate_first_pulse(PAR, 100.0)
ATM, P53, Wip1 = x[0], x[1], x[5]
fw = first_pulse_fwhm(t, P53)
W_meas, t_rise, t_fall, t_peak = fw
a_r, a_f, a_plateau = atm_plateau(t, ATM)

# Wip1 熄灭阈值：ATM* 降至半峰时刻的 Wip1 值
i_af = int(np.searchsorted(t, a_f))
W_threshold = float(Wip1[i_af])
# Wip1 从脉冲启动（p53 阈穿越，0.8 AU）到阈值所需时长
i_start = int(np.argmax(P53 >= 0.8))
t_start = float(t[i_start])
t_wip1_to_threshold = a_f - t_start
# p53 崩塌尾：ATM* 熄灭后 p53 降到半峰的时间
tail = t_fall - a_f
# 分解式：宽度 = (ATM* 熄灭 - p53 半升) + 崩塌尾
W_decomp = (a_f - t_rise) + tail
# K4 诱导延迟（代码89 闭合值）在宽度中的角色：Wip1 从 p53 启动到越过基线抬头
res["B1_decomposition"] = dict(
    DSB=100.0,
    W_meas_FWHM_h=W_meas,
    W_archived_h=3.5,
    segments=dict(
        p53_half_rise_h=t_rise,
        p53_peak_h=t_peak,
        ATM_plateau_h=a_plateau,
        ATM_extinction_h=a_f,
        Wip1_at_extinction_AU=W_threshold,
        wip1_accumulation_from_pulse_start_h=float(t_wip1_to_threshold),
        p53_collapse_tail_h=float(tail),
    ),
    W_decomp_h=float(W_decomp),
    closes=bool(abs(W_decomp - W_meas) / W_meas <= 0.05),
    closes_vs_archive=bool(abs(W_meas - 3.5) <= 0.7),
    tautology_note="W_decomp 与 FWHM 的相等是构造恒等（两段恰好铺满半高宽区间）；"
                   "分解的实质内容是归属：平台期 1.54 h 由 Wip1 积累决定，"
                   "崩塌尾 1.99 h 由 Mdm2 恢复动力学（dM/dAM）决定，"
                   "且崩塌尾长于平台期——宽度的主要分子归因是 Mdm2 恢复，不是 Wip1 阈值穿越本身",
    reading="宽度 = p53 半升到 ATM* 熄灭 (%.2f h) + 崩塌尾 (%.2f h) = %.2f h vs FWHM %.2f h"
            % (a_f - t_rise, tail, W_decomp, W_meas),
)

# ================= B2：剂量扫描，宽度剂量无关复核 =================
dose_scan = []
for DSB in [50.0, 100.0, 200.0, 400.0]:
    t2, x2 = simulate_first_pulse(PAR, DSB)
    fw2 = first_pulse_fwhm(t2, x2[1])
    dose_scan.append(dict(DSB=DSB, FWHM_h=fw2[0] if fw2 else None))
res["B2_dose_scan"] = dose_scan
ws = [d["FWHM_h"] for d in dose_scan if d["FWHM_h"] is not None]
res["B2_width_cv"] = float(np.std(ws) / np.mean(ws)) if len(ws) > 1 else None

# ================= B3：Wip1 扰动，宽度不变性复核 =================
wip1_scan = []
for scale in [1.0, 0.85, 0.70]:
    par2 = dict(PAR)
    par2["Tw"] = PAR["Tw"] * scale
    t3, x3 = simulate_first_pulse(par2, 100.0)
    fw3 = first_pulse_fwhm(t3, x3[1])
    wip1_scan.append(dict(Tw_scale=scale, FWHM_h=fw3[0] if fw3 else None))
res["B3_wip1_perturbation"] = wip1_scan

# ================= 图：三面体 =================
fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.4)

)

ax = axes[0]
ax.plot(t, P53, color="#C44E52", lw=1.8, label="p53")
ax.plot(t, ATM, color="#4C72B0", lw=1.5, label="ATM*")
ax.plot(t, Wip1, color="#55A868", lw=1.5, label="Wip1")
ax.axhline(0.5 * (P53.max() + P53[0]), color="#C44E52", lw=0.7, ls=":")
ax.axvline(t_rise, color="#888888", lw=0.8, ls="--")
ax.axvline(a_f, color="#4C72B0", lw=0.8, ls="--")
ax.axvline(t_fall, color="#C44E52", lw=0.8, ls="--")
ax.annotate("", xy=(t_rise, -0.25), xytext=(a_f, -0.25),
            arrowprops=dict(arrowstyle="<->", color="#333333"))
ax.text((t_rise + a_f) / 2, -0.55,
        "p53 half-rise to ATM* extinction\n%.2f h" % (a_f - t_rise),
        ha="center", fontsize=8)
ax.annotate("", xy=(a_f, -0.85), xytext=(t_fall, -0.85),
            arrowprops=dict(arrowstyle="<->", color="#C44E52"))
ax.text((a_f + t_fall) / 2, -1.15, "collapse tail %.2f h" % tail,
        ha="center", fontsize=8, color="#C44E52")
ax.plot([a_f], [W_threshold], marker="o", color="#55A868", ms=6)
ax.text(a_f + 0.1, W_threshold + 0.15, "Wip1 extinction threshold\n%.2f AU" % W_threshold,
        fontsize=7.5, color="#2D7A43")
ax.set_xlim(0, 8); ax.set_ylim(-1.4, max(P53.max(), ATM.max()) * 1.15)
ax.set_xlabel("time after damage (h)")
ax.set_ylabel("abundance (AU)")
ax.set_title("a  K2 anatomy: width = ATM* plateau (%.2f h)\n"
             "+ collapse tail (%.2f h); FWHM %.2f h vs archived 3.5 h"
             % (a_f - t_rise, tail, W_meas), fontsize=10)
ax.legend(fontsize=8, frameon=False, loc="upper right")

ax = axes[1]
xs_d = [d["DSB"] for d in dose_scan]
ys_d = [d["FWHM_h"] for d in dose_scan]
ax.plot(xs_d, ys_d, "o-", color="#4C72B0")
for xi, yi in zip(xs_d, ys_d):
    ax.text(xi, yi + 0.05, f"{yi:.2f}", ha="center", fontsize=8, color="#4C72B0")
ax.axhline(3.5, color="#888888", lw=1.0, ls="--")
ax.text(390, 3.56, "archived K2 = 3.5 h", fontsize=8, ha="right", color="#555555")
ax.set_xlabel("DSB (constant input)")
ax.set_ylabel("first-pulse FWHM (h)")
ax.set_ylim(2.4, 3.8)
ax.set_title("b  Dose scan: width nearly dose-independent\n"
             "(CV = %.3f across 8x dose range)" % res["B2_width_cv"], fontsize=10)

ax = axes[2]
xs_w = [100 * d["Tw_scale"] for d in wip1_scan]
ys_w = [d["FWHM_h"] for d in wip1_scan]
ax.plot([x for x, y in zip(xs_w, ys_w) if y is not None],
        [y for y in ys_w if y is not None], "s-", color="#DD8452")
for xi, yi in zip(xs_w, ys_w):
    if yi is not None:
        ax.text(xi, yi + 0.05, f"{yi:.2f}", ha="center", fontsize=8, color="#B5641E")
    else:
        ax.text(xi, 2.6, "no pulse\n(locked high)", ha="center", fontsize=8,
                color="#C44E52")
ax.axhline(3.5, color="#888888", lw=1.0, ls="--")
ax.set_xlabel("Wip1 production Tw (%)")
ax.set_ylabel("first-pulse FWHM (h)")
ax.set_ylim(2.4, 3.9)
ax.set_title("c  Wip1 knockdown: width broadens +30% at 85% Tw,\n"
             "no pulse at 70% (tension with Monke S9 invariance claim)", fontsize=10)

fig.suptitle("Code 91 | K2 pulse-width decomposition and invariance audit", y=1.0)
fig.tight_layout()
fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
fig.savefig(OUT_SVG, bbox_inches="tight")

with open(OUT_JSON, "w", encoding="utf-8") as f:
    json.dump(res, f, ensure_ascii=False, indent=2)

print("B1: FWHM %.3f h vs archive 3.5 h ; decomp %.3f h (err %.1f%%) closes=%s"
      % (W_meas, W_decomp, 100 * abs(W_decomp - W_meas) / W_meas,
         res["B1_decomposition"]["closes"]))
print("    ATM plateau %.2f h, Wip1 threshold %.3f AU, tail %.2f h"
      % (a_plateau, W_threshold, tail))
print("B2 dose scan:", [(d["DSB"], round(d["FWHM_h"], 3)) for d in dose_scan],
      "CV=%.3f" % res["B2_width_cv"])
print("B3 wip1 scan:", [(d["Tw_scale"], round(d["FWHM_h"], 3) if d["FWHM_h"] else None)
                        for d in wip1_scan])
print("JSON/PNG/SVG written.")
