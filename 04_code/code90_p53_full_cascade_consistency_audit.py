# -*- coding: utf-8 -*-
"""
代码90：p53 全链路一致性审计
命题：剂量 -> DSB -> γH2AX -> ATM -> p53 脉冲引擎 -> 计数/首脉冲 -> 下游解码，
      每一环的常数必须与前后环对得上账；对不上的显式判"张力"或"缺口"。
方法：
  A 静态链路审计：14 条链路检查项，取自代码81-89 判词卡与常数登记 v01/v02 的归档值；
  B 动态 onset 分解：Mönke 2025 勘误版从基态定点阶跃到 DSB=100（log 型），
    数值提取 ATM* 半升时间、p53 阈值穿越、首峰时刻，分解 t1 的链上各段；
  C 矛盾登记汇总。
纪律：确定性，种子标记 20260925；JSON + PNG/SVG + 判词卡。
"""
import json
import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

SEED_TAG = 20260925
OUT_JSON = r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911\03_细胞线3\结果\代码90_全链路审计_结果.json"
OUT_PNG = r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911\03_细胞线3\代码与图\代码90_全链路审计_两面体.png"
OUT_SVG = OUT_PNG.replace(".png", ".svg")

plt.rcParams.update({"svg.fonttype": "none", "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False})

# ---------------- Mönke 2025 勘误版（同代码89） ----------------
PAR = dict(A=30.5, P=22.0, C=1.4, g=2.5, dAM=20.0,
           Tm=1.2, TM=4.0, Tw=1.2, TW=1.0,
           dA=0.16, dP=0.1, dm=1.0, dM=2.0, dw=1.3, dW=2.3,
           kA=0.5, kWA=0.14, kMP=0.15, kPm=1.0, kPw=1.0,
           R=2.0, Smax=0.2, gam=9.0)


def rhs_vec(t, x, S):
    p = PAR
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


res = {"seed": SEED_TAG}

# ================= A. 静态链路审计（14 环） =================
# 状态：closed / tension / gap
links = [
    dict(id="L0", en="dose->DSB", name="dose -> DSB count", constant="35 DSB/Gy",
         sources="代码81 锚定；代码82 双录4；代码87/88 沿用",
         status="closed",
         note="跨代码一致，但本身是假设口径（未实测换算），假设旗标保留"),
    dict(id="L1", en="DSB->gH2AX", name="DSB -> γH2AX", constant="半峰1-3 min，峰 30 min（K5）",
         sources="常数登记 v01 K5（多来源一致）",
         status="closed", note="分钟级，远小于下游所有时标"),
    dict(id="L2", en="gH2AX->ATM", name="γH2AX -> ATM 激活", constant="缺剂量分级单细胞时序",
         sources="v01 G1 registered gap",
         status="gap", note="ATM-SPARK/ATOMIC 探针存在但无公开剂量分级数据"),
    dict(id="L3", en="ATM->p53 cat.", name="ATM -> p53（催化）", constant="kcat/Km 未实测",
         sources="v02 G5（用户点名催化速率，恰在缺口）",
         status="gap", note="模型用唯象保护因子 R=2 代替"),
    dict(id="L4", en="p53 engine", name="p53 脉冲引擎（K1/K3/K4）", constant="5.5 h / 2.0 h / 1.25 h",
         sources="代码89：三常数全部分解闭合",
         status="closed", note="周期 5.95 vs 5.48（8.5% 尾巴登记）"),
    dict(id="L5", en="counting N(D)", name="剂量 -> 脉冲计数 N(D)", constant="计数律 D_c：归档 0.261 vs 代码82 0.99 Gy",
         sources="代码82 验证靶标1：3.8 倍差",
         status="tension", note="方向对（饱和形状 R2=0.995），绝对刻度差源于 L0 换算假设"),
    dict(id="L6", en="first pulse t1", name="首脉冲 t1", constant="分布（首峰 2.5-3.0 h），剂量无关",
         sources="D1；Lahav 2004；代码87 平坦 2.5-3.0 h 全剂量",
         status="closed", note="分布本体是感知层签名，不是缺陷"),
    dict(id="L7", en="dispersion 5.8", name="弥散比 t1/IPI", constant="归档 5.8 = (240/100)^2，5 Gy 口径",
         sources="代码86 口径错误已撤回；代码87/88 得 6.65-7.32",
         sources_detail="更正登记 2026-09-25",
         status="closed", note="同量级首合账；86 撤回已登记"),
    dict(id="L8", en="clock dose-indep.", name="时钟剂量无关架构", constant="剂量只走门控+终止",
         sources="代码87/88 闭环",
         status="closed", note="0.3 Gy 低剂量持续振荡见 L9"),
    dict(id="L9", en="low-dose osc.", name="低剂量振荡分数", constant="文献 0.3 Gy 约 35% 振荡",
         sources="代码87/88 判词卡",
         status="tension", note="损伤驱动终止造不出低剂量振荡；文献内部张力，已登记"),
    dict(id="L10", en="p21 tracking", name="下游 p21 转录追踪", constant="延迟精确值待提取",
         sources="v01 G3（Hafner 2020）",
         status="gap", note="定性同期振荡已证，定量延迟待从原文图提取"),
    dict(id="L11", en="info budget", name="信息账 I(D;fate)", constant="0.55-0.63 bits（修正后）",
         sources="代码81 的 1.325 因依赖剂量相关时序被 C1 证伪",
         status="closed", note="撤回与修正均已登记；论文从未引用 1.325"),
    dict(id="L12", en="Wip1 RNAi", name="Wip1 RNAi 表型", constant="幅度方向：模型升 vs Batchelor 降 45%",
         sources="代码82 验证靶标5",
         status="tension", note="失控与变宽预言复现；幅度方向冲突登记为模型-实验不符项"),
    dict(id="L13", en="width K2", name="脉冲宽度 K2 = 3.5 h", constant="分子层未分解",
         sources="v02 登记表 D 节",
         status="gap", note="由 ATM* 平台期决定，非纯延迟和，待分解"),
]
res["A_static_links"] = links
n_closed = sum(1 for l in links if l["status"] == "closed")
n_tension = sum(1 for l in links if l["status"] == "tension")
n_gap = sum(1 for l in links if l["status"] == "gap")
res["A_summary"] = dict(total=len(links), closed=n_closed, tension=n_tension, gap=n_gap)

# ================= B. 动态 onset 分解（t1 链上各段） =================
DSB_STEP = 100.0
S_basal = PAR["Smax"] * np.log(2.0 / PAR["gam"] + 1)   # 背景约 2 DSB
S_step = PAR["Smax"] * np.log(DSB_STEP / PAR["gam"] + 1)

x_basal = fsolve(lambda x: rhs_vec(0, x, S_basal),
                 [0.05, 0.3, 0.3, 0.3, 0.3, 0.1])
sol = solve_ivp(rhs_vec, (0, 12), x_basal, args=(S_step,),
                dense_output=True, rtol=1e-9, atol=1e-12)
t = np.linspace(0, 12, 12001)
x = sol.sol(t)
ATM, P53 = x[0], x[1]

# γH2AX 段：外部常数，峰值 0.5 h（K5）
t_gamma = 0.5
# ATM* 半升（相对其首峰最大值）
atm_peak = ATM.max()
atm_half = 0.5 * (atm_peak + ATM[0])
i_atm = np.argmax(ATM >= atm_half) if np.any(ATM >= atm_half) else len(t) - 1
t_atm = t[i_atm]
# p53 阈值穿越（检测阈 0.8 AU，代码82 口径）
i_thr = np.argmax(P53 >= 0.8) if np.any(P53 >= 0.8) else len(t) - 1
t_p53_thr = t[i_thr]
# p53 首峰
i_peak = np.argmax(P53[:np.searchsorted(t, 10)])
t_peak = t[i_peak]

t1_pred_onset = t_gamma + t_p53_thr   # γH2AX 段 + ATM/p53 上升段（onset 口径）
res["B_onset_decomposition"] = dict(
    S_basal=float(S_basal), S_step=float(S_step),
    basal_fixed_point=[float(v) for v in x_basal],
    gammaH2AX_segment_h=t_gamma,
    ATM_half_rise_h=float(t_atm),
    p53_threshold_cross_h=float(t_p53_thr),
    p53_first_peak_h=float(t_peak),
    t1_pred_onset_h=float(t1_pred_onset),
    t1_measured_band_h=[2.5, 3.0],
    t1_definition="文献 t1 = 首峰时刻（D1：MCF7 2-3 h）；onset 是阈穿越，不等于 t1",
    finding="感知段（γH2AX 0.5 + ATM 上升 %.2f + p53 阈值 %.2f）合计仅 %.2f h；"
            "t1 的大头是 p53 积累爬坡到峰（%.2f h），不是感知延迟"
            % (t_atm - t_gamma, t_p53_thr - t_atm, t_p53_thr, t_peak - t_p53_thr),
    peak_closes=bool(t_peak <= 3.0 + 0.5),
)

# 链上时序排序检查：γH2AX(0.5) < ATM 半升 < p53 阈值 < 首峰 < 周期 5.5
ordering_ok = bool(t_gamma <= t_atm <= t_p53_thr <= t_peak <= 5.5)
res["B_ordering_ok"] = ordering_ok

# ================= C. 矛盾登记汇总 =================
res["C_tensions_register"] = [
    "L5: 计数律 D_c 3.8 倍差（换算假设）",
    "L9: 0.3 Gy 低剂量振荡（文献内部张力）",
    "L12: Wip1 RNAi 幅度方向（模型-实验不符）",
    "B: 模拟首峰 %.2f h vs 实测 t1 带 2.5-3.0 h（阶跃口径差）" % t_peak,
]

# ================= 图：两面体 =================
fig = plt.figure(figsize=(12.5, 7.6))
gs = fig.add_gridspec(2, 1, height_ratios=[1.15, 1.0], hspace=0.32)

# (a) 链路图
ax = fig.add_subplot(gs[0])
status_color = {"closed": "#55A868", "tension": "#DD8452", "gap": "#C44E52"}
n = len(links)
cols = 7
for k, lk in enumerate(links):
    row = k // cols
    col = k % cols
    if row == 1:
        col = cols - 1 - col  # 第二行反向，蛇形
    xc = 0.06 + col * (0.88 / (cols - 1))
    yc = 0.72 if row == 0 else 0.18
    ax.annotate("", xy=(0, 0))  # noop keep autoscale off
    ax.add_patch(plt.Rectangle((xc - 0.055, yc - 0.13), 0.11, 0.26,
                               facecolor=status_color[lk["status"]],
                               edgecolor="none", alpha=0.9,
                               transform=ax.transAxes))
    ax.text(xc, yc + 0.045, lk["id"], transform=ax.transAxes,
            ha="center", fontsize=8.5, fontweight="bold", color="white")
    ax.text(xc, yc - 0.045, lk["en"], transform=ax.transAxes,
            ha="center", fontsize=6.5, color="white", wrap=True)
    lk["_pos"] = (xc, yc, row)
# 箭头
order_idx = list(range(cols)) + list(range(cols, n))
for a, b in zip(order_idx[:-1], order_idx[1:]):
    xa, ya, ra = links[a]["_pos"]
    xb, yb, rb = links[b]["_pos"]
    if ra == rb:
        ar = FancyArrowPatch((xa + (0.055 if xb > xa else -0.055), ya),
                             (xb + (-0.055 if xb > xa else 0.055), yb),
                             transform=ax.transAxes, arrowstyle="-|>",
                             mutation_scale=10, color="#555555", lw=1.0)
    else:
        ar = FancyArrowPatch((xa, ya - 0.13), (xb, yb + 0.13),
                             transform=ax.transAxes, arrowstyle="-|>",
                             mutation_scale=10, color="#555555", lw=1.0,
                             connectionstyle="arc3,rad=0.0")
    ax.add_patch(ar)
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.axis("off")
handles = [plt.Rectangle((0, 0), 1, 1, color=status_color[s]) for s in
           ["closed", "tension", "gap"]]
ax.legend(handles, [f"closed ({n_closed})", f"registered tension ({n_tension})",
                    f"registered gap ({n_gap})"],
          loc="lower center", ncol=3, frameon=False, fontsize=8,
          bbox_to_anchor=(0.5, -0.14))
ax.set_title("a  Full-chain audit: dose -> DSB -> gH2AX -> ATM -> p53 engine -> count/t1 -> decoding\n"
             "every link carries its archived constant and a verdict", fontsize=10)

# (b) onset 分解
ax = fig.add_subplot(gs[1])
segs = [("gH2AX formation\n(K5: 0.5 h)", t_gamma, "#4C72B0"),
        ("ATM* rise to half-max\n(model: %.2f h)" % t_atm, max(t_atm - t_gamma, 0.0), "#55A868"),
        ("p53 accumulation to 0.8 AU\n(model: %.2f h)" % (t_p53_thr - t_atm),
         max(t_p53_thr - t_atm, 0.0), "#8172B3"),
        ("rise to first peak\n(model: %.2f h)" % (t_peak - t_p53_thr),
         max(t_peak - t_p53_thr, 0.0), "#C44E52")]
left = 0
for lb, w, c in segs:
    ax.barh([0], [w], left=left, color=c, height=0.42, label=lb)
    if w > 0.25:
        ax.text(left + w / 2, 0, f"{w:.2f} h", ha="center", va="center",
                fontsize=8, color="white")
    left += w
ax.axvspan(2.5, 3.0, color="#DD8452", alpha=0.25)
ax.text(2.75, 0.32, "measured first-peak band 2.5-3.0 h\n(Lahav 2004 / code 87)",
        ha="center", fontsize=8, color="#B5641E")
ax.text(t_peak, 0.32, f"model first peak {t_peak:.2f} h",
        ha="left", fontsize=8, color="#C44E52")
ax.set_xlim(0, 6); ax.set_ylim(-0.55, 0.62)
ax.set_yticks([])
ax.set_xlabel("time after damage (h)")
ax.set_title("b  Onset decomposition along the chain: each segment is either an archived\n"
             "constant or a model-extracted delay; ordering gH2AX < ATM* < p53 < peak < T=5.5 h",
             fontsize=10)
ax.legend(fontsize=6.5, loc="upper right", frameon=False)

fig.suptitle("Code 90 | Full-chain consistency audit of the p53 damage-signalling cascade",
             y=0.985)
fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
fig.savefig(OUT_SVG, bbox_inches="tight")

with open(OUT_JSON, "w", encoding="utf-8") as f:
    json.dump(res, f, ensure_ascii=False, indent=2)

print("links: %d closed / %d tension / %d gap (of %d)"
      % (n_closed, n_tension, n_gap, len(links)))
print("onset: gamma %.2f -> ATM half %.2f -> p53 thr %.2f -> peak %.2f h"
      % (t_gamma, t_atm, t_p53_thr, t_peak))
print("first peak %.2f h vs band 2.5-3.0 (+0.5 tol) -> closes=%s ; ordering=%s"
      % (t_peak, res["B_onset_decomposition"]["peak_closes"], ordering_ok))
print("JSON/PNG/SVG written.")
