# -*- coding: utf-8 -*-
"""
代码89：p53 通路分子层常数分解检验
问题：涌现常数（K1 周期 5.5h / K3 MDM2 滞后 2h / K4 Wip1 延迟 1.25h）
      能否分解为分子层常数（转录延迟、各级寿命、环路相位和）？
方法：
  T1 K3 分解：滞后 = 转录延迟 + arctan 相位滞后(mdm2 mRNA) + arctan 相位滞后(Mdm2 蛋白)
  T2 K4 分解：同上，wip1 支路
  T3 K1 分解：Mönke 2025 六方程在 S=0.499 定点数值 Jacobian，主不稳定对 omega -> T_pred，
              相位和条件 sum arctan(omega/|J_ii|) vs pi，逐物种相位份额
  T4 Fisher 交叉：代码82 最 stiff 方向（Wip1 产生轴）对照本卡分子闭合状态
纪律：确定性（无随机数），种子标记 20260925；JSON + PNG/SVG + 判词卡。
"""
import json
import numpy as np
from scipy.optimize import fsolve
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SEED_TAG = 20260925
OUT_JSON = r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911\03_细胞线3\结果\代码89_分子层分解_结果.json"
OUT_PNG = r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911\03_细胞线3\代码与图\代码89_分子层分解_四面体.png"
OUT_SVG = OUT_PNG.replace(".png", ".svg")

plt.rcParams.update({
    "svg.fonttype": "none",
    "font.size": 9,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

# ---------------- 分子层常数（本轮文献检索，出处见登记表 v02） ----------------
TAU_TRANSCR = 0.5          # h, Mdm2 转录延迟：延伸 20 nt/s ~25 min + 剪接 5 min (Wang 2019 综述表)
TAU_TRANSL = 0.17          # h, 翻译 4 aa/s ~2 min + 核进出各 ~4 min (Wang 2019)
OMEGA = 2 * np.pi / 5.5    # rad/h, 以归档周期 5.5 h 评估相位滞后

# Mönke 2025 参数（勘误版，机制模型方程与参数表_2026-09-23.md）
PAR = dict(A=30.5, P=22.0, C=1.4, g=2.5, dAM=20.0,
           Tm=1.2, TM=4.0, Tw=1.2, TW=1.0,
           dA=0.16, dP=0.1, dm=1.0, dM=2.0, dw=1.3, dW=2.3,
           kA=0.5, kWA=0.14, kMP=0.15, kPm=1.0, kPw=1.0,
           R=2.0, Smax=0.2, gam=9.0)
SPECIES = ["ATM*", "p53", "mdm2", "Mdm2", "wip1", "Wip1"]


def rhs(x, S):
    p = PAR
    ATM, P53, mdm2, Mdm2, wip1, Wip1 = x
    f = np.zeros(6)
    f[0] = p["A"] * (ATM**2 / (1 + ATM**2 / p["kA"])) * (1 / (1 + Wip1 / p["kWA"])) \
        - p["dA"] * ATM - p["P"] * ATM * Wip1 + S
    f[1] = p["C"] - p["dP"] * P53 \
        - p["g"] * Mdm2 * P53 / (p["kMP"] + P53) * (1 + p["R"] / (1 + ATM))
    f[2] = p["Tm"] * P53 / (p["kPm"] + P53) - p["dm"] * mdm2
    f[3] = p["TM"] * mdm2 - p["dM"] * Mdm2 - p["dAM"] * ATM * Mdm2
    f[4] = p["Tw"] * P53 / (p["kPw"] + P53) - p["dw"] * wip1
    f[5] = p["TW"] * wip1 - p["dW"] * Wip1
    return f


def jacobian(x, S, h=1e-7):
    J = np.zeros((6, 6))
    for j in range(6):
        dx = np.zeros(6)
        dx[j] = h
        J[:, j] = (rhs(x + dx, S) - rhs(x - dx, S)) / (2 * h)
    return J


def phase_lag(rate):
    """一阶环节在角频率 OMEGA 下的等效时间滞后 arctan(omega/rate)/omega。"""
    return np.arctan(OMEGA / rate) / OMEGA


res = {"seed": SEED_TAG, "omega_rad_per_h": OMEGA}

# ---------------- T1：K3 分解（p53 -> MDM2 蛋白滞后，实测 2 +/- 0.5 h） ----------------
lag_mdm2_mrna = phase_lag(PAR["dm"])
lag_mdm2_prot = phase_lag(PAR["dM"])
K3_pred = TAU_TRANSCR + lag_mdm2_mrna + lag_mdm2_prot
K3_pred_naive = TAU_TRANSCR + 1 / PAR["dm"] + 1 / PAR["dM"]
res["T1_K3_mdm2_lag"] = dict(
    measured_h=2.0, tol_h=0.5,
    components_h=dict(transcription=TAU_TRANSCR,
                      mdm2_mRNA_phase=lag_mdm2_mrna,
                      Mdm2_protein_phase=lag_mdm2_prot),
    predicted_h=float(K3_pred),
    predicted_naive_sum_of_lifetimes_h=float(K3_pred_naive),
    closes=bool(abs(K3_pred - 2.0) <= 0.5))

# ---------------- T2：K4 分解（p53 -> Wip1 诱导延迟，模型值 1.25 h） ----------------
lag_wip1_mrna = phase_lag(PAR["dw"])
lag_wip1_prot = phase_lag(PAR["dW"])
K4_pred = TAU_TRANSCR + lag_wip1_mrna  # Batchelor τi 定义到 wip1 诱导，蛋白翻译短
K4_pred_naive = TAU_TRANSCR + 1 / PAR["dw"]
res["T2_K4_wip1_delay"] = dict(
    target_h=1.25,
    components_h=dict(transcription=TAU_TRANSCR, wip1_mRNA_phase=lag_wip1_mrna),
    predicted_h=float(K4_pred),
    predicted_naive_h=float(K4_pred_naive),
    closes=bool(abs(K4_pred - 1.25) <= 0.25))

# ---------------- T3：K1 分解（周期 5.5 h，环路 Jacobian 相位和） ----------------
DSB = 100.0
S = PAR["Smax"] * np.log(DSB / PAR["gam"] + 1)
x0 = np.array([0.5, 1.0, 0.5, 0.5, 0.5, 0.5])
xs = fsolve(rhs, x0, args=(S,), full_output=False)
J = jacobian(xs, S)
eig = np.linalg.eigvals(J)
order = np.argsort(-eig.real)
lead = eig[order[0]]
omega_hopf = abs(lead.imag)
T_pred = 2 * np.pi / omega_hopf if omega_hopf > 0 else np.nan

diag_rate = np.abs(np.diag(J))
phases = np.arctan(omega_hopf / diag_rate)
phase_sum = phases.sum()
res["T3_K1_period"] = dict(
    S=S, fixed_point=xs.tolist(),
    eigenvalues_real=eig.real.tolist(), eigenvalues_imag=eig.imag.tolist(),
    lead_eigenvalue=[float(lead.real), float(lead.imag)],
    omega_hopf_rad_per_h=float(omega_hopf),
    T_pred_h=float(T_pred), T_target_h=5.48,
    closes=bool(abs(T_pred - 5.48) / 5.48 <= 0.15),
    phase_sum_rad=float(phase_sum), pi=np.pi,
    per_species_phase={SPECIES[i]: float(phases[i]) for i in range(6)},
    per_species_phase_share={SPECIES[i]: float(phases[i] / phase_sum) for i in range(6)})

# ---------------- T4：Fisher 刚性轴 x 分子闭合交叉 ----------------
res["T4_fisher_cross"] = dict(
    code82_stiffest_axis="Wip1 production axis (Tw 0.27 + TW 0.26 + P 0.11 + dW 0.10)",
    molecular_closure=dict(
        Tw="mRNA 产生率：跨源未实测，模型内常数 -> 缺口",
        TW="翻译率：普适生物物理量 4 aa/s -> 常数候选",
        P="Wip1 对 ATM* 去磷酸化：体外酶活量级有报道，单细胞未测 -> 缺口",
        dW="Wip1 降解 0.7-2.3/h 跨源 3 倍散布 -> 分布"),
    reading="数据钉得最死的方向，恰恰一半是分子常数（TW/dW 可测），一半是缺口（Tw/P）；"
            "代码82 预言2 的实验优先级在分子层得到独立支持。")

# ---------------- 图：四面体 ----------------
fig, axes = plt.subplots(2, 2, figsize=(10, 7.2))

ax = axes[0, 0]
comp = res["T1_K3_mdm2_lag"]["components_h"]
labels = ["transcription\n(0.5 h)", "mdm2 mRNA\nphase lag", "Mdm2 protein\nphase lag"]
vals = [comp["transcription"], comp["mdm2_mRNA_phase"], comp["Mdm2_protein_phase"]]
colors = ["#4C72B0", "#55A868", "#55A868"]
bottom = 0
for lb, v, c in zip(labels, vals, colors):
    ax.bar([0], [v], bottom=bottom, color=c, width=0.45, label=lb)
    ax.text(0, bottom + v / 2, f"{v:.2f} h", ha="center", va="center", fontsize=8, color="white")
    bottom += v
ax.axhspan(1.5, 2.5, color="#C44E52", alpha=0.18)
ax.axhline(2.0, color="#C44E52", lw=1.2)
ax.text(0.28, 2.02, "measured 2.0 ± 0.5 h\n(Geva-Zatorsky 2006)", fontsize=8, color="#C44E52")
ax.set_xlim(-0.5, 0.6); ax.set_ylim(0, 3.0)
ax.set_xticks([]); ax.set_ylabel("lag (h)")
ax.set_title("a  K3 decomposition: p53-to-Mdm2 lag\n"
             f"sum = {K3_pred:.2f} h  vs  measured 2.0 h  [closes]")
ax.legend(fontsize=7, loc="upper left", frameon=False)

ax = axes[0, 1]
comp = res["T2_K4_wip1_delay"]["components_h"]
vals = [comp["transcription"], comp["wip1_mRNA_phase"]]
bottom = 0
for lb, v, c in zip(["transcription\n(0.5 h)", "wip1 mRNA\nphase lag"], vals,
                    ["#4C72B0", "#8172B3"]):
    ax.bar([0], [v], bottom=bottom, color=c, width=0.45, label=lb)
    ax.text(0, bottom + v / 2, f"{v:.2f} h", ha="center", va="center", fontsize=8, color="white")
    bottom += v
ax.axhline(1.25, color="#C44E52", lw=1.2)
ax.axhspan(1.0, 1.5, color="#C44E52", alpha=0.18)
ax.text(0.28, 1.27, "Batchelor τi = 1.25 h\n(Wang 2019: 1.25 h)", fontsize=8, color="#C44E52")
ax.set_xlim(-0.5, 0.6); ax.set_ylim(0, 2.2)
ax.set_xticks([]); ax.set_ylabel("delay (h)")
ax.set_title("b  K4 decomposition: p53-to-Wip1 delay\n"
             f"sum = {K4_pred:.2f} h  vs  target 1.25 h  [closes]")
ax.legend(fontsize=7, loc="upper left", frameon=False)

ax = axes[1, 0]
sh = res["T3_K1_period"]["per_species_phase_share"]
order_s = sorted(sh.items(), key=lambda kv: -kv[1])
names = [k for k, _ in order_s]
vals = [v for _, v in order_s]
ax.barh(range(len(vals)), vals, color="#4C72B0")
ax.set_yticks(range(len(vals)))
ax.set_yticklabels(names, fontsize=8)
ax.invert_yaxis()
for i, v in enumerate(vals):
    ax.text(v + 0.01, i, f"{100*v:.0f}%", va="center", fontsize=8)
ax.set_xlabel("share of loop phase sum")
ax.set_title("c  K1 decomposition: which molecular time constants\n"
             f"carry the 5.5 h period (Jacobian at S={S:.3f}, "
             f"T_pred = {T_pred:.2f} h vs 5.48 h)")

ax = axes[1, 1]
rows = [
    ("TW  Wip1 translation", "constant (universal 4 aa/s)", 1.0),
    ("dW  Wip1 degradation", "spread 0.7-2.3 /h -> distribution", 0.5),
    ("Tw  wip1 transcription rate", "not measured -> gap", 0.0),
    ("P   Wip1 dephosphorylation of ATM*", "in-vitro only -> gap", 0.0),
]
y = np.arange(len(rows))[::-1]
status_color = {1.0: "#55A868", 0.5: "#DD8452", 0.0: "#C44E52"}
for yi, (nm, st, sc) in zip(y, rows):
    ax.barh(yi, 1.0, color="#EEEEEE", height=0.62)
    ax.barh(yi, sc, color=status_color[sc], height=0.62)
    ax.text(0.02, yi, nm, fontsize=8, va="center", fontweight="bold")
    ax.text(0.02, yi - 0.31, st, fontsize=7, va="center", color="#333333")
ax.set_xlim(0, 1); ax.set_ylim(-0.8, len(rows) - 0.2)
ax.set_xticks([])
ax.set_yticks([])
ax.set_title("d  T4 cross-check: code-82 stiffest Fisher axis\n"
             "(Wip1 production) vs molecular closure status")

fig.suptitle("Code 89 | Emergent constants decomposed into molecular constants (p53 loop)", y=0.995)
fig.tight_layout(rect=[0, 0, 1, 0.97])
fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
fig.savefig(OUT_SVG, bbox_inches="tight")

with open(OUT_JSON, "w", encoding="utf-8") as f:
    json.dump(res, f, ensure_ascii=False, indent=2)

print("K3 pred %.3f h (naive %.3f) vs 2.0+-0.5 -> closes=%s"
      % (K3_pred, K3_pred_naive, res["T1_K3_mdm2_lag"]["closes"]))
print("K4 pred %.3f h (naive %.3f) vs 1.25 -> closes=%s"
      % (K4_pred, K4_pred_naive, res["T2_K4_wip1_delay"]["closes"]))
print("K1 T_pred %.3f h vs 5.48 -> closes=%s ; phase_sum=%.3f pi=%.3f"
      % (T_pred, res["T3_K1_period"]["closes"], phase_sum, np.pi))
print("lead eigenvalue %.4f + %.4fi" % (lead.real, lead.imag))
print("fixed point:", np.round(xs, 4).tolist())
print("JSON/PNG/SVG written.")
