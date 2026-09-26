# -*- coding: utf-8 -*-
"""
代码50 · p53 级联的框架自检仿真  v1.2.0 (2026-08-15) 双臂版
============================================================
目的：级联框架（v4.2）三腿 + 信息损失恒等式，在 p53 通路最小随机模型上的数值自检。

谱系：v1.0 单臂（数字臂）→ v1.2 双臂（外部审计修订）。
外部审计指出：v1.0 幅度从不携带剂量信息（恒定幅度），"幅度死于简并"演示力度弱。
修订为双臂：
  臂A 数字臂：幅度恒定，剂量→脉冲数（Lahav 2004 定性事实简化，p53 忠实臂）；
  臂B 模拟臂：幅度随剂量线性增加 a(D)=1+0.25(D-1)，剂量双编码（通用模拟系统臂）。
两臂回答不同问题：A=忠实通路的行为；B=腿1 的强检验（幅度有真信息可死）。

审计的双录（如实记录，双方不动数字）：
  审计预言"fold 在臂B 保持（比值消增益）"——实测否定：fold I=0.401→0.178 衰减。
  仿射 b≠0 下 fold=(ka+b)/(km+b)，k 消不掉 → fold 泄漏，
  与框架接口代数表"fold 仅 b=0 幸存"一行一致。

其余假设：
  A2: 仿射观测/接口：y = k·p53 + b；k ~ LogNormal(0, σ_k) 逐细胞不可辨识；b=0.3。
  A3: 解码器：幅度 / fold / 计数（突出度/MAD 阈值，对仿射不变）。
      差分比解码器在脉冲串模型中定义不稳（脉冲位置与定点不对齐），留待连续响应模型。

纪律：过往项目资产（线1 p53 审计，已封存）仅为背景，不进入本模型。
输出：控制台打印 + 图 结果/仿真_p53级联_三腿核验.png。种子 20260815。纯 CPU。
"""
import numpy as np
from scipy.signal import find_peaks

SEED = 20260815
rng = np.random.default_rng(SEED)

D2 = np.array([1., 2., 3., 4., 5.])     # 损伤水平（映射脉冲数）
AMP_MU, AMP_SD = 1.0, 0.15
AMP_DOSE_COEF = 0.25                    # 臂B：幅度-剂量系数
B_BASAL = 0.3
T_WIN, DT = 24.0, 0.1
SIG_LIST = [0.0, 0.2, 0.4, 0.6]


# ---------- 模型（双臂） ----------
def gen_trace3(D, k, b, arm):
    N = int(round(D)) + (1 if rng.random() < 0.2 else 0)
    ts = np.arange(0, T_WIN, DT)
    y = np.zeros_like(ts)
    for j in range(N):
        pt = 2.0 + j * 5.5 + rng.normal(0, 0.3)
        if pt < T_WIN - 1:
            mu = AMP_MU if arm == "A" else AMP_MU + AMP_DOSE_COEF * (D - 1)
            a = max(rng.normal(mu, AMP_SD * mu), 0.05)
            y += a * np.exp(-0.5 * ((ts - pt) / 0.6) ** 2)
    base = b + 0.05 * np.abs(rng.standard_normal(len(ts)))
    return ts, k * y + base, N


# ---------- 解码器 ----------
def dec_amp(ts, y):
    return y.max()

def dec_fold(ts, y):
    return y.max() / max(np.median(y[:30]), 1e-9)

def dec_count2(ts, y):
    mad = np.median(np.abs(y - np.median(y))) + 1e-9
    pk, _ = find_peaks(y, prominence=3.0 * mad, distance=int(4.0 / DT))
    return len(pk)

FNS = {"幅度": dec_amp, "fold": dec_fold, "计数": dec_count2}


# ---------- 工具 ----------
def mi_plugin(D, Y, n_bins=12):
    yq = np.quantile(Y, np.linspace(0, 1, n_bins + 1))
    yq[0] -= 1e-9; yq[-1] += 1e-9
    Yb = np.clip(np.digitize(Y, yq) - 1, 0, n_bins - 1)
    H = lambda c: -np.sum((c / c.sum()) * np.log2(c / c.sum() + 1e-12))
    HY = H(np.bincount(Yb, minlength=n_bins))
    HYgD = 0.0
    for d in np.unique(D):
        m = D == d
        HYgD += m.mean() * H(np.bincount(Yb[m], minlength=n_bins))
    return max(HY - HYgD, 0.0)

def lowpass(y, tau):
    a = DT / tau
    out = np.zeros_like(y)
    for i in range(1, len(y)):
        out[i] = out[i - 1] + a * (y[i] - out[i - 1])
    return out


# ---------- 实验 1：腿 1 双臂 ----------
def run_arm(arm, sk, ncell=4000):
    Ds = rng.choice(len(D2), ncell)
    out = {n: [] for n in FNS}
    for d in Ds:
        k = np.exp(rng.normal(0, sk))
        ts, y, N = gen_trace3(D2[d], k, B_BASAL, arm)
        for n, fn in FNS.items():
            out[n].append(fn(ts, y))
    return {n: mi_plugin(Ds, np.array(out[n], float)) for n in FNS}

print("=" * 70)
print("实验1 · 腿1 双臂：仿射接口下谁幸存（σ_k 扫描，4000 细胞/点）")
print("=" * 70)
leg1 = {}
for arm, label in [("A", "数字臂（幅度恒定，p53忠实）"), ("B", "模拟臂（幅度载剂量）")]:
    print(f"--- 臂{arm} {label} ---")
    for sk in SIG_LIST:
        leg1[(arm, sk)] = run_arm(arm, sk)
        r = leg1[(arm, sk)]
        print(f"  σ_k={sk}:  幅度 I={r['幅度']:.3f}  fold I={r['fold']:.3f}  计数 I={r['计数']:.3f}")

# ---------- 实验 2：腿 3（臂A） ----------
print("\n" + "=" * 70)
print("实验2 · 腿3：下游积分器带宽错配杀计数（臂A，σ_k=0.4）")
print("=" * 70)
leg3 = {}
for tau in [0.5, 1.0, 2.0, 3.0, 5.0, 8.0]:
    Ds = rng.choice(len(D2), 2000); Ys = []
    for d in Ds:
        k = np.exp(rng.normal(0, 0.4))
        ts, y, N = gen_trace3(D2[d], k, B_BASAL, "A")
        Ys.append(dec_count2(ts, lowpass(y, tau)))
    leg3[tau] = mi_plugin(Ds, np.array(Ys, float))
    print(f"积分器 τ={tau}h:  I(D;计数)={leg3[tau]:.3f}")
print("（脉冲间距 5.5h：τ 逼近间距时计数信息崩溃）")

# ---------- 实验 3：信息损失恒等式（臂A） ----------
print("\n" + "=" * 70)
print("实验3 · 信息损失恒等式（臂A，σ_k=0.4）")
print("=" * 70)
Ds = rng.choice(len(D2), 6000); Ns = []; Yc = []
for d in Ds:
    k = np.exp(rng.normal(0, 0.4))
    ts, y, N = gen_trace3(D2[d], k, B_BASAL, "A")
    Ns.append(N); Yc.append(dec_count2(ts, y))
HD = np.log2(len(D2))
IDN = mi_plugin(Ds, np.array(Ns, float), n_bins=7)
IDY = mi_plugin(Ds, np.array(Yc, float))
L1, L2 = HD - IDN, IDN - IDY
print(f"H(D)={HD:.3f}  I(D;N)={IDN:.3f}  I(D;Y_count)={IDY:.3f}")
print(f"L1={L1:.3f}  L2={L2:.3f}  恒等式: {IDY:.3f}+{L1:.3f}+{L2:.3f}={IDY+L1+L2:.3f} ≈ H(D)={HD:.3f}")

# ---------- 图 ----------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
fm.fontManager.addfont("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
_avail = {f.name for f in fm.fontManager.ttflist}
for _cand in ["Noto Sans CJK SC", "Noto Sans CJK JP", "Noto Sans CJK HK", "Noto Serif CJK SC"]:
    if _cand in _avail:
        matplotlib.rcParams["font.family"] = [_cand]
        break
matplotlib.rcParams["axes.unicode_minus"] = False

fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
ax = axes[0]
for n, mk in [("幅度", "o"), ("fold", "s"), ("计数", "^")]:
    ax.plot(SIG_LIST, [leg1[("B", sk)][n] for sk in SIG_LIST], mk + "-", label=n + "(臂B)")
    ax.plot(SIG_LIST, [leg1[("A", sk)][n] for sk in SIG_LIST], mk + "--", alpha=0.45, label=n + "(臂A)")
ax.set_xlabel("增益涨落 σ_k"); ax.set_ylabel("I(D; 读出) / bits")
ax.legend(fontsize=7, ncol=2); ax.set_title("腿1: 双臂对照——幅度死(臂B实锤)、fold泄漏、计数活")

ax = axes[1]
taus = sorted(leg3)
ax.plot(taus, [leg3[t] for t in taus], 'o-', c='darkred')
ax.axvline(5.5, ls=':', c='gray'); ax.text(5.6, 0.9, "脉冲间距5.5h", fontsize=8)
ax.set_xlabel("下游积分器时间常数 τ (h)"); ax.set_ylabel("I(D; 计数) / bits")
ax.set_title("腿3: 带宽错配杀计数")

axes[2].bar(["L1\n损伤→脉冲串", "I(D;Y)\n端到端幸存", "L2\n脉冲串→读出"],
            [L1, IDY, L2], color=["#c0392b", "#27ae60", "#e67e22"])
axes[2].set_ylabel("bits"); axes[2].set_title("信息损失分解 (H(D)=2.32)")

plt.tight_layout()
plt.savefig("/mnt/agents/output/03_细胞线3/结果/仿真_p53级联_三腿核验.png", dpi=150)
print("\n图已保存：03_细胞线3/结果/仿真_p53级联_三腿核验.png")
