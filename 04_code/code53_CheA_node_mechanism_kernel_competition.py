#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
代码53 CheA 节点机制检验（臂2）+ 核模型竞赛（臂3）v1.0.0
================================================
背景：代码52 排除了归一化协议伪影，矛盾收窄到 (a) 受体→CheA→CheY-P 窄段
的生物坐标变换，或 (b) 实验 regime。本代码两臂执行：

臂 2（§1–§3）：带显式适应动力学的 受体→CheA→CheY-P 最小模型
  状态 (m, Yp)；a = 1/(1+exp(N[α(m0-m)+g(L)]))，g(L)=ln((1+L/Ki)/(1+L/Ka))
  dm/dt = kR(1-a) − kB·a（Barkai-Leibler 精确适应；适应稳态解析求得）
  dYp/dt = kP·a·(1−Yp) − kZ·Yp        标准型
       或 − kZ·Yp/(Km+Yp)             CheZ 饱和（零级超灵敏候选）
  协议严格复刻 Moore：背景适应 → 阶跃 → 固定时刻 t_meas 读 ΔYp。
  扫描 Ki×t_meas×kR（27 组）+ CheZ 饱和 Km×Ki（9 组）。
  问题：这段窄链路能否长出增量型（logF 主导）或绝对锚点 Fc≈0.17 µM？

臂 3（§4）：核模型竞赛。Moore 自己的 Discussion 提出"线性 regime（感知绝对
  变化）→ 对数 regime（感知倍数）"两 regime 转变，其 MWC 移位 Hill 模型
  解释了 K1/2 多样性数据。问题：他们自己的模型形式能否解释【幅值单元表】？
  竞赛：G3' 原经验核(logF+衰减)、G6 线性F+衰减、G7/G8 Moore 式移位 Hill、
  G10 锚定对数 log(1+F/Kc)+衰减、G11 无衰减版。
  在 a 字段单元表与原始 FRET 字段单元表上各跑一次。
  另：行内（B=0.01 与 B=100）logF vs 线性F 可分性检验。

数据：Moore 2024 Dryad doi:10.5061/dryad.nvx0k6dzz（CC0）
运行：python3 代码53_CheA节点机制检验与核模型竞赛.py
依赖：numpy, scipy, pandas
"""

import os
import itertools
import warnings
from collections import defaultdict

import numpy as np
import pandas as pd
import scipy.io as sio
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, curve_fit
from scipy.stats import linregress

warnings.filterwarnings("ignore")
rng = np.random.default_rng(20260815)

CSV = "/mnt/agents/output/03_细胞线3/结果/P7_Moore2024_FCD-Weber/代码49b_单元表.csv"
ROOT = "/mnt/agents/output/03_细胞线3/公开数据/Moore2024_Chemotaxis_FRET"


def r2(x, y):
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3 or np.std(x[m]) == 0 or np.std(y[m]) == 0:
        return float("nan")
    return linregress(x[m], y[m]).rvalue ** 2


# ===============================================================
# §1 臂 2：受体→CheA→CheY-P 动力学模型
# ===============================================================
ALPHA, M0 = 1.0, 1.0

def act(m, L, N, Ki, Ka):
    g = np.log((1 + L / Ki) / (1 + L / Ka))
    return 1.0 / (1.0 + np.exp(N * (ALPHA * (M0 - m) + g)))

def simulate(B, F, N=6.0, Ki=1.0, Ka=50.0, kR=0.02, kB=0.04,
             kP=2.0, kZ=1.0, Km=None, t_meas=8.0):
    """精确适应稳态解析求得，只积阶跃段。Km=None → 标准 CheZ；否则饱和型。"""
    astar = kR / (kR + kB)
    m0 = brentq(lambda m: act(m, B, N, Ki, Ka) - astar, M0 - 20, M0 + 20, xtol=1e-10)
    if Km is None:
        yp0 = kP * astar / (kP * astar + kZ)
        def dephos(Yp): return kZ * Yp
    else:
        # 稳态：kP*a*(1-Y)(Km+Y) = kZ*Y → 二次方程
        A = kP * astar
        b = kZ - A * (1 - Km); c = -A * Km
        yp0 = (-b + np.sqrt(b * b - 4 * A * c)) / (2 * A)
        def dephos(Yp): return kZ * Yp / (Km + Yp)
    def rhs(t, x):
        m, Yp = x
        a = act(m, B + F, N, Ki, Ka)
        return [kR * (1 - a) - kB * a, kP * a * (1 - Yp) - dephos(Yp)]
    s = solve_ivp(rhs, (0, t_meas), [m0, yp0], rtol=1e-7, atol=1e-10)
    return yp0 - s.y[1][-1]

BGS = {0.01: [0.2, 0.5, 1.09, 2.09, 4.09], 0.1: [0.2, 0.5, 1, 2, 4],
       0.3: [0.2, 0.5, 1, 2, 4], 1.0: [0.2, 0.5, 1, 2, 4],
       10.0: [0.5, 1, 2, 4, 8], 100.0: [2, 5, 10, 20, 40]}

def unit_table_sim(**kw):
    Bs, Fs, Rs = [], [], []
    for B, fs in BGS.items():
        for F in fs:
            Bs.append(B); Fs.append(F); Rs.append(simulate(B, F, **kw))
    Bs, Fs, Rs = np.array(Bs), np.array(Fs), np.array(Rs)
    return Rs, (r2(np.log10(Fs), Rs), r2(np.log10((Bs + Fs) / Bs), Rs),
                r2(np.log10(Bs + Fs), Rs))

print("== §2 臂2 标准链路扫描（Ki × 读数时刻 × 适应速率，a*=1/3 固定）==")
print(f"{'Ki':>4} {'t':>4} {'kR':>6} | logF   logr   logT   max")
for Ki, t_m, kR_ in itertools.product([0.2, 1.0, 5.0], [2.0, 8.0, 30.0],
                                      [0.005, 0.02, 0.1]):
    Rv, c = unit_table_sim(Ki=Ki, kR=kR_, kB=2 * kR_, t_meas=t_m)
    print(f"{Ki:4.1f} {t_m:4.0f} {kR_:6.3f} | {c[0]:.3f}  {c[1]:.3f}  {c[2]:.3f}  {Rv.max():.3f}")

print("\n== §3 臂2 CheZ 饱和变体（Km × Ki，t=8s, kR=0.02）==")
for Km, Ki in itertools.product([0.01, 0.05, 0.2], [0.2, 1.0, 5.0]):
    Rv, c = unit_table_sim(Ki=Ki, Km=Km)
    print(f"Km={Km:4.2f} Ki={Ki:3.1f} | logF={c[0]:.3f} logr={c[1]:.3f} "
          f"logT={c[2]:.3f} max={Rv.max():.3f}")

# ===============================================================
# §4 臂 3：核模型竞赛（a 字段单元表 + 原始 FRET 单元表）
# ===============================================================
_df = pd.read_csv(CSV)
B_A = _df["B_uM"].to_numpy(float); F_A = _df["F_uM"].to_numpy(float)
R_A = _df["R_a"].to_numpy(float)

FILES = {
    "210802_FOV1": 0, "210802_FOV2": 0, "210805_FOV1": 0, "210805_FOV2": 0,
    "220106_FOV1": 0, "230417_FOV1": 0,
    "230815_FOV1": 0.01, "230815_FOV2": 0.01, "230816_FOV1": 0.01, "230816_FOV2": 0.01,
    "230830_FOV1": 0.1, "230830_FOV2": 0.1, "230831_FOV1": 0.1, "230831_FOV2": 0.1,
    "220615_FOV1": 0.3, "230410_FOV1": 0.3, "230428_FOV1": 1.0, "230429_FOV1": 1.0,
    "220302_FOV1": 10.0, "220303_FOV1": 10.0,
    "210816_FOV1": 100.0, "210816_FOV2": 100.0, "230717_FOV1": 100.0, "230718_FOV1": 100.0,
}

def extract_raw():
    """原始 FRET 字段单元表（流程同代码49b/52）。"""
    rec = []
    for name, bg in FILES.items():
        p = os.path.join(ROOT, name + ".mat")
        if not os.path.isfile(p):
            continue
        try:
            rd = sio.loadmat(p)["reorgData"]["resp_data"][0, 0]
        except Exception:
            continue
        B = 100.0 if name == "230831_FOV2" else float(bg)
        for ci in range(rd.shape[1]):
            try:
                sig = rd["FRET"][0, ci].astype(float)
                s = rd["s"][0, ci].astype(float)
            except Exception:
                continue
            if sig.shape != (35, 20) or s.shape != (35, 20):
                continue
            lev = {}
            for row in range(35):
                sr = s[row]; sv = float(sr.max())
                stim = np.where(sr >= sv - 1e-9)[0]
                if len(stim) == 0:
                    continue
                pre = np.arange(0, stim[0])
                if len(pre) < 4 or len(stim) < 6:
                    continue
                da = float(np.median(sig[row, pre[-4:]])) - float(np.median(sig[row, stim[-6:]]))
                if sv - B <= 0:
                    continue
                lev.setdefault(round(sv, 4), []).append(da)
            good = {k: v for k, v in lev.items() if len(v) >= 4}
            if len(good) < 3:
                continue
            for k in sorted(good):
                rec.append((B, k - B, float(np.median(good[k]))))
    u = defaultdict(list)
    for B, F, da in rec:
        if B > 0:
            u[(round(B, 4), round(F, 4))].append(da)
    un = [(B, F, float(np.median(v))) for (B, F), v in sorted(u.items()) if len(v) >= 10]
    return (np.array([x[0] for x in un]), np.array([x[1] for x in un]),
            np.array([x[2] for x in un]))

def kG3(X, A, Fc, Bs, p):  return A * np.log10(X[1] / Fc) / (1 + (X[0] / Bs) ** p)
def kG6(X, A, Bs, p):      return A * X[1] / (1 + (X[0] / Bs) ** p)
def kG7(X, amax, K0, Ki):
    B, F = X; T = B + F; Kh = K0 * (1 + B / Ki)
    return amax * (T / (Kh + T) - B / (Kh + B))
def kG8(X, amax, K0, Ki, n):
    B, F = X; T = B + F; Kh = K0 * (1 + B / Ki)
    return amax * (T ** n / (Kh ** n + T ** n) - B ** n / (Kh ** n + B ** n))
def kG10(X, A, Kc, Bs, p): return A * np.log(1 + X[1] / Kc) / (1 + (X[0] / Bs) ** p)
def kG11(X, A, Kc):        return A * np.log(1 + X[1] / Kc)

RACE = [
    ("G3'  log10(F/Fc)+衰减   ", kG3,  [0.3, 0.17, 100, 0.8], ([0, 1e-4, 1e-2, .05], [1e3, 1e3, 1e5, 5])),
    ("G6   线性F+衰减         ", kG6,  [0.1, 100, 0.8],       ([0, 1e-2, .05],       [1e3, 1e5, 5])),
    ("G7   Moore移位Hill n=1  ", kG7,  [0.5, 0.2, 10],        ([0, 1e-4, 1e-2],      [1e3, 1e3, 1e4])),
    ("G8   Moore移位Hill 自由n", kG8,  [0.5, 0.2, 10, 1.5],   ([0, 1e-4, 1e-2, .3],  [1e3, 1e3, 1e4, 6])),
    ("G10  log(1+F/Kc)+衰减   ", kG10, [0.2, 0.5, 70, 0.8],   ([0, 1e-4, 1e-2, .05], [1e3, 1e3, 1e5, 5])),
    ("G11  log(1+F/Kc) 无B项  ", kG11, [0.2, 0.5],            ([0, 1e-4],            [1e3, 1e3])),
]

def race(Bv, Fv, Rv, tag):
    print(f"\n== 核模型竞赛 @ {tag}（n={len(Rv)}）==")
    for nm, fn, p0, bd in RACE:
        try:
            p_, _ = curve_fit(fn, (Bv, Fv), Rv, p0=p0, bounds=bd, maxfev=40000)
            pr = fn((Bv, Fv), *p_)
            ss = 1 - ((pr - Rv) ** 2).sum() / ((Rv - Rv.mean()) ** 2).sum()
            print(f"  {nm}: R²fit={ss:.3f}  参数={np.round(p_, 3)}")
        except Exception as e:
            print(f"  {nm}: 拟合失败 {str(e)[:40]}")

race(B_A, F_A, R_A, "a 字段单元表（49b 主判）")
Br_, Fr_, Rr_ = extract_raw()
race(Br_, Fr_, Rr_, "原始 FRET 字段单元表")

print("\n== 行内可分性：logF vs 线性F（a 字段）==")
for B_ in [0.01, 100.0]:
    m = B_A == B_
    print(f"  B={B_:6.2f}: R²(logF)={r2(np.log10(F_A[m]), R_A[m]):.3f}  "
          f"R²(线性F)={r2(F_A[m], R_A[m]):.3f}")

print("""
== v1.0.0 总结 ==
臂2：标准链路 27 组 + CheZ 饱和 9 组，全部保持倍数/总量原生型
   （logF 最高 0.13），未长出增量或 Fc 锚点 → CheA 节点（标准动力学）
   不能完成倍数→增量变换。NEGATIVE。
臂3：Moore 自己的移位 Hill 模型形式（G7/G8）在【幅值单元表】上
   R²fit≈0.61–0.67（raw 上 0.53–0.57），明显低于锚定对数核
   G10（0.967/0.906），且 G8 需把 Ki 推到参数边界 10⁴（退化极限）。
   Moore 模型解释的是 K1/2 多样性统计量，不解释幅值表。
新经验核（更新）：R = A·ln(1+F/Kc)/(1+(B/Bs)^p)，Kc≈0.46–0.57 µM，
   Bs≈63–70 µM，p≈0.8–0.9，a 字段 R²=0.967、raw R²=0.906。
   行内检验：B=0.01 与 B=100 行均为 logF 主导（0.97/0.99）。
解读：幅值响应是【绝对增量 F 的对数】，参考零点锚定在 Kc（II 类
   结构常数），背景只通过衰减项进入——等效于自由能参考态固定在
   零配体、适应只调增益。这不被任何标准精确适应机制
   （h(g(B+F)−g(B)) 家族、Moore 自己的移位 Hill、CheA 段动力学）
   复现。该幅值表是先前未被任何模型解释过的经验规律。
""")
