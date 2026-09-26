#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
代码70 混合配体验证臂 P-1/P-2 判决 v1.0.0
================================================
预注册文件：04_细胞线4/结果/预注册_混合配体验证臂_P1P2.md
（2026-08-16 封存于数据拆封之前；本代码严格按其中管线与判决线执行）

验证臂（Moore 2024 Dryad，本轮经浏览器过 Anubis 验证后下载）：
  V1: 230822_FOV1, 230823_FOV1/2   背景 10µM Lasp+100µM meAsp，前景 meAsp
  V2: 231017_FOV1/2, 231018_FOV1/2 背景 100µM meAsp+0.1µM Lasp，前景 Lasp
有效背景折算（§18.3/§19.4，ρ=23）：V1 B_eff=330，V2 B_eff=102.3 µM meAsp 当量。

预测（预注册写定）：
  P-1a 结构：logF 坐标 R² − logr 坐标 R² ≥ +0.3
  P-1b 数值：V1 da(F)=0.0677·ln(1+F/0.570)；V2 da(F_L)=0.1682·ln(1+40.35·F_L)
        容差 |实测中位−预测| ≤ 0.08（全体共有 F 档）
  P-2  K1/2（总有效轴）：V1=388 µM（±0.15 dex；若 F 上限不足则预言高删失）；
        V2=122 µM 当量=前景 0.86 µM Lasp（±0.15 dex，可估需 n≥20）
管线与代码58 完全相同（a 字段，da=刺激前4点中位−刺激末6点中位；
逐细胞 K1/2：log10 轴线性插值过 0.5，D.max()<0.5 删失）。
运行：python3 代码70_混合配体验证臂_P1P2判决.py
================================================
"""
import os, warnings
from collections import defaultdict
import numpy as np
import scipy.io as sio

warnings.filterwarnings("ignore")
ROOT = "/mnt/agents/output/03_细胞线3/公开数据/Moore2024_Chemotaxis_FRET"

ARMS = {
    # s_sub: 从 s 字段档位得到前景 F 需减去的量。
    # V1 前景 meAsp 与背景主成分同配体 → s 记录总 meAsp（105=100+5），F=sv−100；
    # V2 前景 Lasp 跨配体 → s 记录前景（§18 数据事实），但背景含 0.1 Lasp，
    #   存在"总 Lasp"歧义 → 两种读法都算并比较拟合优度（见输出 §1b 歧义核验）。
    "V1": {"files": ["230822_FOV1", "230823_FOV1", "230823_FOV2"],
           "bg_txt": "10uM Lasp + 100uM meAsp", "B_eff": 330.0, "rho_fg": 1.0, "s_sub": 100.0},
    "V2": {"files": ["231017_FOV1", "231017_FOV2", "231018_FOV1", "231018_FOV2"],
           "bg_txt": "100uM meAsp + 0.1uM Lasp", "B_eff": 102.3, "rho_fg": 23.0, "s_sub": 0.0},
    "V2_alt": {"files": ["231017_FOV1", "231017_FOV2", "231018_FOV1", "231018_FOV2"],
           "bg_txt": "100uM meAsp + 0.1uM Lasp（s 按总 Lasp 读）", "B_eff": 102.3, "rho_fg": 23.0, "s_sub": 0.1},
}
A_, Kc_, Bs_, p_ = 0.2268, 0.570, 69.7, 0.778   # G10（P7 重拟合，预注册数值）
RHO = 23.0
TOL_AMP, TOL_K12 = 0.08, 0.15

def g10(F_eff, B_eff):
    return A_ * np.log(1 + F_eff / Kc_) / (1 + (B_eff / Bs_) ** p_)

def unit_table(files, s_sub=0.0):
    lev_all, percell = defaultdict(list), []
    for name in files:
        p = os.path.join(ROOT, name + ".mat")
        d = sio.loadmat(p)["reorgData"]["resp_data"][0, 0]
        for ci in range(d.shape[1]):
            try:
                a = d["a"][0, ci].astype(float); s = d["s"][0, ci].astype(float)
            except Exception:
                continue
            if a.shape != (35, 20) or s.shape != (35, 20):
                continue
            lev = defaultdict(list)
            for row in range(35):
                sr = s[row]; sv = float(sr.max())
                stim = np.where(sr >= sv - 1e-9)[0]
                pre = np.arange(0, stim[0]) if len(stim) else None
                if pre is None or len(pre) < 4 or len(stim) < 6:
                    continue
                da = float(np.median(a[row, pre[-4:]])) - float(np.median(a[row, stim[-6:]]))
                lev[round(sv - s_sub, 4)].append(da)
            good = {k: float(np.median(v)) for k, v in lev.items() if len(v) >= 4}
            if len(good) >= 3:
                percell.append(good)
                for k, v in good.items():
                    lev_all[k].append(v)
    return lev_all, percell

def percell_k12(percell, B_eff, rho_fg):
    """总有效轴 K1/2；返回 (可估值数组, 细胞总数)"""
    out = []
    for good in percell:
        ks = sorted(good)
        T = B_eff + rho_fg * np.array(ks)
        D = np.array([good[k] for k in ks])
        if D.max() < 0.5:
            continue
        i = int(np.argmax(D >= 0.5))
        if i == 0:
            continue
        x0, x1 = np.log10(T[i - 1]), np.log10(T[i])
        y0, y1 = D[i - 1], D[i]
        out.append(10 ** (x0 + (0.5 - y0) * (x1 - x0) / (y1 - y0)))
    return np.array(out), len(percell)

print("== 0. 文件身份核验（s 字段档位 vs BackgroundList 条件）==")
for arm, cfg in ARMS.items():
    for name in cfg["files"]:
        d = sio.loadmat(os.path.join(ROOT, name + ".mat"))["reorgData"]["resp_data"][0, 0]
        s0 = d["s"][0, 0].astype(float)
        smm = d["s_min_max"][0, 0].astype(float)
        print(f"  {name}: s 行档位={sorted(set(np.round(s0.max(axis=1),3)))[:8]}... "
              f"s_min_max 值域=[{smm.min():.3g},{smm.max():.3g}]  (背景: {cfg['bg_txt']})")

tables = {arm: unit_table(cfg["files"], cfg.get("s_sub", 0.0)) for arm, cfg in ARMS.items()}

print("\n== 1. 幅值单元表（前景轴）与 P-1b 零拟合对账 ==")
verdict_p1b = {}
for arm, cfg in ARMS.items():
    lev_all, percell = tables[arm]
    tag = "（歧义对照，不进判决）" if arm.endswith("_alt") else ""
    print(f"\n[{arm}] {len(percell)} 细胞；B_eff={cfg['B_eff']} µM 当量 {tag}")
    devs = []
    for k in sorted(lev_all):
        obs = float(np.median(lev_all[k]))
        pred = g10(cfg["rho_fg"] * k, cfg["B_eff"])
        devs.append(abs(obs - pred))
        print(f"   F={k:8.3f}: 实测 da={obs:+.3f} (n={len(lev_all[k])})  预测={pred:+.3f}  偏差={obs-pred:+.3f}")
    mx = max(devs)
    verdict_p1b[arm] = mx <= TOL_AMP
    print(f"   → 最大偏差 {mx:.3f}，容差 {TOL_AMP}：{'命中' if verdict_p1b[arm] else '未命中'}")

print("\n== 2. P-1a 结构判决（坐标检验）==")
for arm, cfg in ARMS.items():
    lev_all, _ = tables[arm]
    ks = sorted(lev_all)
    F = np.array(ks); R = np.array([np.median(lev_all[k]) for k in ks])
    F_eff = cfg["rho_fg"] * F
    def r2(xv, yv):
        A = np.vstack([xv, np.ones_like(xv)]).T
        c = np.linalg.lstsq(A, yv, rcond=None)[0]
        r = yv - A @ c
        return 1 - np.sum(r**2) / np.sum((yv - yv.mean())**2)
    r2_logF = r2(np.log10(F_eff), R)
    r2_logr = r2(np.log10((cfg["B_eff"] + F_eff) / cfg["B_eff"]), R)
    ok = (r2_logF - r2_logr) >= 0.3
    print(f"[{arm}] R²(logF)={r2_logF:.3f}  R²(logr)={r2_logr:.3f}  差={r2_logF-r2_logr:+.3f}  "
          f"{'命中(≥+0.3)' if ok else '未命中'}")

print("\n== 3. P-2 逐细胞 K1/2（总有效轴）==")
for arm, cfg in ARMS.items():
    lev_all, percell = tables[arm]
    v, ntot = percell_k12(percell, cfg["B_eff"], cfg["rho_fg"])
    fmax = max(lev_all)
    pred_total = 1.17 * (1.95 + cfg["B_eff"])
    print(f"\n[{arm}] 可估 {len(v)}/{ntot} 细胞（删失 {100*(1-len(v)/ntot):.0f}%）；"
          f"前景 F 上限={fmax}（实测轴）")
    print(f"   预注册预测：总轴 K1/2={pred_total:.1f} µM 当量 ±{TOL_K12} dex")
    if len(v) >= 20:
        med = float(np.median(v))
        dd = abs(np.log10(med / pred_total))
        print(f"   实测中位={med:.1f} µM 当量，IQR=[{np.percentile(v,25):.1f},{np.percentile(v,75):.1f}]")
        print(f"   偏差={np.log10(med/pred_total):+.3f} dex → {'命中' if dd<=TOL_K12 else '未命中'}")
        if cfg["rho_fg"] != 1:
            fg = (med - cfg["B_eff"]) / cfg["rho_fg"]
            print(f"   前景轴读数：F*={fg:.3f} µM Lasp（预注册点估计 0.86）")
    else:
        print(f"   可估细胞 <20 → 按预注册规则：无法判决（若为高删失则与 V1 删失预言定性一致："
              f"{'是' if len(v)/ntot<0.5 else '否'}）")

print("\n== 4. 附带观察项（断裂侧，非判决）：幅值坐标结构在混合背景下是否保持 ==")
print("（见第2节 R² 差；断裂形式化检验留待联合拟合，本代码不设判决）")
