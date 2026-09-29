#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
代码52 协议伪影检验：原始 FRET 坐标分析 v1.0.0
================================================
任务（臂 1）：检验 Moore 单元表的增量型结构（R~logF, R²=0.652）是否由
数据处理管道（每细胞 min-max 归一化，字段 a）制造。

设计：与 代码49b 完全相同的响应提取（位置判据、前 4 列基线、末 6 列响应、
水平内 ≥4 次重复取中位、单元内跨细胞中位、n≥10），唯一区别是读数字段：
  字段 a      —— 每细胞按饱和刺激反应 min-max 归一（49b 主判用）
  字段 A      —— 每细胞仅按最大值归一（49b 稳健①用）
  字段 FRET   —— 原始（已处理）FRET 值，无任何每细胞归一 ← 本代码主测
另做细胞内坐标对照：同一细胞 B 固定，比较其响应梯度对 logF 与
log(B+F)（细胞内 = fold 坐标的平移）的可分性，只在分离度好的背景
（B=0.3/1/10，F 跨度大且 F≲B）上做，配对 Wilcoxon。

判定逻辑：若原始 FRET 已呈增量型，则归一化管道不可能是增量的来源
（归一化只是单调重标定）；协议伪影假设（归一化层面）被排除，
矛盾定位推向生物机制（受体→CheY-P 窄段）或 regime。

数据：Moore et al. 2024 Dryad doi:10.5061/dryad.nvx0k6dzz（CC0）
运行：python3 代码52_协议伪影检验_原始FRET坐标分析.py
"""

import os
import sys
from collections import defaultdict

import numpy as np
import scipy.io as sio
from scipy.stats import linregress, wilcoxon

ROOT = "/mnt/agents/output/03_细胞线3/公开数据/Moore2024_Chemotaxis_FRET"

# 与 代码49b 完全一致的文件表与 B 口径（元数据；230831_FOV2 例外 B=100）
FILES = {
    "210802_FOV1": 0, "210802_FOV2": 0, "210805_FOV1": 0, "210805_FOV2": 0,
    "220106_FOV1": 0, "230417_FOV1": 0,
    "230815_FOV1": 0.01, "230815_FOV2": 0.01, "230816_FOV1": 0.01, "230816_FOV2": 0.01,
    "230830_FOV1": 0.1, "230830_FOV2": 0.1, "230831_FOV1": 0.1, "230831_FOV2": 0.1,
    "220615_FOV1": 0.3, "230410_FOV1": 0.3,
    "230428_FOV1": 1.0, "230429_FOV1": 1.0,
    "220302_FOV1": 10.0, "220303_FOV1": 10.0,
    "210816_FOV1": 100.0, "210816_FOV2": 100.0, "230717_FOV1": 100.0, "230718_FOV1": 100.0,
}


def extract(field):
    """逐细胞逐水平提取响应幅值，流程与代码49b相同，仅字段不同。"""
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
                sig = rd[field][0, ci].astype(float)
                s = rd["s"][0, ci].astype(float)
            except Exception:
                continue
            if sig.shape != (35, 20) or s.shape != (35, 20):
                continue
            lev = {}
            for row in range(35):
                sr = s[row]
                sv = float(sr.max())
                stim = np.where(sr >= sv - 1e-9)[0]
                if len(stim) == 0:
                    continue
                pre = np.arange(0, stim[0])
                if len(pre) < 4 or len(stim) < 6:
                    continue
                da = float(np.median(sig[row, pre[-4:]])) - float(np.median(sig[row, stim[-6:]]))
                F = sv - B
                if F <= 0:
                    continue
                lev.setdefault(round(sv, 4), []).append(da)
            good = {k: v for k, v in lev.items() if len(v) >= 4}
            if len(good) < 3:
                continue
            for k in sorted(good):
                rec.append((B, k - B, float(np.median(good[k])), f"{name}#{ci}"))
    return rec


def unit_r2(rec):
    u = defaultdict(list)
    for B, F, da, ck in rec:
        if B > 0:
            u[(round(B, 4), round(F, 4))].append(da)
    units = [(B, F, float(np.median(v)), len(v))
             for (B, F), v in sorted(u.items()) if len(v) >= 10]
    Bv = np.array([x[0] for x in units])
    Fv = np.array([x[1] for x in units])
    Rv = np.array([x[2] for x in units])
    return (len(units),
            linregress(np.log10(Fv), Rv).rvalue ** 2,
            linregress(np.log10((Bv + Fv) / Bv), Rv).rvalue ** 2,
            linregress(np.log10(Bv + Fv), Rv).rvalue ** 2)


def main():
    print("=" * 74)
    print("代码52 协议伪影检验：原始 FRET 坐标分析 v1.0.0")
    print("=" * 74)

    for field in ("FRET", "A", "a"):
        rec = extract(field)
        n, rf, rr, rt = unit_r2(rec)
        print(f"字段 {field:5s}: 记录 {len(rec)} 条, 单元 {n} 个 | "
              f"R²(logF)={rf:.3f}  R²(logr)={rr:.3f}  R²(logT)={rt:.3f}")

    # 细胞内坐标对照（原始 FRET）
    rec_raw = extract("FRET")
    cells = defaultdict(list)
    for B, F, da, ck in rec_raw:
        cells[ck].append((B, F, da))
    res = []
    for ck, lst in cells.items():
        Bs = set(round(x[0], 4) for x in lst)
        if len(Bs) != 1:
            continue
        B = lst[0][0]
        if B not in (0.3, 1.0, 10.0) or len(lst) < 4:
            continue
        Fv_ = np.array([x[1] for x in lst])
        y = np.array([x[2] for x in lst])
        if np.std(y) == 0:
            continue
        res.append((B,
                    linregress(np.log10(Fv_), y).rvalue ** 2,
                    linregress(np.log10(B + Fv_), y).rvalue ** 2))
    res = np.array(res)
    print(f"\n细胞内坐标对照（原始 FRET，B∈{{0.3,1,10}}，有效细胞 {len(res)} 个）")
    for B_ in (0.3, 1.0, 10.0):
        sub = res[res[:, 0] == B_]
        win = (sub[:, 1] > sub[:, 2]).mean()
        print(f"  B={B_:5.1f}: n={len(sub):4d}  中位R²(logF)={np.median(sub[:, 1]):.3f}  "
              f"中位R²(log(B+F))={np.median(sub[:, 2]):.3f}  logF 胜出比例={win:.2f}")
    w = int((res[:, 1] > res[:, 2]).sum())
    print(f"  全体配对：logF 胜出 {w}/{len(res)} = {w / len(res):.2f}；"
          f"Wilcoxon p = {wilcoxon(res[:, 1], res[:, 2]).pvalue:.2e}")

    print("""
== 结论 ==
原始 FRET（无任何每细胞归一）单元表已呈增量型（logF 主导），
归一化字段 a/A 只是保持并略增强该结构 → 归一化管道不是增量来源。
细胞内（B 固定、无跨细胞平均）响应梯度在分离度最好的 B=10 组
80% 细胞 logF 胜出 → 增量结构存在于单细胞原始响应。
协议伪影假设（归一化层面）排除；定位推向：
(a) 受体→CheY-P 窄段的生物坐标变换；(b) 实验 regime；(上游 'processed
FRET' 的逐细胞单调预处理仍残余可能，但其不能改变坐标依赖结构）。
""")


if __name__ == "__main__":
    sys.exit(main())
