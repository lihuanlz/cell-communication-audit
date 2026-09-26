#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
代码64 D2R 上游速率排序检验（原始时程深挖）v1.0.0 —— 数据质量不足，如实记负
================================================
任务：检验"通路动力学是否配体无关"。单活性态+平衡级联预言：同一通路
内各配体的归一化响应波形应当重合（速率排序是通路属性）。若配体给
通路留下动力学指纹，则存在更深层的态分辨。

数据：Klein Herenbrink 2016 Supplementary Data 1（6 通路的单浓度
饱和刺激时程，n=3–4 次重复/配体）。

结果（v1.0.0 实跑）：
  直接 t50 提取给出配体间巨大差异（如 Gαi1: dopamine 0.9 min vs
  bifeprunox 17.9 min）——看似"配体指纹"。
  ⚠️ 但质量控制暴露三个致命问题：
  1. 单调度普遍 0.3–0.7（轨迹漂移/噪声大，非干净单相上升）；
  2. 平台归一化振幅出现物理不合理值（pERK bifeprunox=多巴胺的572%、
     Gαi1 aripiprazole=804%），与 SI Table 7 的 Emax 排序矛盾——
     说明基线/平台估计被漂移污染；
  3. pERK 时程窗口从 30 min 才开始（捕获的是衰减尾），其 t50 无定义。
判定：**该批原始时程不足以支持动力学排序检验**——t50 的表观配体
差异不可信，既不能作为断裂证据，也不能作为无断裂证据。如实记负。
动力学断裂的可用证据保持为代码63 检验3（SI T7 的 pEC50 时程，
那是逐时点完整浓度-响应拟合，统计量级完全不同）。

教训（双录）：单浓度时程 + 低重复 + 基线漂移的组合不能做波形归一化
比较；此类检验需要逐时点浓度-响应曲面（time × dose 矩阵）。

运行：python3 代码64_D2R上游速率排序检验.py
依赖：numpy, pandas, openpyxl
"""

import numpy as np
import pandas as pd

XLSX = "/mnt/agents/output/04_细胞线4/公开数据/D2R_KleinHerenbrink2016/SupplementaryData1.xlsx"

def parse_sheet(sh):
    df = pd.read_excel(XLSX, sheet_name=sh, header=None)
    t = pd.to_numeric(df.iloc[3:, 0], errors="coerce")
    curves = {}
    c = 1
    while c < df.shape[1]:
        name = df.iloc[1, c]
        if pd.isna(name):
            c += 1; continue
        cols = [cc for cc in range(c, min(c + 4, df.shape[1]))
                if str(df.iloc[2, cc]).startswith("n")]
        if not cols:
            c += 1; continue
        Y = df.iloc[3:, cols].apply(pd.to_numeric, errors="coerce")
        ok = ~(t.isna() | Y.mean(axis=1).isna())
        curves[str(name)] = (t[ok].to_numpy(float),
                             Y[ok].mean(axis=1).to_numpy(float))
        c = max(cols) + 1
    return curves

def t50(t, y):
    base = np.mean(y[:3]); plat = np.mean(y[int(len(y) * 0.8):])
    yn = (y - base) / (plat - base + 1e-12)
    for i in range(1, len(yn)):
        if (yn[i-1] - 0.5) * (yn[i] - 0.5) <= 0 and yn[i] != yn[i-1]:
            return t[i-1] + (0.5 - yn[i-1]) * (t[i]-t[i-1]) / (yn[i]-yn[i-1])
    return np.nan

def main():
    sheets = {sh: parse_sheet(sh)
              for sh in ["CAMYEL cAMP", "pERK12", "B-arrestin2", "Gai1", "Gao"]}
    print("== 表观 t50（min）与质量控制 ==")
    for sh, cu in sheets.items():
        spans = {lig: np.mean(y[int(len(y)*0.8):]) - np.mean(y[:3])
                 for lig, (t, y) in cu.items()}
        d0 = abs(spans.get("Dopamine", 1)) or 1
        print(f"\n{sh}")
        for lig, (t, y) in cu.items():
            base = np.mean(y[:3]); plat = np.mean(y[int(len(y)*0.8):])
            yn = (y - base) / (plat - base + 1e-12)
            mono = np.mean(np.diff(yn) >= -0.05)
            print(f"  {lig:13s}: t50={t50(t,y):6.1f} min, 幅={100*abs(spans[lig])/d0:6.1f}% DA, "
                  f"单调度={mono:.2f}, 窗口=[{t.min():.0f},{t.max():.0f}] min")
    print("""
== 判定 ==
单调度 0.3–0.7、归一化振幅物理不合理（pERK 572%、Gαi1 804%）、
pERK 窗口只含衰减尾 → 表观 t50 配体差异不可信。
该批单浓度时程不足以做动力学排序检验（记负，双录）。
动力学断裂的有效证据保持为代码63 检验3（SI T7 pEC50 时程）。
""")

if __name__ == "__main__":
    main()
