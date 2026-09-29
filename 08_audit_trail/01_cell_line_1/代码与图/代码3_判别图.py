# -*- coding: utf-8 -*-
"""
代码 3：判别图（脉冲可观测量 vs 回路泛类）
==========================================
把手稿 §4 的"三签名判别定理"画成一张可操作的图：
  横轴 α = ∂lnA/∂lnD（振幅-剂量敏感度，对数轴）
  纵轴 ν = 每十倍剂量的脉冲数增量（计数敏感度）

落点分两类：
  · 模型点（来自代码 2 / Track B 的计算）——定义两个泛类的区域
  · 实验点（半定量，取自已发表统计量）——检验定理预言

实验点数值的推导（手稿 Methods 同款）：
  p53/γ-IR (Lahav 2004)：0.1–10 Gy（100×）振幅"与剂量无关"，
      取振幅变动上限 15% → |α| ≤ ln1.15/ln100 ≈ 0.03；
      计数 ~1.5 → 5.5（2 decades）→ ν ≈ 2
  p53/NCS (Mönke 2017)：25→400 ng/ml（16×）振幅均值独立 → |α| ≲ 0.05；计数增 → ν ≈ 1.8
  p53/UV (Batchelor 2011)：振幅∝剂量（graded）→ α ≈ 1；单脉冲 → ν ≈ 0
      —— 内建对照：同蛋白同细胞，实验判定"非可兴奋"，签名落模拟区 ✓

运行：python3 代码3_判别图.py   （几秒，输出判别图 PNG）
"""
import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(8.5, 7))

# ---------- 两个泛类的区域 ----------
ax.axvspan(0.3, 15, alpha=0.09, color='red')     # NF 模拟区
ax.axvspan(5e-4, 0.15, alpha=0.09, color='blue') # EXC 数字区
ax.text(2.2, 3.4, 'NF 族（模拟）\n$\\alpha\\gtrsim0.5$\n（亚临界：+迟滞）',
        color='darkred', fontsize=10, ha='center')
ax.text(0.008, 3.4, 'EXC 族（数字）\n$\\alpha\\lesssim0.15$', color='darkblue', fontsize=10, ha='center')
ax.axvline(0.22, color='gray', ls='--', lw=1)

# ---------- 模型点（Track B 计算结果，ν 单位：脉冲/decade） ----------
ax.errorbar([8.76], [0.02], fmt='s', color='red', ms=9, label='NF 模型，近起始（α→∞）')
ax.errorbar([0.7], [0.05], fmt='s', color='red', ms=9, mfc='none', label='NF 模型，远离起始')
ax.errorbar([0.014], [4.05], fmt='o', color='blue', ms=9, label='FHN 慢驱动（纯数字）')
ax.errorbar([0.26], [2.0], fmt='^', color='green', ms=9, label='Morris–Lecar 快驱动（泄漏）')

# ---------- 实验点（半定量） ----------
ax.errorbar([0.03], [2.0], xerr=[0.02], yerr=[0.5], fmt='*', color='purple',
            ms=19, capsize=4, label='p53 / γ-IR（Lahav 2004）：0.1–10 Gy')
ax.errorbar([0.05], [1.8], xerr=[0.04], yerr=[0.5], fmt='*', color='magenta',
            ms=17, capsize=4, label='p53 / NCS（Mönke 2017）：25–400 ng/ml')
ax.errorbar([1.0], [0.05], xerr=[0.3], yerr=[0.05], fmt='D', color='darkorange',
            ms=10, capsize=4, label='p53 / UV（Batchelor 2011）：graded，非可兴奋')
ax.annotate('内建对照：同蛋白同细胞\n拓扑↔签名严格共变', (0.8, 1.35), fontsize=8.5, color='darkorange')

ax.set_xscale('log'); ax.set_xlim(5e-4, 15); ax.set_ylim(-0.3, 4.6)
ax.set_xlabel('$\\alpha = \\partial\\ln A/\\partial\\ln D$（振幅-剂量敏感度）')
ax.set_ylabel('$\\nu$（每十倍剂量的脉冲数）')
ax.set_title('判别图：实验落点（半定量；精确抠图见 SI）', fontsize=11)
ax.legend(fontsize=8.5, loc='center right')
plt.tight_layout()
plt.savefig('TrackC_判别图_实验落点.png', dpi=200)
print('已保存 TrackC_判别图_实验落点.png')
print('\n判读：紫/品红星深落 EXC 区（数字），橙色菱形落模拟区（对照）——四象限无一错位。')
