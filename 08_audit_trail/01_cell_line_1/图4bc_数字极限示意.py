# -*- coding: utf-8 -*-
"""图 4(b,c) | 共享数字极限示意 + 决策元件 vs 测量元件（手稿图 4 补全）"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig = plt.figure(figsize=(14, 4.8))

# ---------- (b) 共享数字极限：两条路径，一个极限 ----------
axb = fig.add_subplot(1, 2, 1)
axb.set_xlim(0, 10); axb.set_ylim(0, 6); axb.axis('off')
axb.add_patch(FancyBboxPatch((0.2, 3.6), 3.2, 1.7, boxstyle='round,pad=0.08', fc='#eaf1fb', ec='#1f77b4', lw=1.8))
axb.text(1.8, 4.45, '测量 / 共价循环\n化学饱和 κ→0\n（标度因子从似然消去）', fontsize=9.5, ha='center', va='center')
axb.add_patch(FancyBboxPatch((0.2, 0.7), 3.2, 1.7, boxstyle='round,pad=0.08', fc='#fdeaea', ec='#d62728', lw=1.8))
axb.text(1.8, 1.55, '脉冲编码器\n拓扑解耦 α→0\n（振幅与剂量脱钩）', fontsize=9.5, ha='center', va='center')
axb.add_patch(FancyBboxPatch((5.6, 2.1), 4.2, 1.8, boxstyle='round,pad=0.1', fc='#f4f8e8', ec='#6b8e23', lw=2.2))
axb.text(7.7, 3.0, '同一个数字极限\n似然退化为计数统计\nFisher 谱：单一非零方向沿计数\n（条件定理，SI S2）', fontsize=10, ha='center', va='center', fontweight='bold', color='#4a6512')
axb.add_patch(FancyArrowPatch((3.5, 4.45), (5.5, 3.6), arrowstyle='-|>', mutation_scale=20, color='#1f77b4', lw=2))
axb.add_patch(FancyArrowPatch((3.5, 1.55), (5.5, 2.4), arrowstyle='-|>', mutation_scale=20, color='#d62728', lw=2))
axb.set_title('(b) 共享数字极限：κ→0 ≡ α→0', fontsize=12)

# ---------- (c) 模拟–数字轴上的决策元件与测量元件 ----------
axc = fig.add_subplot(1, 2, 2)
axc.set_xlim(-0.5, 10.5); axc.set_ylim(0, 6); axc.axis('off')
# 轴线
axc.add_patch(FancyArrowPatch((0.3, 1.0), (10.2, 1.0), arrowstyle='-|>', mutation_scale=22, color='#333', lw=2))
axc.text(0.3, 0.45, '模拟区', fontsize=11, ha='center', fontweight='bold')
axc.text(10.2, 0.45, '数字极限', fontsize=11, ha='center', fontweight='bold')
axc.text(5.2, 0.45, '校准需求 → 信息丰富但校准饥渴（内点最优 κ*）', fontsize=8.5, ha='center', color='#666')
# 决策元件（静态数字）
axc.add_patch(FancyBboxPatch((7.6, 3.6), 2.7, 1.5, boxstyle='round,pad=0.08', fc='#fff4e0', ec='#e69500', lw=1.8))
axc.text(8.95, 4.35, '决策元件\nGK 零级开关（静态数字）\n≈1 bit/快照：阈值穿越报告', fontsize=8.8, ha='center', va='center')
axc.add_patch(FancyArrowPatch((8.95, 3.55), (9.4, 1.15), arrowstyle='-|>', mutation_scale=14, color='#e69500'))
# 测量元件（动态数字）
axc.add_patch(FancyBboxPatch((7.3, 1.7), 3.1, 1.4, boxstyle='round,pad=0.08', fc='#e8f6ef', ec='#2ca02c', lw=1.8))
axc.text(8.85, 2.4, '测量元件\np53 脉冲串（动态数字）\n对数范围计数编码：N ∝ ln D', fontsize=8.8, ha='center', va='center')
axc.add_patch(FancyArrowPatch((8.85, 1.65), (9.8, 1.15), arrowstyle='-|>', mutation_scale=14, color='#2ca02c'))
# 模拟区标记
axc.add_patch(FancyBboxPatch((0.7, 2.2), 2.9, 1.6, boxstyle='round,pad=0.08', fc='#eef0f3', ec='#666', lw=1.5))
axc.text(2.15, 3.0, '模拟区\nPdPC 中段 κ、NF 振子\n信息富有 ∝ 校准饥渴', fontsize=8.8, ha='center', va='center')
axc.add_patch(FancyArrowPatch((2.15, 2.15), (2.15, 1.15), arrowstyle='-|>', mutation_scale=14, color='#666'))
axc.set_title('(c) 决策元件 vs 测量元件：两种数字，两种功能', fontsize=12)

plt.tight_layout()
plt.savefig('图4bc_数字极限示意.png', dpi=200)
print('已保存 图4bc_数字极限示意.png')
