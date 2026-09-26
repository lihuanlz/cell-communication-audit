# -*- coding: utf-8 -*-
"""
代码 11：判别图 v0.9（模型侧 + 实验侧落点）
==========================================
代码 3 的 v0.9 更新：实验点从半定量升级为公开数据审计定量值
（Batchelor 2011 源图数字化；Albeck 2013 源图数字化；Tay 2010 原文定量句 + SI 读数），
并新增 ERK / NF-κB 实验落点与 NF-κB 模型-实验落差箭头。

落点数值来源：
  实验侧（星/菱形标记）：
    p53/γ-IR (Lahav 2004)   α≈0.03（报告统计量上限）, ν≈2
    p53/NCS (Mönke 2017)    α≈0.05（报告统计量上限）, ν≈1.8
    p53/NCS (Batchelor 2011, 源图) α=0.00±0.28（对数轴取 0.008 绘制）, ν≈2（计数随 DSB 增加）
    p53/UV  (Batchelor 2011, 源图) α=+0.99±0.32, ν≈0（单脉冲）
    ERK     (Albeck 2013, 源图)    |α|≲0.05（取 0.02 绘制）, ν≈+10/decade（10→50 pg/ml 区间，
                                    示例轨迹 8→15 脉冲/24h；200 pg/ml 饱和）——计数主信道
    NF-κB   (Tay 2010, 原文句+SI)  α=0.151, ν≈+0.6/decade（1.7→4.0 次 / 0.01→100 ng/ml）
  模型侧（方/圆/三角标记）：
    NF 模型近起始 (8.76, 0.02)；NF 模型远离起始 (0.7, 0.05)
    NF-κB Krishna 软模式 (0.714, 0.05)——与 Tay 实验点间箭头=模型-实验落差（待决）
    FHN 慢驱动 (0.014, 4.05)；Morris–Lecar 快驱动 (0.26, 2.0)
    PdPC 正反馈开关（bistable 静态决策元件，ν=0）置于右缘

运行：python3 代码11_判别图_v09.py
"""
import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(9.5, 7.5))

# ---------- 区域 ----------
ax.axvspan(0.3, 15, alpha=0.08, color='red')      # NF 模拟区
ax.axvspan(5e-4, 0.15, alpha=0.08, color='blue')  # EXC 数字区
ax.axvspan(0.15, 0.3, alpha=0.12, color='green')  # 漏数字带
ax.text(2.4, 9.6, 'NF 族（模拟）\n$\\alpha\\gtrsim0.5$\n（亚临界：+迟滞）',
        color='darkred', fontsize=10, ha='center')
ax.text(0.006, 9.6, 'EXC 族（数字）\n$\\alpha\\lesssim0.15$', color='darkblue', fontsize=10, ha='center')
ax.text(0.212, 9.6, '漏数字带', color='darkgreen', fontsize=8.5, ha='center')
ax.axvline(0.22, color='gray', ls='--', lw=1)

# ---------- 模型点 ----------
ax.plot([8.76], [0.02], 's', color='red', ms=9, label='模型：NF 近起始（α→∞）')
ax.plot([0.7], [0.05], 's', color='red', ms=9, mfc='none', label='模型：NF 远离起始')
ax.plot([0.714], [0.05], 's', color='darkred', ms=10, mfc='none', mew=2,
        label='模型：NF-κB Krishna 软模式（α=0.714）')
ax.plot([0.014], [4.05], 'o', color='blue', ms=9, label='模型：FHN 慢驱动（纯数字）')
ax.plot([0.26], [2.0], '^', color='green', ms=9, label='模型：Morris–Lecar 快驱动（泄漏）')
ax.plot([11], [0.02], 'p', color='darkred', ms=10, label='模型：PdPC 正反馈开关（bistable）')

# ---------- 实验点（v0.9 公开数据审计） ----------
ax.errorbar([0.03], [2.0], xerr=[0.02], yerr=[0.5], fmt='*', color='purple',
            ms=17, capsize=4, label='实验：p53/γ-IR（Lahav 2004，报告统计量）')
ax.errorbar([0.05], [1.8], xerr=[0.04], yerr=[0.5], fmt='*', color='magenta',
            ms=15, capsize=4, label='实验：p53/NCS（Mönke 2017，报告统计量）')
ax.errorbar([0.008], [2.2], yerr=[0.5], fmt='*', color='indigo', ms=19, capsize=4,
            label='实验：p53/NCS（Batchelor 2011 源图，α=0.00±0.28）')
ax.errorbar([0.99], [0.05], xerr=[0.32], yerr=[0.05], fmt='D', color='darkorange',
            ms=10, capsize=4, label='实验：p53/UV（Batchelor 2011 源图，α=+0.99±0.32）')
ax.errorbar([0.02], [10.0], yerr=[2.5], fmt='*', color='royalblue', ms=19, capsize=4,
            label='实验：ERK（Albeck 2013 源图，|α|≲0.05，ν≈+10/dec 后饱和）')
ax.errorbar([0.151], [0.6], xerr=[0.05], yerr=[0.3], fmt='*', color='teal',
            ms=17, capsize=4, label='实验：NF-κB（Tay 2010 原文句+SI，α=0.151，ν≈+0.6/dec）')

# ---------- 模型-实验落差箭头 ----------
ax.annotate('', xy=(0.20, 0.55), xytext=(0.66, 0.08),
            arrowprops=dict(arrowstyle='->', color='darkred', lw=1.8, ls='-'))
ax.text(0.42, 0.9, 'NF-κB 模型-实验落差\n（结构性：参数不可达——代码12）', color='darkred', fontsize=8.5, ha='center')

ax.annotate('内建对照：同蛋白同细胞系\nα 从 0 到 1 严格共变', (0.35, 1.5), fontsize=8.5, color='darkorange')
ax.annotate('v0.8 预言 α≳0.5 被证伪：\nERK 实为定幅+频率调制（数字侧）', (0.06, 7.6), fontsize=8.5, color='royalblue')

ax.set_xscale('log'); ax.set_xlim(5e-4, 15); ax.set_ylim(-0.4, 11.5)
ax.set_xlabel('$\\alpha = \\partial\\ln A/\\partial\\ln D$（振幅-剂量敏感度）')
ax.set_ylabel('$\\nu$（每十倍剂量的脉冲数增量）')
ax.set_title('判别图 v0.9：模型侧落点 + 实验侧落点（公开数据源图/原文定量）', fontsize=11)
ax.legend(fontsize=7.8, loc='center right', framealpha=0.95)
plt.tight_layout()
plt.savefig('判别图_v09_实验落点.png', dpi=200)
print('已保存 判别图_v09_实验落点.png')
