# -*- coding: utf-8 -*-
"""图 1 | 审计流程示意（手稿图 1 实现）：一回路一列，五步机器一行。"""
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

fig, ax = plt.subplots(figsize=(13.5, 8))
ax.set_xlim(0, 13.5); ax.set_ylim(0, 8); ax.axis('off')

steps = ['回路', '可观测量', '简并群 G', 'Fisher 分层', '报告协议']
left = ['PdPC 共价循环\n（静态传感器）',
        '磷酸化分数 u(ξ)\nGK 曲线族 κ: $10^{-3}$–1',
        '乘性三维群（定理 1）\n$k_1k_2$ / $K_mS_T$ / 联合重标定',
        '形状 χ≈30（可辨识）\nvs 绝对 ~$10^{17}$–$10^{19}$（壁垒）',
        '深饱和单边报界\n$n_H$ = 1+1/(2κ) 闭式']
right = ['p53 脉冲发生器\n（动态编码器，FHN 可兴奋）',
         '脉冲三元组\n振幅 A、计数 N、间隔 T',
         '振幅标度群\n产生率 × 成熟 × 增益',
         '零方向 = G 轨道切向\n（代码13，余弦 1.000000）',
         '三签名判别（定理 4）\n起始 × 迟滞 × α']

xs, w, h = [2.1, 8.6], 3.6, 1.05
ys = [6.6, 5.25, 3.9, 2.55, 1.2]
for i, s in enumerate(steps):
    ax.text(0.35, ys[i]+h/2, f'步骤{i+1}\n{s}' if i else s, fontsize=10,
            ha='left', va='center', fontweight='bold', color='#333')
for x, col, ec in [(xs[0], left, '#1f77b4'), (xs[1], right, '#d62728')]:
    for y, txt in zip(ys, col):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0.06',
                                    fc='white', ec=ec, lw=1.8))
        ax.text(x+w/2, y+h/2, txt, fontsize=9.5, ha='center', va='center')
    for j in range(len(ys)-1):
        ax.add_patch(FancyArrowPatch((x+w/2, ys[j]-0.02), (x+w/2, ys[j+1]+h+0.02),
                                     arrowstyle='-|>', mutation_scale=16, color='#555'))

# 底部共享数字极限
ax.add_patch(FancyBboxPatch((3.2, 0.05), 7.1, 0.72, boxstyle='round,pad=0.08',
                            fc='#f4f8e8', ec='#6b8e23', lw=2))
ax.text(6.75, 0.41, '共享数字极限：κ→0（化学饱和）≡ α→0（拓扑解耦）→ 似然退化为计数统计（§5 条件定理，SI S2）',
        fontsize=10.5, ha='center', va='center', fontweight='bold', color='#4a6512')
for x in xs:
    ax.add_patch(FancyArrowPatch((x+w/2, ys[-1]-0.02), (6.75, 0.82),
                                 arrowstyle='-|>', mutation_scale=16, color='#6b8e23', ls='--'))

# 列标题
ax.text(xs[0]+w/2, 7.55, '静态路径', fontsize=13, ha='center', fontweight='bold', color='#1f77b4')
ax.text(xs[1]+w/2, 7.55, '动态路径', fontsize=13, ha='center', fontweight='bold', color='#d62728')

plt.tight_layout()
plt.savefig('图1_审计流程示意.png', dpi=200)
print('已保存 图1_审计流程示意.png')
