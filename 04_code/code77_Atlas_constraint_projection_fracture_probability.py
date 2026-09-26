# -*- coding: utf-8 -*-
"""
代码77_Atlas约束投影与断裂概率_BpC.py
版本 v1.0.0 · 2026-08-16

「从审计师到建筑师」第一步（按数据现实修正后的方案）：
GPCR Atlas 无原始剂量-反应曲线，DeepSeek 原案的"约束嵌入拟合"无原料；
修正为两个不需要原始数据的建设性算子：

  B′｜约束投影（最小代价修复）：
     组内恒等式要求 q_i ≡ logRA_i − tc_i 全组相等（= −log E_sys 实现值）。
     投影 = 把每条目的 (tc, logRA) 最小平方移动到恒等式流形 q_i = q̄ 上：
       等分修复下每参数移动 (q_i−q̄)/2 dex；组总代价 = Σ(q_i−q̄)²/2。
     输出：每组修复代价、每参数移动量、最大缺口配体（"为一致性买单最多的是谁"）。
     立场声明：投影输出是"框架为真假设下最省事的解释"，不是"修复后的真值"。

  C｜断裂概率混合模型：
     dev_entry 两组分零均值正态混合（健康 N(0,σ_h) + 断裂 N(0,σ_b)），EM 估计，
     无阈值输出每条目 P(断裂)。σ_h 同时是"健康文献噪声"的实测刻度。

  交叉验证：组级修复代价 vs 组内平均 P(断裂) 的 Spearman ρ。

输入：atlas_audit_entries.csv（代码71 产物）
输出：atlas_repair_projection.csv（组级）、atlas_fracture_probability.csv（条目级）

判决（2026-08-16 运行）：
  B′：269 组中 24.2% 需要 >0.3 dex 级修复——与代码71（26%）、代码73（28.2% 保守口径）
      三方法收敛。
  C：EM 收敛健康组占比 74.5%（σ_h=0.089 dex, σ_b=0.68 dex）——无阈值独立复证
      "约 3/4 健康"结论；P>0.9 条目占 15.0% ≈ |dev|>0.3 的 14.5%。
      σ_h=0.089 dex 实测了健康文献的报告噪声刻度（大于代码72 M2 的合法机制散布
      0.004–0.05，差额即报告/数字化噪声）。
  交叉：ρ=0.962，两算子高度互证。
"""
import numpy as np, pandas as pd

ENT = '/mnt/agents/output/04_细胞线4/结果/atlas_audit/atlas_audit_entries.csv'
OUT = '/mnt/agents/output/04_细胞线4/结果/atlas_audit'

def main():
    E = pd.read_csv(ENT).dropna(subset=['dev_entry'])
    E['grp'] = (E['doi'].astype(str)+'|'+E['receptor'].astype(str)+'|'+
                E['Measured process'].astype(str)+'|'+E['Pathway level'].astype(str)+'|'+
                E['Cell line'].astype(str)+'|'+E['Primary effector subtype'].astype(str))
    E['q'] = E['logRA'] - E['tc']

    # ---- B′ 约束投影 ----
    recs = []
    for gname, d in E.groupby('grp'):
        q = d['q'].to_numpy(); n = len(q)
        if n < 2: continue
        qstar = q.mean(); move = q - qstar
        i_w = np.argmax(np.abs(move))
        recs.append({'grp':gname,'n':n,'cost_total_dex2':(move**2).sum()/2,
                     'cost_rms_dex':np.sqrt(((move/2)**2).mean()),
                     'worst_ligand':d['ligand'].iloc[i_w],'worst_move_dex':move[i_w]})
    BJ = pd.DataFrame(recs)
    E['repair_move'] = E['q'] - E['grp'].map(E.groupby('grp')['q'].mean())
    E['repair_per_param'] = E['repair_move']/2

    # ---- C 混合模型 EM ----
    x = E['dev_entry'].to_numpy()
    pi, s_h, s_b = 0.75, 0.05, 0.8
    for it in range(2000):
        lh = pi*np.exp(-x**2/(2*s_h**2))/s_h
        lb = (1-pi)*np.exp(-x**2/(2*s_b**2))/s_b
        pb = lb/(lh+lb)
        pi_n = 1-pb.mean()
        sh_n = np.sqrt(((1-pb)*x**2).sum()/(1-pb).sum())
        sb_n = np.sqrt((pb*x**2).sum()/pb.sum())
        if max(abs(pi_n-pi), abs(sh_n-s_h), abs(sb_n-s_b)) < 1e-10:
            pi,s_h,s_b = pi_n,sh_n,sb_n; break
        pi,s_h,s_b = pi_n,sh_n,sb_n
    E['p_broken'] = pb

    BJ.to_csv(f'{OUT}/atlas_repair_projection.csv', index=False)
    E.to_csv(f'{OUT}/atlas_fracture_probability.csv', index=False)

    from scipy.stats import spearmanr
    chk = E.groupby('grp').agg(pmean=('p_broken','mean')).reset_index().merge(BJ, on='grp')
    rho = spearmanr(chk.cost_rms_dex, chk.pmean).statistic

    print(f"[B′] 组数 {len(BJ)}；修复代价 RMS 中位 {BJ['cost_rms_dex'].median():.3f} dex；"
          f"需 >0.3dex 级修复的组 {(BJ['cost_rms_dex']>0.15).mean():.1%}")
    print(f"[C] EM 收敛 @{it} 步：π_健康={pi:.3f}，σ_h={s_h:.4f}，σ_b={s_b:.3f} dex；"
          f"P>0.9 条目 {(pb>0.9).mean():.1%}（对照 |dev|>0.3：{(np.abs(x)>0.3).mean():.1%}）")
    print(f"[交叉] 修复代价 vs 断裂概率 Spearman ρ = {rho:.3f}")
    print("[B′ 代价最高 5 组]")
    print(BJ.nlargest(5,'cost_rms_dex')[['grp','n','cost_rms_dex','worst_ligand','worst_move_dex']].to_string(index=False))

if __name__ == '__main__':
    main()
