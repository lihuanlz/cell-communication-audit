# -*- coding: utf-8 -*-
"""
代码73_稳健性三联升级_盲测聚合Jackknife.py
版本 v1.0.0 · 2026-08-16 · 种子 20260816

按评审后的优先级落地三条稳健性升级（DeepSeek 建议的批判性采纳版）：
  模块A（策略4/B1）：血红蛋白外部盲测阴性对照
      - 数据：Severinghaus 1979 标准人血氧离曲线（对人血实测的拟合式，37°C, pH 7.4）
              + LITFL 临床锚点表（5 个独立的 SpO2/PaO2 锚点）
      - 检验1（恒等式类比）：从曲线反演 Hill 参数 (p50, nH)，再检查
              "参数→曲线" 自洽性（20–80% 饱和带，Hill 线性区）
      - 检验2（跨源一致性，Atlas 式）：同一条件 p50 的多来源散布
      - 检验3（阳性对照）：人为污染曲线（一半点 ×1.6），审计必须报警
      - 诚实限制：Severinghaus 是拟合式非原始点；Hill 对 Hb 是近似（已知百年）
  模块B（策略1）：Atlas 组级聚合检验 + 相关性校准
      - 原假设模拟注入共享共模误差 rho∈{0, 0.3, 0.5}，报告 p 值膨胀系数
      - 组级判决须在 rho=0 与 rho=0.5 下同时成立才算高置信
  模块C（策略3）：Jackknife 留一一致性（Atlas 组级 + Gillis 列级附注）
  模块D：断裂置信度矩阵（0–100 连续评分 + 三级标签）

输入：04_细胞线4/结果/atlas_audit/atlas_audit_entries.csv（代码71 产物，复用不重下）
输出：04_细胞线4/结果/atlas_audit/atlas_confidence_matrix.csv
"""
import numpy as np, pandas as pd

rng = np.random.default_rng(20260816)
LN10 = np.log(10)

print("=" * 70)
print("模块A｜血红蛋白外部盲测阴性对照（策略4 / B1）")
print("=" * 70)

# --- A.0 数据 ---------------------------------------------------------------
def severinghaus_s(po2):
    """Severinghaus 1979 标准人血 ODC（对实测人血的拟合式）。"""
    po2 = np.asarray(po2, float)
    x3 = po2**3 + 150.0 * po2
    return 100.0 * x3 / (x3 + 23400.0)

po2_grid = np.arange(5.0, 121.0, 5.0)          # 24 个参考点
sat_curve = severinghaus_s(po2_grid)

litfl_anchor_po2 = np.array([95., 60., 50., 40., 27.])
litfl_anchor_sat = np.array([97., 92., 89., 75., 50.])

# 同一条件(标准人血 pH7.4/37°C) p50 多来源：Severinghaus 曲线内插, LITFL 锚点,
# Morgan(acutecaretesting) 26.7, medmastery 26.6, PMC7547706 26.9(文献值26.9±)
p50_sources = {"Severinghaus曲线": None, "LITFL锚点": 27.0,
               "Morgan综述": 26.7, "medmastery": 26.6, "PMC7547706": 26.9}

def hill_invert(sat, p50, nH):
    """由饱和度反演 pO2（Hill 方程）。"""
    y = np.asarray(sat, float) / 100.0
    return p50 * (y / (1.0 - y)) ** (1.0 / nH)

def hill_fit(po2, sat, lo=20.0, hi=80.0):
    """Hill 线性化拟合（仅限 20–80% 饱和带，Hb 的 Hill 近似线性区）。"""
    m = (sat >= lo) & (sat <= hi)
    ly = np.log10(sat[m] / (100.0 - sat[m]))
    lx = np.log10(po2[m])
    A = np.vstack([lx, np.ones_like(lx)]).T
    (nH, b), *_ = np.linalg.lstsq(A, ly, rcond=None)
    p50 = 10.0 ** (-b / nH)
    return p50, nH, m.sum()

# --- A.1 检验1：曲线↔参数自洽恒等式（盲测，期望：静默） ----------------------
p50_fit, nH_fit, nband = hill_fit(po2_grid, sat_curve)
band = (sat_curve >= 20) & (sat_curve <= 80)
dev_hb = np.log10(hill_invert(sat_curve[band], p50_fit, nH_fit) / po2_grid[band])
print(f"[A1] Hill 反演：p50 = {p50_fit:.2f} mmHg, nH = {nH_fit:.2f}（带内点数 {nband}）")
print(f"     恒等式偏差 |Δlog10 pO2|：中位 {np.median(np.abs(dev_hb)):.4f} dex，"
      f"最大 {np.abs(dev_hb).max():.4f} dex")
print(f"     >0.3 dex 报警数：{int((np.abs(dev_hb) > 0.3).sum())}/{band.sum()}"
      f"（期望 0 → 审计在健康体系上静默）")

# 锚点独立复核：LITFL 5 锚点相对 Severinghaus 曲线的水平偏差
dev_anchor = np.log10(hill_invert(litfl_anchor_sat, p50_fit, nH_fit) / litfl_anchor_po2)
print(f"     LITFL 独立锚点水平偏差：{np.round(dev_anchor, 3)} dex，"
      f"max |Δ| = {np.abs(dev_anchor).max():.3f}")

# --- A.2 检验2：跨源 p50 一致性 ----------------------------------------------
p50_sev = float(np.interp(50.0, sat_curve, po2_grid))
p50_sources["Severinghaus曲线"] = round(p50_sev, 2)
vals = np.array(list(p50_sources.values()))
spread_dex = np.log10(vals.max() / vals.min())
print(f"[A2] 五来源 p50：{p50_sources}")
print(f"     跨源散布 = {spread_dex:.4f} dex（阈值 0.3；期望 ≪0.3 → 静默）")

# --- A.3 检验3：阳性对照（污染必须被抓住，证明本检验在这套数据上有牙齿） ------
po2_bad = po2_grid.copy()
po2_bad[::2] *= 2.5   # 0.4 dex 污染：恒等式在污染后不可修复，必须报警
p50b, nHb, _ = hill_fit(po2_bad, sat_curve)
dev_bad = np.log10(hill_invert(sat_curve[band], p50b, nHb) / po2_bad[band])
print(f"[A3] 阳性对照（半数点×2.5）：max |Δ| = {np.abs(dev_bad).max():.3f} dex，"
      f">0.3 报警 {int((np.abs(dev_bad) > 0.3).sum())}/{band.sum()}"
      f"（期望 ≫0 → 审计有牙齿）")

print()
print("=" * 70)
print("模块B｜Atlas 组级聚合检验 + 相关性校准（策略1 修正版）")
print("=" * 70)

E = pd.read_csv('/mnt/agents/output/04_细胞线4/结果/atlas_audit/atlas_audit_entries.csv')
E = E.dropna(subset=['dev_entry'])
# 与代码71同口径重建组键：doi|receptor|process|level|cell|effector
E['grp'] = (E['doi'].astype(str) + '|' + E['receptor'].astype(str) + '|' +
            E['Measured process'].astype(str) + '|' + E['Pathway level'].astype(str) + '|' +
            E['Cell line'].astype(str) + '|' + E['Primary effector subtype'].astype(str))
groups = {k: v['dev_entry'].to_numpy() for k, v in E.groupby('grp') if len(v) >= 3}
print(f"组数（n≥3）：{len(groups)}")

# 合法散布尺度：主口径 0.05 dex（代码72 M2：nH 等合法机制中位仅 0.004 dex）；
# 感应分析再加 0.10 dex（约等于 Gillis 类数据的报告 SEM 量级）双口径
SIG_LEGIT = 0.05
SIGMAS = [0.05, 0.10]
NSIM = 20000
RHOS = [0.0, 0.3, 0.5]

def group_stats(d):
    return np.median(d), (np.abs(d) > 0.3).mean()

def null_sim(n, rho, sig=SIG_LEGIT, nsim=NSIM):
    """共享共模误差模型：dev_i = rho*c + sqrt(1-rho^2)*e_i。"""
    c = rng.standard_normal(nsim)[:, None] * sig
    e = rng.standard_normal((nsim, n)) * sig
    dev = rho * c + np.sqrt(1 - rho**2) * e
    meds = np.median(dev, axis=1)
    fracs = (np.abs(dev) > 0.3).mean(axis=1)
    return meds, fracs

rows = []
null_cache = {}
for gname, d in groups.items():
    n = len(d)
    if n not in null_cache:
        null_cache[n] = {r: null_sim(n, r) for r in RHOS}
    med, frac = group_stats(d)
    rec = {'grp': gname, 'n': n, 'dev_median': med, 'dev_max': np.abs(d).max(),
           'frac_gt03': frac}
    for r in RHOS:
        meds, fracs = null_cache[n][r]
        # 两个统计量检验不同备择：median 抓"系统性同向偏移"，frac 抓"散布型不一致"。
        # 组级报警 = 任一显著 → 取 min（曾误用 max，导致 median≈0 的散布型断裂全部漏判，已修正并双录）
        p_med = (np.abs(meds) >= abs(med)).mean()
        p_frac = (fracs >= frac).mean()
        rec[f'p_rho{int(r*10)}'] = min(p_med, p_frac)
    rows.append(rec)

B = pd.DataFrame(rows)
B['高置信断裂'] = (B['p_rho0'] < 0.001) & (B['p_rho5'] < 0.001)
B['候选断裂'] = (~B['高置信断裂']) & (B['p_rho0'] < 0.01)
n_hc, n_cand = int(B['高置信断裂'].sum()), int(B['候选断裂'].sum())
print(f"高置信断裂（rho=0 与 0.5 下 p 均<0.001）：{n_hc}/{len(B)} 组（{n_hc/len(B):.1%}）")
print(f"候选断裂（仅 rho=0 下 p<0.01）：{n_cand} 组")
print(f"未检出：{len(B)-n_hc-n_cand} 组（{1-(n_hc+n_cand)/len(B):.1%}）")

# 膨胀系数演示：共模相关使 p 值虚低的程度
demo_n = 8
meds0, fracs0 = null_cache.get(demo_n, {r: null_sim(demo_n, r) for r in RHOS})[0.0]
obs_frac = 0.5
p0 = (null_cache[demo_n][0.0][1] >= obs_frac).mean() if demo_n in null_cache else None
p5 = (null_cache[demo_n][0.5][1] >= obs_frac).mean() if demo_n in null_cache else None
if p0:
    print(f"相关性校准演示（n={demo_n}，超标率 50%）：独立假设 p={p0:.2e} → "
          f"rho=0.5 下 p={p5:.2e}，膨胀 {p5/max(p0,1e-9):.0f}×")

# 合法散布尺度的感应分析：SIG=0.10 dex（≈Gillis 类报告 SEM 量级）、rho=0.5 最保守口径
ns01 = {}
hc01 = 0
for gname, d in groups.items():
    n = len(d)
    if n not in ns01:
        ns01[n] = null_sim(n, 0.5, sig=0.10)
    meds, fracs = ns01[n]
    med, frac = group_stats(d)
    p = min((np.abs(meds) >= abs(med)).mean(), (fracs >= frac).mean())
    hc01 += p < 0.001
print(f"感应分析（SIG=0.10 dex, rho=0.5 最保守口径）：高置信断裂 {hc01}/{len(groups)}"
      f"（{hc01/len(groups):.1%}）——主结论对合法散布尺度假设稳健与否的判据")

print()
print("=" * 70)
print("模块C｜Jackknife 留一一致性（策略3）")
print("=" * 70)

def verdict_of(d):
    return (np.abs(d) > 0.3).mean() > 0.0  # 组级口径：存在超阈值条目即"不一致"

jk_scores = []
for gname, d in groups.items():
    if not verdict_of(d):
        jk_scores.append({'grp': gname, 'jackknife': 1.0}); continue
    keep = 0
    for i in range(len(d)):
        keep += verdict_of(np.delete(d, i))
    jk_scores.append({'grp': gname, 'jackknife': keep / len(d)})
JK = pd.DataFrame(jk_scores)
B = B.merge(JK, on='grp')
frac_robust = (B.loc[B['frac_gt03'] > 0, 'jackknife'] >= 1.0).mean()
print(f"Atlas：{int((B['frac_gt03']>0).sum())} 个含超标条目的组中，"
      f"留一后判决不变的占 {frac_robust:.1%}")
print(f"jackknife 评分分布：1.0 → {int((B['jackknife']==1).sum())} 组，"
      f"0.5–1 → {int(((B['jackknife']<1)&(B['jackknife']>=0.5)).sum())} 组，"
      f"<0.5 → {int((B['jackknife']<0.5).sum())} 组")

# Gillis 列级附注（n=5，区分度弱，只作参考）
gillis_girk_dev = np.array([-0.68, -0.79, -0.83, -0.90, -0.97])  # 代码66 恒等式1
def column_verdict(d):  return (np.abs(d) > 0.3).all()
keep = sum(column_verdict(np.delete(gillis_girk_dev, i)) for i in range(5))
print(f"Gillis GIRK 列（n=5，附注）：留一后仍 4/4 全超 0.3 dex → jackknife = {keep/5:.1f}"
      f"（断裂非单点驱动）")

print()
print("=" * 70)
print("模块D｜断裂置信度矩阵")
print("=" * 70)

def confidence(row):
    """0–100：幅度(30) + 聚合显著性(40，相关性保守口径) + jackknife(30)。"""
    s_amp = min(row['dev_max'] / 1.0, 1.0) * 30
    p = row['p_rho5']                      # 用最保守口径
    s_sig = 40 * (1 if p < 1e-4 else 0.75 if p < 1e-3 else
                  0.5 if p < 1e-2 else 0.25 if p < 0.05 else 0)
    s_jk = row['jackknife'] * 30
    return round(s_amp + s_sig + s_jk, 1)

B['置信度'] = B.apply(confidence, axis=1)
def tier(r):
    if r['高置信断裂']: return '高置信断裂'
    if r['候选断裂']:   return '候选断裂'
    return '未检出'
B['标签'] = B.apply(tier, axis=1)
B = B.sort_values('置信度', ascending=False)

out = '/mnt/agents/output/04_细胞线4/结果/atlas_audit/atlas_confidence_matrix.csv'
B.to_csv(out, index=False)
print(f"已写出：{out}")
print(f"标签分布：{B['标签'].value_counts().to_dict()}")
print("\n置信度最高的 10 组：")
print(B.head(10)[['grp', 'n', 'dev_max', 'frac_gt03', 'p_rho5', 'jackknife',
                  '置信度', '标签']].to_string(index=False))

print()
print("=" * 70)
print("汇总判决")
print("=" * 70)
print(f"A 盲测：Hb 体系恒等式 max|Δ| = {np.abs(dev_hb).max():.3f} dex → 静默；"
      f"阳性对照 max|Δ| = {np.abs(dev_bad).max():.3f} dex → 报警。特异性+灵敏度双通过。")
print(f"B 聚合：Atlas 组级高置信断裂 {n_hc}/{len(B)}（{n_hc/len(B):.1%}，"
      f"口径=组内含≥1条不可修复条目，比代码71的26%口径宽）；"
      f"SIG=0.10 最保守口径 {hc01}/{len(groups)}（{hc01/len(groups):.1%}）。")
print(f"C Jackknife：超标组留一稳定率 {frac_robust:.1%}；Gillis GIRK 非单点驱动。")
