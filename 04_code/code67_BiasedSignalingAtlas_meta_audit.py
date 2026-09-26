# -*- coding: utf-8 -*-
"""
代码67_BiasedSignalingAtlas元审计.py  v1.0.0
============================================================
细胞线4 · GPCR偏向性断裂审计 · 元审计：Biased Signaling Atlas 全库

数据：protwis/gpcrdb_data GitHub 仓库
  ligand_data/bias_data/Biased_ligand_single_pathway_data.xlsx（11.7 MB，9041 行，
  212 篇文献的 ligand-receptor-pathway 活动记录；本卷存档副本于
  04_细胞线4/公开数据/BiasedSignalingAtlas/）

设计（三条检验，全部无需拟合）：
  检验A（群体恒等式）：Black–Leff nH=1 下，同一（文献×受体×测定×细胞系）组内
    任意两个配体满足  Δlog(τ/KA) = Δlog(Emax%/EC50)  —— 严格恒等（推导：
    Emax%∝τ/(1+τ)、EC50=KA/(1+τ)，比值∝τ/KA，组内比例常数相消）。
    用 Atlas 自存的 log(τ/KA) 与 Emax、pEC50 两两做差对账。
    偏差≠错误（nH≠1 的合法拟合、文献直接报 τ/KA 都会偏），但分布尾部
    与极端组值得标记；这是"群体层面不一致地图"。
  检验B（Emax 质控）：Emax%>150 或负值的行数（以完全激动剂为参考时不可能，
    以部分激动剂为参考时合法——只标记分布，不逐行定罪）。
  检验C（传递性闭合）：Δlog(τ/KA) 同文献同配体跨通路的三角闭合，
    bias(p1,p2)+bias(p2,p3) ?= bias(p1,p3)，应严格=0（纯算术检验）。
  交叉核对：定位 Gillis 2020（aaz3140）条目，检验其在检验A中的表现，
    并核查 Atlas 转录的是哪一列 GIRK 数据。
============================================================
"""
import numpy as np, pandas as pd, itertools, os, sys

PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                    "公开数据", "BiasedSignalingAtlas", "Biased_ligand_single_pathway_data.xlsx")
if not os.path.exists(PATH):  # 兜底：/tmp
    PATH = "/tmp/bias/atlas.xlsx"
df = pd.read_excel(PATH, sheet_name="Data")
df.columns = [c.replace("\n", " ") for c in df.columns]

tc  = pd.to_numeric(df["Transduction Coefficient [log(τ/KA)]"], errors="coerce")
rtc = pd.to_numeric(df["Relative Transduction Coefficient [Δlog(τ/KA)]"], errors="coerce")
em  = pd.to_numeric(df["Alt 1) Quantitative efficacy"], errors="coerce")
act = pd.to_numeric(df["Alt 1) Quantitative activity"], errors="coerce")
mt, un = df["Measure type"], df["Unit"]
pEC = np.where(mt == "pEC50", act,
        np.where(mt == "EC50",
                 -np.log10(act*np.where(un=="nM",1e-9,np.where(un=="µM",1e-6,
                              np.where(un=="M",1.0,np.where(un=="pM",1e-12,np.nan))))),
                 np.nan))
df["pEC50u"], df["tc"], df["em"], df["rtc"] = pEC, tc, em, rtc

print(f"总行数 {len(df)}；有 log(τ/KA) {tc.notna().sum()}；有 Δlog(τ/KA) {rtc.notna().sum()}；"
      f"有 Emax {em.notna().sum()}；有效 pEC50/EC50 {np.isfinite(pEC).sum()}")

# ---------- 检验A ----------
sub = df[tc.notna() & em.notna() & np.isfinite(pEC) & (em > 0)].copy()
sub["logRA"] = np.log10(sub["em"]) + sub["pEC50u"]
sub["grp"] = (sub["Reference DOI or PMID"].astype(str) + "|" +
              sub["Receptor UniProt entry name or code"].astype(str) + "|" +
              sub["Measured process"].astype(str) + "|" + sub["Pathway level"].astype(str) + "|" +
              sub["Cell line"].astype(str) + "|" + sub["Primary effector subtype"].astype(str))
pairs = []
for k, v in sub.groupby("grp"):
    d = v.drop_duplicates("Ligand tested for bias or func. Sel. Name")
    if len(d) < 2: continue
    for i, j in itertools.combinations(d.index, 2):
        dtc = d.loc[i, "tc"] - d.loc[j, "tc"]
        dra = d.loc[i, "logRA"] - d.loc[j, "logRA"]
        pairs.append((d.loc[i, "Reference DOI or PMID"], abs(dtc - dra)))
pairs = pd.DataFrame(pairs, columns=["doi", "dev"])
print("\n=== 检验A：|Δlog(τ/KA) − Δlog(Emax/EC50)|（恒等式应≈0）===")
print(f"多配体组两两对数: {len(pairs)}")
for q in [50, 75, 90, 95, 99]:
    print(f"  P{q}: {np.percentile(pairs['dev'], q):.2f} dex")
print(f"  >0.3 dex 占比: {(pairs['dev']>0.3).mean():.1%}；>1 dex: {(pairs['dev']>1).mean():.1%}；max {pairs['dev'].max():.2f}")
per = pairs.groupby("doi")["dev"].median()
print(f"  涉及文献 {len(per)} 篇；中位偏差>0.3 dex 的文献 {(per>0.3).sum()} 篇（{(per>0.3).mean():.0%}）")
print("  中位偏差最大的5篇:")
for doi, v in per.sort_values(ascending=False).head(5).items():
    print(f"    {doi}: {v:.2f} dex")

# ---------- 检验B ----------
print("\n=== 检验B：Emax 质控 ===")
print(f"  Emax>110%: {(df['em']>110).sum()} 行；>150%: {(df['em']>150).sum()} 行；>300%: {(df['em']>300).sum()} 行；负值: {(df['em']<0).sum()} 行")
print("  （以部分激动剂为参考时>100%合法——只标记，不定罪）")

# ---------- 检验C ----------
sub2 = df[rtc.notna()].copy()
sub2["pwy"] = (sub2["Pathway level"].astype(str) + "|" + sub2["Measured process"].astype(str)
               + "|" + sub2["Primary effector subtype"].astype(str))
sub2["grp2"] = (sub2["Reference DOI or PMID"].astype(str) + "|" +
                sub2["Receptor UniProt entry name or code"].astype(str) + "|" +
                sub2["Ligand tested for bias or func. Sel. Name"].astype(str) + "|" +
                sub2["Cell line"].astype(str))
viol = []
for k, v in sub2.groupby("grp2"):
    d = v.drop_duplicates("pwy")
    if len(d) < 3: continue
    rs = d["rtc"].tolist()
    for i, j, kk in itertools.combinations(range(len(rs)), 3):
        viol.append(abs((rs[i]-rs[j]) + (rs[j]-rs[kk]) - (rs[i]-rs[kk])))
viol = np.array(viol)
print("\n=== 检验C：Δlog(τ/KA) 跨通路三角闭合 ===")
print(f"  三角数 {len(viol)}；中位 {np.median(viol):.4f}；>0.01 占比 {(viol>0.01).mean():.3f}；max {viol.max():.4f}")

# ---------- 交叉核对 Gillis 2020 ----------
print("\n=== 交叉核对：Gillis 2020（scisignal.aaz3140）在 Atlas 中的转录 ===")
gil = df[df["Reference DOI or PMID"] == "10.1126/scisignal.aaz3140"]
print(f"  Atlas 收录 {len(gil)} 行 / {gil['Ligand tested for bias or func. Sel. Name'].nunique()} 配体")
girk = gil[gil["Measured process"] == "Activation"][["Ligand tested for bias or func. Sel. Name", "tc", "em", "pEC50u"]]
print("  GIRK（Activation）行——注意 Emax 为 β-CNA 处理列（DAMGO=75），未处理列未收录:")
print(girk.to_string(index=False))
gilp = pairs[pairs["doi"] == "10.1126/scisignal.aaz3140"]
print(f"  Gillis 组在检验A中: {len(gilp)} 对, 中位偏差 {gilp['dev'].median():.2f} dex, "
      f">0.3 占比 {(gilp['dev']>0.3).mean():.0%}, max {gilp['dev'].max():.2f}")
