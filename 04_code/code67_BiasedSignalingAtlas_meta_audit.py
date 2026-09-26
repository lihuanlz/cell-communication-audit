# -*- coding: utf-8 -*-
"""
代码67_BiasedSignalingAtlas元审计.py  v1.0.0
============================================================
Cell line 4 · GPCR biased-signaling fracture audit · meta-audit: the full Biased Signaling Atlas

Data: protwis/gpcrdb_data GitHub repository
  ligand_data/bias_data/Biased_ligand_single_pathway_data.xlsx (11.7 MB, 9041 rows,
  ligand-receptor-pathway activity records from 212 publications; archived copy of this volume at
  04_细胞线4/公开数据/BiasedSignalingAtlas/)

Design (three tests, all fit-free):
  Test A (population identity): under Black–Leff nH=1, within the same (publication × receptor × assay × cell line) group
    any two ligands satisfy  Δlog(τ/KA) = Δlog(Emax%/EC50)  — a strict identity (derivation:
    Emax%∝τ/(1+τ), EC50=KA/(1+τ), ratio∝τ/KA; the proportionality constant cancels within a group).
    Reconcile pairwise differences of the Atlas's own log(τ/KA) against Emax and pEC50.
    Deviation≠error (legitimate fits with nH≠1, or literature reporting τ/KA directly, both deviate), but the distribution tail
    and extreme groups deserve flagging; this is a "population-level inconsistency map".
  Test B (Emax QC): number of rows with Emax%>150 or negative values (impossible when referenced to a full agonist,
    legitimate when referenced to a partial agonist — only flag the distribution, do not convict row by row).
  Test C (transitivity closure): triangular closure of Δlog(τ/KA) across pathways for the same publication and ligand,
    bias(p1,p2)+bias(p2,p3) ?= bias(p1,p3), which should be strictly =0 (pure arithmetic test).
  Cross-check: locate the Gillis 2020 (aaz3140) entries, examine their behavior in test A,
    and verify which GIRK data column the Atlas transcribed.
============================================================
"""
import numpy as np, pandas as pd, itertools, os, sys

PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                    "公开数据", "BiasedSignalingAtlas", "Biased_ligand_single_pathway_data.xlsx")
if not os.path.exists(PATH):  # fallback: /tmp
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

print(f"total rows {len(df)}; with log(τ/KA) {tc.notna().sum()}; with Δlog(τ/KA) {rtc.notna().sum()};"
      f"with Emax {em.notna().sum()}; valid pEC50/EC50 {np.isfinite(pEC).sum()}")

# ---------- Test A ----------
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
print("\n=== Test A: |Δlog(τ/KA) − Δlog(Emax/EC50)| (identity should ≈0) ===")
print(f"pairwise comparisons in multi-ligand groups: {len(pairs)}")
for q in [50, 75, 90, 95, 99]:
    print(f"  P{q}: {np.percentile(pairs['dev'], q):.2f} dex")
print(f"  fraction >0.3 dex: {(pairs['dev']>0.3).mean():.1%}; >1 dex: {(pairs['dev']>1).mean():.1%}; max {pairs['dev'].max():.2f}")
per = pairs.groupby("doi")["dev"].median()
print(f"  publications involved {len(per)}; publications with median deviation >0.3 dex {(per>0.3).sum()} ({(per>0.3).mean():.0%})")
print("  5 publications with the largest median deviation:")
for doi, v in per.sort_values(ascending=False).head(5).items():
    print(f"    {doi}: {v:.2f} dex")

# ---------- Test B ----------
print("\n=== Test B: Emax QC ===")
print(f"  Emax>110%: {(df['em']>110).sum()} rows; >150%: {(df['em']>150).sum()} rows; >300%: {(df['em']>300).sum()} rows; negative: {(df['em']<0).sum()} rows")
print("  (>100% is legitimate when referenced to a partial agonist — flag only, no conviction)")

# ---------- Test C ----------
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
print("\n=== Test C: Δlog(τ/KA) cross-pathway triangular closure ===")
print(f"  triangles {len(viol)}; median {np.median(viol):.4f}; fraction >0.01 {(viol>0.01).mean():.3f}; max {viol.max():.4f}")

# ---------- cross-check Gillis 2020 ----------
print("\n=== Cross-check: transcription of Gillis 2020 (scisignal.aaz3140) in the Atlas ===")
gil = df[df["Reference DOI or PMID"] == "10.1126/scisignal.aaz3140"]
print(f"  Atlas contains {len(gil)} rows / {gil['Ligand tested for bias or func. Sel. Name'].nunique()} ligands")
girk = gil[gil["Measured process"] == "Activation"][["Ligand tested for bias or func. Sel. Name", "tc", "em", "pEC50u"]]
print("  GIRK (Activation) rows — note Emax is the β-CNA-treated column (DAMGO=75); the untreated column was not transcribed:")
print(girk.to_string(index=False))
gilp = pairs[pairs["doi"] == "10.1126/scisignal.aaz3140"]
print(f"  Gillis group in test A: {len(gilp)} pairs, median deviation {gilp['dev'].median():.2f} dex, "
      f" fraction >0.3 {(gilp['dev']>0.3).mean():.0%}, max {gilp['dev'].max():.2f}")
