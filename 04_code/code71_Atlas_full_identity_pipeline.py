# -*- coding: utf-8 -*-
"""
代码71_Atlas全库恒等式流水线.py  v1.0.0
============================================================
Cell line 4 · A1: Biased Signaling Atlas full-database audit pipeline (full-scale version of Code 67)

Screen all 9041 Atlas rows one by one, producing a machine-readable "inconsistency map" data layer:
  (1) Identity A (within-group ligand pairs): Δlog(τ/KA) ?= Δlog(Emax/EC50), strict under Black–Leff n=1
  (2) Emax QC flags (>110/150/300%, negative values)
  (3) Δlog(τ/KA) cross-pathway triangle closure (in-database arithmetic check)
Aggregation dimensions: paper × receptor × pathway × cell line; each record gives deviation median/max,
sample counts, QC flag counts. Outputs:
  atlas_audit_entries.csv   entry level (with logRA and tc differences)
  atlas_audit_groups.csv    group-level aggregation
  atlas_audit_papers.csv    paper-level aggregation
  atlas_audit_summary.json  full-database summary + extreme-group list
Interpretation discipline (as in Code 67): legitimate nH≠1 fits and paper-reported τ/KA both produce deviations;
this pipeline produces an "inconsistency map", not an "error list".
Run: python3 代码71_Atlas全库恒等式流水线.py
============================================================
"""
import numpy as np, pandas as pd, itertools, os, json

BASE = os.path.dirname(os.path.abspath(__file__))
PATH = os.path.join(BASE, "..", "公开数据", "BiasedSignalingAtlas",
                    "Biased_ligand_single_pathway_data.xlsx")
OUTDIR = os.path.join(BASE, "..", "结果", "atlas_audit")
os.makedirs(OUTDIR, exist_ok=True)

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
df["doi"] = df["Reference DOI or PMID"].astype(str)
df["receptor"] = df["Receptor UniProt entry name or code"].astype(str)
df["ligand"] = df["Ligand tested for bias or func. Sel. Name"].astype(str)

# ---------- Entry level ----------
mask = tc.notna() & em.notna() & np.isfinite(pEC) & (em > 0)
ent = df[mask].copy()
ent["logRA"] = np.log10(ent["em"]) + ent["pEC50u"]
ent["grp"] = (ent["doi"] + "|" + ent["receptor"] + "|" + ent["Measured process"].astype(str)
              + "|" + ent["Pathway level"].astype(str) + "|" + ent["Cell line"].astype(str)
              + "|" + ent["Primary effector subtype"].astype(str))
# Within-group reference: group median of logRA − tc as the group constant (≈ realized value of −log E_sys)
gmed = ent.groupby("grp").apply(lambda x: (x["logRA"]-x["tc"]).median(), include_groups=False)
ent["dev_entry"] = (ent["logRA"] - ent["tc"]) - ent["grp"].map(gmed)
ent["doi_"] = ent["doi"]; ent["receptor_"] = ent["receptor"]
ent_out = ent[["doi_","receptor_","ligand","Measured process","Pathway level","Cell line",
               "Primary effector subtype","em","pEC50u","tc","logRA","dev_entry"]].rename(
               columns={"doi_":"doi","receptor_":"receptor"})
ent_out.to_csv(os.path.join(OUTDIR,"atlas_audit_entries.csv"), index=False)

# ---------- Group level ----------
rows=[]
for k, v in ent.groupby("grp"):
    d = v.drop_duplicates("ligand")
    if len(d) < 2: continue
    devs=[]
    for i,j in itertools.combinations(d.index,2):
        devs.append(abs((d.loc[i,"tc"]-d.loc[j,"tc"]) - (d.loc[i,"logRA"]-d.loc[j,"logRA"])))
    devs=np.array(devs)
    parts=k.split("|")
    rows.append({"grp":k,"doi":parts[0],"receptor":parts[1],"process":parts[2],
                 "level":parts[3],"cell":parts[4],"effector":parts[5],
                 "n_ligands":len(d),"n_pairs":len(devs),
                 "dev_median":float(np.median(devs)),"dev_max":float(devs.max())})
grp_df=pd.DataFrame(rows).sort_values("dev_median",ascending=False)
grp_df.to_csv(os.path.join(OUTDIR,"atlas_audit_groups.csv"), index=False)

# ---------- Paper level ----------
pap = grp_df.groupby("doi").agg(
    n_groups=("grp","size"), n_pairs=("n_pairs","sum"),
    dev_median=("dev_median","median"), dev_max=("dev_max","max"),
    receptors=("receptor", lambda s: "|".join(sorted(set(s))[:5]))).reset_index()
pap=pap.sort_values("dev_median",ascending=False)
pap.to_csv(os.path.join(OUTDIR,"atlas_audit_papers.csv"), index=False)

# ---------- Emax QC ----------
qc = {"emax_gt110": int((df["em"]>110).sum()), "emax_gt150": int((df["em"]>150).sum()),
      "emax_gt300": int((df["em"]>300).sum()), "emax_neg": int((df["em"]<0).sum())}

# ---------- Triangle closure ----------
sub2 = df[rtc.notna()].copy()
sub2["pwy"] = (sub2["Pathway level"].astype(str)+"|"+sub2["Measured process"].astype(str)
               +"|"+sub2["Primary effector subtype"].astype(str))
sub2["grp2"] = sub2["doi"]+"|"+sub2["receptor"]+"|"+sub2["ligand"]+"|"+sub2["Cell line"].astype(str)
viol=[]
for k,v in sub2.groupby("grp2"):
    d=v.drop_duplicates("pwy")
    if len(d)<3: continue
    rs=d["rtc"].tolist()
    for i,j,kk in itertools.combinations(range(len(rs)),3):
        viol.append(abs((rs[i]-rs[j])+(rs[j]-rs[kk])-(rs[i]-rs[kk])))
viol=np.array(viol)

summary = {
  "rows_total": int(len(df)),
  "entries_with_tc_em_pec": int(mask.sum()),
  "groups_multi_ligand": int(len(grp_df)),
  "pairs_total": int(grp_df["n_pairs"].sum()),
  "pair_dev_quantiles_dex": {f"P{q}": float(np.percentile(
        np.concatenate([np.full(int(r.n_pairs), r.dev_median) for r in grp_df.itertuples()]), q))
        for q in [50,75,90,95,99]} if len(grp_df) else {},
  "groups_dev_median_gt_0.3": int((grp_df["dev_median"]>0.3).sum()),
  "papers_total": int(len(pap)),
  "papers_dev_median_gt_0.3": int((pap["dev_median"]>0.3).sum()),
  "emax_qc": qc,
  "closure_triangles": int(len(viol)),
  "closure_violations_gt_0.01": int((viol>0.01).sum()),
  "worst_groups": grp_df.head(20).to_dict("records"),
  "worst_papers": pap.head(20).to_dict("records"),
}
with open(os.path.join(OUTDIR,"atlas_audit_summary.json"),"w") as f:
    json.dump(summary, f, ensure_ascii=False, indent=1)

print(f"Entry level: {len(ent_out)} rows → atlas_audit_entries.csv")
print(f"Group level: {len(grp_df)} groups, pairs {int(grp_df['n_pairs'].sum())} → atlas_audit_groups.csv")
print(f"Paper level: {len(pap)} papers → atlas_audit_papers.csv")
print(f"Groups dev median>0.3: {(grp_df['dev_median']>0.3).sum()}/{len(grp_df)}; "
      f"papers dev median>0.3: {(pap['dev_median']>0.3).sum()}/{len(pap)}")
print(f"Triangle closures: {len(viol)}, violations {(viol>0.01).sum()}")
print(f"Emax QC: {qc}")
print("→ atlas_audit_summary.json")
