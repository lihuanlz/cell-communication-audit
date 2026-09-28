# -*- coding: utf-8 -*-
"""
P14 audit: do downstream decoders read p53 pulse vs sustained shape differently
at matched total exposure (AUC)?
Data: Jimenez-Asins et al. 2022, Mol Syst Biol 18:e10588 (Lahav lab),
GitHub albajimenezasins/Proteomics_MSB_2022 (RNA-seq TPM + TMT proteomics).
Pre-registered: 预注册卡_P14_p53下游误读_2026-09-28.md (amendment A1).
"""
import os, json, warnings
import numpy as np
import pandas as pd
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, 'data', 'msb2022_repo', 'Proteomics_MSB_2022-main',
                   'Figures_1_2_4_5_6', 'Raw_Data')

# ---------- sample map ----------
sm = pd.read_excel(os.path.join(RAW, 'samples_summary.xlsx'), header=None)
sm.columns = ['S', 'name', 'hour', 'rep', 'cond']  # cond: p / s
# untreated t0 samples named A0/B0/C0 with cond 'p' but shared by both arms
mapping = {r['S']: (r['name'], float(r['hour']), r['cond']) for _, r in sm.iterrows()}

tpm = pd.read_csv(os.path.join(RAW, 'tximport-tpm.csv'), index_col=0)
print('tpm', tpm.shape, 'index sample:', tpm.index[:3].tolist())

# gene symbol map
gm = pd.read_excel(os.path.join(RAW, 'GRCh38_93_genemap.xlsx'))
print('genemap cols', gm.columns.tolist())
idcol = [c for c in gm.columns if 'gene_id' in c.lower() or 'ensembl' in c.lower()]
symcol = [c for c in gm.columns if 'name' in c.lower() or 'symbol' in c.lower()]
id2sym = dict(zip(gm[idcol[0]].astype(str).str.split('.').str[0], gm[symcol[0]].astype(str)))
idx = tpm.index.astype(str).str.split('.').str[0]
tpm['sym'] = [id2sym.get(g, g) for g in idx]
tpm = tpm[tpm['sym'] != 'nan']
tpm = tpm.groupby('sym').max(numeric_only=True)  # collapse duplicates

# ---------- organize into condition x hour matrices ----------
def get_matrix(cond):
    cols0 = [s for s, (n, h, c) in mapping.items() if h == 0]  # A0 B0 C0 shared
    hours = sorted({h for _, (n, h, c) in mapping.items() if c == cond and h > 0})
    mat = {}
    base = tpm[cols0].mean(axis=1)
    for h in hours:
        cols = [s for s, (n, hh, c) in mapping.items() if c == cond and hh == h]
        mat[h] = tpm[cols].mean(axis=1)
    return base, pd.DataFrame(mat)

base_p, P = get_matrix('p')   # pulsatile (10 Gy X-ray)
base_s, S = get_matrix('s')   # sustained (Nutlin-3)
print('pulsatile hours', P.columns.tolist(), 'sustained hours', S.columns.tolist())

FCp = (P.T / base_p.replace(0, np.nan)).T
FCs = (S.T / base_s.replace(0, np.nan)).T

# ---------- verify AUC matching from TMT proteomics TP53 ----------
pq = pd.read_excel(os.path.join(RAW, 'DL_p53_Dynamics_Protein_Quant_sum-values[1].xlsx'),
                   sheet_name='protein_quant_19467', header=0)
tp53 = pq[pq['Gene Symbol'].astype(str).str.upper() == 'TP53']
print('TP53 rows:', len(tp53))
auc_report = {}
if len(tp53):
    row = tp53.iloc[0]
    hrs = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    for arm, reps in (('Pulse', ['A', 'B']), ('Sustain', ['C', 'D'])):
        for rep in reps:
            try:
                u = float(row['Untreated_%s~rq_126_sn sum' % rep])
                vals = []
                for h in hrs:
                    k1 = '%s_%dh_%s~' % (arm, h, rep)
                    col = [c for c in pq.columns if str(c).startswith(k1)]
                    vals.append(float(row[col[0]]) / u if col else np.nan)
                auc = np.nansum(vals)  # 1 h spacing
                auc_report['%s_%s' % (arm, rep)] = round(auc, 2)
            except Exception as e:
                auc_report['%s_%s' % (arm, rep)] = 'ERR %s' % e
print('p53 protein AUC (0-9h, fold-to-untreated x h):', auc_report)

# ---------- main test ----------
common = [h for h in [6, 7, 8, 9] if h in FCp.columns and h in FCs.columns]
# induced set: max FC >= 2 in either arm over matched hours
ind = ((FCp[common].max(axis=1) >= 2) | (FCs[common].max(axis=1) >= 2))
genes = FCp.index[ind]
print('induced genes (maxFC>=2 either arm, t6-9h):', len(genes))

rows = []
for g in genes:
    rp = FCp.loc[g, common].values
    rs = FCs.loc[g, common].values
    with np.errstate(divide='ignore', invalid='ignore'):
        log2R = np.log2(rs / rp)
    rows.append({'gene': g, 'log2R_mean': np.nanmean(log2R),
                 'log2R_max_abs': np.nanmax(np.abs(log2R)),
                 'FCp_max': np.nanmax(FCp.loc[g, common]),
                 'FCs_max': np.nanmax(FCs.loc[g, common])})
R = pd.DataFrame(rows).set_index('gene')
frac_up_s = (R['log2R_mean'] > 1).mean()
frac_up_p = (R['log2R_mean'] < -1).mean()
frac_mid = ((R['log2R_mean'].abs()) <= 1).mean()
print('fraction sustained>>pulsed (>2x): %.3f' % frac_up_s)
print('fraction pulsed>>sustained (>2x): %.3f' % frac_up_p)
print('fraction within 2x: %.3f' % frac_mid)

# replicate-noise null: use t0 replicate scatter as FC noise scale
t0cols = [s for s, (n, h, c) in mapping.items() if h == 0]
t0 = tpm.loc[genes, t0cols]
t0cv = (t0.std(axis=1) / t0.mean(axis=1)).replace([np.inf], np.nan).dropna()
print('t0 replicate CV median: %.3f (null FC noise ~sqrt(2)*CV)' % t0cv.median())

# ---------- known fate genes ----------
fate = ['CDKN1A', 'SFN', 'BBC3', 'BAX', 'FAS', 'PML', 'MDM2', 'GADD45A',
        'DDB2', 'XPC', 'POLH', 'RRM2B', 'SESN1', 'SESN2', 'TIGAR', 'PTEN',
        'TP53I3', 'FDXR', 'BCL2', 'MCL1', 'NOXA', 'PMAIP1', 'ZMAT3', 'BTG2']
print('\n%-10s %8s %8s %10s' % ('gene', 'FCp9h', 'FCs9h', 'log2R(s/p)'))
fate_rows = []
for g in fate:
    if g in FCp.index:
        fp = float(FCp.loc[g, 9]) if 9 in FCp.columns else np.nan
        fs = float(FCs.loc[g, 9]) if 9 in FCs.columns else np.nan
        lr = np.log2(fs / fp) if fp > 0 and fs > 0 else np.nan
        fate_rows.append({'gene': g, 'FCp_9h': fp, 'FCs_9h': fs, 'log2R_9h': lr})
        print('%-10s %8.2f %8.2f %+10.2f' % (g, fp, fs, lr))

# ---------- persist ----------
out = {
    'auc_p53_protein': auc_report,
    'n_induced': int(len(genes)),
    'frac_sustained_higher_2x': float(frac_up_s),
    'frac_pulsed_higher_2x': float(frac_up_p),
    'frac_within_2x': float(frac_mid),
    't0_cv_median': float(t0cv.median()),
    'fate_genes': fate_rows,
}
R.to_csv(os.path.join(HERE, 'p14_pergene_log2R.csv'))
with open(os.path.join(HERE, 'p14_results.json'), 'w') as f:
    json.dump(out, f, indent=1, ensure_ascii=False)
print('\nsaved p14_pergene_log2R.csv + p14_results.json')
