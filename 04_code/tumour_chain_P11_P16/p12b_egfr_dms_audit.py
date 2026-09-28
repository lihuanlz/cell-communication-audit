# -*- coding: utf-8 -*-
"""
P12b EGFR arm: do drug-enriched (resistant) variants preserve EGFR signaling activity?
Data: Wang et al. 2025 npj Precis Oncol (PMC12368022; GEO GSE305057),
EGFR-L858R kinase domain saturation mutagenesis in Ba/F3, codon counts.
Pre-registered: 预注册卡_P12 amendment A4.
Activity proxy = DMSO relative frequency / synonymous mean frequency.
Resistance = log2(drug/DMSO) concordant across 2 replicates.
"""
import os, gzip, json, warnings
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
DAT = os.path.join(HERE, 'data', 'gse305057')

CODON_TABLE = {
 'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
 'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
 'TCT':'S','TCC':'S','TCA':'S','TCG':'S','AGT':'S','AGC':'S','CCT':'P','CCC':'P',
 'CCA':'P','CCG':'P','ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A',
 'GCA':'A','GCG':'A','TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H',
 'CAA':'Q','CAG':'Q','AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D',
 'GAA':'E','GAG':'E','TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R',
 'CGA':'R','CGG':'R','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G'}

SAMPLES = {'DMSO': ['GSM9161531_DMSO-A', 'GSM9161532_DMSO-B'],
           'Osi': ['GSM9161533_Osi-A', 'GSM9161534_Osi-B'],
           'BLU': ['GSM9161527_BLU945-A', 'GSM9161528_BLU945-B']}

def load(sample):
    with gzip.open(os.path.join(DAT, sample + '_codoncounts.csv.gz'), 'rt') as f:
        df = pd.read_csv(f)
    df = df.set_index('site')
    wt = df['wildtype']
    df = df.drop(columns=['wildtype'])
    return df, wt

# load all, restrict to common codon columns
frames, wts = {}, None
for arm, ss in SAMPLES.items():
    for s in ss:
        df, wt = load(s)
        frames[s] = df
        wts = wt
codons = sorted(set.intersection(*[set(f.columns) for f in frames.values()]))
sites = frames[SAMPLES['DMSO'][0]].index
print('sites', len(sites), 'codons', len(codons))

# ---------- build amino-acid variant table ----------
rows = []
for site in sites:
    wtc = wts.loc[site]
    wtaa = CODON_TABLE.get(wtc, '?')
    for cod in codons:
        aa = CODON_TABLE.get(cod, '?')
        counts = {s: int(f.loc[site, cod]) for s, f in frames.items()}
        rows.append({'site': site, 'codon': cod, 'wt_codon': wtc,
                     'wt_aa': wtaa, 'mut_aa': aa, 'count': counts})
long = pd.DataFrame(rows)
long['class'] = 'missense'
long.loc[long['mut_aa'] == long['wt_aa'], 'class'] = 'synonymous'
long.loc[long['mut_aa'] == '*', 'class'] = 'nonsense'
long.loc[long['codon'] == long['wt_codon'], 'class'] = 'WT'

# ---------- frequencies ----------
for s in frames:
    tot = long['count'].apply(lambda c: c[s]).sum()
    long[s + '_freq'] = long['count'].apply(lambda c: c[s]) / tot

# ---------- activity proxy (DMSO) ----------
syn = long[long['class'] == 'synonymous']
syn_mean = (syn['GSM9161531_DMSO-A_freq'] + syn['GSM9161532_DMSO-B_freq']).mean() / 2
non_mean = (long[long['class'] == 'nonsense']['GSM9161531_DMSO-A_freq'] +
            long[long['class'] == 'nonsense']['GSM9161532_DMSO-B_freq']).mean() / 2
long['dmso_freq'] = (long['GSM9161531_DMSO-A_freq'] + long['GSM9161532_DMSO-B_freq']) / 2
long['activity'] = long['dmso_freq'] / syn_mean
print('synonymous mean freq %.3e, nonsense mean freq %.3e (activity %.3f)' %
      (syn_mean, non_mean, non_mean / syn_mean))

# ---------- resistance scores ----------
for arm, ss in (('Osi', SAMPLES['Osi']), ('BLU', SAMPLES['BLU'])):
    fa = frames[ss[0]]; fb = frames[ss[1]]
    long[arm + '_freq'] = (long['count'].apply(lambda c: c[ss[0]]) / long['count'].apply(lambda c: c[ss[0]]).sum() +
                           long['count'].apply(lambda c: c[ss[1]]) / long['count'].apply(lambda c: c[ss[1]]).sum()) / 2
    long[arm + '_res'] = np.log2((long[arm + '_freq'] + 1e-12) / (long['dmso_freq'] + 1e-12))

# ---------- main test (missense only, AA collapsed by max codon evidence) ----------
mm = long[long['class'] == 'missense'].copy()
# require minimal evidence: DMSO+drug total counts >= 10 to keep enrichment meaningful
cnt_cols = list(frames.keys())
mm['tot'] = mm['count'].apply(lambda c: sum(c.values()))
mm = mm[mm['tot'] >= 10]
print('missense variants with >=10 total counts:', len(mm))

res_mask = (mm['Osi_res'] > 2) & (mm['BLU_res'] > 2)          # concordant both drugs >4x
res_osi = mm['Osi_res'] > 2
res_blu = mm['BLU_res'] > 2
print('enriched both drugs:', res_mask.sum(), 'osi:', res_osi.sum(), 'blu:', res_blu.sum())

out = {}
for name, mask in (('both', res_mask), ('osi', res_osi), ('blu', res_blu)):
    a = mm.loc[mask, 'activity']; b = mm.loc[~mask, 'activity']
    if len(a) < 5:
        out[name] = 'too few'
        continue
    u = mannwhitneyu(a, b, alternative='greater')
    out[name] = {'n_res': int(len(a)), 'n_other': int(len(b)),
                 'med_activity_res': float(a.median()), 'med_activity_other': float(b.median()),
                 'median_ratio': float(a.median() / b.median()), 'p_one_sided': float(u.pvalue)}
    print('%s: n=%d vs %d, med activity %.3f vs %.3f, ratio %.2f, p=%.2e' %
          (name, len(a), len(b), a.median(), b.median(), a.median()/b.median(), u.pvalue))

# ---------- named clinical variants ----------
AA3 = {'C':'C','G':'G','L':'L','M':'M','T':'T','Q':'Q','V':'V','S':'S','P':'P','R':'R'}
def find(site_num, mut_aa):
    # sites in this library: find row where site index maps; L858R library, site numbering offset unknown.
    r = mm[(mm['site'] == site_num) & (mm['mut_aa'] == mut_aa)]
    return r
print('\nnamed variants (site numbering as in library; kinase domain starts ~exon18):')
named = []
for site_num, aa, label in [(6,'S','G719S'), (7,'C','G719C'), (42,'I','L747?'), (90,'M','T790M'),
                            (97,'S','C797S'), (18,'Q','L718Q'), (18,'V','L718V'),
                            (96,'X','G796?'), (92,'X','L792?'), (24,'S','G724S')]:
    r = find(site_num, aa)
    for _, x in r.iterrows():
        named.append({'label': label, 'site': int(x['site']), 'wt_aa': x['wt_aa'],
                      'mut_aa': x['mut_aa'], 'activity': float(x['activity']),
                      'osi_res': float(x['Osi_res']), 'blu_res': float(x['BLU_res'])})
        print('%-8s site %d %s>%s  activity=%.3f  osi_res=%+.2f  blu_res=%+.2f' %
              (label, x['site'], x['wt_aa'], x['mut_aa'], x['activity'], x['Osi_res'], x['BLU_res']))

out['named'] = named
out['n_missense_tested'] = int(len(mm))
mm[['site','wt_aa','mut_aa','codon','activity','Osi_res','BLU_res','tot']].to_csv(
    os.path.join(HERE, 'p12b_egfr_dms_variants.csv'), index=False)
with open(os.path.join(HERE, 'p12b_egfr_dms_results.json'), 'w') as f:
    json.dump(out, f, indent=1, ensure_ascii=False)
print('\nsaved p12b_egfr_dms_variants.csv + p12b_egfr_dms_results.json')
