# -*- coding: utf-8 -*-
"""
P15 audit: does oncogenic Ras/Raf mutation decouple pERK from upstream EGFR
in patient tumors (TCGA RPPA, cBioPortal pancan atlas)?
Pre-registered: 预注册卡_P15_TCGA信道解耦_2026-09-28.md
"""
import os, json, time
import numpy as np
import pandas as pd
import requests
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'data')
API = 'https://www.cbioportal.org/api'
STUDIES = ['coadread_tcga_pan_can_atlas_2018', 'luad_tcga_pan_can_atlas_2018',
           'skcm_tcga_pan_can_atlas_2018', 'paad_tcga_pan_can_atlas_2018']
PERK = [-29, -30]           # MAPK1/3_PT202_Y204
EGFR_TOT, EGFR_P = 1956, -19
RASRAF = {'KRAS': 3845, 'BRAF': 673, 'HRAS': 3265, 'NRAS': 4893}

def post(path, body, tries=4):
    for a in range(tries):
        try:
            r = requests.post(API + path, json=body, timeout=300)
            if r.status_code == 200:
                return r.json()
        except Exception as e:
            print('retry', a, e)
        time.sleep(5 * (a + 1))
    raise RuntimeError('failed: ' + path)

def get_rppa(study):
    return post('/molecular-profiles/%s_rppa/molecular-data/fetch?projection=SUMMARY&pageSize=5000000'
                % study,
                {'entrezGeneIds': PERK + [EGFR_TOT, EGFR_P],
                 'sampleListId': study + '_rppa'})

def get_mut(study):
    return post('/molecular-profiles/%s_mutations/mutations/fetch?projection=DETAILED&pageSize=5000000'
                % study,
                {'entrezGeneIds': list(RASRAF.values()),
                 'sampleListId': study + '_rppa'})

def analyze(study):
    rp = pd.DataFrame(get_rppa(study))
    mu = pd.DataFrame(get_mut(study))
    if rp.empty:
        return None
    wide = rp.pivot_table(index='sampleId', columns='entrezGeneId',
                          values='value', aggfunc='mean')
    wide.columns = ['g%s' % c for c in wide.columns]
    wide['pERK'] = wide[['g-29', 'g-30']].mean(axis=1)
    wide = wide.rename(columns={'g1956': 'EGFR', 'g-19': 'pEGFR'})
    # mutation classification
    mu = mu[mu['mutationType'] != 'Silent'] if 'mutationType' in mu else mu
    cls = {}
    for _, r in mu.iterrows():
        sid = r.get('sampleId'); gene = r.get('gene', {}).get('hugoGeneSymbol', '')
        pc = str(r.get('proteinChange', ''))
        cur = cls.get(sid, set())
        tag = None
        if gene == 'BRAF' and 'V600E' in pc: tag = 'BRAF_V600E'
        elif gene == 'KRAS' and pc.startswith('G12'): tag = 'KRAS_G12'
        elif gene == 'KRAS' and pc.startswith('Q61'): tag = 'KRAS_Q61'
        elif gene == 'NRAS' and pc.startswith('Q61'): tag = 'NRAS_Q61'
        elif gene == 'BRAF': tag = 'BRAF_other'
        elif gene in ('KRAS', 'HRAS', 'NRAS'): tag = 'RAS_other'
        if tag: cur.add(tag)
        cls[sid] = cur
    wide['group'] = 'WT'
    strong = {'BRAF_V600E', 'KRAS_Q61', 'NRAS_Q61'}
    for sid, tags in cls.items():
        if sid not in wide.index: continue
        if tags & strong: wide.loc[sid, 'group'] = 'strong'
        elif tags & {'KRAS_G12'}: wide.loc[sid, 'group'] = 'KRAS_G12'
        elif tags: wide.loc[sid, 'group'] = 'other_mut'
    wide = wide.dropna(subset=['pERK', 'EGFR'])
    res = {'study': study, 'n': int(len(wide)),
           'counts': wide['group'].value_counts().to_dict()}
    # slope pERK ~ EGFR per group
    slopes = {}
    for g in ['WT', 'KRAS_G12', 'strong', 'other_mut']:
        sub = wide[wide['group'] == g]
        if len(sub) >= 30:
            sl, ic, rv, pv, se = stats.linregress(sub['EGFR'], sub['pERK'])
            cv = float(sub['pERK'].std() / abs(sub['pERK'].mean())) if sub['pERK'].mean() != 0 else np.nan
            slopes[g] = {'n': int(len(sub)), 'slope': float(sl), 'slope_se': float(se),
                         'r': float(rv), 'p': float(pv), 'cv_perk': cv,
                         'mean_perk': float(sub['pERK'].mean())}
    res['slopes'] = slopes
    if 'WT' in slopes and 'strong' in slopes:
        res['slope_ratio_strong_vs_WT'] = slopes['strong']['slope'] / slopes['WT']['slope'] if slopes['WT']['slope'] else np.nan
    if 'WT' in slopes and 'KRAS_G12' in slopes:
        res['slope_ratio_G12_vs_WT'] = slopes['KRAS_G12']['slope'] / slopes['WT']['slope'] if slopes['WT']['slope'] else np.nan
    # interaction test on strong+WT
    sub = wide[wide['group'].isin(['WT', 'strong'])].copy()
    if len(sub) >= 60 and (sub['group'] == 'strong').sum() >= 30:
        sub['mut'] = (sub['group'] == 'strong').astype(float)
        X = np.column_stack([np.ones(len(sub)), sub['EGFR'], sub['mut'], sub['EGFR'] * sub['mut']])
        y = sub['pERK'].values
        beta, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        resid = y - X @ beta
        dof = len(sub) - 4
        s2 = resid @ resid / dof
        se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
        tint = beta[3] / se[3]
        res['interaction'] = {'beta': float(beta[3]), 't': float(tint),
                              'p': float(2 * stats.t.sf(abs(tint), dof))}
    wide.to_csv(os.path.join(OUT, 'p15_%s_persample.csv' % study.split('_')[0]))
    return res

def main():
    allres = []
    for st in STUDIES:
        print('===', st, flush=True)
        r = analyze(st)
        if r:
            allres.append(r)
            print(json.dumps(r, indent=1, ensure_ascii=False)[:1200])
    with open(os.path.join(HERE, 'p15_tcga_results.json'), 'w') as f:
        json.dump(allres, f, indent=1)
    print('saved p15_tcga_results.json')

if __name__ == '__main__':
    main()
