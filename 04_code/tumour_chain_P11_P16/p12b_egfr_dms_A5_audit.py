# -*- coding: utf-8 -*-
"""P12b EGFR arm, amendment A5: curated resistance lists vs all missense activity.
Saves results to p12b_egfr_dms_A5_results.json (archival fix)."""
import os, json
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
v = pd.read_csv(os.path.join(HERE, 'p12b_egfr_dms_variants.csv'))
v['residue'] = v['site'] + 711

curated = [(718,'Q'),(718,'V'),(714,'R'),(716,'T'),(725,'M'),(728,'E'),(754,'E'),
           (754,'N'),(771,'S'),(771,'T'),(783,'I'),(791,'L'),(791,'K'),(863,'S'),
           (895,'N'),(929,'I'),(971,'L')]
clinical = [(790,'M'),(797,'S'),(724,'S'),(792,'H'),(796,'S'),(796,'C'),(796,'R')]

def get(residue, aa):
    r = v[(v['residue'] == residue) & (v['mut_aa'] == aa)]
    return float(r['activity'].max()) if len(r) else np.nan

out = {}
for name, lst in (('curated', curated), ('clinical', clinical)):
    acts = np.array([get(r_, a_) for r_, a_ in lst], float)
    ok = ~np.isnan(acts)
    a = acts[ok]
    u = mannwhitneyu(a, v['activity'], alternative='greater')
    out[name] = {
        'n': int(ok.sum()), 'activities': [round(float(x), 4) for x in a],
        'median': float(np.median(a)),
        'all_missense_median': float(v['activity'].median()),
        'median_ratio': float(np.median(a) / v['activity'].median()),
        'mw_one_sided_p': float(u.pvalue),
        'percentiles': [round(float((v['activity'] < x).mean()), 3) for x in a],
    }
    print(name, json.dumps(out[name], indent=1)[:500])
json.dump(out, open(os.path.join(HERE, 'p12b_egfr_dms_A5_results.json'), 'w'), indent=1)
print('saved p12b_egfr_dms_A5_results.json')
