# -*- coding: utf-8 -*-
"""
P13 robustness recheck: bin-count sensitivity + KSG continuous estimator.
Same data and response definition as p13_mutual_info_audit.py (locked).
Adds: nbins in {4,6,8,12,16} and KSG (sklearn mutual_info_regression, discrete X),
for the EGF main pairs. Bootstrap 300 for KSG deltas.
"""
import os, json, warnings
import numpy as np
from scipy.io import loadmat
from sklearn.feature_selection import mutual_info_regression
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
MAT = os.path.join(HERE, 'data', 'gillies2020_suppl', 'Gillies2020_SourceData_Fig2.mat')
m = loadmat(MAT, squeeze_me=False, struct_as_record=False)
d = m['MEF_alldata'][0, 0]
t = m['new_tseries'][0]
PRE = t < 0
POST = (t > 0) & (t <= 1800)
CONDS = ['noGFNaN', 'EGF0_1', 'EGF1', 'EGF10']
PAIRS = ['G12C', 'G12D', 'G12V', 'Q61R', 'BRafV600E']

def cell_responses(line, cond):
    arr = getattr(getattr(d, line)[0, 0], cond)
    outs = []
    for rep in arr.flat:
        if isinstance(rep, np.ndarray) and rep.ndim == 2:
            base = np.nanmedian(rep[:, PRE], axis=1)
            peak = np.nanmax(rep[:, POST], axis=1)
            r = peak - base
            outs.append(r[np.isfinite(r)])
    return np.concatenate(outs)

DATA = {L: [cell_responses(L, c) for c in CONDS] for L in ['KrasWT'] + PAIRS}

def mi_mm(x, y, nbins):
    n = len(x)
    qs = np.quantile(y, np.linspace(0, 1, nbins + 1))
    qs[0] -= 1e-12; qs[-1] += 1e-12
    yb = np.clip(np.digitize(y, qs) - 1, 0, nbins - 1)
    Kx = len(np.unique(x))
    pxy = np.zeros((Kx, nbins))
    for xi, yi in zip(x, yb):
        pxy[xi, yi] += 1
    pxy /= n
    px = pxy.sum(1, keepdims=True); py = pxy.sum(0, keepdims=True)
    nz = pxy > 0
    I = np.sum(pxy[nz] * np.log2(pxy[nz] / (px @ py)[nz]))
    return I - (np.count_nonzero(pxy) - Kx - nbins + 1) / (2 * n * np.log(2))

def mi_ksg(x, y):
    # sklearn: I in nats; X = discrete dose column, y = continuous response
    mi = mutual_info_regression(x.reshape(-1, 1).astype(float), y, discrete_features=[True],
                                n_neighbors=4, random_state=0)
    return float(mi[0]) / np.log(2)  # to bits

def xy(per_dose, idx=None):
    if idx is None:
        xs = np.concatenate([np.full(len(r), k) for k, r in enumerate(per_dose)])
        ys = np.concatenate(per_dose)
    else:
        xs = np.concatenate([np.full(len(r), k) for k, r in enumerate(per_dose)])
        ys = np.concatenate([r[idx[k]] for k, r in enumerate(per_dose)])
    return xs, ys

rng = np.random.default_rng(7)
res = {'bins': {}, 'ksg': {}}
for nb in [4, 6, 8, 12, 16]:
    row = {}
    Iwt = mi_mm(*xy(DATA['KrasWT']), nb)
    row['KrasWT'] = round(Iwt, 3)
    for L in PAIRS:
        Im = mi_mm(*xy(DATA[L]), nb)
        row[L] = round(Im - Iwt, 3)
    res['bins'][str(nb)] = row
    print('nbins=%2d  WT=%.3f  ' % (nb, Iwt) +
          '  '.join('%s=%+.3f' % (L, row[L]) for L in PAIRS))

print('\nKSG (bits):')
Iwt_k = mi_ksg(*xy(DATA['KrasWT']))
res['ksg']['KrasWT'] = round(Iwt_k, 3)
print('KrasWT=%.3f' % Iwt_k)
for L in PAIRS:
    Im = mi_ksg(*xy(DATA[L]))
    # bootstrap delta CI
    deltas = np.empty(300)
    pdose_m, pdose_w = DATA[L], DATA['KrasWT']
    for b in range(300):
        im_b = mi_ksg(*xy(pdose_m, [rng.integers(0, len(r), len(r)) for r in pdose_m]))
        iw_b = mi_ksg(*xy(pdose_w, [rng.integers(0, len(r), len(r)) for r in pdose_w]))
        deltas[b] = im_b - iw_b
    ci = np.percentile(deltas, [2.5, 97.5])
    res['ksg'][L] = {'I': round(Im, 3), 'delta': round(Im - Iwt_k, 3),
                     'CI95': [round(float(ci[0]), 3), round(float(ci[1]), 3)]}
    print('%-10s I=%.3f delta=%+.3f CI[%+.3f,%+.3f]' % (L, Im, Im - Iwt_k, ci[0], ci[1]))

with open(os.path.join(HERE, 'p13_robustness.json'), 'w') as f:
    json.dump(res, f, indent=1)
print('\nsaved p13_robustness.json')
