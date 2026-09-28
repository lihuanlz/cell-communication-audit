# -*- coding: utf-8 -*-
"""
P13 audit: does oncogenic Ras/Raf mutation reduce mutual information
I(ligand dose; single-cell ERK output)?
Data: Gillies et al. 2020, Mol Syst Biol 16:e9518 (PMC7569415), Source Data Fig2.
Pre-registered: 预注册卡_P13_信道容量互信息_2026-09-28.md (amendment A1).
Output metric (locked): per-cell response = max(EKAR3, 0<t<=30min) - median(EKAR3, t<0).
MI: equiprobable 8-bin discretization per cell line, plug-in + Miller-Madow.
Uncertainty: cell-level bootstrap (2000), paired for delta-I vs KrasWT.
"""
import os, sys, json, warnings
import numpy as np
from scipy.io import loadmat
warnings.filterwarnings('ignore')

HERE = os.path.dirname(os.path.abspath(__file__))
MAT = os.path.join(HERE, 'data', 'gillies2020_suppl', 'Gillies2020_SourceData_Fig2.mat')

m = loadmat(MAT, squeeze_me=False, struct_as_record=False)
d = m['MEF_alldata'][0, 0]
t = m['new_tseries'][0]

# EGF dose ladder: noGF (=0), 0.1, 1, 10 ng/ml
DOSES = {'noGFNaN': 0.0, 'EGF0_1': 0.1, 'EGF1': 1.0, 'EGF10': 10.0}
LIGANDS = {
    'EGF':  ['noGFNaN', 'EGF0_1', 'EGF1', 'EGF10'],
    'AREG': ['noGFNaN', 'AREG1', 'AREG10', 'AREG100'],
    'HGF':  ['noGFNaN', 'HGF1', 'HGF10', 'HGF100'],
    'PDGF': ['noGFNaN', 'PDGF1', 'PDGF10', 'PDGF100'],
    'FGF':  ['noGFNaN', 'FGF0_1', 'FGF1', 'FGF10'],
    'IGF':  ['noGFNaN', 'IGF1', 'IGF10', 'IGF100'],
}
PRE = t < 0
POST = (t > 0) & (t <= 1800)  # 0-30 min
LINES = ['KrasWT', 'HRasWT', 'NrasWT', 'G12C', 'G12D', 'G12V', 'Q61R', 'BRafV600E']


def cell_responses(line, cond):
    """Pool replicates; return per-cell response vector."""
    arr = getattr(getattr(d, line)[0, 0], cond)
    outs = []
    for rep in arr.flat:
        if not isinstance(rep, np.ndarray) or rep.ndim != 2:
            continue
        base = np.nanmedian(rep[:, PRE], axis=1)
        peak = np.nanmax(rep[:, POST], axis=1)
        r = peak - base
        r = r[np.isfinite(r)]
        outs.append(r)
    return np.concatenate(outs) if outs else np.array([])


def mi_mm(x, y, nbins=8):
    """Plug-in MI with Miller-Madow correction. x discrete labels, y continuous."""
    n = len(x)
    qs = np.quantile(y, np.linspace(0, 1, nbins + 1))
    qs[0] -= 1e-12; qs[-1] += 1e-12
    yb = np.clip(np.digitize(y, qs) - 1, 0, nbins - 1)
    Kx = len(np.unique(x)); Ky = nbins
    # joint histogram
    Kxy = Kx * Ky
    pxy = np.zeros((Kx, Ky))
    for xi, yi in zip(x, yb):
        pxy[xi, yi] += 1
    pxy /= n
    px = pxy.sum(1, keepdims=True); py = pxy.sum(0, keepdims=True)
    nz = pxy > 0
    I = np.sum(pxy[nz] * np.log2(pxy[nz] / (px @ py)[nz]))
    Kxy_eff = np.count_nonzero(pxy)
    I_mm = I - (Kxy_eff - Kx - Ky + 1) / (2 * n * np.log(2))
    return I, I_mm


def build_xy(line, conds):
    xs, ys = [], []
    for k, c in enumerate(conds):
        r = cell_responses(line, c)
        xs.append(np.full(len(r), k))
        ys.append(r)
    return np.concatenate(xs), np.concatenate(ys)


def boot_delta(line_mut, conds, rng, nboot=2000):
    """Paired bootstrap: resample cells per dose in both lines with same draw indices."""
    data = {}
    for line in (line_mut, 'KrasWT'):
        per_dose = [cell_responses(line, c) for c in conds]
        data[line] = per_dose
    def I_of(per_dose):
        x = np.concatenate([np.full(len(r), k) for k, r in enumerate(per_dose)])
        y = np.concatenate(per_dose)
        return mi_mm(x, y)[1]
    I0_mut = I_of(data[line_mut]); I0_wt = I_of(data['KrasWT'])
    deltas = np.empty(nboot)
    for b in range(nboot):
        ds = []
        for line in (line_mut, 'KrasWT'):
            resampled = []
            for r in data[line]:
                idx = rng.integers(0, len(r), len(r))
                resampled.append(r[idx])
            ds.append(I_of(resampled))
        deltas[b] = ds[0] - ds[1]
    return I0_mut, I0_wt, np.percentile(deltas, [2.5, 97.5]), deltas.mean()


def main():
    rng = np.random.default_rng(20260928)
    print('time range: %.0f to %.0f s, %d points' % (t[0], t[-1], len(t)))
    out = {'ligands': {}, 'egf_main': {}}
    for lig, conds in LIGANDS.items():
        # check all lines have all conds
        ok = all(hasattr(getattr(d, L)[0, 0], c) for L in LINES for c in conds)
        if not ok:
            out['ligands'][lig] = 'conditions incomplete, skipped'
            print(lig, 'conditions incomplete, skipped')
            continue
        res = {}
        I_wt_all = {}
        # per-line MI first
        for L in LINES:
            x, y = build_xy(L, conds)
            I, Imm = mi_mm(x, y)
            res[L] = {'I_plugin': float(I), 'I_mm': float(Imm), 'n': int(len(y))}
        # bootstrap deltas vs KrasWT
        for L in ['G12C', 'G12D', 'G12V', 'Q61R', 'BRafV600E']:
            Im, Iw, ci, dm = boot_delta(L, conds, rng)
            res[L]['I_mm_boot'] = float(Im)
            res[L]['delta_vs_KrasWT'] = float(Im - Iw)
            res[L]['delta_CI95'] = [float(ci[0]), float(ci[1])]
        # KrasWT self CI via bootstrap
        per_dose = [cell_responses('KrasWT', c) for c in conds]
        boots = np.empty(500)
        for b in range(500):
            rs = [r[rng.integers(0, len(r), len(r))] for r in per_dose]
            x = np.concatenate([np.full(len(r), k) for k, r in enumerate(rs)])
            y = np.concatenate(rs)
            boots[b] = mi_mm(x, y)[1]
        res['KrasWT']['I_CI95'] = [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))]
        out['ligands'][lig] = res
        print('\n=== %s ===' % lig)
        for L in LINES:
            r = res[L]
            line = '%-10s I_mm=%.3f n=%d' % (L, r['I_mm'], r['n'])
            if 'delta_vs_KrasWT' in r:
                line += '  delta=%+.3f CI[%+.3f,%+.3f]' % (
                    r['delta_vs_KrasWT'], r['delta_CI95'][0], r['delta_CI95'][1])
            if 'I_CI95' in r:
                line += '  I_CI=[%.3f,%.3f]' % tuple(r['I_CI95'])
            print(line)
    with open(os.path.join(HERE, 'p13_mi_results.json'), 'w') as f:
        json.dump(out, f, indent=1)
    print('\nsaved p13_mi_results.json')


if __name__ == '__main__':
    main()
