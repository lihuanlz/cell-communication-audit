# P9 audit: Bcd profile absolute-anchor test (frozen prereg 2026-09-27)
# Data: Nikolic et al. 2024 repo (Liu 2013 Bcd profiles; Petkova 2019 gap genes)
import numpy as np
from scipy.io import loadmat
from scipy import stats
import json, os

SEED = 20260927
rng = np.random.default_rng(SEED)
HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, 'data_nikolic2024', 'Scale-invariance-in-early-embryonic-development-main')

def bootstrap_ci(x, y, n=2000):
    slopes = []
    N = len(x)
    for _ in range(n):
        idx = rng.integers(0, N, N)
        if np.std(x[idx]) < 1e-12:
            continue
        s, *_ = stats.linregress(x[idx], y[idx])
        slopes.append(s)
    slopes = np.array(slopes)
    return np.percentile(slopes, [2.5, 97.5]), np.std(slopes)

def theil_sen(x, y):
    # simple Theil-Sen on log-log
    N = len(x)
    sl = []
    idx = rng.integers(0, N, size=(min(20000, N*(N-1)//2), 2))
    for a, b in idx:
        if a != b and abs(x[a]-x[b]) > 1e-12:
            sl.append((y[a]-y[b])/(x[a]-x[b]))
    return float(np.median(sl))

# ---------- P9-1: per-embryo lambda in absolute units ----------
m = loadmat(os.path.join(D, 'Fig5_rawProfiles_Bcd.mat'))
print('keys:', [k for k in m if not k.startswith('__')])
L = np.asarray(m['L']).ravel().astype(float)
xs = np.asarray(m['xs']).ravel().astype(float)
P = np.asarray(m['profiles_bcd']).astype(float)   # 90 x 582 (rows=pos, cols=embryo)
print('embryos:', P.shape, 'L range:', L.min(), L.max())

win = (xs >= 0.10) & (xs <= 0.50)      # frozen window (relative)
lams, Ls, r2s, B0s, mono = [], [], [], [], []
for j in range(P.shape[1]):
    prof = P[:, j]
    xw = xs[win] * L[j]                # ABSOLUTE microns
    yw = prof[win]
    ok = np.isfinite(yw) & (yw > 0)
    if ok.sum() < 8:
        continue
    xw, yw = xw[ok], yw[ok]
    ly = np.log(yw)
    s, b, r, _, _ = stats.linregress(xw, ly)
    if s >= 0:
        continue
    # monotonicity fraction within window
    d = np.diff(yw)
    mfrac = np.mean(d <= 0)
    if r*r < 0.8 or mfrac < 0.6:       # frozen filter
        continue
    lams.append(-1.0/s); Ls.append(L[j]); r2s.append(r*r)
    b0 = prof[(xs >= 0.05) & (xs <= 0.10)]
    B0s.append(np.nanmean(b0)); mono.append(mfrac)

lams = np.array(lams); Ls = np.array(Ls); r2s = np.array(r2s); B0s = np.array(B0s)
print(f'P9-1 valid embryos: {len(lams)}/{P.shape[1]}  median R2={np.median(r2s):.4f}')

x = np.log(Ls); y = np.log(lams)
alpha, a0, r_a, p_a, se_a = stats.linregress(x, y)
(ci_lo, ci_hi), sd_boot = bootstrap_ci(x, y)
alpha_ts = theil_sen(x, y)
res = dict(P9_1=dict(n=len(lams), alpha=float(alpha), se=float(se_a),
                     ci95=[float(ci_lo), float(ci_hi)], theil_sen=float(alpha_ts),
                     median_lambda_um=float(np.median(lams)),
                     L_range=[float(Ls.min()), float(Ls.max())]))
print(res['P9_1'])

# ---------- P9-3 exploratory: B0 amplitude vs L ----------
g, g0, r_g, _, _ = stats.linregress(np.log(Ls), np.log(B0s))
(ci_lo_g, ci_hi_g), sd_g = bootstrap_ci(np.log(Ls), np.log(B0s))
res['P9_3'] = dict(gamma_B0=float(g), ci95=[float(ci_lo_g), float(ci_hi_g)])
print(res['P9_3'])

# ---------- P9-2: Hb boundary scaling (same repo, gap genes) ----------
try:
    gmat = loadmat(os.path.join(D, 'rawProfiles_gapGenes_Hb_Gt_Kni_Kr.mat'))
    d = gmat['data']                      # (1,609) struct array, per-embryo fields
    n_emb = d.shape[1]
    xb, Lb = [], []
    for j in range(n_emb):
        Lg = float(np.asarray(d[0, j]['L']).ravel()[0])
        xsok = np.asarray(d[0, j]['xs']).ravel().astype(float)
        pok = np.asarray(d[0, j]['Hb']).ravel().astype(float)
        ok = np.isfinite(pok)
        if ok.sum() < 20:
            continue
        xsok, pok = xsok[ok], pok[ok]
        half = 0.5*np.nanmax(pok)
        # anterior boundary: first crossing from above going posterior
        above = pok >= half
        idx = np.where(np.diff(above.astype(int)) == -1)[0]
        if len(idx) == 0:
            continue
        i0 = idx[0]
        x0, x1 = xsok[i0], xsok[i0+1]
        y0, y1 = pok[i0], pok[i0+1]
        if y1 == y0:
            continue
        frac = (half - y0)/(y1 - y0)
        xb_rel = x0 + frac*(x1-x0)
        xb.append(xb_rel*Lg); Lb.append(Lg)   # ABSOLUTE boundary position
    xb = np.array(xb); Lb = np.array(Lb)
    ab, _, _, _, se_b = stats.linregress(np.log(Lb), np.log(xb))
    (cb_lo, cb_hi), _ = bootstrap_ci(np.log(Lb), np.log(xb))
    res['P9_2'] = dict(n=len(xb), alpha_boundary=float(ab), se=float(se_b),
                       ci95=[float(cb_lo), float(cb_hi)])
    print(res['P9_2'])
except Exception as e:
    res['P9_2'] = dict(error=str(e))
    print('P9-2 failed:', e)

out = os.path.join(HERE, 'p9_bcd_anchor_results.json')
json.dump(res, open(out, 'w'), indent=2)
print('saved', out)
