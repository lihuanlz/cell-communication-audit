# Cascade audit: GPCR (muscarinic ACh) -> Ca2+ single-cell peaks
# Keshelava 2018 SD3: 27 experiments, 433 cells, 7 doses x 5 reps per cell.
# Row->dose mapping: rows are descending dose (population mean peak declines monotonically 0.215->0.018 across rows; non-responding lowest dose at last row, matching the paper's statement that cells do not respond to the lowest ACh). Absolute ladder from Methods: [10,3,1.5,0.75,0.5,0.25,0.1] uM; the SD2 MATLAB script carries a legacy ladder [4,2,1,0.5,0.25,0.13,0.06], registered as provenance discrepancy (all statistics here are monotone-invariant except absolute EC50, which uses the Methods ladder).
import pickle, numpy as np, json
from scipy import stats

out = pickle.load(open('keshelava_parsed.pkl','rb'))
CONCS = np.array([10,3,1.5,0.75,0.5,0.25,0.1])  # uM, Methods ladder (Keshelava 2018 Methods); rows descending, row0 = 10 uM

cells = []  # (experiment, cell_id, matrix 7x5)
for name, doses in out.items():
    for cid, rows_ in doses.items():
        a = np.array([r for r in rows_ if len(r)==5])
        if a.shape == (7,5):
            cells.append((name, cid, a))
print('usable cells:', len(cells))

# ---------- A. population dose-response per experiment: EC50/n dispersion ----------
def hill(c, top, ec50, n):
    return top * c**n / (ec50**n + c**n)
from scipy.optimize import curve_fit
exp_params = {}
for name in sorted(set(c[0] for c in cells)):
    mats = [c[2] for c in cells if c[0]==name]
    pop = np.array([m.mean(axis=1) for m in mats])      # cells x 7
    mu = pop.mean(axis=0)                                # 7 doses
    if mu.max() <= 0: continue
    try:
        p,_ = curve_fit(hill, CONCS, mu, p0=[mu.max(), 0.3, 1.0],
                        bounds=([0,0.01,0.1],[5*mu.max(),50,10]), maxfev=20000)
        exp_params[name] = dict(top=float(p[0]), ec50=float(p[1]), n=float(p[2]),
                                ncells=len(mats))
    except Exception as e:
        exp_params[name] = dict(error=str(e))
tops = [v['top'] for v in exp_params.values() if 'top' in v]
ec50 = [v['ec50'] for v in exp_params.values() if 'ec50' in v]
ns   = [v['n'] for v in exp_params.values() if 'n' in v]
print('A. per-experiment Hill fits:', len(ec50))
print('   top   median %.3f  CV %.2f' % (np.median(tops), np.std(tops)/np.mean(tops)))
print('   EC50  median %.3f uM  CV %.2f  range %.3f-%.3f' % (np.median(ec50), np.std(ec50)/np.mean(ec50), min(ec50), max(ec50)))
print('   n     median %.2f  CV %.2f' % (np.median(ns), np.std(ns)/np.mean(ns)))

# ---------- B. per-cell repeatability vs across-cell dispersion ----------
rat_cv_within, cv_across = [], []
for d in range(7):
    pools = []
    for name,cid,a in cells:
        pools.append(a[d])
    pools = np.array(pools)   # cells x 5 reps at dose d
    per_cell_cv = pools.std(axis=1) / np.clip(pools.mean(axis=1),1e-9,None)
    rat_cv_within.append(np.median(per_cell_cv[np.isfinite(per_cell_cv)]))
    cv_across.append(pools.mean(axis=1).std() / pools.mean(axis=1).mean())
print('B. median within-cell CV (5 reps) per dose:', [round(x,3) for x in rat_cv_within])
print('   across-cell CV of mean per dose:      ', [round(x,3) for x in cv_across])

# rank preservation: Spearman of per-cell mean response between adjacent doses
rhos = []
for d in range(6):
    x = np.array([c[2][d].mean() for c in cells])
    y = np.array([c[2][d+1].mean() for c in cells])
    rhos.append(stats.spearmanr(x,y).statistic)
print('   Spearman rank corr adjacent doses:', [round(r,2) for r in rhos])

# ---------- C. adjacent-dose discrimination: raw amplitude vs self-normalised ----------
def auc(x, y):
    # P(x>y) via rank
    r = stats.rankdata(np.concatenate([x,y]))
    return (r[:len(x)].sum() - len(x)*(len(x)+1)/2) / (len(x)*len(y))
auc_raw, auc_ratio = [], []
for d in range(6):
    hi = np.array([c[2][d].mean() for c in cells])
    lo = np.array([c[2][d+1].mean() for c in cells])
    auc_raw.append(auc(hi, lo))
    # self-normalised: response at dose d+1 relative to cell's own max-dose response
    ref = np.array([c[2][0].mean() for c in cells])
    hi_r = hi/np.clip(ref,1e-9,None); lo_r = lo/np.clip(ref,1e-9,None)
    auc_ratio.append(auc(hi_r, lo_r))
print('C. AUC adjacent doses, raw amplitude:   ', [round(a,3) for a in auc_raw])
print('   AUC adjacent doses, self-normalised:', [round(a,3) for a in auc_ratio])

# ---------- D. responding fraction (peak > threshold) vs dose ----------
thr = 0.02
fracs = []
for d in range(7):
    m = np.array([c[2][d].mean() for c in cells])
    fracs.append(float((m>thr).mean()))
print('D. fraction of cells with mean peak > 0.02 per dose:', [round(f,2) for f in fracs])

json.dump({'exp_params':exp_params,
           'within_cell_CV':rat_cv_within, 'across_cell_CV':cv_across,
           'spearman_adjacent':rhos, 'auc_raw':auc_raw, 'auc_ratio':auc_ratio,
           'responding_fraction':fracs, 'n_cells':len(cells)},
          open('cascade_audit_results.json','w'), indent=1)
print('saved cascade_audit_results.json')
