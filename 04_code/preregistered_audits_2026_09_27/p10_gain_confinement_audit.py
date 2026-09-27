# P10 audit: gain-distribution confinement across systems (frozen prereg 2026-09-27)
# ERK (Chavez-Abiega 2022) + NF-kB (Wang 2022), identical pipeline, SEED=20260927
import numpy as np, pandas as pd, os, json
from scipy.io import loadmat
from scipy.signal import find_peaks
from scipy.stats import spearmanr

SEED = 20260927
rng = np.random.default_rng(SEED)
ROOT = r'D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911'
PROM, MINDIST, MINP = 0.10, 8, 3      # frozen

def cell_pulse_amps(trace):
    base = np.nanpercentile(trace, 20)
    pk, _ = find_peaks(trace, prominence=PROM, distance=MINDIST)
    amps = trace[pk] - base
    amps = amps[np.isfinite(amps) & (amps > 0)]
    return amps

def stats_block(amps_list):
    """amps_list: list of per-cell amplitude arrays (>=MINP pulses)."""
    wc_i, m_i, fh, sh = [], [], [], []
    for a in amps_list:
        m = a.mean()
        if m <= 0: continue
        wc_i.append(a.std(ddof=1)/m); m_i.append(m)
        h = len(a)//2
        if h >= 1 and len(a)-h >= 1:
            fh.append(a[:h].mean()); sh.append(a[h:].mean())
    wc_i = np.array(wc_i); m_i = np.array(m_i)
    across_cv = m_i.std(ddof=1)/m_i.mean()
    R = across_cv/np.median(wc_i)
    rho = spearmanr(fh, sh).statistic
    N = len(wc_i)
    idx = rng.integers(0, N, size=(2000, N))
    bs = m_i[idx]; bw = wc_i[idx]
    Rs = (bs.std(axis=1, ddof=1)/bs.mean(axis=1))/np.median(bw, axis=1)
    ci = np.percentile(Rs, [2.5, 97.5])
    return dict(n_cells=N, within_cv_median=float(np.median(wc_i)),
                across_cv=float(across_cv), R=float(R), R_ci95=[float(ci[0]),float(ci[1])],
                rho_splithalf=float(rho))

res = {}

# ---------- ERK (Chavez-Abiega) ----------
erk_dir = os.path.join(ROOT, '01_细胞线', '公开数据', 'ERK_Akt_GPCR_ChavezAbiega2022', 'Figure_1', 'Data')
group_Rs = []
all_report = {}
for f in ['All_Ex_Histamine_DMSO.csv', 'All_Ex_UK_DMSO.csv']:
    df = pd.read_csv(os.path.join(erk_dir, f))
    lig = f.split('_')[2]
    per_cond_R = []
    for cond, g0 in df.groupby('Condition'):
        amps_list = []
        for obj, g1 in g0.groupby('Unique_Object'):
            g1 = g1.sort_values('Time_in_min')
            tr = g1['CN_ERK'].values.astype(float)
            if len(tr) < 30: continue
            a = cell_pulse_amps(tr)
            if len(a) >= MINP: amps_list.append(a)
        if len(amps_list) >= 30:
            s = stats_block(amps_list)
            per_cond_R.append(s['R'])
            all_report[f'{lig}@{cond}'] = s
    if per_cond_R:
        group_Rs.append(float(np.median(per_cond_R)))
res['ERK'] = dict(per_condition=all_report, median_R_across_conditions=group_Rs)
print('ERK conditions:', len(all_report), 'median R:', group_Rs)

# ---------- NF-kB (Wang sequential) ----------
m = loadmat(os.path.join(ROOT, '03_细胞线3', '公开数据', 'Wang2022_Sequential_NFkB', 'scmat_sequentialstim.mat'))
sc = np.asarray(m['scmatcomb_norm'], dtype=float)   # 11267 x 172
amps_list = []
for i in range(sc.shape[0]):
    tr = sc[i]
    if np.isfinite(tr).mean() < 0.9: continue
    a = cell_pulse_amps(tr)
    if len(a) >= MINP: amps_list.append(a)
res['NFkB'] = stats_block(amps_list)
print('NFkB:', res['NFkB'])

out = os.path.join(ROOT, '05_主线纲领与设计', '级联线', 'p10_gain_confinement_results.json')
json.dump(res, open(out, 'w'), indent=1)
print('saved', out)
