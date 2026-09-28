# P12a digitization: ERK-KTR baseline (DMSO) violin outlines, Fig 1B, JBC 2022 (PMC9358475).
# Reconstruct per-cell C/N distribution from violin width profile; compute CV.
import numpy as np
from PIL import Image
import os, csv

HERE = os.path.dirname(os.path.abspath(__file__))
img = np.array(Image.open(os.path.join(HERE, 'data', 'erk_fig1.jpg')).convert('L')).astype(float)

# (name, axis_x, bottom_tick_y=0.4, px_per_0.2, violin0 window)
panels = {
 'H1299 (NRAS-Q61K)':  dict(bot=233, step=27.6, x0=548, x1=602, ytop=88,  ybot=240),
 'HCT-116 (KRAS-G13D)':dict(bot=562, step=23.7, x0=548, x1=602, ytop=410, ybot=568),
 'SH-SY5Y (WT Ras/Raf)':dict(bot=894, step=23.7, x0=548, x1=602, ytop=720, ybot=900),
}

def digitize(p):
    ys, ws, vals = [], [], []
    for y in range(p['ytop'], p['ybot']):
        row = img[y, p['x0']:p['x1']]
        dark = np.where(row < 180)[0]
        if len(dark) >= 2:
            w = dark.max() - dark.min()
            if w >= 2:
                v = 0.4 + (p['bot'] - y) / p['step'] * 0.2
                ys.append(y); ws.append(w); vals.append(v)
    vals = np.array(vals); ws = np.array(ws, dtype=float)
    mean = np.average(vals, weights=ws)
    var = np.average((vals-mean)**2, weights=ws)
    std = np.sqrt(var)
    return vals, ws, mean, std, std/mean

res = {}
print('violin-width reconstruction of DMSO (baseline) per-cell C/N distribution:')
for name, p in panels.items():
    vals, ws, mean, std, cv = digitize(p)
    res[name] = (mean, std, cv)
    print('  %-20s rows=%3d  mean=%.3f  std=%.3f  CV=%.3f' % (name, len(vals), mean, std, cv))

names = list(res)
print()
print('pairwise |dCV|:')
for i in range(3):
    for j in range(i+1,3):
        d = abs(res[names[i]][2]-res[names[j]][2])
        print('  %s vs %s: %.3f (criterion >0.15)' % (names[i], names[j], d))

# permutation test on |dCV| is not meaningful with n=3 per group of one;
# card fallback: bootstrap over digitization uncertainty instead.
# digitization uncertainty: vary width threshold +/-40 gray levels and window +/-3 px
print()
print('robustness over thresholds/windows:')
cvs = {n: [] for n in names}
for thr in (140, 180, 220):
    for dx in (-3, 0, 3):
        for name, p in panels.items():
            q = dict(p); q['x0'] = p['x0']+dx; q['x1'] = p['x1']+dx
            ys, ws, vals = [], [], []
            for y in range(q['ytop'], q['ybot']):
                row = img[y, q['x0']:q['x1']]
                dark = np.where(row < thr)[0]
                if len(dark) >= 2:
                    w = dark.max()-dark.min()
                    if w >= 2:
                        ys.append(y); ws.append(w); vals.append(0.4+(q['bot']-y)/q['step']*0.2)
            vals=np.array(vals); ws=np.array(ws,float)
            m=np.average(vals,weights=ws); s=np.sqrt(np.average((vals-m)**2,weights=ws))
            cvs[name].append(s/m)
for n in names:
    arr=np.array(cvs[n])
    print('  %-20s CV=%.3f +/- %.3f (range %.3f-%.3f)' % (n, arr.mean(), arr.std(), arr.min(), arr.max()))

with open(os.path.join(HERE,'data','P12a_erk_baseline_cv.csv'),'w',newline='') as f:
    w=csv.writer(f); w.writerow(['cell_line','mutation','mean_CN','std_CN','CV'])
    mut={'H1299 (NRAS-Q61K)':'NRAS-Q61K','HCT-116 (KRAS-G13D)':'KRAS-G13D','SH-SY5Y (WT Ras/Raf)':'none (ALK-F1174L)'}
    for n,(m,s,c) in res.items():
        w.writerow([n.split(' ')[0], mut[n], round(m,4), round(s,4), round(c,4)])
print()
print('CSV saved: data/P12a_erk_baseline_cv.csv')
