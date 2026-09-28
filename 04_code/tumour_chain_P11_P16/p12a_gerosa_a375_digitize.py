# P12a extension (Amendment A3): digitize Gerosa 2020 Fig 1H DMSO violin (A375, BRAF-V600E)
# as 4th cell line. Same method as p12a_erk_violin_digitize.py.
import numpy as np
from PIL import Image
import os, csv

HERE = os.path.dirname(os.path.abspath(__file__))
img = np.array(Image.open(os.path.join(HERE, 'data', 'gerosa_fig1.jpg')).convert('L')).astype(float)

# locate axis: vertical dark line near x 540-560, y 480-640
best_x, best_len = None, 0
for x in range(535, 570):
    dark = np.sum(img[480:640, x] < 100)
    if dark > best_len: best_len, best_x = dark, x
ax = best_x
tick_ys = []
for y in range(480, 645):
    if np.mean(img[y, ax-6:ax-1] < 100) > 0.8:
        tick_ys.append(y)
ticks = []
for y in tick_ys:
    if not ticks or y - ticks[-1][-1] > 2: ticks.append([y])
    else: ticks[-1].append(y)
ticks = [int(np.mean(t)) for t in ticks]
print('axis_x =', ax, 'ticks =', ticks)
# expected ticks: 0.5, 1.0, 1.5, 2.0 (bottom to top)
if len(ticks) == 4:
    vals = [0.5, 1.0, 1.5, 2.0]  # bottom-up
    py = np.polyfit(ticks, vals, 1)  # value = py[0]*y + py[1]
else:
    raise SystemExit('unexpected tick count, inspect manually')
print('calibration: value = %.5f*y + %.3f' % (py[0], py[1]))

def digitize(x0, x1, thr=210):
    vs, ws = [], []
    for y in range(480, 645):
        row = img[y, x0:x1]
        dark = np.where(row < thr)[0]
        if len(dark) >= 2:
            w = dark.max() - dark.min()
            if w >= 3:
                vs.append(py[0]*y + py[1]); ws.append(w)
    vs = np.array(vs); ws = np.array(ws, float)
    m = np.average(vs, weights=ws)
    s = np.sqrt(np.average((vs-m)**2, weights=ws))
    return m, s, s/m

# DMSO violin: middle, grey. scan candidate windows
for x0, x1 in [(560,580),(565,585),(570,590),(558,592)]:
    m, s, cv = digitize(x0, x1)
    print('window %d-%d: mean=%.3f std=%.3f CV=%.3f' % (x0, x1, m, s, cv))
