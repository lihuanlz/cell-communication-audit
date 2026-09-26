# -*- coding: utf-8 -*-
"""
Code 11: discrimination map v0.9 (model-side + experiment-side placements)
==========================================
v0.9 update of code 3: experimental points upgraded from semi-quantitative to quantitative values audited from public data
(Batchelor 2011 source figure digitized; Albeck 2013 source figure digitized; Tay 2010 main-text quantitative sentence + SI readouts),
with new ERK / NF-κB experimental placements and the NF-κB model-experiment gap arrow.

Sources of the placement values:
  Experiment side (star/diamond markers):
    p53/γ-IR (Lahav 2004)   α≈0.03 (upper limit of the reported statistic), ν≈2
    p53/NCS (Mönke 2017)    α≈0.05 (upper limit of the reported statistic), ν≈1.8
    p53/NCS (Batchelor 2011, source figure) α=0.00±0.28 (plotted at 0.008 on the log axis), ν≈2 (count increases with DSB)
    p53/UV  (Batchelor 2011, source figure) α=+0.99±0.32, ν≈0 (single pulse)
    ERK     (Albeck 2013, source figure)    |α|≲0.05 (plotted at 0.02), ν≈+10/decade (over the 10→50 pg/ml range,
                                    example trajectory 8→15 pulses/24h; saturation at 200 pg/ml) — count is the main channel capacity
    NF-κB   (Tay 2010, main-text sentence+SI)  α=0.151, ν≈+0.6/decade (1.7→4.0 events / 0.01→100 ng/ml)
  Model side (square/circle/triangle markers):
    NF model near onset (8.76, 0.02); NF model far from onset (0.7, 0.05)
    NF-κB Krishna soft mode (0.714, 0.05) — arrow to the Tay experimental point = model-experiment gap (pending)
    FHN slow drive (0.014, 4.05); Morris–Lecar fast drive (0.26, 2.0)
    PdPC positive-feedback switch (bistable static decision element, ν=0) placed at the right edge

Run: python3 代码11_判别图_v09.py
"""
import numpy as np
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(9.5, 7.5))

# ---------- regions ----------
ax.axvspan(0.3, 15, alpha=0.08, color='red')      # NF simulated region
ax.axvspan(5e-4, 0.15, alpha=0.08, color='blue')  # EXC digital region
ax.axvspan(0.15, 0.3, alpha=0.12, color='green')  # leaky-digital band
ax.text(2.4, 9.6, 'NF family (simulated)\n$\\alpha\\gtrsim0.5$\n(subcritical: +hysteresis)',
        color='darkred', fontsize=10, ha='center')
ax.text(0.006, 9.6, 'EXC family (digital)\n$\\alpha\\lesssim0.15$', color='darkblue', fontsize=10, ha='center')
ax.text(0.212, 9.6, 'leaky-digital band', color='darkgreen', fontsize=8.5, ha='center')
ax.axvline(0.22, color='gray', ls='--', lw=1)

# ---------- model points ----------
ax.plot([8.76], [0.02], 's', color='red', ms=9, label='model: NF near onset (α→∞)')
ax.plot([0.7], [0.05], 's', color='red', ms=9, mfc='none', label='model: NF far from onset')
ax.plot([0.714], [0.05], 's', color='darkred', ms=10, mfc='none', mew=2,
        label='model: NF-κB Krishna soft mode (α=0.714)')
ax.plot([0.014], [4.05], 'o', color='blue', ms=9, label='model: FHN slow drive (purely digital)')
ax.plot([0.26], [2.0], '^', color='green', ms=9, label='model: Morris–Lecar fast drive (leaky)')
ax.plot([11], [0.02], 'p', color='darkred', ms=10, label='model: PdPC positive-feedback switch (bistable)')

# ---------- experimental points (v0.9 public-data audit) ----------
ax.errorbar([0.03], [2.0], xerr=[0.02], yerr=[0.5], fmt='*', color='purple',
            ms=17, capsize=4, label='experiment: p53/γ-IR (Lahav 2004, reported statistic)')
ax.errorbar([0.05], [1.8], xerr=[0.04], yerr=[0.5], fmt='*', color='magenta',
            ms=15, capsize=4, label='experiment: p53/NCS (Mönke 2017, reported statistic)')
ax.errorbar([0.008], [2.2], yerr=[0.5], fmt='*', color='indigo', ms=19, capsize=4,
            label='experiment: p53/NCS (Batchelor 2011 source figure, α=0.00±0.28)')
ax.errorbar([0.99], [0.05], xerr=[0.32], yerr=[0.05], fmt='D', color='darkorange',
            ms=10, capsize=4, label='experiment: p53/UV (Batchelor 2011 source figure, α=+0.99±0.32)')
ax.errorbar([0.02], [10.0], yerr=[2.5], fmt='*', color='royalblue', ms=19, capsize=4,
            label='experiment: ERK (Albeck 2013 source figure, |α|≲0.05, ν≈+10/dec then saturation)')
ax.errorbar([0.151], [0.6], xerr=[0.05], yerr=[0.3], fmt='*', color='teal',
            ms=17, capsize=4, label='experiment: NF-κB (Tay 2010 main-text sentence+SI, α=0.151, ν≈+0.6/dec)')

# ---------- model-experiment gap arrow ----------
ax.annotate('', xy=(0.20, 0.55), xytext=(0.66, 0.08),
            arrowprops=dict(arrowstyle='->', color='darkred', lw=1.8, ls='-'))
ax.text(0.42, 0.9, 'NF-κB model-experiment gap\n(structural: parameter unreachable — code 12)', color='darkred', fontsize=8.5, ha='center')

ax.annotate('built-in control: same protein, same cell line\nα covaries strictly from 0 to 1', (0.35, 1.5), fontsize=8.5, color='darkorange')
ax.annotate('v0.8 prediction α≳0.5 falsified:\nERK is actually fixed-amplitude + frequency-modulated (digital side)', (0.06, 7.6), fontsize=8.5, color='royalblue')

ax.set_xscale('log'); ax.set_xlim(5e-4, 15); ax.set_ylim(-0.4, 11.5)
ax.set_xlabel('$\\alpha = \\partial\\ln A/\\partial\\ln D$ (amplitude-dose sensitivity)')
ax.set_ylabel('$\\nu$ (pulse-count increment per tenfold dose)')
ax.set_title('Discrimination map v0.9: model-side + experiment-side placements (public source figures / main-text quantification)', fontsize=11)
ax.legend(fontsize=7.8, loc='center right', framealpha=0.95)
plt.tight_layout()
plt.savefig('判别图_v09_实验落点.png', dpi=200)
print('Saved 判别图_v09_实验落点.png')
