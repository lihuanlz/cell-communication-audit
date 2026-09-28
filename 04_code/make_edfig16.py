# -*- coding: utf-8 -*-
"""EDFig16 v2: the cancer chain evidence map (P11-P16). Fixed spacing."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(sys.executable).parent.parent.parent))
from daimon_runtime import setup_plot
setup_plot()
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

ROOT = Path(r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911\05_主线纲领与设计\论文_细胞通讯审计_2026-09-23")

C_OK = '#2e7d32'; C_PART = '#f9a825'; C_FAIL = '#c62828'
C_BOUND = '#546e7a'; C_FROZEN = '#6a1b9a'

fig, ax = plt.subplots(figsize=(14.5, 7.6))
ax.set_xlim(0, 14.5); ax.set_ylim(0, 7.6); ax.axis('off')

def box(x, y, w, h, color, title, lines, sub=''):
    ax.add_patch(mpatches.FancyBboxPatch((x, y), w, h,
                 boxstyle='round,pad=0.08', fc='white', ec=color, lw=2.2))
    ax.text(x + w/2, y + h - 0.30, title, ha='center', va='top',
            fontsize=10.5, fontweight='bold', color=color)
    yy = y + h - 0.82
    for ln in lines:
        ax.text(x + w/2, yy, ln, ha='center', va='top', fontsize=8.0, color='#222')
        yy -= 0.38
    if sub:
        ax.text(x + w/2, y + 0.22, sub, ha='center', va='bottom',
                fontsize=7.5, style='italic', color=color)

def arrow(x1, y1, x2, y2):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='-|>', lw=1.8, color='#555'))

Y1, H1, Y2, H2 = 4.4, 2.75, 0.75, 2.55
box(0.3, Y1, 2.5, H1, C_FAIL, 'P11 | regime switch',
    ['prediction reversed: slow repair', 'does not count faster; it switches',
     'p53 pulses -> sustained regime', 'rho = +0.857, p = 0.0068 (7 lines)'],
    'falsified; claim contracted')
box(3.15, Y1, 2.5, H1, C_OK, 'P14 | downstream misread',
    ['two decoder classes coexist:', 'arrest genes amplified by sustained',
     'shape (CDKN1A +2.3 bits, AUC-', 'corrected); repair genes count AUC'],
    'confirmed (registered correction)')
box(6.0, Y1, 2.5, H1, C_PART, 'P13 | channel capacity',
    ['I(dose; ERK): 0.43 bit WT falls to', '0.09 (Q61R) and 0.01 (BRAF-V600E);',
     'graded with allele strength;', 'G12V unaffected; 5 bin sets + KSG'],
    'partial; graded by allele')
box(8.85, Y1, 2.5, H1, C_OK, 'P12b | resistance walks',
    ['ABL: clinical mutants 2.94x closer', 'to catalytic zero direction (p = 0.033);',
     'EGFR DMS: 2.2x (p = 0.002);', 'L718X counterexample branch found'],
    'confirmed; two branches')
box(11.7, Y1, 2.5, H1, C_BOUND, 'P15 | tissue edge',
    ['bulk RPPA: no coupling even in', 'wild type (|r| <= 0.18 across four',
     'cancer types, ~1,000 tumours);', 'mean pERK up in KRAS-G12 LUAD'],
    'falsified at bulk resolution')
for x0 in (2.80, 5.65, 8.50, 11.35):
    arrow(x0, Y1 + H1/2, x0 + 0.35, Y1 + H1/2)

box(0.3, Y2, 3.2, H2, C_FROZEN, 'P16 | frozen prediction',
    ['KRAS-G13D: delta I in [-0.25, -0.05]', 'bits vs WT; ordering G12V > G13D',
     '> G12D; frozen 2026-09-28'],
    'awaiting data')
box(3.8, Y2, 3.2, H2, C_BOUND, 'P12a | width unchanged',
    ['ERK gain distribution width', 'untouched by mutation:',
     'CV 0.339 / 0.353 / 0.353', '|dCV| <= 0.014'],
    'falsified; informative null')
box(7.3, Y2, 3.2, H2, C_PART, 'Fig. 6 side evidence',
    ['ERK deactivation half-time:', 'WT 6.2 min -> 2.8 min (Q61R);',
     'feedback/phosphatase', 'renormalisation direction'],
    'exploratory (6 reps per line)')
box(10.8, Y2, 3.4, H2, C_OK, 'Thesis',
    ['cancer is not a cut wire:', 'a cheaper language, read off',
     'the old price list by', 'downstream decoders'],
    'chain summary')

ax.text(7.25, 7.4, 'Extended Data Fig. 16 | The cancer chain: hijacked, not interrupted (P11-P16)',
        ha='center', fontsize=12.5, fontweight='bold')
ax.text(7.25, 0.18, 'Every link pre-registered before data access; verdicts quoted from frozen cards (SI S19). Colours: green = confirmed, amber = partial/graded,\n'
        'red = falsified (claim contracted), grey = boundary/null, purple = frozen prospective.',
        ha='center', va='bottom', fontsize=7.8, color='#555')

for d in (ROOT / 'submission_package' / '02_figures', ROOT / 'figures_svg'):
    d.mkdir(exist_ok=True)
    fig.savefig(d / 'EDFig16_cancer_chain.png', dpi=300, bbox_inches='tight')
    fig.savefig(d / 'EDFig16_cancer_chain.svg', bbox_inches='tight')
print('saved EDFig16 v2 png+svg')
