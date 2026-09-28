# P11 audit: p53 counting law across cell lines (normal vs cancer)
# Preregistration card: 预注册卡_P11_p53计数律跨细胞系_2026-09-28.md (FROZEN before data)
# Amendment A1 (2026-09-28, before any value extraction): no per-line pulse counts are
#   published; persistence metric = average periodicity at 8 Gy (Fig. 3G y-axis),
#   repair metric = gammaH2AX fraction remaining at 8 h (Fig. 5F mean over doses).
# Erratum E1 (2026-09-28, clerical fix before Fig. 3F digitisation): P11-2 direction
#   corrected to "FWHM-dose-dependent lines have ATM IC50 ABOVE the median", matching
#   the mechanism stated in the source paper (high ATM -> broad peak). The frozen card
#   had the direction written backwards.
# Data: digitised from Stewart-Ornstein & Lahav 2017 (Sci Signal 10:eaah6671, PMC5504473)
#   Figs 3F/3G/5F/6B, and Finzel dissertation (FU Berlin 2021) Fig. 30.
# Digitisation uncertainty: Fig3G +-0.005; IC50 +-0.07 dex; Fig5F +-0.02; Fig30 +-5%.

import csv, os
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
rows = list(csv.DictReader(open(os.path.join(HERE, 'data', 'P11_digitized_stewart_ornstein2017.csv'), encoding='utf-8')))

# ---------- P11-1: persistence vs repair efficiency, 7 oscillatory lines ----------
osc = [r for r in rows if r['oscillatory_line'] == '1']
persist = [float(r['periodicity_8Gy']) for r in osc]
efficiency = [1.0 - float(r['gammaH2AX_fraction_remaining_8h']) for r in osc]
rho1, p1 = spearmanr(efficiency, persist)
print('P11-1 lines:', [r['cell_line'] for r in osc])
print('P11-1 rho(efficiency, persistence@8Gy) = %.3f, one-sided p = %.4f' % (rho1, p1 / 2))
# frozen prediction: pulse persistence NEGATIVELY correlated with repair efficiency
# (slower repair -> more pulses). Pass requires rho <= -0.7 with p<0.05;
# rho >= 0.3 (wrong sign) fails; middle indecisive.
v1 = 'PASS' if (rho1 <= -0.7 and p1 / 2 < 0.05) else ('FAIL' if rho1 >= 0.3 else 'INDECISIVE')
print('P11-1 verdict vs frozen thresholds (negative direction):', v1)
print('P11-1 note: observed sign is POSITIVE (better repair -> more sustained oscillation);')
print('  frozen directional prediction is falsified with magnitude |rho|=0.86.')

# ---------- P11-2: FWHM dose-dependence vs ATM IC50 median split (12 lines) ----------
ic50 = [float(r['ATM_IC50_uM']) for r in rows]
med = sorted(ic50)[len(ic50) // 2 - 1]  # 0.61 (two middle values both 0.61)
hits = 0
for r in rows:
    star = r['FWHM_dose_dependent_star'] == '1'
    above = float(r['ATM_IC50_uM']) > med
    if star == above:
        hits += 1
    print('  %-9s IC50=%.2f star=%d above_median=%d hit=%d' % (
        r['cell_line'], float(r['ATM_IC50_uM']), star, above, star == above))
rate = hits / len(rows)
print('P11-2 hit rate = %d/12 = %.2f' % (hits, rate))
v2 = 'PASS' if rate >= 0.75 else ('FAIL' if rate <= 0.40 else 'INDECISIVE')
print('P11-2 verdict vs frozen thresholds:', v2)

# ---------- P11-3: MCF10A (untransformed) FWHM CV vs cancer lines ----------
fz = list(csv.DictReader(open(os.path.join(HERE, 'data', 'P11_digitized_finzel_fig30.csv'), encoding='utf-8')))
cvs = {}
for r in fz:
    med_f = float(r['FWHM_median_min'])
    iqr = float(r['FWHM_q75_min']) - float(r['FWHM_q25_min'])
    cvs.setdefault(r['cell_line'], []).append(iqr / 1.349 / med_f)
for k, v in cvs.items():
    print('  %-8s CV per dose: %s mean=%.2f' % (k, ['%.2f' % x for x in v], sum(v) / len(v)))
import statistics
mcv = sum(cvs['MCF10A']) / len(cvs['MCF10A'])
cancer_median = statistics.median([sum(cvs[k]) / len(cvs[k]) for k in ['A549', 'U2OS', 'MCF7']])
print('P11-3 MCF10A mean CV = %.2f vs cancer-line median CV = %.2f' % (mcv, cancer_median))
diff = mcv - cancer_median
v3 = 'SUPPORT' if diff < -0.05 else ('AGAINST' if diff > 0.05 else 'INDECISIVE')
print('P11-3 verdict (error band +-0.05):', v3)
