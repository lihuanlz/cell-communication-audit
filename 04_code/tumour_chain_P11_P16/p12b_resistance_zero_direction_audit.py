# P12b verdict: do clinical BCR-ABL resistance mutations preserve catalytic efficacy
# (output-preserving direction) better than predicted-but-not-clinical mutations?
# Data: EVER paper (Liu/Pei/Lai, Commun Biol 2019, PMC6952392) Supplementary Data 1,
# sheet 'KM&Kcat' (absolute) and 'IC50_normalization' (WT=1).
# Pre-registered card: 预注册卡_P12_增益分布与耐药零方向_2026-09-28.md (FROZEN before values taken).
import math
from itertools import combinations
import csv, os

HERE = os.path.dirname(os.path.abspath(__file__))

# (mutation, kcat s-1, KM mM) from sheet KM&Kcat
kdata = {
 'WT':   (10.8, 0.0218),
 'L248R':(7.08, 0.288),
 'Y253H':(17.9, 0.0147),
 'E255V':(24.8, 0.029),
 'T315I':(46.9, 0.0186),
 'T315V':(86.9, 0.0658),
 'L248M':(205,  0.0353),
 'Y253C':(54.4, 0.0277),
 'V299L':(95.5, 0.0545),
 'V299M':(35.7, 0.0111),
 'I314V':(4.78, 0.0128),
 'M318I':(35.4, 0.00298),
 'Y320D':(119,  0.0761),
 'F382Y':(56,   0.00675),
}

# classification per pre-registered card
pc_main  = ['Y253H', 'E255V']                 # clinically observed, main test (T315I gatekeeper excluded)
gatekeeper = ['T315I']                        # steric gatekeeper, registered separately
pnc      = ['Y253C','V299M','I314V','F382Y','M318I','Y320D','L248M','V299L']  # predicted, not clinical (the 8 tested)
ambiguous = ['L248R','T315V']                 # rare/clinical-status ambiguous, sensitivity only

wt_eff = kdata['WT'][0]/kdata['WT'][1]
eff = {m: v[0]/v[1] for m, v in kdata.items()}
dlog = {m: abs(math.log10(eff[m]/wt_eff)) for m in eff if m != 'WT'}

def mwu_one_sided(group_small, group_large):
    """one-sided Mann-Whitney: H1 group_small < group_large. exact via rank enumeration."""
    allv = [('a', v) for v in group_small] + [('b', v) for v in group_large]
    n1, n2 = len(group_small), len(group_large)
    vals = sorted([v for v in group_small] + [v for v in group_large])
    # average ranks
    ranks = {}
    i = 0
    while i < len(vals):
        j = i
        while j+1 < len(vals) and vals[j+1] == vals[i]:
            j += 1
        r = (i + j) / 2 + 1
        for k in range(i, j+1):
            ranks[(vals[k], k)] = r
        i = j + 1
    rs = sorted(group_small)
    ri = sorted(group_large)
    # assign ranks by order
    seq = sorted([(v, 'a') for v in group_small] + [(v, 'b') for v in group_large])
    rk = []
    i = 0
    while i < len(seq):
        j = i
        while j+1 < len(seq) and seq[j+1][0] == seq[i][0]:
            j += 1
        r = (i + j)/2 + 1
        for k in range(i, j+1):
            rk.append((seq[k][1], r))
        i = j+1
    R1 = sum(r for g, r in rk if g == 'a')
    U1 = R1 - n1*(n1+1)/2
    # exact one-sided p: enumerate all choices of n1 ranks among N
    N = n1 + n2
    allranks = [r for _, r in rk]
    cnt = 0
    tot = 0
    for combo in combinations(range(N), n1):
        s = sum(allranks[i] for i in combo)
        tot += 1
        if s <= R1 + 1e-9:
            cnt += 1
    return U1, cnt/tot

pc_vals  = [dlog[m] for m in pc_main]
pnc_vals = [dlog[m] for m in pnc]

print('catalytic efficacy kcat/KM (WT=%.1f s-1 mM-1)' % wt_eff)
for m in ['Y253H','E255V','T315I','T315V','L248R','Y253C','V299M','I314V','F382Y','M318I','Y320D','L248M','V299L']:
    fold = eff[m]/wt_eff
    print('  %-6s eff=%8.1f  fold=%6.2f  |dlog10|=%.3f' % (m, eff[m], fold, dlog[m]))

pc_vals.sort(); pnc_vals.sort()
med = lambda v: (v[len(v)//2] if len(v)%2 else (v[len(v)//2-1]+v[len(v)//2])/2)
med_pc, med_pnc = med(pc_vals), med(pnc_vals)
fold_ratio = 10**(med_pnc - med_pc)
U, p_main = mwu_one_sided(pc_vals, pnc_vals)
print()
print('MAIN TEST (card rule 1): p-c median |dlog|=%.3f  p-nc median=%.3f  fold ratio=%.2f (criterion >2)' % (med_pc, med_pnc, fold_ratio))
print('  one-sided Mann-Whitney U=%.0f  exact p=%.4f (criterion p<0.05)' % (U, p_main))
verdict_main = (p_main < 0.05) and (fold_ratio > 2)
print('  => joint criterion:', 'PASS' if verdict_main else 'FAIL')

# sensitivity 1: include T315I in clinical group
pc2 = sorted(pc_vals + [dlog['T315I']])
U2, p2 = mwu_one_sided(pc2, pnc_vals)
print()
print('SENSITIVITY A (T315I included in clinical, n=3): medians %.3f vs %.3f, fold=%.2f, U=%.0f, p=%.4f'
      % (med(pc2), med_pnc, 10**(med_pnc-med(pc2)), U2, p2))

# sensitivity 2: also add ambiguous rare-clinical T315V (L248R excluded: 20-fold LOSS of efficacy, status disputed)
pc3 = sorted(pc2 + [dlog['T315V']])
U3, p3 = mwu_one_sided(pc3, pnc_vals)
print('SENSITIVITY B (+T315V, n=4): medians %.3f vs %.3f, fold=%.2f, U=%.0f, p=%.4f'
      % (med(pc3), med_pnc, 10**(med_pnc-med(pc3)), U3, p3))

# descriptive: fraction of p-nc above 6-fold efficacy change (paper's own statement)
big = [m for m in pnc if eff[m]/wt_eff > 6]
print()
print('p-nc mutants with >6-fold efficacy change:', big, '(paper claims 4 of 8)')
print('clinical group max fold change: %.2f (paper claims <4 for Y253H/E255V/T315I; T315I measured here = %.2f)'
      % (max(eff[m]/wt_eff for m in ['Y253H','E255V']), eff['T315I']/wt_eff))

# save CSV
with open(os.path.join(HERE, 'data', 'P12b_abl_catalytic_efficacy.csv'), 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['mutation','kcat_s-1','KM_mM','kcat_over_KM','fold_vs_WT','abs_dlog10','class'])
    for m, (kc, km) in kdata.items():
        if m == 'WT':
            cls = 'WT'
        elif m in pc_main: cls = 'clinical_main'
        elif m in gatekeeper: cls = 'clinical_gatekeeper_excluded'
        elif m in pnc: cls = 'predicted_not_clinical'
        else: cls = 'ambiguous_sensitivity'
        w.writerow([m, kc, km, round(kc/km,1), round(kc/km/wt_eff,3) if m!='WT' else 1.0,
                    round(dlog.get(m,0),4), cls])
print()
print('CSV saved: data/P12b_abl_catalytic_efficacy.csv')
