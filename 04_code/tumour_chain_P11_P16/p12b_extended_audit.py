# P12b extended (Amendment A1): expanded clinical group + EGFR second system.
# Sources: EVER all_data sheet (Commun Biol 2019, PMC6952392, canonical);
#          Yun et al. 2008 PNAS (PMC2538882) Table 2 for EGFR.
import math
from itertools import combinations
import csv, os

HERE = os.path.dirname(os.path.abspath(__file__))

# ABL from EVER all_data (kcat s-1, KM uM), canonical
abl = {
 'WT':   (10.8, 21.8),
 'Y253H':(19.5, 19.2),   # p-c
 'E255V':(24.8, 29.0),   # p-c
 'T315I':(46.9, 28.6),   # p-c gatekeeper (steric), excluded
 'Y253C':(54.4, 27.7),   # p-nc
 'V299M':(35.7, 11.1),   # p-nc
 'I314V':(4.78, 12.8),   # p-nc
 'F382Y':(56.0, 6.75),   # p-nc
 'M318I':(35.4, 2.98),   # p-nc
 'Y320D':(119,  76.1),   # p-nc
 'L248M':(205,  35.3),   # p-nc
 'V299L':(95.5, 54.5),   # p-nc
 'I315M':(17.4, 24.2),   # p-c (ponatinib, from T315I background) - Amendment A1
}
# EGFR from Yun 2008 Table 2 (kcat s-1, KM uM, eff in uM-1 s-1 given)
egfr = {
 'WT':          (0.026, 5.2),
 'T790M':       (0.137, 5.9),    # clinical resistance, ATP-affinity mechanism - Amendment A1
 'L858R':       (1.484, 148.0),  # driver, NOT resistance - excluded from groups
 'L858R/T790M': (0.456, 8.4),    # clinical resistance on driver background
}

pc_abl  = ['Y253H', 'E255V', 'I315M']
pnc     = ['Y253C','V299M','I314V','F382Y','M318I','Y320D','L248M','V299L']
gate    = ['T315I']

wt_abl = abl['WT'][0]/abl['WT'][1]
dlog_abl = {m: abs(math.log10((v[0]/v[1])/wt_abl)) for m, v in abl.items() if m != 'WT'}

wt_egfr = egfr['WT'][0]/egfr['WT'][1]
eff_egfr = {m: v[0]/v[1] for m, v in egfr.items()}
dlog_egfr_vsWT = {m: abs(math.log10(eff_egfr[m]/wt_egfr)) for m in egfr if m != 'WT'}
dlog_t790m_vs_driver = abs(math.log10(eff_egfr['L858R/T790M']/eff_egfr['L858R']))

def mwu_exact(g1, g2):
    n1, n2 = len(g1), len(g2)
    seq = sorted([(v,'a') for v in g1]+[(v,'b') for v in g2])
    rk=[]; i=0
    while i < len(seq):
        j=i
        while j+1<len(seq) and seq[j+1][0]==seq[i][0]: j+=1
        r=(i+j)/2+1
        for k in range(i,j+1): rk.append((seq[k][1],r))
        i=j+1
    R1=sum(r for g,r in rk if g=='a'); U1=R1-n1*(n1+1)/2
    allr=[r for _,r in rk]; cnt=tot=0
    for combo in combinations(range(n1+n2), n1):
        tot+=1
        if sum(allr[i] for i in combo) <= R1+1e-9: cnt+=1
    return U1, cnt/tot

med = lambda v: (sorted(v)[len(v)//2] if len(v)%2 else (sorted(v)[len(v)//2-1]+sorted(v)[len(v)//2])/2)

print('=== ABL (all_data canonical) |dlog10 kcat/KM| ===')
for m in pc_abl+gate+pnc:
    fold=(abl[m][0]/abl[m][1])/wt_abl
    print('  %-6s fold=%6.2f  |dlog|=%.3f' % (m, fold, dlog_abl[m]))

print()
print('=== EGFR (Yun 2008) ===')
for m in ['T790M','L858R','L858R/T790M']:
    print('  %-12s eff=%.5f  fold vs WT=%.2f  |dlog|=%.3f' % (m, eff_egfr[m], eff_egfr[m]/wt_egfr, dlog_egfr_vsWT[m]))
print('  L858R/T790M vs L858R driver background: fold=%.2f |dlog|=%.3f' % (eff_egfr['L858R/T790M']/eff_egfr['L858R'], dlog_t790m_vs_driver))

# main expanded test: ABL clinical (n=3, gatekeeper excl) + EGFR T790M (vs WT caliber)
pc_ext = [dlog_abl[m] for m in pc_abl] + [dlog_egfr_vsWT['T790M']]
pnc_v  = [dlog_abl[m] for m in pnc]
U,p = mwu_exact(pc_ext, pnc_v)
print()
print('EXPANDED MAIN TEST (Amendment A1): clinical n=%d %s' % (len(pc_ext), [round(v,3) for v in pc_ext]))
print('  medians: clinical %.3f vs p-nc %.3f ; fold ratio %.2f' % (med(pc_ext), med(pnc_v), 10**(med(pnc_v)-med(pc_ext))))
print('  one-sided MW U=%.0f exact p=%.4f (criteria p<0.05 and fold>2)' % (U,p))
print('  => %s' % ('PASS' if (p<0.05 and 10**(med(pnc_v)-med(pc_ext))>2) else 'FAIL'))

# sensitivity: EGFR vs driver background instead
pc_ext2 = [dlog_abl[m] for m in pc_abl] + [dlog_t790m_vs_driver]
U2,p2 = mwu_exact(pc_ext2, pnc_v)
print('SENSITIVITY (T790M vs driver background, |dlog|=%.3f): U=%.0f p=%.4f fold=%.2f'
      % (dlog_t790m_vs_driver, U2, p2, 10**(med(pnc_v)-med(pc_ext2))))

# sensitivity: all_data vs KM&Kcat sheet discrepancy for Y253H/T315I
wt2 = 10.8/0.0218/1000.0  # KM&Kcat sheet: KM in mM -> convert to uM-1s-1 basis consistent? use fold only
alt = {'Y253H': abs(math.log10((17.9/14.7)/(10.8/21.8))),
       'T315I': abs(math.log10((46.9/18.6)/(10.8/21.8)))}
print('DATA-QC sensitivity (KM&Kcat sheet): Y253H |dlog|=%.3f (all_data %.3f), T315I %.3f (%.3f)'
      % (alt['Y253H'], dlog_abl['Y253H'], alt['T315I'], dlog_abl['T315I']))

with open(os.path.join(HERE,'data','P12b_extended_clinical_expansion.csv'),'w',newline='') as f:
    w=csv.writer(f); w.writerow(['system','mutation','kcat','KM_uM','fold_vs_WT','abs_dlog10','class','source'])
    for m,v in abl.items():
        if m=='WT': continue
        cls = 'clinical' if m in pc_abl else ('gatekeeper_excluded' if m in gate else 'predicted_not_clinical')
        w.writerow(['ABL',m,v[0],v[1],round((v[0]/v[1])/wt_abl,3),round(dlog_abl[m],4),cls,'EVER all_data PMC6952392'])
    for m in ['T790M','L858R/T790M']:
        w.writerow(['EGFR',m,egfr[m][0],egfr[m][1],round(eff_egfr[m]/wt_egfr,3),round(dlog_egfr_vsWT[m],4),'clinical','Yun2008 PMC2538882 Table2'])
print()
print('CSV saved: data/P12b_extended_clinical_expansion.csv')
