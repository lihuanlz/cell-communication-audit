# P12b final (Amendments A1+A2): clinical n=5 vs predicted-not-clinical n=8.
# ABL: EVER all_data (PMC6952392). EGFR T790M: vs driver background, median of 3 measurements
#   (Yun 2008 PMC2538882: 5.41x vs L858R; Kasuga rs-6573462 Table 1: 4.50x vs ex19del, 2.10x vs L858R).
# EGFR C797S: median of 3 background-relative measurements (1.69x, 1.15x, 1/0.77=1.30x) = 1.30x? -> see code.
import math
from itertools import combinations
import csv, os

HERE = os.path.dirname(os.path.abspath(__file__))

abl = {
 'WT':   (10.8, 21.8),
 'Y253H':(19.5, 19.2), 'E255V':(24.8, 29.0), 'T315I':(46.9, 28.6),
 'Y253C':(54.4, 27.7), 'V299M':(35.7, 11.1), 'I314V':(4.78, 12.8), 'F382Y':(56.0, 6.75),
 'M318I':(35.4, 2.98), 'Y320D':(119, 76.1), 'L248M':(205, 35.3), 'V299L':(95.5, 54.5),
 'I315M':(17.4, 24.2),
}
wt_abl = abl['WT'][0]/abl['WT'][1]
dlog_abl = {m: abs(math.log10((v[0]/v[1])/wt_abl)) for m, v in abl.items() if m != 'WT'}

t790m_folds = [5.41, 4.50, 2.10]          # vs driver background, 3 measurements
c797s_folds = [1.69, 1.15, 1/0.77]        # vs driver background, 3 measurements (0.77 inverted)
t790m = sorted(t790m_folds)[1]; c797s = sorted(c797s_folds)[1]
d_t790m = math.log10(t790m); d_c797s = math.log10(c797s)
print('T790M vs driver: folds %s -> median %.2f, |dlog|=%.3f' % (t790m_folds, t790m, d_t790m))
print('C797S vs driver: folds %s -> median %.2f, |dlog|=%.3f' % ([round(f,2) for f in c797s_folds], c797s, d_c797s))

pc  = [('ABL Y253H', dlog_abl['Y253H']), ('ABL E255V', dlog_abl['E255V']),
       ('ABL I315M', dlog_abl['I315M']), ('EGFR T790M', d_t790m), ('EGFR C797S', d_c797s)]
pnc = [('ABL '+m, dlog_abl[m]) for m in ['Y253C','V299M','I314V','F382Y','M318I','Y320D','L248M','V299L']]

def mwu_exact(g1, g2):
    n1,n2=len(g1),len(g2)
    seq=sorted([(v,'a') for v in g1]+[(v,'b') for v in g2])
    rk=[];i=0
    while i<len(seq):
        j=i
        while j+1<len(seq) and seq[j+1][0]==seq[i][0]: j+=1
        r=(i+j)/2+1
        for k in range(i,j+1): rk.append((seq[k][1],r))
        i=j+1
    R1=sum(r for g,r in rk if g=='a'); U1=R1-n1*(n1+1)/2
    allr=[r for _,r in rk]; cnt=tot=0
    for combo in combinations(range(n1+n2),n1):
        tot+=1
        if sum(allr[i] for i in combo)<=R1+1e-9: cnt+=1
    return U1,cnt/tot

med=lambda v:(sorted(v)[len(v)//2] if len(v)%2 else (sorted(v)[len(v)//2-1]+sorted(v)[len(v)//2])/2)
pc_v=[v for _,v in pc]; pnc_v=[v for _,v in pnc]
U,p=mwu_exact(pc_v,pnc_v)
fr=10**(med(pnc_v)-med(pc_v))
print()
print('FINAL TEST (A1+A2): clinical n=5 %s' % [round(v,3) for v in sorted(pc_v)])
print('  p-nc n=8 %s' % [round(v,3) for v in sorted(pnc_v)])
print('  medians %.3f vs %.3f ; fold ratio %.2f' % (med(pc_v),med(pnc_v),fr))
print('  one-sided MW U=%.0f exact p=%.4f' % (U,p))
print('  => %s' % ('PASS (p<0.05 and fold>2)' if (p<0.05 and fr>2) else 'FAIL'))

with open(os.path.join(HERE,'data','P12b_final_expansion.csv'),'w',newline='') as f:
    w=csv.writer(f); w.writerow(['mutation','abs_dlog10_catEff_vs_ref','class','note'])
    for m,v in pc: w.writerow([m,round(v,4),'clinical','EGFR entries vs driver background'])
    for m,v in pnc: w.writerow([m,round(v,4),'predicted_not_clinical','EVER all_data'])
    w.writerow(['ABL T315I',round(dlog_abl['T315I'],4),'gatekeeper_excluded','steric'])
print('CSV saved: data/P12b_final_expansion.csv')
