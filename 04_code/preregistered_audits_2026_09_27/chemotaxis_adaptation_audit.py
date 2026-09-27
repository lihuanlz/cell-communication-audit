# Cascade stage-2 audit v2: baseline restoration in Moore 2024 FRET archive
# Data structure (empirically read): resp_data fields are (35 events x 20 samples);
# each event = one 40 s pulse onset, 10 samples (5 s) pre-onset OFF baseline +
# 10 samples (5 s) post-onset response window. Adaptation plateau inside the pulse
# is NOT exported; what is testable is baseline restoration across dose history:
# does the pre-pulse baseline depend on the PREVIOUS pulse's dose?
import scipy.io as sio, numpy as np, glob, json

files = sorted(glob.glob(r'03_细胞线3/公开数据/Moore2024_Chemotaxis_FRET/*.mat'))
events = []
ncells = 0
for f in files:
    try:
        rd = sio.loadmat(f)['reorgData']['resp_data'][0,0]
    except Exception:
        continue
    for ci in range(rd.shape[1]):
        c = rd[0,ci]
        try:
            s = c['s'].astype(int); a = c['a']
        except Exception:
            continue
        if s.shape != (35, 20): continue
        ncells += 1
        for i in range(35):
            dose_pre = int(round(np.median(s[i, :10])))   # should be 0 (OFF)
            dose = int(round(np.median(s[i, 10:])))        # pulse dose
            b = np.nanmedian(a[i, :10])
            post = a[i, 10:] - b
            r = post[np.nanargmax(np.abs(post))]
            prev_dose = int(round(np.median(s[i-1, 10:]))) if i > 0 else None
            events.append(dict(ev=i, dose=dose, prev=prev_dose,
                               baseline=float(b), resp=float(r)))
print('cells:', ncells, 'events:', len(events))

# 1. baseline by previous dose (history dependence = incomplete adaptation)
import numpy as np
by_prev = {}
for e in events:
    if e['prev'] is None: continue
    by_prev.setdefault(e['prev'], []).append(e['baseline'])
print('\nbaseline (pre-pulse, OFF window) by PREVIOUS dose:')
for d in sorted(by_prev):
    v = np.array(by_prev[d])
    print(f'  prev {d:3d} uM: n={len(v):5d}  baseline median {np.median(v):+.4f}  IQR [{np.percentile(v,25):+.4f},{np.percentile(v,75):+.4f}]')

# 2. response amplitude by current dose (sanity: monotone dose-response)
by_dose = {}
for e in events:
    if e['dose'] == 0: continue
    by_dose.setdefault(e['dose'], []).append(e['resp'])
print('\nresponse amplitude by CURRENT dose:')
for d in sorted(by_dose):
    v = np.array(by_dose[d])
    print(f'  dose {d:3d} uM: n={len(v):5d}  resp median {np.median(v):+.4f}')

# 3. quantitative: baseline shift after high (40/80) vs low (2/10) previous dose
hi = np.array([e['baseline'] for e in events if e['prev'] in (40, 80)])
lo = np.array([e['baseline'] for e in events if e['prev'] in (2, 10)])
resp_all = np.abs([e['resp'] for e in events if e['dose'] > 0])
print(f'\nbaseline after high-dose pulse: {np.median(hi):+.4f} (n={len(hi)})')
print(f'baseline after low-dose pulse:  {np.median(lo):+.4f} (n={len(lo)})')
print(f'difference: {np.median(hi)-np.median(lo):+.4f}  vs median |resp| {np.median(resp_all):.4f}')
json.dump(dict(ncells=ncells, events=len(events),
               baseline_by_prev={str(d): [float(np.median(by_prev[d])), len(by_prev[d])] for d in by_prev},
               baseline_after_hi=float(np.median(hi)), baseline_after_lo=float(np.median(lo)),
               median_abs_resp=float(np.median(resp_all))),
          open(r'05_主线纲领与设计/级联线/chemotaxis_adaptation_audit.json', 'w'), indent=1)
print('saved json')
