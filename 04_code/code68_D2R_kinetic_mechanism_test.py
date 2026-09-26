# -*- coding: utf-8 -*-
"""
代码68_D2R动力学机制检验.py  v1.0.0
============================================================
Cell line 4 · deepening round: kinetic-mechanism adjudication of the D2R breakage

Question: Code 63 ruled that the three high-affinity partial agonists bifeprunox/aripiprazole/cariprazine
show a −1.2…−3.5 dex breakage between binding state (pKi≈8.9–10.4) and functional state (pKA_func≈6.9–7.8),
and that the bias factor reverses over time (bifeprunox vs ropinirole −0.42 → +1.93 dex).
Klein Herenbrink 2016 itself proposed a qualitative "kinetic context" explanation (slow-dissociating ligands not at equilibrium).
This code makes the quantitative adjudication: using the kon/koff measured in that paper's SI Tables 5/6, construct the non-equilibrium occupancy
  ρ(L,t) = [kon·L/(kon·L+koff)]·(1−e^{−(kon·L+koff)·t})
substitute into the operational model (n=1, E = Em·τρ/(1+ρ(τ−1))), with τ inverted from the 90' Emax (Esys=100,
ropinirole full-agonist reference), generate apparent concentration–response curves and extract apparent pEC50(t),
and reconcile against the measured pEC50 trajectories at 8 time points of SI Table 7. **Zero free parameters** (Kd not fitted,
kinetic measurements used directly).

Ruling logic: non-equilibrium can only press apparent potency temporarily below its equilibrium value; as t→∞ the model floor is
pKd+log10(1+τ) ≥ pKd. If the observed pEC50 lies more than 1 dex below that floor at any time,
the kinetic explanation fails quantitatively and the breakage survives.

Data (Klein Herenbrink 2016, Nat Commun 7:10842, SI):
  Table 5: [3H]spiperone tracer competition binding kinetics (measurable only for 3 slow ligands)
  Table 6: PPHT-red fluorescent tracer Tag-lite binding kinetics (7 ligands)
  Table 7: pEC50/Emax for cAMP/Gαo/CI × 2–90 min × 7 ligands
  — note: Table 5 vs Table 6 differ 82-fold in kon for cariprazine (pKd 7.50 vs 9.56);
    the paper's own binding data are internally inconsistent; recorded as-is with a two-tracer sensitivity analysis.
============================================================
"""
import numpy as np
from scipy.interpolate import interp1d

LIGS7 = ["Ropinirole","Dopamine","Aripiprazole","Cariprazine","Bifeprunox","Pardoprunox","S-3PPP"]
# SI Table 6 (PPHT)
KON6 = {"Ropinirole":1.46e6,"Dopamine":3.14e5,"Aripiprazole":1.01e9,"Cariprazine":1.27e9,
        "Bifeprunox":1.84e8,"Pardoprunox":1.25e8,"S-3PPP":3.25e6}
KOFF6= {"Ropinirole":2.60,"Dopamine":2.00,"Aripiprazole":0.21,"Cariprazine":0.35,
        "Bifeprunox":0.01,"Pardoprunox":2.28,"S-3PPP":1.51}
# SI Table 5 (spiperone, three slow ligands only)
KON5 = {"Aripiprazole":1.31e8,"Cariprazine":1.55e7,"Bifeprunox":1.07e8}
KOFF5= {"Aripiprazole":0.14,"Cariprazine":0.49,"Bifeprunox":0.01}
PKI  = {"Aripiprazole":9.43,"Cariprazine":8.90,"Bifeprunox":10.36}  # Table 5 pKi

TIMES=[2,5,10,15,30,45,60,75,90]
T7={
"cAMP":{
 "Ropinirole":([8.02,7.65,7.38,7.37,7.46,7.01,6.73,6.97,6.76],[100]*9),
 "Dopamine":([8.40,8.07,7.71,7.65,7.78,7.29,6.97,7.30,7.03],[98,98,98,99,102,101,105,98,102]),
 "Aripiprazole":([7.26,7.47,7.58,7.65,7.76,7.56,7.45,7.71,7.53],[75,76,75,72,71,64,54,71,65]),
 "Cariprazine":([8.30,8.38,8.38,8.38,8.49,8.37,8.22,8.40,8.16],[77,77,75,72,72,66,55,72,67]),
 "Bifeprunox":([7.60,7.91,8.13,8.27,8.62,8.56,8.51,8.80,8.69],[89,92,91,89,89,88,81,90,88]),
 "Pardoprunox":([8.98,8.69,8.39,8.33,8.46,8.03,7.82,8.04,7.83],[87,87,86,87,85,82,74,84,81]),
 "S-3PPP":([7.10,6.78,6.43,6.37,6.47,5.99,5.82,6.04,5.82],[73,73,70,70,71,64,53,66,60])},
"Gao":{
 "Ropinirole":([6.97,7.02,7.17,7.33,7.51,7.38,7.29,7.40,7.36],[100]*9),
 "Dopamine":([7.26,7.31,7.40,7.49,7.30,7.35,7.35,7.47,7.33],[93,91,97,95,101,100,101,99,95]),
 "Aripiprazole":([6.91,7.07,7.36,7.60,7.64,7.72,7.68,7.39,7.49],[80,71,62,63,68,60,60,87,69]),
 "Cariprazine":([7.97,8.05,8.14,8.25,8.36,8.55,8.58,8.71,8.71],[71,70,67,65,67,61,60,52,63]),
 "Bifeprunox":([7.28,7.57,7.89,8.12,8.39,8.52,8.49,8.66,8.67],[94,80,77,79,87,88,80,94,90]),
 "Pardoprunox":([7.78,7.77,8.00,8.21,8.70,8.71,8.64,8.34,8.58],[82,79,76,81,84,85,86,80,87]),
 "S-3PPP":([6.46,6.54,6.57,6.67,7.02,6.78,6.79,6.62,6.77],[64,68,73,76,74,72,72,68,77])},
"CI":{
 "Ropinirole":([8.33,8.29,8.10,7.97,7.77,7.71,7.62,7.55,7.46],[100]*9),
 "Dopamine":([8.55,8.58,8.36,8.25,8.05,8.03,7.97,7.91,7.86],[95,96,101,103,105,106,106,106,106]),
 "Aripiprazole":([7.65,8.34,8.62,8.80,8.74,9.09,9.17,9.28,9.37],[27,30,25,24,22,21,19,18,17]),
 "Cariprazine":([8.94,9.28,9.38,9.45,9.50,9.52,9.55,9.57,9.60],[33,39,34,35,33,33,32,31,30]),
 "Bifeprunox":([7.79,8.46,8.80,9.04,9.26,9.47,9.58,9.64,9.70],[61,67,57,58,53,51,49,46,44]),
 "Pardoprunox":([9.30,9.41,9.22,9.20,9.06,9.02,8.95,8.87,8.80],[75,76,70,67,62,59,57,54,52]),
 "S-3PPP":([7.66,7.67,7.58,7.58,7.50,7.41,7.39,7.40,7.41],[46,46,36,32,26,23,21,19,17])}}

Lgrid = np.logspace(-11, -4.5, 400)
xgrid = np.log10(Lgrid)

def apparent(lig, t, tau, Em, kon, koff):
    kobs = kon*Lgrid + koff
    req  = Lgrid*kon/kobs
    r    = req*(1-np.exp(-kobs*t))
    y    = Em*tau*r/(1+r*(tau-1))
    em   = y[-1]; half = 0.5*em
    if not np.all(np.isfinite(y)) or y.max() < half or y[0] > half:
        return np.nan
    return -float(np.interp(half, y, xgrid))   # pEC50 = -log10(EC50)

def run(KON, KOFF, tag, ligands):
    print(f"\n########## tracer {tag} ##########")
    for a in ["cAMP","Gao","CI"]:
        print(f"--- {a} ---  (residual = model − obs, dex)")
        for lig in ligands:
            em90 = T7[a][lig][1][8]
            tau  = em90/(100-em90) if em90 < 99.5 else 200.0
            pred = [apparent(lig, t, tau, 100, KON[lig], KOFF[lig]) for t in TIMES]
            obs  = T7[a][lig][0]
            pKd  = -np.log10(KOFF[lig]/KON[lig])
            dm   = pred[8]-pred[0] if np.isfinite(pred[0]) and np.isfinite(pred[8]) else np.nan
            print(f"{lig:14s} pKd={pKd:5.2f} τ={tau:6.1f}  drift model={dm:+.2f}/obs={obs[8]-obs[0]:+.2f}"
                  f"  90' residual={(pred[8]-obs[8]):+.2f}")
            print("    residuals 2'→90': " + " ".join(f"{(p-o):+.2f}" if np.isfinite(p) else "  nan"
                                                 for p,o in zip(pred,obs)))

print("Table 5 vs Table 6 kinetic-affinity self-comparison:",
      {l: f"pKd5={-np.log10(KOFF5[l]/KON5[l]):.2f} / pKd6={-np.log10(KOFF6[l]/KON6[l]):.2f} / pKi={PKI[l]}"
       for l in KON5})
run(KON6, KOFF6, "Table 6 (PPHT-red)", LIGS7)
run(KON5, KOFF5, "Table 5 ([3H]spiperone)", list(KON5))
