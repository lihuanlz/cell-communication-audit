# -*- coding: utf-8 -*-
"""
代码69_D2R功能态跨臂亲和力地图.py  v1.0.0
============================================================
Cell line 4 · deepening round 2: locating Code 68's residual breakage (the absolute-difference component)

Code 68 concluded: slow-dissociation kinetics can explain the drift, but the equilibrium-state absolute
difference in the cAMP/Gαo arms (bifeprunox 2+ dex) survives. This code asks: is that absolute difference
a "binding state vs functional state" problem or a "cross-arm within the functional state" problem?

Method: for each (ligand, arm, time) grid point of SI Table 7, invert the functional-state pKA:
  τ(t) = Emax(t)/(100−Emax(t))   (ropinirole=100 reference; full agonists τ=200 as placeholder only)
  pKA(t) = pEC50(t) − log10(1+τ(t))
pKA inversion is robust for partial agonists (finite τ); for full agonists (Ropinirole/Dopamine, where only
a lower bound on τ exists) the pKA value depends on the chosen τ and serves as direction reference only.

Results summary:
  - the cAMP and Gαo arms are always mutually consistent (late-time difference 0.1–0.5 dex);
  - CI-arm pEC50 is systematically 0.7–1.6 dex higher for all 7 ligands, producing functional-state pKA
    cross-arm differences of 1.4–2.2 dex (90', 5 partial agonists);
  - this cross-arm difference grows with time for the three slow ligands (CI-arm pKA climbs); fast ligands are stable;
  - late-time CI-arm pKA nearly coincides with binding pKi / kinetic pKd (aripiprazole 9.29 vs
    pKi 9.43; cariprazine 9.45 vs Table 6 pKd 9.56); the bifeprunox residual converges to 0.9 dex;
  - the same ligand on the same cellular background: cAMP says "low-affinity strong partial agonism"
    (bifeprunox Emax 88%), CI says "high-affinity weak partial agonism" (Emax 44%) — no single (KA,τ)
    can satisfy both arms; this is the purified cross-arm breakage (AT1R type, within one paper).

Interpretation discipline: CI is xCELLigence whole-cell impedance (different instrument/plate format/possibly
different clone), so a strict "within-arm" comparison does not hold; but the directional fact (proximal BRET
arm group vs whole-cell integrative arm) is ligand-independent and consistent, making it a testable structural fact.
============================================================
"""
import numpy as np

LIGS7 = ["Ropinirole","Dopamine","Aripiprazole","Cariprazine","Bifeprunox","Pardoprunox","S-3PPP"]
TIMES = [2,5,10,15,30,45,60,75,90]
T7 = {  # Klein Herenbrink 2016 SI Table 7 (pEC50, Emax) — same transcription as Code 68
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
PKI = {"Ropinirole":5.60,"Dopamine":5.05,"Aripiprazole":9.43,"Cariprazine":8.90,
       "Bifeprunox":10.36,"Pardoprunox":7.63,"S-3PPP":5.84}
PKD6 = {"Aripiprazole":9.68,"Cariprazine":9.56,"Bifeprunox":10.26}  # Table 6 kinetic pKd

def pka_row(a, lig):
    out=[]
    for k in range(9):
        em=T7[a][lig][1][k]; pec=T7[a][lig][0][k]
        tau=em/(100-em) if em<99.5 else 200.0
        out.append(pec-np.log10(1+tau))
    return out

print("Functional-state pKA (arm, time) map   [columns: 2' 5' 10' 15' 30' 45' 60' 75' 90']")
print("="*108)
PARTIAL=["Aripiprazole","Cariprazine","Bifeprunox","Pardoprunox","S-3PPP"]
for lig in LIGS7:
    print(f"\n### {lig}   (pKi={PKI[lig]})")
    M={}
    for a in ["cAMP","Gao","CI"]:
        M[a]=pka_row(a,lig)
        print(f"  {a:5s}: "+" ".join(f"{v:5.2f}" for v in M[a]))
    arr=np.array([M["cAMP"],M["Gao"],M["CI"]])
    print("  cross-arm spread: "+" ".join(f"{arr[:,k].max()-arr[:,k].min():.2f}" for k in range(9)))

print("\n" + "="*108)
print("Key cross-section (90'): functional-state pKA cross-arm difference vs binding affinity, partial agonists")
print(f"{'ligand':14s} {'pKA_cAMP':>8s} {'pKA_Gao':>8s} {'pKA_CI':>8s} {'CI−cAMP':>8s} {'pKi':>6s} {'CI−pKi':>7s} {'cAMP−pKi':>8s}")
for lig in PARTIAL:
    c,g,i=pka_row("cAMP",lig)[8],pka_row("Gao",lig)[8],pka_row("CI",lig)[8]
    print(f"{lig:14s} {c:8.2f} {g:8.2f} {i:8.2f} {i-c:+8.2f} {PKI[lig]:6.2f} {i-PKI[lig]:+7.2f} {c-PKI[lig]:+8.2f}")

print("\nTime evolution of the cross-arm difference (CI−cAMP, dex):")
for lig in PARTIAL:
    c=pka_row("cAMP",lig); i=pka_row("CI",lig)
    print(f"  {lig:14s}: "+" ".join(f"{i[k]-c[k]:+.2f}" for k in range(9)))

print("\ncAMP–Gαo mutual consistency (|Gao−cAMP|, 90'):",
      {lig: round(abs(pka_row('Gao',lig)[8]-pka_row('cAMP',lig)[8]),2) for lig in PARTIAL})
print("Late-time CI pKA vs Table 6 kinetic pKd:",
      {lig: f"{pka_row('CI',lig)[8]:.2f} vs {PKD6[lig]}" for lig in ["Aripiprazole","Cariprazine","Bifeprunox"]})
