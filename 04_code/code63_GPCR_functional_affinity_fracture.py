#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code 63 GPCR functional-affinity fracture test: D2R multi-pathway operational model audit v1.1.0
================================================
Cell line 4, first cut. Task: port the "fracture audit" from chemotaxis to GPCR.
Data: Klein Herenbrink et al. 2016, Nat Commun 7:10842 (D2L receptor, CHO
cells, same lab same system), SI Table 1 (7 ligands × 6 pathways log(τ/KA) ± SEM,
5 min), SI Table 5 ([3H]spiperone competition binding pKi ± SEM), SI Table 6 (fluorescent
PPHT tracer kinetics pKd, independent replication), SI Table 7 (cAMP/GαoB pEC50+Emax time course).

Audit design (internal constraints of the operational model; breaking any one is a fracture):
  Test 1 (binding anchor): pKi anchors τ → predicts Emax vs measured.
  Test 2 (functional-affinity inversion): pKA_func=logR−log τ; cross-pathway conservation +
    consistency with pKi, SEM fully propagated.
  Test 3 (time invariance): SI T7 all-ligand pEC50 time course (2–90 min);
    the equilibrium model requires ΔlogR not to drift.
  Self-audit ①: zero-model criterion pKi−pEC50>0 (under BL, pEC50=pKA+log(1+τ)≥pKA,
    under the single-KA assumption it must not lie significantly below pKi) — no fitting assumptions.
  Self-audit ②: second tracer (PPHT kinetics pKd) independently replicates the binding anchor.
  Self-audit ③: sensitivity to ±10% perturbation of system maximum E_sys.
  Self-audit ④: full error propagation (σ_Δ includes logR/pKi/Emax).

Results (v1.1.0 actual run):
  Test 2: Δ=pKA_func−pKi: bifeprunox −3.5 (21–23σ), aripiprazole
    −2.6/−2.9 (18–26σ), cariprazine −1.2/−1.4 (9–11σ); conventional agonists
    |Δ|≤1.3 (dopamine +1.2/1.3 in the direction expected for two-state shift). pKA_func conserved
    across pathways (between-arm difference 0.08–0.31 dex).
  Self-audit ①: zero-model reproduction: pKi−pEC50 = +1.96/+0.52/+2.45 dex (the three high-affinity
    partial agonists), negative for the other ligands (spare-receptor direction, allowed under single KA).
  Self-audit ②: PPHT replication: same three ligands pKd−pEC50 = +2.19/+1.15/+2.39 dex;
    the two tracers differ by only 0.06–0.63 dex — the anchor itself is solid.
  Self-audit ③: with E_sys 100→110, Δ shifts only 0.1–0.35 dex; classification unchanged.
  Test 3: all full-agonist potencies shift down −1.15…−1.37 dex by 90 min (desensitization); the three partial
    agonists stable or shifted up (bifeprunox +1.09); bifeprunox vs ropinirole
    bias factor −0.42 at 2' → +1.93 dex at 90', a 2.35 dex swing with sign reversal
    (SEM 0.02–0.03) → equilibrium bias quantification gives opposite answers depending on readout time.

Verdict: the single-affinity equilibrium operational model fractures on the D2R data, localized to the high-affinity
partial agonists — exactly the class to which all significant-bias reports belong (SI T3 asterisks cluster here).
All four method self-audits pass; the fracture is robust.

⚠️ dual-record: PTH1R (Sachdev 2024, Source Data archived) lacks a binding anchor;
on normalized curves Black–Leff τ/KA is unidentifiable when nH≠1 (boundary-hitting artifacts
+5.77/+9.58 dex, no conclusion drawn). Lesson: without an independent affinity anchor, marginal
operational statistics do not constitute a falsifiable test.

run: python3 代码63_GPCR功能亲和力_断裂检验.py
dependencies: numpy
"""

import numpy as np

LIGS = ["Dopamine", "Ropinirole", "Aripiprazole", "Cariprazine",
        "Bifeprunox", "Pardoprunox", "S-3PPP"]
logR = {
"Dopamine":    {"cAMP": (8.06,.05), "pERK": (8.56,.05), "Gai1": (7.17,.09), "GaoB": (7.24,.08), "CI": (8.55,.05), "Barr": (6.33,.13)},
"Ropinirole":  {"cAMP": (7.67,.04), "pERK": (8.44,.05), "Gai1": (6.75,.09), "GaoB": (7.08,.08), "CI": (8.32,.05), "Barr": (6.48,.15)},
"Aripiprazole":{"cAMP": (7.36,.06), "pERK": (6.48,.54), "Gai1": (7.90,.16), "GaoB": (6.94,.11), "CI": (7.83,.15), "Barr": (6.17,.58)},
"Cariprazine": {"cAMP": (8.27,.06), "pERK": (7.45,.47), "Gai1": (8.10,.14), "GaoB": (7.91,.11), "CI": (8.88,.11), "Barr": None},
"Bifeprunox":  {"cAMP": (7.88,.05), "pERK": (7.32,.10), "Gai1": (8.60,.11), "GaoB": (7.50,.10), "CI": (8.30,.07), "Barr": (5.75,.16)},
"Pardoprunox": {"cAMP": (8.64,.05), "pERK": (9.06,.08), "Gai1": (8.10,.13), "GaoB": (7.70,.10), "CI": (9.29,.06), "Barr": (6.80,.29)},
"S-3PPP":      {"cAMP": (6.65,.06), "pERK": (7.39,.14), "Gai1": (6.60,.12), "GaoB": (6.40,.11), "CI": (7.34,.10), "Barr": None},
}
pKi = {"Dopamine": (5.05,.06), "Ropinirole": (5.60,.08), "Aripiprazole": (9.43,.06),
       "Cariprazine": (8.90,.07), "Bifeprunox": (10.36,.09), "Pardoprunox": (7.63,.06),
       "S-3PPP": (5.84,.06)}
pKd2 = {"Dopamine": 5.18, "Ropinirole": 5.73, "Aripiprazole": 9.66, "Cariprazine": 9.53,
        "Bifeprunox": 10.30, "Pardoprunox": 7.75, "S-3PPP": 6.11}
EMAX5 = {"cAMP": {"Dopamine": (98,2.0), "Ropinirole": (100,0), "Aripiprazole": (76,2.0),
                  "Cariprazine": (77,2.4), "Bifeprunox": (92,2.0), "Pardoprunox": (87,1.5), "S-3PPP": (73,2.0)},
         "GaoB": {"Dopamine": (91,2.0), "Ropinirole": (100,0), "Aripiprazole": (71,5.0),
                  "Cariprazine": (70,3.7), "Bifeprunox": (80,3.6), "Pardoprunox": (79,4.3), "S-3PPP": (68,3.0)}}
PEC50_5 = {"Dopamine": 8.07, "Ropinirole": 7.65, "Aripiprazole": 7.47, "Cariprazine": 8.38,
           "Bifeprunox": 7.91, "Pardoprunox": 8.69, "S-3PPP": 6.78}
T = [2, 5, 10, 15, 30, 45, 60, 75, 90]
PEC50_T = {
"Ropinirole":  [8.02,7.65,7.38,7.37,7.46,7.01,6.73,6.97,6.76],
"Dopamine":    [8.40,8.07,7.71,7.65,7.78,7.29,6.97,7.30,7.03],
"Aripiprazole":[7.26,7.47,7.58,7.65,7.76,7.56,7.45,7.71,7.53],
"Cariprazine": [8.30,8.38,8.38,8.38,8.49,8.37,8.22,8.40,8.16],
"Bifeprunox":  [7.60,7.91,8.13,8.27,8.62,8.56,8.51,8.80,8.69],
"Pardoprunox": [8.98,8.69,8.39,8.33,8.46,8.03,7.82,8.04,7.83],
"S-3PPP":      [7.10,6.78,6.43,6.37,6.47,5.99,5.82,6.04,5.82],
}

def main():
    print("== Test 2: functional KA inversion, full error propagation ==")
    print(f"{'ligand':13s} {'path':5s} {'pKA_func':>8s} {'pKi':>6s} {'Δ(dex)':>8s} {'σ(Δ)':>6s} {'signif.':>6s}")
    for lig in LIGS:
        for path in ["cAMP", "GaoB"]:
            em, eme = EMAX5[path][lig]
            if em >= 99.5:
                print(f"{lig:13s} {path:5s}   τ→∞ (full agonism, KA indistinguishable)")
                continue
            tau = em / (100 - em)
            lr, lre = logR[lig][path]; pk, pke = pKi[lig]
            pka = lr - np.log10(tau)
            s_logtau = (1 / np.log(10)) * 100 * eme / (em * (100 - em))
            sd = np.sqrt(lre**2 + pke**2 + s_logtau**2)
            print(f"{lig:13s} {path:5s} {pka:8.2f} {pk:6.2f} {pka-pk:+8.2f} {sd:6.2f} {abs(pka-pk)/sd:5.1f}σ")

    print("\n== Self-audit ①: zero-model criterion pKi−pEC50 (>0.3 = fracture, no fitting assumptions) ==")
    for lig in LIGS:
        d = pKi[lig][0] - PEC50_5[lig]
        tag = "⚠fracture" if d > 0.3 else ("spare-receptor direction (expected for full agonism)" if d < -0.3 else "consistent")
        print(f"  {lig:13s} pKi={pKi[lig][0]:5.2f} pEC50={PEC50_5[lig]:5.2f} diff={d:+.2f}  {tag}")

    print("\n== Self-audit ②: second tracer PPHT kinetics pKd replication ==")
    for lig in LIGS:
        print(f"  {lig:13s} pKd_ppht={pKd2[lig]:5.2f} two-tracer diff={pKd2[lig]-pKi[lig][0]:+.2f} pKd−pEC50={pKd2[lig]-PEC50_5[lig]:+.2f}")

    print("\n== Self-audit ③: system-maximum perturbation (E_sys=100/105/110, cAMP arm Δ) ==")
    for esys in [100, 105, 110]:
        out = []
        for lig in LIGS:
            em = min(EMAX5["cAMP"][lig][0], esys - 1)
            tau = em / (esys - em)
            out.append(f"{lig[:4]}={logR[lig]['cAMP'][0]-np.log10(tau)-pKi[lig][0]:+.2f}")
        print(f"  E_sys={esys}: " + ", ".join(out))

    print("\n== Test 3: pEC50 time drift (2'→90') and bias-factor swing ==")
    for lig in LIGS:
        v = PEC50_T[lig]
        print(f"  {lig:13s}: {v[0]:.2f}→{v[-1]:.2f}  Δ={v[-1]-v[0]:+.2f} dex")
    for lig in ["Bifeprunox", "Aripiprazole", "Cariprazine", "S-3PPP"]:
        d = [a - b for a, b in zip(PEC50_T[lig], PEC50_T["Ropinirole"])]
        print(f"  {lig:13s} vs Ropi: 2'={d[0]:+.2f} → 90'={d[-1]:+.2f}  swing={max(d)-min(d):.2f} dex")

if __name__ == "__main__":
    main()
