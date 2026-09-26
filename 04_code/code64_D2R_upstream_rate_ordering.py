#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code 64 D2R upstream rate-ordering test (deep dig into raw time courses) v1.0.0 — data quality insufficient, negative recorded as-is
================================================
Task: test "whether pathway kinetics are ligand-independent". Single-active-state + equilibrium cascade
predicts: within one pathway all ligands' normalized response waveforms should coincide (rate ordering
is a pathway property). If ligands leave kinetic fingerprints on the pathway, deeper state resolution exists.

Data: Klein Herenbrink 2016 Supplementary Data 1 (single-concentration saturating-stimulus time courses
for 6 pathways, n=3–4 repeats/ligand).

Results (v1.0.0 actual run):
  Direct t50 extraction gives huge ligand-to-ligand differences (e.g. Gαi1: dopamine 0.9 min vs
  bifeprunox 17.9 min) — an apparent "ligand fingerprint".
  ⚠️ But quality control exposes three fatal problems:
  1. Monotonicity generally 0.3–0.7 (trajectory drift / high noise, not a clean monophasic rise);
  2. Plateau-normalized amplitudes take physically unreasonable values (pERK bifeprunox = 572% of dopamine,
     Gαi1 aripiprazole = 804%), contradicting the Emax ordering of SI Table 7 —
     meaning baseline/plateau estimates are contaminated by drift;
  3. The pERK time-course window starts only at 30 min (capturing the decay tail); its t50 is undefined.
Verdict: **this batch of raw time courses is insufficient to support a kinetic-ordering test** — the apparent
ligand differences in t50 are not credible, usable neither as breakage evidence nor as no-breakage evidence. Negative recorded as-is.
The usable evidence for kinetic breakage remains Code 63 Test 3 (the SI T7 pEC50 time course,
which is a full concentration–response fit at each time point, a completely different statistical scale).

Lesson (double-recorded): single-concentration time courses + low replication + baseline drift cannot support
waveform-normalization comparisons; such tests need per-time-point concentration–response surfaces (time × dose matrix).

Run: python3 代码64_D2R上游速率排序检验.py
Dependencies: numpy, pandas, openpyxl
"""

import numpy as np
import pandas as pd

XLSX = "/mnt/agents/output/04_细胞线4/公开数据/D2R_KleinHerenbrink2016/SupplementaryData1.xlsx"

def parse_sheet(sh):
    df = pd.read_excel(XLSX, sheet_name=sh, header=None)
    t = pd.to_numeric(df.iloc[3:, 0], errors="coerce")
    curves = {}
    c = 1
    while c < df.shape[1]:
        name = df.iloc[1, c]
        if pd.isna(name):
            c += 1; continue
        cols = [cc for cc in range(c, min(c + 4, df.shape[1]))
                if str(df.iloc[2, cc]).startswith("n")]
        if not cols:
            c += 1; continue
        Y = df.iloc[3:, cols].apply(pd.to_numeric, errors="coerce")
        ok = ~(t.isna() | Y.mean(axis=1).isna())
        curves[str(name)] = (t[ok].to_numpy(float),
                             Y[ok].mean(axis=1).to_numpy(float))
        c = max(cols) + 1
    return curves

def t50(t, y):
    base = np.mean(y[:3]); plat = np.mean(y[int(len(y) * 0.8):])
    yn = (y - base) / (plat - base + 1e-12)
    for i in range(1, len(yn)):
        if (yn[i-1] - 0.5) * (yn[i] - 0.5) <= 0 and yn[i] != yn[i-1]:
            return t[i-1] + (0.5 - yn[i-1]) * (t[i]-t[i-1]) / (yn[i]-yn[i-1])
    return np.nan

def main():
    sheets = {sh: parse_sheet(sh)
              for sh in ["CAMYEL cAMP", "pERK12", "B-arrestin2", "Gai1", "Gao"]}
    print("== Apparent t50 (min) and quality control ==")
    for sh, cu in sheets.items():
        spans = {lig: np.mean(y[int(len(y)*0.8):]) - np.mean(y[:3])
                 for lig, (t, y) in cu.items()}
        d0 = abs(spans.get("Dopamine", 1)) or 1
        print(f"\n{sh}")
        for lig, (t, y) in cu.items():
            base = np.mean(y[:3]); plat = np.mean(y[int(len(y)*0.8):])
            yn = (y - base) / (plat - base + 1e-12)
            mono = np.mean(np.diff(yn) >= -0.05)
            print(f"  {lig:13s}: t50={t50(t,y):6.1f} min, amp={100*abs(spans[lig])/d0:6.1f}% DA, "
                  f"monotonicity={mono:.2f}, window=[{t.min():.0f},{t.max():.0f}] min")
    print("""
== Verdict ==
Monotonicity 0.3–0.7, physically unreasonable normalized amplitudes (pERK 572%, Gαi1 804%),
pERK window contains only the decay tail -> apparent t50 ligand differences are not credible.
This batch of single-concentration time courses is insufficient for a kinetic-ordering test (negative, double-recorded).
Valid evidence for kinetic breakage remains Code 63 Test 3 (SI T7 pEC50 time course).
""")

if __name__ == "__main__":
    main()
