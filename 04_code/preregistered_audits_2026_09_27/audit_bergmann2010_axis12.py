# -*- coding: utf-8 -*-
"""
Morphogen-line audit axis 1+2: Bergmann et al. 2010 (PMC2858443) Datasets S1/S2.

Axis 1 (relative vs absolute readout):
  Dataset S1 gives only aggregate sigma(x/L) in % EL; without per-embryo lengths
  we cannot convert to um. Decisive evidence is therefore Dataset S2 scaling
  coefficient S = beta * Lbar / xbar from regression x = alpha + beta*L.
  S = 1  -> perfect scaling (relative readout)
  S > 1  -> hyper-scaling (boundary shifts more than embryo length)
  S < 1  -> hypo-scaling (absolute-anchor component)

Axis 2 (mechanism discrimination, per Bergmann's four-model comparison):
  Nuclear-density degradation model predicts hyper-scaling anterior and
  near-perfect scaling mid/posterior at 1x bcd, with S dose-dependent.
  We independently recompute statistics from their Table S2.

Pre-registered checks:
  P1: at 1x bcd, anterior boundaries (gt, eve stripes 1-2, <35% EL) have S>1
      with 68% CI excluding 1 (hyper-scaling confirmed).
  P2: mid/posterior boundaries (40-70% EL) have S consistent with 1
      (|S-1| <= CI) in the majority of cases.
  P3: S decreases with bcd dosage for the same relative position
      (anterior hyper-scaling weakens at 4x) -> position is read off the Bcd
      concentration profile, not off pure geometry.
  P4: precision sigma(x/L) ~ 1-4% EL across all genes/doses (no abrupt
      degradation posteriorly) from Dataset S1.
"""
import csv, math, pathlib, statistics

HERE = pathlib.Path(__file__).parent
DATA = HERE / "data"

def load(name):
    with open(DATA / name, encoding="utf-8") as f:
        rows = list(csv.reader(f))
    header = rows[0]
    out = []
    for r in rows[1:]:
        out.append({
            "x_pct": float(r[0]),
            "v1": float(r[1]),   # sigma (S1) or S (S2)
            "v2": float(r[2]),   # error (S1) or 68% CI (S2)
            "gene": r[3].strip(),
            "costain": r[4].strip(),
            "dose": int(r[5]),
        })
    return header, out

h1, s1 = load("bergmann2010_datasetS1_precision.csv")
h2, s2 = load("bergmann2010_datasetS2_scaling.csv")
assert len(s1) == len(s2) == 84
# rows must correspond one-to-one
for a, b in zip(s1, s2):
    assert abs(a["x_pct"] - b["x_pct"]) < 1e-9 and a["gene"] == b["gene"] and a["dose"] == b["dose"]

def domain(x):
    if x < 35: return "anterior(<35%EL)"
    if x <= 70: return "mid(35-70%EL)"
    return "posterior(>70%EL)"

print("=== Axis 2: scaling coefficient S by domain and dose ===")
print(f"{'domain':<20}{'dose':<5}{'n':<4}{'S median':<9}{'S range':<18}{'|S-1|<=CI frac':<14}{'S>1 & CI excl 1':<6}")
for dom in ["anterior(<35%EL)", "mid(35-70%EL)", "posterior(>70%EL)"]:
    for dose in (1, 2, 4):
        sub = [r for r in s2 if domain(r["x_pct"]) == dom and r["dose"] == dose]
        if not sub: continue
        Ss = [r["v1"] for r in sub]
        ok1 = sum(1 for r in sub if abs(r["v1"] - 1) <= r["v2"]) / len(sub)
        hyper = sum(1 for r in sub if r["v1"] - 1 > r["v2"])   # S>1 with CI excluding 1
        print(f"{dom:<20}{dose:<5}{len(sub):<4}{statistics.median(Ss):<9.3f}"
              f"{f'{min(Ss):.2f}-{max(Ss):.2f}':<18}{ok1:<14.2f}{hyper}/{len(sub)}")

print()
print("=== P1: anterior hyper-scaling at 1x bcd (S>1, CI excludes 1) ===")
ant1 = [r for r in s2 if r["dose"] == 1 and r["x_pct"] < 35]
n_h = sum(1 for r in ant1 if r["v1"] - 1 > r["v2"])
print(f"anterior 1x boundaries: {len(ant1)}, hyper-scaling: {n_h}, "
      f"S median {statistics.median(r['v1'] for r in ant1):.3f}")
for r in ant1[:6]:
    print(f"  {r['gene']} x={r['x_pct']:.1f}% S={r['v1']:.2f} CI={r['v2']:.2f}")

print()
print("=== P3: dose dependence of S at matched boundaries (gene, co-stain, ~same x) ===")
# group by (gene, x rounded to nearest int) across doses
from collections import defaultdict
grp = defaultdict(dict)
for r in s2:
    grp[(r["gene"], round(r["x_pct"]))][r["dose"]] = r["v1"]
drops, flat, rises = 0, 0, 0
examples = []
for k, d in grp.items():
    if 1 in d and 4 in d:
        if d[4] < d[1] - 0.2: drops += 1; examples.append((k, d[1], d[4]))
        elif d[4] > d[1] + 0.2: rises += 1
        else: flat += 1
print(f"matched 1x vs 4x pairs: S drops {drops}, flat {flat}, rises {rises}")
for k, a, b in examples[:8]:
    print(f"  {k[0]} x~{k[1]}%  S(1x)={a:.2f} -> S(4x)={b:.2f}")

print()
print("=== P4: precision sigma(x/L) from Dataset S1 ===")
sig = [r["v1"] for r in s1]
print(f"all 84 boundaries: sigma median {statistics.median(sig):.2f}% EL, "
      f"range {min(sig):.2f}-{max(sig):.2f}%")
for dose in (1, 2, 4):
    sub = [r["v1"] for r in s1 if r["dose"] == dose]
    print(f"  dose {dose}x: n={len(sub)}, median {statistics.median(sub):.2f}%, "
          f"range {min(sub):.2f}-{max(sub):.2f}%")

print()
print("=== Axis 2 global: weighted mean S by dose (all genes) ===")
for dose in (1, 2, 4):
    sub = [r for r in s2 if r["dose"] == dose]
    w = [1 / r["v2"] ** 2 for r in sub if r["v2"] > 0]
    S = [r["v1"] for r in sub if r["v2"] > 0]
    wm = sum(a * b for a, b in zip(S, w)) / sum(w)
    print(f"  dose {dose}x: weighted mean S = {wm:.3f} (n={len(S)})")

print()
print("=== posterior tail check (>56% EL, eve stripes 5-7 region) ===")
post = [r for r in s2 if r["x_pct"] > 56]
ok = sum(1 for r in post if abs(r["v1"] - 1) <= r["v2"])
print(f"n={len(post)}, consistent with S=1: {ok} ({ok/len(post):.0%}), "
      f"S median {statistics.median(r['v1'] for r in post):.3f}")
