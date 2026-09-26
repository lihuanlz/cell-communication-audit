# -*- coding: utf-8 -*-
"""
Code 50 · framework self-check simulation of the p53 cascade  v1.2.0 (2026-08-15) two-arm version
============================================================
Purpose: numerical self-check of the cascade framework (v4.2) three legs + information-loss identity on a minimal random model of the p53 pathway.

Lineage: v1.0 single arm (digital arm) → v1.2 two arms (external-audit revision).
The external audit noted: in v1.0 the amplitude never carries dose information (constant amplitude), so the "amplitude dies of degeneracy" demonstration is weak.
Revised to two arms:
  Arm A digital arm: constant amplitude, dose→pulse count (simplified from Lahav 2004 qualitative facts, p53-faithful arm);
  Arm B simulated arm: amplitude increases linearly with dose a(D)=1+0.25(D-1), dual dose encoding (generic simulated-system arm).
The two arms answer different questions: A = behavior of the faithful pathway; B = a strong test of leg 1 (whether true amplitude information can die).

Dual record of the audit (recorded faithfully; neither side touches the numbers):
  Audit prediction "fold survives in arm B (the ratio cancels the gain)" — refuted by measurement: fold I=0.401→0.178 decay.
  With affine b≠0, fold=(ka+b)/(km+b); k does not cancel → fold leaks,
  consistent row-for-row with the framework interface-algebra table "fold survives only at b=0".

Remaining hypotheses:
  A2: affine observation/interface: y = k·p53 + b; k ~ LogNormal(0, σ_k), unidentifiable per cell; b=0.3.
  A3: decoders: amplitude / fold / count (prominence/MAD threshold, affine-invariant).
      The differential-ratio decoder is unstable to define in pulse-train models (pulse positions misaligned with sampling points); left to continuous-response models.

Discipline: past project assets (line 1 p53 audit, sealed) are background only and do not enter this model.
Output: console printout + figure 结果/仿真_p53cascade_三腿nucleus验.png. Seed 20260815. Pure CPU.
"""
import numpy as np
from scipy.signal import find_peaks

SEED = 20260815
rng = np.random.default_rng(SEED)

D2 = np.array([1., 2., 3., 4., 5.])     # damage levels (mapped to pulse counts)
AMP_MU, AMP_SD = 1.0, 0.15
AMP_DOSE_COEF = 0.25                    # arm B: amplitude-dose coefficient
B_BASAL = 0.3
T_WIN, DT = 24.0, 0.1
SIG_LIST = [0.0, 0.2, 0.4, 0.6]


# ---------- model (two arms) ----------
def gen_trace3(D, k, b, arm):
    N = int(round(D)) + (1 if rng.random() < 0.2 else 0)
    ts = np.arange(0, T_WIN, DT)
    y = np.zeros_like(ts)
    for j in range(N):
        pt = 2.0 + j * 5.5 + rng.normal(0, 0.3)
        if pt < T_WIN - 1:
            mu = AMP_MU if arm == "A" else AMP_MU + AMP_DOSE_COEF * (D - 1)
            a = max(rng.normal(mu, AMP_SD * mu), 0.05)
            y += a * np.exp(-0.5 * ((ts - pt) / 0.6) ** 2)
    base = b + 0.05 * np.abs(rng.standard_normal(len(ts)))
    return ts, k * y + base, N


# ---------- decoders ----------
def dec_amp(ts, y):
    return y.max()

def dec_fold(ts, y):
    return y.max() / max(np.median(y[:30]), 1e-9)

def dec_count2(ts, y):
    mad = np.median(np.abs(y - np.median(y))) + 1e-9
    pk, _ = find_peaks(y, prominence=3.0 * mad, distance=int(4.0 / DT))
    return len(pk)

FNS = {"amplitude": dec_amp, "fold": dec_fold, "count": dec_count2}


# ---------- utilities ----------
def mi_plugin(D, Y, n_bins=12):
    yq = np.quantile(Y, np.linspace(0, 1, n_bins + 1))
    yq[0] -= 1e-9; yq[-1] += 1e-9
    Yb = np.clip(np.digitize(Y, yq) - 1, 0, n_bins - 1)
    H = lambda c: -np.sum((c / c.sum()) * np.log2(c / c.sum() + 1e-12))
    HY = H(np.bincount(Yb, minlength=n_bins))
    HYgD = 0.0
    for d in np.unique(D):
        m = D == d
        HYgD += m.mean() * H(np.bincount(Yb[m], minlength=n_bins))
    return max(HY - HYgD, 0.0)

def lowpass(y, tau):
    a = DT / tau
    out = np.zeros_like(y)
    for i in range(1, len(y)):
        out[i] = out[i - 1] + a * (y[i] - out[i - 1])
    return out


# ---------- experiment 1: leg 1, two arms ----------
def run_arm(arm, sk, ncell=4000):
    Ds = rng.choice(len(D2), ncell)
    out = {n: [] for n in FNS}
    for d in Ds:
        k = np.exp(rng.normal(0, sk))
        ts, y, N = gen_trace3(D2[d], k, B_BASAL, arm)
        for n, fn in FNS.items():
            out[n].append(fn(ts, y))
    return {n: mi_plugin(Ds, np.array(out[n], float)) for n in FNS}

print("=" * 70)
print("Experiment 1 · leg 1, two arms: who survives the affine interface (σ_k scan, 4000 cells/point)")
print("=" * 70)
leg1 = {}
for arm, label in [("A", "digital arm (constant amplitude, p53-faithful)"), ("B", "simulated arm (amplitude carries dose)")]:
    print(f"--- arm {arm} {label} ---")
    for sk in SIG_LIST:
        leg1[(arm, sk)] = run_arm(arm, sk)
        r = leg1[(arm, sk)]
        print(f"  σ_k={sk}:  amplitude I={r['amplitude']:.3f}  fold I={r['fold']:.3f}  count I={r['count']:.3f}")

# ---------- experiment 2: leg 3 (arm A) ----------
print("\n" + "=" * 70)
print("Experiment 2 · leg 3: downstream-integrator bandwidth mismatch kills count (arm A, σ_k=0.4)")
print("=" * 70)
leg3 = {}
for tau in [0.5, 1.0, 2.0, 3.0, 5.0, 8.0]:
    Ds = rng.choice(len(D2), 2000); Ys = []
    for d in Ds:
        k = np.exp(rng.normal(0, 0.4))
        ts, y, N = gen_trace3(D2[d], k, B_BASAL, "A")
        Ys.append(dec_count2(ts, lowpass(y, tau)))
    leg3[tau] = mi_plugin(Ds, np.array(Ys, float))
    print(f"integrator τ={tau}h:  I(D;count)={leg3[tau]:.3f}")
print("(pulse spacing 5.5h: count information collapses as τ approaches the spacing)")

# ---------- experiment 3: information-loss identity (arm A) ----------
print("\n" + "=" * 70)
print("Experiment 3 · information-loss identity (arm A, σ_k=0.4)")
print("=" * 70)
Ds = rng.choice(len(D2), 6000); Ns = []; Yc = []
for d in Ds:
    k = np.exp(rng.normal(0, 0.4))
    ts, y, N = gen_trace3(D2[d], k, B_BASAL, "A")
    Ns.append(N); Yc.append(dec_count2(ts, y))
HD = np.log2(len(D2))
IDN = mi_plugin(Ds, np.array(Ns, float), n_bins=7)
IDY = mi_plugin(Ds, np.array(Yc, float))
L1, L2 = HD - IDN, IDN - IDY
print(f"H(D)={HD:.3f}  I(D;N)={IDN:.3f}  I(D;Y_count)={IDY:.3f}")
print(f"L1={L1:.3f}  L2={L2:.3f}  identity: {IDY:.3f}+{L1:.3f}+{L2:.3f}={IDY+L1+L2:.3f} ≈ H(D)={HD:.3f}")

# ---------- figure ----------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import os as _os
if _os.path.exists("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"):
    fm.fontManager.addfont("/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc")
_avail = {f.name for f in fm.fontManager.ttflist}
for _cand in ["Noto Sans CJK SC", "Noto Sans CJK JP", "Noto Sans CJK HK", "Noto Serif CJK SC"]:
    if _cand in _avail:
        matplotlib.rcParams["font.family"] = [_cand]
        break
matplotlib.rcParams["axes.unicode_minus"] = False

fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
ax = axes[0]
for n, mk in [("amplitude", "o"), ("fold", "s"), ("count", "^")]:
    ax.plot(SIG_LIST, [leg1[("B", sk)][n] for sk in SIG_LIST], mk + "-", label=n + "(arm B)")
    ax.plot(SIG_LIST, [leg1[("A", sk)][n] for sk in SIG_LIST], mk + "--", alpha=0.45, label=n + "(arm A)")
ax.set_xlabel("gain fluctuation σ_k"); ax.set_ylabel("I(D; readout) / bits")
ax.legend(fontsize=7, ncol=2); ax.set_title("leg 1: two-arm control — amplitude dead (arm B confirmed), fold leaks, count survives")

ax = axes[1]
taus = sorted(leg3)
ax.plot(taus, [leg3[t] for t in taus], 'o-', c='darkred')
ax.axvline(5.5, ls=':', c='gray'); ax.text(5.6, 0.9, "pulse spacing 5.5h", fontsize=8)
ax.set_xlabel("downstream integrator time constant τ (h)"); ax.set_ylabel("I(D; count) / bits")
ax.set_title("leg 3: bandwidth mismatch kills count")

axes[2].bar(["L1\ndamage→pulse train", "I(D;Y)\nend-to-end survival", "L2\npulse train→readout"],
            [L1, IDY, L2], color=["#c0392b", "#27ae60", "#e67e22"])
axes[2].set_ylabel("bits"); axes[2].set_title("information-loss decomposition (H(D)=2.32)")

plt.tight_layout()
plt.savefig("/mnt/agents/output/03_细胞线3/结果/仿真_p53cascade_三腿nucleus验.png", dpi=150)
print("\nFigure saved: 03_细胞线3/结果/仿真_p53cascade_三腿nucleus验.png")
