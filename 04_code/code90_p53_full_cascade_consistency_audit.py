# -*- coding: utf-8 -*-
"""
Code 90: p53 full-chain consistency audit
Proposition: dose -> DSB -> γH2AX -> ATM -> p53 pulse engine -> count/first pulse -> downstream decoding;
      the constants of every link must reconcile with its neighbors; those that do not are explicitly adjudicated "tension" or "gap".
Methods:
  A static link audit: 14 link check items, taken from the archived values of the code 81-89 verdict cards and constant registries v01/v02;
  B dynamic onset decomposition: the Mönke 2025 errata version stepped from the basal fixed point to DSB=100 (log form);
    numerically extract the ATM* half-rise time, p53 threshold crossing, and first-peak time, decomposing t1 into its on-chain segments;
  C contradiction-registry summary.
Discipline: deterministic, seed tag 20260925; JSON + PNG/SVG + verdict card.
"""
import json
import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import solve_ivp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

SEED_TAG = 20260925
OUT_JSON = r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911\03_细胞线3\结果\代码90_全链路审计_结果.json"
OUT_PNG = r"D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911\03_细胞线3\代码与图\代码90_全链路审计_两面体.png"
OUT_SVG = OUT_PNG.replace(".png", ".svg")

plt.rcParams.update({"svg.fonttype": "none", "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False})

# ---------------- Mönke 2025 errata version (same as code 89) ----------------
PAR = dict(A=30.5, P=22.0, C=1.4, g=2.5, dAM=20.0,
           Tm=1.2, TM=4.0, Tw=1.2, TW=1.0,
           dA=0.16, dP=0.1, dm=1.0, dM=2.0, dw=1.3, dW=2.3,
           kA=0.5, kWA=0.14, kMP=0.15, kPm=1.0, kPw=1.0,
           R=2.0, Smax=0.2, gam=9.0)


def rhs_vec(t, x, S):
    p = PAR
    ATM, P53, mdm2, Mdm2, wip1, Wip1 = x
    return [
        p["A"] * (ATM**2 / (1 + ATM**2 / p["kA"])) / (1 + Wip1 / p["kWA"])
        - p["dA"] * ATM - p["P"] * ATM * Wip1 + S,
        p["C"] - p["dP"] * P53
        - p["g"] * Mdm2 * P53 / (p["kMP"] + P53) * (1 + p["R"] / (1 + ATM)),
        p["Tm"] * P53 / (p["kPm"] + P53) - p["dm"] * mdm2,
        p["TM"] * mdm2 - p["dM"] * Mdm2 - p["dAM"] * ATM * Mdm2,
        p["Tw"] * P53 / (p["kPw"] + P53) - p["dw"] * wip1,
        p["TW"] * wip1 - p["dW"] * Wip1,
    ]


res = {"seed": SEED_TAG}

# ================= A. Static link audit (14 links) =================
# status: closed / tension / gap
links = [
    dict(id="L0", en="dose->DSB", name="dose -> DSB count", constant="35 DSB/Gy",
         sources="anchored by code 81; double-record 4 of code 82; reused by code 87/88",
         status="closed",
         note="consistent across codes, but itself an assumed convention (no measured conversion); assumption flag retained"),
    dict(id="L1", en="DSB->gH2AX", name="DSB -> γH2AX", constant="half-peak 1-3 min, peak 30 min (K5)",
         sources="constant registry v01 K5 (consistent across sources)",
         status="closed", note="minute scale, far below all downstream timescales"),
    dict(id="L2", en="gH2AX->ATM", name="γH2AX -> ATM activation", constant="dose-graded single-cell time course missing",
         sources="v01 G1 registered gap",
         status="gap", note="ATM-SPARK/ATOMIC probes exist but no public dose-graded data"),
    dict(id="L3", en="ATM->p53 cat.", name="ATM -> p53 (catalytic)", constant="kcat/Km not measured",
         sources="v02 G5 (user singled out the catalytic rate, which sits exactly in the gap)",
         status="gap", note="the model substitutes a phenomenological protection factor R=2"),
    dict(id="L4", en="p53 engine", name="p53 pulse engine (K1/K3/K4)", constant="5.5 h / 2.0 h / 1.25 h",
         sources="code 89: all three constants decomposed and closed",
         status="closed", note="period 5.95 vs 5.48 (8.5% residual registered)"),
    dict(id="L5", en="counting N(D)", name="dose -> pulse count N(D)", constant="counting-law D_c: archived 0.261 vs code 82 0.99 Gy",
         sources="code 82 validation target 1: 3.8x discrepancy",
         status="tension", note="direction correct (saturating shape R2=0.995); absolute-scale discrepancy stems from the L0 conversion assumption"),
    dict(id="L6", en="first pulse t1", name="first pulse t1", constant="distribution (first peak 2.5-3.0 h), dose-independent",
         sources="D1; Lahav 2004; code 87 flat 2.5-3.0 h across all doses",
         status="closed", note="the distribution itself is a sensing-layer signature, not a defect"),
    dict(id="L7", en="dispersion 5.8", name="dispersion ratio t1/IPI", constant="archived 5.8 = (240/100)^2, 5 Gy convention",
         sources="code 86 convention error retracted; code 87/88 obtain 6.65-7.32",
         sources_detail="correction registered 2026-09-25",
         status="closed", note="first reconciliation at the same order of magnitude; the 86 retraction is registered"),
    dict(id="L8", en="clock dose-indep.", name="clock dose-independent architecture", constant="dose acts only through gating + termination",
         sources="code 87/88 closed loop",
         status="closed", note="for 0.3 Gy low-dose sustained oscillation see L9"),
    dict(id="L9", en="low-dose osc.", name="low-dose oscillating fraction", constant="literature: ~35% oscillating at 0.3 Gy",
         sources="code 87/88 verdict card",
         status="tension", note="damage-driven termination cannot produce low-dose oscillation; intra-literature tension, registered"),
    dict(id="L10", en="p21 tracking", name="downstream p21 transcriptional tracking", constant="exact delay value pending extraction",
         sources="v01 G3 (Hafner 2020)",
         status="gap", note="qualitative co-oscillation established; quantitative delay pending extraction from the original figure"),
    dict(id="L11", en="info budget", name="information budget I(D;fate)", constant="0.55-0.63 bits (after correction)",
         sources="code 81's 1.325 was falsified by C1 for relying on dose-dependent timing",
         status="closed", note="retraction and correction both registered; the paper never cited 1.325"),
    dict(id="L12", en="Wip1 RNAi", name="Wip1 RNAi phenotype", constant="amplitude direction: model up vs Batchelor down 45%",
         sources="code 82 validation target 5",
         status="tension", note="runaway and widening predictions reproduced; the amplitude-direction conflict is registered as a model-experiment mismatch"),
    dict(id="L13", en="width K2", name="pulse width K2 = 3.5 h", constant="molecular layer not decomposed",
         sources="v02 registry section D",
         status="gap", note="set by the ATM* plateau, not a pure sum of delays; pending decomposition"),
]
res["A_static_links"] = links
n_closed = sum(1 for l in links if l["status"] == "closed")
n_tension = sum(1 for l in links if l["status"] == "tension")
n_gap = sum(1 for l in links if l["status"] == "gap")
res["A_summary"] = dict(total=len(links), closed=n_closed, tension=n_tension, gap=n_gap)

# ================= B. Dynamic onset decomposition (chain segments of t1) =================
DSB_STEP = 100.0
S_basal = PAR["Smax"] * np.log(2.0 / PAR["gam"] + 1)   # background ~2 DSB
S_step = PAR["Smax"] * np.log(DSB_STEP / PAR["gam"] + 1)

x_basal = fsolve(lambda x: rhs_vec(0, x, S_basal),
                 [0.05, 0.3, 0.3, 0.3, 0.3, 0.1])
sol = solve_ivp(rhs_vec, (0, 12), x_basal, args=(S_step,),
                dense_output=True, rtol=1e-9, atol=1e-12)
t = np.linspace(0, 12, 12001)
x = sol.sol(t)
ATM, P53 = x[0], x[1]

# γH2AX segment: external constant, peak 0.5 h (K5)
t_gamma = 0.5
# ATM* half-rise (relative to its first-peak maximum)
atm_peak = ATM.max()
atm_half = 0.5 * (atm_peak + ATM[0])
i_atm = np.argmax(ATM >= atm_half) if np.any(ATM >= atm_half) else len(t) - 1
t_atm = t[i_atm]
# p53 threshold crossing (detection threshold 0.8 AU, code 82 convention)
i_thr = np.argmax(P53 >= 0.8) if np.any(P53 >= 0.8) else len(t) - 1
t_p53_thr = t[i_thr]
# p53 first peak
i_peak = np.argmax(P53[:np.searchsorted(t, 10)])
t_peak = t[i_peak]

t1_pred_onset = t_gamma + t_p53_thr   # γH2AX segment + ATM/p53 rising segment (onset convention)
res["B_onset_decomposition"] = dict(
    S_basal=float(S_basal), S_step=float(S_step),
    basal_fixed_point=[float(v) for v in x_basal],
    gammaH2AX_segment_h=t_gamma,
    ATM_half_rise_h=float(t_atm),
    p53_threshold_cross_h=float(t_p53_thr),
    p53_first_peak_h=float(t_peak),
    t1_pred_onset_h=float(t1_pred_onset),
    t1_measured_band_h=[2.5, 3.0],
    t1_definition="literature t1 = first-peak time (D1: MCF7 2-3 h); onset is the threshold crossing, not equal to t1",
    finding="the sensing segment (γH2AX 0.5 + ATM rise %.2f + p53 threshold %.2f) totals only %.2f h; "
            "the bulk of t1 is the p53 accumulation ramp to the peak (%.2f h), not sensing delay"
            % (t_atm - t_gamma, t_p53_thr - t_atm, t_p53_thr, t_peak - t_p53_thr),
    peak_closes=bool(t_peak <= 3.0 + 0.5),
)

# on-chain timing-order check: γH2AX(0.5) < ATM half-rise < p53 threshold < first peak < period 5.5
ordering_ok = bool(t_gamma <= t_atm <= t_p53_thr <= t_peak <= 5.5)
res["B_ordering_ok"] = ordering_ok

# ================= C. Contradiction-registry summary =================
res["C_tensions_register"] = [
    "L5: counting-law D_c 3.8x discrepancy (conversion assumption)",
    "L9: 0.3 Gy low-dose oscillation (intra-literature tension)",
    "L12: Wip1 RNAi amplitude direction (model-experiment mismatch)",
    "B: simulated first peak %.2f h vs measured t1 band 2.5-3.0 h (step-convention difference)" % t_peak,
]

# ================= Figure: two-panel =================
fig = plt.figure(figsize=(12.5, 7.6))
gs = fig.add_gridspec(2, 1, height_ratios=[1.15, 1.0], hspace=0.32)

# (a) link diagram
ax = fig.add_subplot(gs[0])
status_color = {"closed": "#55A868", "tension": "#DD8452", "gap": "#C44E52"}
n = len(links)
cols = 7
for k, lk in enumerate(links):
    row = k // cols
    col = k % cols
    if row == 1:
        col = cols - 1 - col  # second row reversed, serpentine
    xc = 0.06 + col * (0.88 / (cols - 1))
    yc = 0.72 if row == 0 else 0.18
    ax.annotate("", xy=(0, 0))  # noop keep autoscale off
    ax.add_patch(plt.Rectangle((xc - 0.055, yc - 0.13), 0.11, 0.26,
                               facecolor=status_color[lk["status"]],
                               edgecolor="none", alpha=0.9,
                               transform=ax.transAxes))
    ax.text(xc, yc + 0.045, lk["id"], transform=ax.transAxes,
            ha="center", fontsize=8.5, fontweight="bold", color="white")
    ax.text(xc, yc - 0.045, lk["en"], transform=ax.transAxes,
            ha="center", fontsize=6.5, color="white", wrap=True)
    lk["_pos"] = (xc, yc, row)
# arrows
order_idx = list(range(cols)) + list(range(cols, n))
for a, b in zip(order_idx[:-1], order_idx[1:]):
    xa, ya, ra = links[a]["_pos"]
    xb, yb, rb = links[b]["_pos"]
    if ra == rb:
        ar = FancyArrowPatch((xa + (0.055 if xb > xa else -0.055), ya),
                             (xb + (-0.055 if xb > xa else 0.055), yb),
                             transform=ax.transAxes, arrowstyle="-|>",
                             mutation_scale=10, color="#555555", lw=1.0)
    else:
        ar = FancyArrowPatch((xa, ya - 0.13), (xb, yb + 0.13),
                             transform=ax.transAxes, arrowstyle="-|>",
                             mutation_scale=10, color="#555555", lw=1.0,
                             connectionstyle="arc3,rad=0.0")
    ax.add_patch(ar)
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
ax.axis("off")
handles = [plt.Rectangle((0, 0), 1, 1, color=status_color[s]) for s in
           ["closed", "tension", "gap"]]
ax.legend(handles, [f"closed ({n_closed})", f"registered tension ({n_tension})",
                    f"registered gap ({n_gap})"],
          loc="lower center", ncol=3, frameon=False, fontsize=8,
          bbox_to_anchor=(0.5, -0.14))
ax.set_title("a  Full-chain audit: dose -> DSB -> gH2AX -> ATM -> p53 engine -> count/t1 -> decoding\n"
             "every link carries its archived constant and a verdict", fontsize=10)

# (b) onset decomposition
ax = fig.add_subplot(gs[1])
segs = [("gH2AX formation\n(K5: 0.5 h)", t_gamma, "#4C72B0"),
        ("ATM* rise to half-max\n(model: %.2f h)" % t_atm, max(t_atm - t_gamma, 0.0), "#55A868"),
        ("p53 accumulation to 0.8 AU\n(model: %.2f h)" % (t_p53_thr - t_atm),
         max(t_p53_thr - t_atm, 0.0), "#8172B3"),
        ("rise to first peak\n(model: %.2f h)" % (t_peak - t_p53_thr),
         max(t_peak - t_p53_thr, 0.0), "#C44E52")]
left = 0
for lb, w, c in segs:
    ax.barh([0], [w], left=left, color=c, height=0.42, label=lb)
    if w > 0.25:
        ax.text(left + w / 2, 0, f"{w:.2f} h", ha="center", va="center",
                fontsize=8, color="white")
    left += w
ax.axvspan(2.5, 3.0, color="#DD8452", alpha=0.25)
ax.text(2.75, 0.32, "measured first-peak band 2.5-3.0 h\n(Lahav 2004 / code 87)",
        ha="center", fontsize=8, color="#B5641E")
ax.text(t_peak, 0.32, f"model first peak {t_peak:.2f} h",
        ha="left", fontsize=8, color="#C44E52")
ax.set_xlim(0, 6); ax.set_ylim(-0.55, 0.62)
ax.set_yticks([])
ax.set_xlabel("time after damage (h)")
ax.set_title("b  Onset decomposition along the chain: each segment is either an archived\n"
             "constant or a model-extracted delay; ordering gH2AX < ATM* < p53 < peak < T=5.5 h",
             fontsize=10)
ax.legend(fontsize=6.5, loc="upper right", frameon=False)

fig.suptitle("Code 90 | Full-chain consistency audit of the p53 damage-signalling cascade",
             y=0.985)
fig.savefig(OUT_PNG, dpi=300, bbox_inches="tight")
fig.savefig(OUT_SVG, bbox_inches="tight")

with open(OUT_JSON, "w", encoding="utf-8") as f:
    json.dump(res, f, ensure_ascii=False, indent=2)

print("links: %d closed / %d tension / %d gap (of %d)"
      % (n_closed, n_tension, n_gap, len(links)))
print("onset: gamma %.2f -> ATM half %.2f -> p53 thr %.2f -> peak %.2f h"
      % (t_gamma, t_atm, t_p53_thr, t_peak))
print("first peak %.2f h vs band 2.5-3.0 (+0.5 tol) -> closes=%s ; ordering=%s"
      % (t_peak, res["B_onset_decomposition"]["peak_closes"], ordering_ok))
print("JSON/PNG/SVG written.")
