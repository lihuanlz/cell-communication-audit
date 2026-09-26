#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
code47b_P5erratum重裁_GPCRcalciumtransient.py
version v1.0.0（2026-08-14）

[origin of the erratum] a post-ruling audit triggered by the user's recheck request: using the authors' SD3 peak table as an independent ground-truth cross-check, a feature-implementation defect of
code 47 v0.2.0 was confirmed - the window peak referenced the pre-experiment mu_0, counting the paper-documented full-trace background drift
(+0.25) and flow-switch micro-traces (synchronous across all cells, 0.01-0.03, sitting on 3*sigma_0~0.011) as "responses":
D1 response rate falsely reported at 87.5% (authors' convention ~6%); the reference channel (mu_0~0.98 flat line) slipped through and was counted as a perfect responder.
The as-run ruling (intermediate) is void-by-erratum and recorded in the ledger; this script re-adjudicates per the pre-registration v0.4 section 11 frozen spec:
local reference (median of 25-5 s before the window) + response threshold max(k*sigma_0, 0.04) + within-experiment mu_0 robust outlier removal (>3xMAD).
The verdict table / units / AUC directions / statistic criteria / timing pipeline (v0.3 section 6) are untouched.

[version history] v0.1.0 first version -> v0.1.1: pulse-timing revision (pre-registration bumped to v0.2, zero-statistic window) -
the original "detect 35 peaks directly on the median trace" was infeasible (background drift drowned low-dose peaks; the low-dose median trace has no peaks at all;
the first run guardrail aborted = empty run #1); changed to de-drift residual peak detection on high doses + anchored 120 s uniform-grid extrapolation.
-> v0.2.0: pre-registration bumped to v0.3 (empty run #2 "detected peaks not all on grid" triggered deep diagnosis) - the uniform-grid hypothesis was falsified
(dose-switch boundaries at irregular 85-240 s intervals, low-dose blocks invisible at population level, terminal ionomycin extra peak);
section 6 fully replaced by a two-stage pipeline "per-cell peak-time clustering -> 35-position chain assembly", including per-experiment timing exclusions (pure metadata criteria,
exclusion reasons logged; <8 experiments admitted or <150 valid cells aborts to an empty run; [150,200) carries a power-reduction qualifier).
Criteria / units / features / clauses / directions remain untouched.

[case] P5: first blind test of the boundary predictor - Keshelava 2018 (Nat Commun 9:876) M3R-GPCR -> Ca2+
transient spikes, 7 ascending ACh doses x 5 pulses each, within-cell design, 27 experiments / 433 cells (paper convention).

[frozen basis] 预注册_P5_边界预测器首测_GPCR钙瞬态_v01.md (entry 89, frozen 2026-08-14).
All criteria, unit definitions, feature definitions and guardrails of this script correspond verbatim to that pre-registration; any revision can only happen within
a zero-statistic window (guardrail abort = empty run) with full trace left (code-46 three-empty-run precedent).

[execution site] the user explicitly instructed on 2026-08-14 that our sandbox executes (the original "user machine executes" constraint was changed by the user and is on file).

[frozen parameter block D1-D8]
D1 data: 27 experiment_XX.dot (rows = time frames at 1 frame/s, columns = trace1..N, no header, ASCII tab-separated)
D2 dose mapping: pulses 1-5=D1(100nM) ... 31-35=D7(10uM); 7 ascending doses; 5 pulses per dose
D3 pulse timing (pre-registration v0.3 section 6): per-cell peak times (sigma_diff convention, prominence>=4*sigma_diff, distance>=90)
   -> 5 s binned population clustering (threshold max(3, ceil(0.08N)), split clusters <60 s merged weighted) -> T_hat = median of regular intervals in [115,125]
   -> anchor = the first ascending cluster passing (tail-window / transient fallback / precursor / ionomycin discrimination) -> walk back 34 steps (strict layer [85,170] -> midpoint vacancy check
   -> relaxed layer (170,250], score 2*count - 0.1*|gap - T_hat|) -> chain validation (no cluster at positions 36/37 +/-30 s, no empty holes in the D7 block,
   holes <=12 and <=3 consecutive, baseline p1-8 >=35 s, matches >=23); onset = peak position - 8 s; timing exclusions registered per experiment
D4 features (k=3; v0.4 erratum convention): baseline = all frames before the first onset; mu_0/sigma_0 = baseline mean/std; pulse window = [onset, onset+30 s];
   base_p = median of 25-5 s before the window; peak_p = window max - base_p; responded_p = (peak_p >= max(k*sigma_0, 0.04));
   tau_p = time of first crossing of base_p + max(k*sigma_0, 0.04) within the window minus onset (no crossing = inf); within-experiment mu_0 outliers (>3xMAD) removed;
   per (cell, dose): count = number responded (0-5); peak = median of the 5 pulse peak_p; tau = median of responded pulse tau_p (all-none = inf)
D5 units: low zone (D1,D2)(D2,D3) + high zone (D5,D6)(D6,D7) = 4 ruling units; mid zone (D3,D4)(D4,D5) descriptive arm
D6 AUC directions: tau-AUC = P(tau_low > tau_high) + 0.5 ties; count/peak-AUC = P(high > low) + 0.5 ties; MW average ranks, inf censoring as ties
D7 criteria: P5-1 hit = 4/4 units tau-AUC in [0.40,0.60), falsified = >=2 units outside [0.40,0.60];
   P5-2 hit = 4/4 count-AUC >= 0.60, falsified = 4/4 < 0.60; P5-3 calibration = >=3/4 peak-AUC >= 0.60 (control arm, not in the overall ruling);
   overall: supported = P5-1 hit AND P5-2 hit; falsified = either side falsified; otherwise intermediate
D8 statistics: seed 20260815; per-cell paired percentile bootstrap x2000; N_MIN=30; engine unit-tested bit-for-bit against scipy before delivery

[cell validity (independent of the authors' config, frozen)] sigma_0 > 0; baseline frames >= 30; all values finite; covering all 35 pulse windows.
[robust arms] k=2 / k=4; desensitization correction (subtracting the global median drift by pulse position within block); largest-experiment purity arm.
[output (main products before the descriptive arm)] 单元table.csv / robust臂.csv / ruling.json / 判定日志.txt / 描述臂.csv
"""

import os
import sys
import json
import numpy as np

ROOT = "/mnt/agents/output/01_细胞线/公开数据/Keshelava2018/SD1/source_data_1"
OUTDIR = "/mnt/agents/output/01_细胞线/结果/Keshelava2018_P5_47berratum重裁"
os.makedirs(OUTDIR, exist_ok=True)
LOG_PATH = os.path.join(OUTDIR, "代码47b_判定日志.txt")
_log = open(LOG_PATH, "w", encoding="utf-8")


def say(msg):
    print(msg)
    _log.write(str(msg) + "\n")
    _log.flush()


def fail(msg):
    say("[ABORT] " + msg)
    _log.close()
    sys.exit(2)


SEED = 20260815
N_BOOT = 2000
N_MIN = 30
K_MAIN = 3
K_ROB = (2, 4)
M_AUC = 0.60
BAND_LO, BAND_HI = 0.40, 0.60
DOSES_NM = [100.0, 250.0, 500.0, 750.0, 1500.0, 3000.0, 10000.0]
PPD = 5                 # pulses per dose
N_DOSE = 7
N_PULSE = 35
PERIOD = 120.0          # s
WIN = 30                # s response window from onset
PEAK_LAG = 8            # onset = detected peak position − 8 s
MIN_BASELINE = 35       # experiment level: baseline before the first onset (s); below -> timing exclusion (pre-registration v0.3 section 6-F)
MIN_EXP = 8             # admitted-experiment lower bound (pre-registration v0.3 section 7 item 6)
MIN_CELLS = 150         # valid-cell lower bound; [150,200) carries the power-reduction qualifier
UNITS = [(0, 1), (1, 2), (4, 5), (5, 6)]      # 0-based; low zone (D1,D2)(D2,D3), high zone (D5,D6)(D6,D7)
MID_UNITS = [(2, 3), (3, 4)]                   # descriptive arm

say("code 47 v0.2.0 - P5 blind ruling (boundary predictor first test, Keshelava 2018 GPCR->Ca2+; pre-registration v0.3)")
say(f"frozen parameters: SEED={SEED} N_BOOT={N_BOOT} k={K_MAIN} sufficiency line={M_AUC} dead band=[{BAND_LO},{BAND_HI})")


# ============ AUC engine (same implementation as codes 44-46: MW average ranks, inf censoring as ties) ============
def auc_mw(x, y):
    """P(X>Y)+0.5 P(X=Y), average ranks. x = the direction-positive group."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    v = np.concatenate([x, y])
    order = np.argsort(v, kind="mergesort")
    ranks = np.empty(len(v), dtype=float)
    i = 0
    while i < len(v):
        j = i
        while j + 1 < len(v) and v[order[j + 1]] == v[order[i]]:
            j += 1
        ranks[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    nx = len(x)
    u = ranks[:nx].sum() - nx * (nx + 1) / 2.0
    return u / (nx * len(y))


def engine_selftest():
    """guardrail 7: bit-for-bit check against scipy + inf-censoring direction + all-censored tie"""
    from scipy.stats import mannwhitneyu
    rng = np.random.default_rng(0)
    for _ in range(5):
        a = rng.normal(0, 1, 40)
        b = rng.normal(0.3, 1, 35)
        u, _ = mannwhitneyu(a, b, alternative="two-sided")
        ref = u / (len(a) * len(b))
        if abs(auc_mw(a, b) - ref) > 1e-12:
            fail("AUC engine inconsistent with scipy (no censoring)")
    a = [1.0, 2.0, np.inf, np.inf]
    b = [0.5, 1.5, 3.0, np.inf]
    # by hand: P(a>b)=4/16 (a3,a4 > b1,b2), tie a1? pair by pair: (1,.5)W (1,1.5)L (1,3)L (1,inf)L
    # (2,.5)W (2,1.5)W (2,3)L (2,inf)L (inf,.5)W (inf,1.5)W (inf,3)W (inf,inf)T ×2
    # W=5? recount: a3=inf > b1,b2,b3 = 3W +T; a4=inf same 3W+T -> W=1+2+3+3=9? a1: 1>.5=1W, other 3L; a2: 2>.5,1.5=2W, 2<3, 2<inf
    # → W=1+2+3+3=9，T=2，AUC=(9+1)/16=0.625
    if abs(auc_mw(a, b) - 0.625) > 1e-12:
        fail(f"AUC engine inf-censoring direction error (got {auc_mw(a, b)})")
    if abs(auc_mw([np.inf, np.inf], [np.inf, np.inf]) - 0.5) > 1e-12:
        fail("AUC engine all-censored tie error")
    say("guardrail 7 AUC engine scipy bit-for-bit unit test passed (no-censoring 5/5 + inf-censoring direction + all-censored tie)")


# ============ data loading + guardrail chain ============
def load_all():
    files = sorted(f for f in os.listdir(ROOT) if f.endswith(".dot"))
    if len(files) != 27:
        fail(f"experiment file count {len(files)} != 27")
    exps = []
    tot_cols = 0
    for f in files:
        arr = np.loadtxt(os.path.join(ROOT, f), delimiter="\t")
        if arr.ndim != 2:
            fail(f"{f} is not a 2D matrix")
        t, c = arr.shape
        if not (4100 <= t <= 4800):
            fail(f"{f} row count {t} out of bounds [4100,4800]")
        if not (5 <= c <= 40):
            fail(f"{f} column count {c} out of bounds [5,40]")
        tot_cols += c
        exps.append((f, arr))
    if not (430 <= tot_cols <= 560):
        fail(f"total trace columns {tot_cols} out of bounds [430,560]")
    say(f"guardrails 1/2/3 passed: 27 experiments, total trace columns {tot_cols}")
    return exps


def pulse_timeline(name, arr):
    """D3 (v0.2.0, pre-registration v0.3 section 6): two-stage pipeline - per-cell peak-time clustering -> 35-position chain assembly.
    Pure timing metadata (peak times + population aggregate trace); produces no per-cell x dose readout statistic.
    Returns (onsets, info); timing exclusion returns (None, reason) - exclusions are pure metadata criteria, logged item by item."""
    from scipy.signal import find_peaks, medfilt
    nrows = arr.shape[0]
    # A. per-cell peak detection (times only)
    allp, ncell = [], 0
    for c in range(arr.shape[1]):
        tr = arr[:, c]
        sd = np.median(np.abs(np.diff(tr))) * 1.4826 / np.sqrt(2)
        if sd <= 0:
            continue
        ncell += 1
        pk, _ = find_peaks(tr, prominence=4 * sd, distance=90)
        allp.extend(pk.tolist())
    if not allp:
        return None, "none可用峰"
    # B. population clustering (5 s bins -> strong bins -> merge clusters <=10 s -> split clusters <60 s merged count-weighted)
    allp = np.sort(allp)
    edges = np.arange(0, nrows + 5, 5)
    h, _ = np.histogram(allp, bins=edges)
    ctr = (edges[:-1] + edges[1:]) / 2
    thr = max(3, int(np.ceil(0.08 * ncell)))
    strong = ctr[h >= thr]
    groups = []
    for c_ in strong:
        if groups and c_ - groups[-1][-1] <= 10:
            groups[-1].append(c_)
        else:
            groups.append([c_])
    raw = [(float(np.mean(g)), float(sum(h[np.searchsorted(ctr, x)] for x in g))) for g in groups]
    merged = []
    for pos, cnt in raw:
        if merged and pos - merged[-1][0] < 60:
            p0, c0 = merged[-1]
            merged[-1] = ((p0 * c0 + pos * cnt) / (c0 + cnt), c0 + cnt)
        else:
            merged.append((pos, cnt))
    if len(merged) < 20:
        return None, f"簇{len(merged)}<20"
    cs = np.array([p for p, _ in merged])
    cnt = np.array([c for _, c in merged])
    med = np.median(arr, axis=1)
    res = med - medfilt(med, 301)

    def transient_ok(A):
        A = int(A)
        if A + 100 >= nrows:
            return False
        during = res[A + 4:A + 28].mean()
        after = res[A + 40:A + 100].mean()
        return after < 0.5 * max(during, 1e-9) or (during < 0.02 and after < 0.02)

    # C. period estimation
    d2 = np.diff(cs)
    reg = d2[(d2 >= 100) & (d2 <= 130)]
    if len(reg) < 5:
        return None, "规则间隔<5"
    That = float(np.median(reg))
    if not (115.0 <= That <= 125.0):
        return None, f"T̂={That:.1f}越界"

    # E. walk-back (strict layer [85,170] -> midpoint vacancy check -> relaxed layer (170,250])
    def walkback(ai):
        p = [0] * N_PULSE
        hole = [False] * N_PULSE
        p[34] = int(cs[ai])
        for i in range(34, 0, -1):
            cand = list(np.where((cs >= p[i] - 170) & (cs <= p[i] - 85))[0])
            if len(cand) == 0:
                if np.any(np.abs(cs - (p[i] - That)) <= 30):
                    return None, None
                cand = list(np.where((cs >= p[i] - 250) & (cs < p[i] - 170))[0])
            best, bs = None, 0.0
            for ci in cand:
                sc = 2.0 * cnt[ci] - 0.1 * abs((p[i] - cs[ci]) - That)
                if sc > bs:
                    best, bs = ci, sc
            if best is None:
                p[i - 1] = p[i] - int(round(That))
                hole[i - 1] = True
            else:
                p[i - 1] = int(cs[best])
        return p, hole

    # D. anchor position (first ascending full pass) + F. chain validation
    for ai in np.argsort(cs):
        A = cs[ai]
        if A + 105 > nrows:
            continue
        if not transient_ok(A):
            continue
        if not np.any((A - cs >= 85) & (A - cs <= 170)):
            continue
        fol = np.where((cs - A > 85) & (cs - A <= 180))[0]
        if any(transient_ok(cs[fi]) for fi in fol):
            continue
        p, hole = walkback(ai)
        if p is None:
            continue
        if np.any(np.abs(cs - (p[0] - That)) <= 30):
            continue
        if np.any(np.abs(cs - (p[0] - 2 * That)) <= 30):
            continue
        if any(hole[30:35]):
            continue
        nh = sum(hole)
        if nh > 12:
            continue
        cons = mx = 0
        for hf in hole:
            cons = cons + 1 if hf else 0
            mx = max(mx, cons)
        if mx > 3:
            continue
        if p[0] - PEAK_LAG < MIN_BASELINE:
            continue
        if N_PULSE - nh < 23:
            continue
        onsets = [int(pi) - PEAK_LAG for pi in p]
        info = {"That": That, "holes": nh, "p1": p[0], "p35": p[34]}
        return onsets, info
    return None, "noneyes效anchor/链"


# ============ featurization (v0.4 section 11 erratum version: local reference + artifact floor + reference-channel interception) ============
H_FLOOR = 0.04   # absolute floor of the response threshold (Fura-2 ratio units, frozen value: between the artifact-cluster ceiling 0.03 and the true-response-cluster floor 0.05)


def featurize(arr, onsets, k):
    """Returns the per-cell feature table: list[dict(count=[], peak=[], tau=[])], indexed by dose.
    v0.4: peak_p = window max - median of 25-5 s before the window (local reference, removes background drift);
    responded = peak_p >= max(k*sigma_0, H_FLOOR); tau = first crossing of base_p+thr within the window minus onset;
    validity adds within-experiment mu_0 robust outlier removal (>3xMAD, intercepting reference channels)."""
    t_end = onsets[-1] + WIN
    base_end = onsets[0]
    # within-experiment mu_0 outlier pre-screen (reference-channel interception)
    mu_list = []
    for c in range(arr.shape[1]):
        v = arr[:, c]
        if np.all(np.isfinite(v[:t_end])) and base_end >= 30:
            mu_list.append(v[:base_end].mean())
        else:
            mu_list.append(np.nan)
    mu_arr = np.array(mu_list)
    mu_med = np.nanmedian(mu_arr)
    mu_mad = np.nanmedian(np.abs(mu_arr - mu_med)) * 1.4826
    outlier = np.abs(mu_arr - mu_med) > 3 * max(mu_mad, 1e-9)
    cells = []
    n_excl = 0
    n_outlier = 0
    for c in range(arr.shape[1]):
        v = arr[:, c]
        if not np.all(np.isfinite(v[:t_end])):
            n_excl += 1
            continue
        b = v[:base_end]
        if len(b) < 30:
            n_excl += 1
            continue
        mu0 = b.mean()
        sd0 = b.std(ddof=1)
        if not (sd0 > 0):
            n_excl += 1
            continue
        if outlier[c]:
            n_outlier += 1
            continue
        thr_abs = max(k * sd0, H_FLOOR)
        cnt = np.zeros(N_DOSE)
        pk = np.zeros(N_DOSE)
        ta = np.full(N_DOSE, np.inf)
        pk_by_dose = [[] for _ in range(N_DOSE)]
        ta_by_dose = [[] for _ in range(N_DOSE)]
        for p, o in enumerate(onsets):
            w = v[o:o + WIN]
            d = p // PPD
            base_p = np.median(v[o - 25:o - 5])
            peak = w.max() - base_p
            pk_by_dose[d].append(peak)
            if peak >= thr_abs:
                cnt[d] += 1
                cross = np.where(w >= base_p + thr_abs)[0]
                ta_by_dose[d].append(float(cross[0]) if len(cross) else np.inf)
        for d in range(N_DOSE):
            pk[d] = np.median(pk_by_dose[d])
            ta[d] = np.median(ta_by_dose[d]) if ta_by_dose[d] else np.inf
        cells.append({"cnt": cnt, "pk": pk, "ta": ta})
    return cells, n_excl + n_outlier


# ============ unit AUC + bootstrap ============
def unit_auc(cells, pair, stat, rng=None, boot=False):
    ia, ib = pair
    if stat == "ta":
        xa = np.array([c["ta"][ia] for c in cells])
        xb = np.array([c["ta"][ib] for c in cells])
        # direction: P(tau_low > tau_high) -> x = low-dose group
        x, y = xa, xb
    else:
        key = "cnt" if stat == "cnt" else "pk"
        x = np.array([c[key][ib] for c in cells])   # high-dose group is the positive direction
        y = np.array([c[key][ia] for c in cells])
    n = len(x)
    out = {"auc": auc_mw(x, y), "n": n}
    if boot:
        aucs = np.empty(N_BOOT)
        idx_all = np.arange(n)
        for b in range(N_BOOT):
            sel = rng.choice(idx_all, size=n, replace=True)
            aucs[b] = auc_mw(x[sel], y[sel])
        out["lo"] = float(np.percentile(aucs, 2.5))
        out["hi"] = float(np.percentile(aucs, 97.5))
    return out


def run_units(cells, pairs, rng, boot):
    res = {}
    for stat in ("ta", "cnt", "pk"):
        res[stat] = []
        for pair in pairs:
            r = unit_auc(cells, pair, stat, rng=rng, boot=boot)
            r["pair"] = f"D{pair[0]+1}↔D{pair[1]+1}"
            res[stat].append(r)
    return res


# ============ main pipeline ============
def main():
    engine_selftest()
    exps = load_all()

    rng = np.random.default_rng(SEED)
    # ---- pulse time tables (pure timing-metadata stage; timing exclusions registered item by item) ----
    timelines = {}
    timing_excl = []
    for name, arr in exps:
        onsets, info = pulse_timeline(name, arr)
        if onsets is None:
            timing_excl.append({"experiment": name, "reason": info})
            say(f"timingexcluded：{name} —— {info}")
            continue
        timelines[name] = onsets
        say(f"timing admitted: {name} (T_hat={info['That']:.1f} holes {info['holes']} p1={info['p1']} p35={info['p35']})")
    say(f"guardrails 4/5 done: admitted {len(timelines)}/27, timing-excluded {len(timing_excl)}")
    if len(timelines) < MIN_EXP:
        fail(f"admitted experiments {len(timelines)} < {MIN_EXP}")

    all_cells = []
    per_exp_cells = []
    tot_excl = 0
    for name, arr in exps:
        if name not in timelines:
            continue
        cells, n_excl = featurize(arr, timelines[name], K_MAIN)
        tot_excl += n_excl
        per_exp_cells.append((name, cells))
        all_cells.extend(cells)
    power_note = ""
    if len(all_cells) < MIN_CELLS:
        fail(f"valid cells {len(all_cells)} < {MIN_CELLS}")
    if len(all_cells) < 200:
        power_note = f"power-reduction qualifier: valid cells {len(all_cells)} in [150,200)"
        say("[QUALIFIER] " + power_note)
    say(f"loading done: valid cells {len(all_cells)} (invalid excluded {tot_excl}; timing-excluded {len(timing_excl)} experiments; paper convention 433, discrepancies from the independent validity rules reported as-is)")

    # ---- main analysis (4 ruling units x 3 statistics, with CI) ----
    main_res = run_units(all_cells, UNITS, rng, boot=True)

    say("\n===== main analysis: 4 ruling units (low zone D1<->D2, D2<->D3; high zone D5<->D6, D6<->D7) =====")
    header = f"{'stat':<8}{'unit':<10}{'AUC':>8}{'CI95':>22}{'n':>6}"
    say(header)
    rows = []
    for stat, label in (("ta", "tau event-time"), ("cnt", "count"), ("pk", "peak (paid)")):
        for r in main_res[stat]:
            say(f"{label:<8}{r['pair']:<10}{r['auc']:>8.4f}   [{r['lo']:.4f},{r['hi']:.4f}]{r['n']:>6}")
            rows.append({"stat": stat, "pair": r["pair"], "auc": r["auc"],
                         "lo": r["lo"], "hi": r["hi"], "n": r["n"]})

    # ---- D7 verdicts (frozen table executed verbatim) ----
    ta_aucs = [r["auc"] for r in main_res["ta"]]
    cnt_aucs = [r["auc"] for r in main_res["cnt"]]
    pk_aucs = [r["auc"] for r in main_res["pk"]]

    n_in = sum(1 for a in ta_aucs if BAND_LO <= a < BAND_HI)
    n_out = 4 - n_in
    if n_in == 4:
        v51 = "hit（死区确认：4/4 单元 τ-AUC∈[0.40,0.60)）"
    elif n_out >= 2:
        v51 = f"证伪（{n_out}/4 单元出 [0.40,0.60] 带）"
    else:
        v51 = f"中间态（{n_in}/4 单元在带内）"

    n_pass = sum(1 for a in cnt_aucs if a >= M_AUC)
    if n_pass == 4:
        v52 = "hit（活区确认：4/4 单元 count-AUC≥0.60）"
    elif n_pass == 0:
        v52 = "证伪（4/4 单元 count-AUC<0.60）"
    else:
        v52 = f"中间态（{n_pass}/4 单元达线）"

    n_cal = sum(1 for a in pk_aucs if a >= M_AUC)
    v53 = f"校准pass（{n_cal}/4 单元 peak-AUC≥0.60）" if n_cal >= 3 else f"数据集级弱decoding限定（仅 {n_cal}/4 单元 peak-AUC≥0.60）"

    hit51 = n_in == 4
    fal51 = n_out >= 2
    hit52 = n_pass == 4
    fal52 = n_pass == 0
    if hit51 and hit52:
        overall = "P5 获support（boundaryprediction器：死区recheckhit ∧ 活区首测hit）"
    elif fal51 or fal52:
        overall = "P5 证伪（boundaryprediction器在this形态上死亡，死亡侧见分款）"
    else:
        overall = "P5 中间态"

    say("\n===== verdicts (frozen table D7 executed verbatim) =====")
    say(f"P5-1 (R1 dead-zone recheck, tau): {v51}  unit values {['%.4f' % a for a in ta_aucs]}")
    say(f"P5-2 (R2 live-zone first test, count): {v52}  unit values {['%.4f' % a for a in cnt_aucs]}")
    say(f"P5-3 (R3 calibration, paid peak, control arm): {v53}  unit values {['%.4f' % a for a in pk_aucs]}")
    say(f"[overall ruling] {overall}")

    # ---- main products saved (before descriptive/robust-arm conclusive content) ----
    import csv
    with open(os.path.join(OUTDIR, "代码47b_单元table.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["stat", "pair", "auc", "lo", "hi", "n"])
        w.writeheader()
        w.writerows(rows)
    verdict = {
        "case": "P5", "dataset": "Keshelava2018_NatCommun9_876", "seed": SEED,
        "n_boot": N_BOOT, "k": K_MAIN, "n_cells": len(all_cells), "n_excluded": tot_excl,
        "experiments_included": len(timelines),
        "timing_exclusions": timing_excl,
        "power_note": power_note,
        "units": rows,
        "verdicts": {"P5-1": v51, "P5-2": v52, "P5-3": v53, "overall": overall},
    }
    with open(os.path.join(OUTDIR, "代码47b_ruling.json"), "w", encoding="utf-8") as f:
        json.dump(verdict, f, ensure_ascii=False, indent=2,
                  default=lambda o: int(o) if isinstance(o, np.integer) else float(o) if isinstance(o, np.floating) else str(o))
    say("main products saved: 单元table.csv / ruling.json / 判定日志.txt")

    # ---- robust arms ----
    say("\n===== robust arms =====")
    rob_rows = []
    for kk in K_ROB:
        cells_k = []
        for name, arr in exps:
            if name not in timelines:
                continue
            ck, _ = featurize(arr, timelines[name], kk)
            cells_k.extend(ck)
        res_k = run_units(cells_k, UNITS, rng, boot=False)
        for stat in ("ta", "cnt"):
            for r in res_k[stat]:
                rob_rows.append({"arm": f"k={kk}", "stat": stat, "pair": r["pair"], "auc": r["auc"], "n": r["n"]})
                say(f"k={kk} {stat} {r['pair']} AUC={r['auc']:.4f}")
    # desensitization-correction arm: subtract the global median peak drift by pulse position within block (same as the paper), recompute cnt/pk
    # collect all (cell, dose, j) peaks
    say("desensitization-correction arm: recompute after subtracting the population median drift by pulse position j in {1..5}")
    # when re-featurizing keep per-pulse peaks - simplified: here only pk/cnt are approximately corrected (tau is unaffected by amplitude drift)
    # global position median
    # (strict implementation: featurize does not return per-pulse values; this arm re-featurizes with k=3 and corrects at the pulse level)
    cells_dc = []
    for name, arr in exps:
        if name not in timelines:
            continue
        onsets = timelines[name]
        t_end = onsets[-1] + WIN
        base_end = onsets[0]
        # within-experiment mu_0 outlier pre-screen (same convention as featurize)
        mu_arr = np.array([v[:base_end].mean() if np.all(np.isfinite(arr[:t_end, ci])) else np.nan
                           for ci, v in enumerate(arr.T)])
        mu_med = np.nanmedian(mu_arr)
        mu_mad = np.nanmedian(np.abs(mu_arr - mu_med)) * 1.4826
        outlier = np.abs(mu_arr - mu_med) > 3 * max(mu_mad, 1e-9)
        # first pass: collect all pulse peaks to estimate the position drift (local reference)
        peaks_j = [[] for _ in range(PPD)]
        valid_cols = []
        for c in range(arr.shape[1]):
            v = arr[:, c]
            if not np.all(np.isfinite(v[:t_end])) or outlier[c]:
                continue
            b = v[:base_end]
            mu0 = b.mean()
            sd0 = b.std(ddof=1)
            if not (sd0 > 0) or len(b) < 30:
                continue
            pps = []
            bases = []
            for p, o in enumerate(onsets):
                bases.append(np.median(v[o - 25:o - 5]))
                pps.append(v[o:o + WIN].max() - bases[-1])
            for p, val in enumerate(pps):
                peaks_j[p % PPD].append(val)
            valid_cols.append((v, mu0, sd0, pps, bases))
        med_j = np.array([np.median(x) for x in peaks_j])
        corr = med_j - med_j.mean()
        for v, mu0, sd0, pps, bases in valid_cols:
            thr_abs = max(K_MAIN * sd0, H_FLOOR)
            cnt = np.zeros(N_DOSE)
            pk = np.zeros(N_DOSE)
            ta = np.full(N_DOSE, np.inf)
            pk_d = [[] for _ in range(N_DOSE)]
            ta_d = [[] for _ in range(N_DOSE)]
            for p, o in enumerate(onsets):
                d = p // PPD
                val = pps[p] - corr[p % PPD]
                pk_d[d].append(val)
                if val >= thr_abs:
                    cnt[d] += 1
                    w = v[o:o + WIN]
                    cross = np.where(w >= bases[p] + thr_abs)[0]
                    ta_d[d].append(float(cross[0]) if len(cross) else np.inf)
            for d in range(N_DOSE):
                pk[d] = np.median(pk_d[d])
                ta[d] = np.median(ta_d[d]) if ta_d[d] else np.inf
            cells_dc.append({"cnt": cnt, "pk": pk, "ta": ta})
    res_dc = run_units(cells_dc, UNITS, rng, boot=False)
    for stat in ("ta", "cnt"):
        for r in res_dc[stat]:
            rob_rows.append({"arm": "脱敏校正", "stat": stat, "pair": r["pair"], "auc": r["auc"], "n": r["n"]})
            say(f"desensitization-corrected {stat} {r['pair']} AUC={r['auc']:.4f}")
    # largest-experiment purity arm
    big = max(per_exp_cells, key=lambda x: len(x[1]))
    res_big = run_units(big[1], UNITS, rng, boot=False)
    for stat in ("ta", "cnt"):
        for r in res_big[stat]:
            rob_rows.append({"arm": f"单实验({big[0]},n={len(big[1])})", "stat": stat, "pair": r["pair"], "auc": r["auc"], "n": r["n"]})
            say(f"single-experiment purity arm {stat} {r['pair']} AUC={r['auc']:.4f}")
    with open(os.path.join(OUTDIR, "代码47b_robust臂.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["arm", "stat", "pair", "auc", "n"])
        w.writeheader()
        w.writerows(rob_rows)
    say("robust arms saved: robust臂.csv")

    # ---- descriptive arm (not adjudicated) ----
    say("\n===== descriptive arm (not adjudicated) =====")
    desc_rows = []
    mid_res = run_units(all_cells, MID_UNITS, rng, boot=True)
    for stat in ("ta", "cnt", "pk"):
        for r in mid_res[stat]:
            desc_rows.append({"section": "中区对", "stat": stat, "key": r["pair"],
                              "auc": r["auc"], "lo": r["lo"], "hi": r["hi"], "n": r["n"]})
            say(f"mid-zone pair {stat} {r['pair']} AUC={r['auc']:.4f} [{r['lo']:.4f},{r['hi']:.4f}]")
    # per-dose response rate (fraction of cells with count>=1) and mean count
    for d in range(N_DOSE):
        rr = np.mean([1.0 if c["cnt"][d] >= 1 else 0.0 for c in all_cells])
        mc = np.mean([c["cnt"][d] for c in all_cells])
        desc_rows.append({"section": "响应率", "stat": "cnt", "key": f"D{d+1}({DOSES_NM[d]}nM)",
                          "auc": rr, "lo": mc, "hi": np.nan, "n": len(all_cells)})
        say(f"D{d+1}({DOSES_NM[d]:.0f}nM) response rate={rr:.3f} mean count={mc:.2f}")
    with open(os.path.join(OUTDIR, "代码47b_描述臂.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["section", "stat", "key", "auc", "lo", "hi", "n"])
        w.writeheader()
        w.writerows(desc_rows)
    say("descriptive arm saved: 描述臂.csv")

    say("\n===== final verdict restatement (frozen table verbatim) =====")
    say(f"P5-1：{v51}")
    say(f"P5-2：{v52}")
    say(f"P5-3：{v53}")
    say(f"[overall ruling] {overall}")
    say("(recorded either way: this log and ruling.json enter the ledger in sync; the frozen ruling is never moved afterwards)")
    _log.close()


if __name__ == "__main__":
    main()
