# -*- coding: utf-8 -*-
"""
Code 44 · P2 blind ruling: amplitude-stratified decay prediction (Msn2, Hansen & Zechner 2021)
====================================================================
Pre-registration document: 05_主线纲领与设计/预注册_P2_幅度分层衰减预言_v01.md (registered 2026-08-14)
Dataset  : Hansen & Zechner 2021 (MSB), Zenodo doi:10.5281/zenodo.2755026

Criteria (locked at registration; no post-hoc modification allowed):
  P2-1  Δ_mol  = AUC_pooled − mean(AUC_within size terciles) ≥ +0.03, ≥6/7 wild-type promoters
  P2-2  Δ_AU   ≤ 0.01 and per-promoter Δ_AU < Δ_mol, ≥6/7 wild-type promoters
  P2-3  |Δ_time| ≤ 0.01 (k=3 main analysis; k∈{2,4} robust arm), ≥6/7 wild-type promoters
  P2-F  any direction reversal → the corresponding clause is recorded as 证伪, double-recorded in the ledger either way

Implementation conventions not detailed in the pre-registration, fixed by this script before unblinding (equally blind-locked, entered in the ledger registry):
  · Promoter-level Δ = arithmetic mean of the Δ of three adjacent dose pairs (100↔275, 275↔690, 690nM↔3µM);
    the 21 adjudication units (3 pairs × 7 promoters) are also reported pair-by-pair in the unit table
  · Tercile boundaries are taken from the pooled sample of the unit's observations (both doses pooled); during
    bootstrap resampling each cell carries its original stratum label (boundaries are not recomputed per resample)
  · The event-time decoder is based on the YFP(AU) channel: per-cell baseline μ₀+kσ₀ (t≤0 frames),
    first threshold-crossing time point with t>0; no crossing = right censoring (recorded as +inf); AUC uses
    Mann–Whitney average ranks (censored cells tie with each other and rank above all crossers)
  · Missingness rule: YFP or cell_size missing in >10% of frames (≥7 of 64 frames) → excluded
  · Bootstrap (2000 per-cell resamples) CIs are reported only for the Δ of the three main adjudication channels
    (mol_peak / au_peak / time_k3) at the main-analysis duration;
    alternative arms (AUC decoder, k∈{2,4}, alternative proxy stratification, other durations) get point estimates
  · Main-analysis duration = 50 min (available for all promoters; deterministic rule: take 50 if present, otherwise the longest)

How to run (user machine, Windows + Spyder):
  1) Unzip HansenZechner_RawData.zip
  2) Point ROOT below to the unzipped 2018_DataForPaper directory
  3) %runfile this script in Spyder; pure CPU, approx. 5–15 min
Dependencies: numpy, scipy. Random seed 20260814 (locked).
Outputs: 代码44_判定日志.txt, 代码44_单元table.csv, 代码44_ruling.json (written next to this script)
"""
import os
import re
import sys
import json
import time as _time
import numpy as np
from scipy.io import loadmat
from scipy.stats import rankdata

# ---------------- user configuration ----------------
ROOT = r"C:/Users/lihua/Desktop/NC/data/HansenZechner2021/2018_DataForPaper"
if len(sys.argv) > 1 and os.path.isdir(sys.argv[1]):
    ROOT = sys.argv[1]
# ---------------------------------------------

SEED = 20260814
N_BOOT = 2000
WT = ["ALD3", "DCS2", "DDR2", "HXK1", "RTN2", "SIP18", "TKL2"]
MUT = ["pSIP18_mut6", "pSIP18_mut21"]
DOSE_TAG = {100: "100nM", 275: "275nM", 690: "690nM", 3000: "3uM"}
PAIRS = [(100, 275), (275, 690), (690, 3000)]
DURS = ["10min", "20min", "30min", "40min", "50min"]
MAIN_DUR = "50min"
MARGIN_P21, MARGIN_P22, MARGIN_P23 = 0.03, 0.01, 0.01
K_MAIN, K_ROB = 3, (2, 4)

LOG_PATH = None
_log_fh = None


def log(msg=""):
    print(msg)
    if _log_fh is not None:
        _log_fh.write(str(msg) + "\n")
        _log_fh.flush()


# ---------------- data loading and per-cell features ----------------
def cond_path(prom, series, dur, dose_tag):
    return os.path.join(ROOT, prom, f"{prom}_{series}_{dur}_{dose_tag}_size.mat")


def extract(fpath):
    """Load one _size.mat file; return a per-cell feature dict + exclusion counts. All rules deterministic."""
    m = loadmat(fpath)
    t = m["time"].ravel().astype(float)
    Y = m["YFP"].astype(float)            # AU: size-normalized concentration
    M = m["YFP_molecules"].astype(float)  # absolute abundance (not normalized)
    S = m["cell_size_pixels"].astype(float)
    n, nF = Y.shape
    b_idx = np.where(t <= 0)[0]
    w_idx = np.where(t >= 0)[0]
    a_idx = np.where(t > 0)[0]
    dt = float(np.median(np.diff(t)))

    # --- exclusion rules (pre-registration §2.2 + script-locked conventions) ---
    miss_y = np.isnan(Y).sum(1)
    miss_s = np.isnan(S).sum(1)
    keep = (miss_y <= int(0.10 * nF)) & (miss_s <= int(0.10 * nF))  # either channel missing >10% of frames → excluded
    dS = np.abs(np.diff(S, axis=1))
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = dS / np.abs(S[:, :-1])
    rel[~np.isfinite(rel)] = np.inf                 # previous frame 0/NaN counts as a jump
    jump = rel.max(1)
    keep &= ~(jump > 0.40)                          # relative size jump between adjacent frames >40%

    mu0 = np.nanmean(Y[:, b_idx], axis=1)
    sd0 = np.nanstd(Y[:, b_idx], axis=1, ddof=1)
    keep &= np.isfinite(mu0) & np.isfinite(sd0) & (mu0 > 0)   # corrected baseline AU ≤ 0 → excluded

    # --- per-cell features ---
    mean_size = np.nanmean(S, axis=1)
    base_mol = np.nanmean(M[:, b_idx], axis=1)      # alternative stratification proxy
    peak_mol = np.nanmax(M[:, w_idx], axis=1)
    peak_au = np.nanmax(Y[:, w_idx], axis=1)
    auc_mol = np.nansum(M[:, w_idx], axis=1) * dt   # alternative amplitude decoder
    auc_au = np.nansum(Y[:, w_idx], axis=1) * dt

    def cross(k):
        thr = (mu0 + k * sd0)[:, None]
        above = Y[:, a_idx] > thr                   # NaN comparison is False, inherently safe
        has = above.any(1)
        first = above.argmax(1)
        return np.where(has, t[a_idx][first], np.inf)

    feat = {
        "mean_size": mean_size, "base_mol": base_mol,
        "mol_peak": peak_mol, "au_peak": peak_au,
        "mol_auc": auc_mol, "au_auc": auc_au,
        "time_k3": cross(K_MAIN),
    }
    for k in K_ROB:
        feat[f"time_k{k}"] = cross(k)
    feat["responded"] = np.isfinite(feat["time_k3"]).astype(float)
    n_excl = int((~keep).sum())
    return {k: v[keep] for k, v in feat.items()}, n, n_excl


# ---------------- statistics: Mann–Whitney AUC and stratification-penalized Δ ----------------
def auc_mw(lo, hi):
    """AUC = P(hi > lo) + 0.5·P(tie); average ranks for ties; NaN dropped, +inf (censored) kept."""
    lo = lo[~np.isnan(lo)]
    hi = hi[~np.isnan(hi)]
    nL, nH = len(lo), len(hi)
    if nL < 2 or nH < 2:
        return np.nan
    r = rankdata(np.concatenate([lo, hi]))
    return float((r[nL:].sum() - nH * (nH + 1) / 2.0) / (nL * nH))


def tercile_labels(x):
    """Split into three strata by pooled-sample 1/3 and 2/3 quantiles; return 0/1/2 labels."""
    q1, q2 = np.quantile(x, [1.0 / 3.0, 2.0 / 3.0])
    return (x > q1).astype(int) + (x > q2).astype(int)


def delta_strat(v_lo, v_hi, l_lo, l_hi):
    """Δ = AUC_pooled − mean(AUC_within strata); skip a stratum if either group has <2 valid cells (<2 usable strata → NaN)."""
    ap = auc_mw(v_lo, v_hi)
    ts = []
    for g in range(3):
        a = auc_mw(v_lo[l_lo == g], v_hi[l_hi == g])
        if not np.isnan(a):
            ts.append(a)
    if np.isnan(ap) or len(ts) < 2:
        return np.nan, ap, np.nan
    return ap - float(np.mean(ts)), ap, float(np.mean(ts))


def boot_delta(v_lo, v_hi, l_lo, l_hi, rng):
    """Per-cell bootstrap, 2000 resamples; cells carry their original stratum labels (locked convention)."""
    nL, nH = len(v_lo), len(v_hi)
    out = np.full(N_BOOT, np.nan)
    for b in range(N_BOOT):
        iL = rng.integers(0, nL, nL)
        iH = rng.integers(0, nH, nH)
        d, _, _ = delta_strat(v_lo[iL], v_hi[iH], l_lo[iL], l_hi[iH])
        out[b] = d
    ok = np.isfinite(out)
    if ok.sum() < 100:
        return np.nan, np.nan
    return tuple(np.percentile(out[ok], [2.5, 97.5]))


# ---------------- main analysis pipeline ----------------
CHANNELS_MAIN = ["mol_peak", "au_peak", "time_k3"]
CHANNELS_ROB = ["time_k2", "time_k4", "mol_auc", "au_auc"]


def run_condition(prom, series, dur, dose):
    p = cond_path(prom, series, dur, DOSE_TAG[dose])
    if not os.path.isfile(p):
        return None
    return extract(p)


def promoter_delta(rows, prom, ch, strat="size", dur=MAIN_DUR):
    """Promoter-level Δ = arithmetic mean of the three dose-pair Δ; any unavailable pair → NaN (counted as fail: conservative, deterministic)."""
    ds = [r["delta"] for r in rows
          if r["promoter"] == prom and r["duration"] == dur
          and r["strat"] == strat and r["channel"] == ch]
    ds = [d for d in ds if np.isfinite(d)]
    return float(np.mean(ds)) if len(ds) == len(PAIRS) else np.nan


def sanitize(o):
    if isinstance(o, dict):
        return {k: sanitize(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [sanitize(v) for v in o]
    if isinstance(o, float) and not np.isfinite(o):
        return None
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    return o


# ---------------- descriptive arm (not adjudicated) ----------------
def fm_features(prom, series, dur):
    """FM/FM4 conditions: threshold-crossing count (μ₀+3σ₀, AU channel), first crossing time, peak."""
    p = cond_path(prom, series, dur, "690nM")
    if not os.path.isfile(p):
        return None
    m = loadmat(p)
    t = m["time"].ravel().astype(float)
    Y = m["YFP"].astype(float)
    S = m["cell_size_pixels"].astype(float)
    nF = Y.shape[1]
    b_idx = np.where(t <= 0)[0]
    a_idx = np.where(t > 0)[0]
    w_idx = np.where(t >= 0)[0]
    miss_y = np.isnan(Y).sum(1)
    miss_s = np.isnan(S).sum(1)
    keep = (miss_y <= int(0.10 * nF)) & (miss_s <= int(0.10 * nF))
    dS = np.abs(np.diff(S, axis=1))
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = dS / np.abs(S[:, :-1])
    rel[~np.isfinite(rel)] = np.inf
    keep &= ~(rel.max(1) > 0.40)
    mu0 = np.nanmean(Y[:, b_idx], axis=1)
    sd0 = np.nanstd(Y[:, b_idx], axis=1, ddof=1)
    keep &= np.isfinite(mu0) & np.isfinite(sd0) & (mu0 > 0)
    thr = (mu0 + 3 * sd0)[:, None]
    above = (Y[:, a_idx] > thr).astype(int)       # NaN → False
    n_cross = (np.diff(above, axis=1) == 1).sum(1).astype(float)
    has = above.any(1)
    first = np.where(has, t[a_idx][above.argmax(1)], np.inf)
    peak = np.nanmax(Y[:, w_idx], axis=1)
    out = {"n_cross": n_cross, "first_time": first, "peak_au": peak}
    return {k: v[keep] for k, v in out.items()}


def fm_descriptive(prom):
    contrasts = [
        ("FM", "8_5min", "FM4", "15minINT", "FM8(8pulse) vs FM4-15min(4pulse)"),
        ("FM", "4_5min", "FM4", "15minINT", "FM4(4pulse) vs FM4-15min(4pulse, same count, different interval)"),
        ("FM", "4_5min", "FM4", "75minINT", "FM4(4pulse) vs FM4-75min(4pulse, same count, different interval)"),
    ]
    for s1, d1, s2, d2, label in contrasts:
        f1 = fm_features(prom, s1, d1)
        f2 = fm_features(prom, s2, d2)
        if f1 is None or f2 is None:
            log(f"    [FM descriptive] {prom} {label}: missing file, skipped")
            continue
        a_cnt = auc_mw(f2["n_cross"], f1["n_cross"])
        a_tim = auc_mw(f1["first_time"], f2["first_time"])
        a_pk = auc_mw(f2["peak_au"], f1["peak_au"])
        log(f"    [FM descriptive] {prom} {label}: "
            f"countAUC={a_cnt:.3f} firstTimeAUC={a_tim:.3f} peakAUC={a_pk:.3f} "
            f"(n={len(f1['n_cross'])}/{len(f2['n_cross'])})")


def cfp_probe(prom):
    """CFP channel as a gain-probe check: correlation of baseline CFP with time-averaged cell size (descriptive)."""
    p = cond_path(prom, "DM", MAIN_DUR, "275nM")
    if not os.path.isfile(p):
        return
    m = loadmat(p)
    t = m["time"].ravel().astype(float)
    C = m["CFP"].astype(float)
    S = m["cell_size_pixels"].astype(float)
    b_idx = np.where(t <= 0)[0]
    cb = np.nanmean(C[:, b_idx], axis=1)
    ms = np.nanmean(S, axis=1)
    ok = np.isfinite(cb) & np.isfinite(ms)
    if ok.sum() < 10:
        log(f"    [CFP probe] {prom}: insufficient valid cells")
        return
    r = float(np.corrcoef(cb[ok], ms[ok])[0, 1])
    log(f"    [CFP probe] {prom}: corr(baselineCFP, size)={r:+.3f} (n={int(ok.sum())})")


# ---------------- main entry ----------------
def main():
    t_start = _time.time()
    global LOG_PATH, _log_fh
    outdir = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
    LOG_PATH = os.path.join(outdir, "代码44_判定日志.txt")
    _log_fh = open(LOG_PATH, "w", encoding="utf-8")

    rng = np.random.default_rng(SEED)
    log("=" * 72)
    log("Code 44 · P2 blind ruling: amplitude-stratified decay prediction")
    log(f"pre-registration: 2026-08-14 v0.1 | seed {SEED} | bootstrap {N_BOOT}")
    log(f"data root directory: {ROOT}")
    log(f"criteria: P2-1 Δ_mol≥+{MARGIN_P21}; P2-2 Δ_AU≤{MARGIN_P22} and <Δ_mol; "
        f"P2-3 |Δ_time|≤{MARGIN_P23} (k={K_MAIN}); each ≥6/7 promoters")
    log("=" * 72)

    rows = []
    for prom in WT + MUT:
        for dur in DURS:
            conds = {}
            ok = True
            for dose in (100, 275, 690, 3000):
                r = run_condition(prom, "DM", dur, dose)
                if r is None:
                    log(f"[warning] missing file: {prom} DM {dur} {DOSE_TAG[dose]}, skipping {prom} {dur}")
                    ok = False
                    break
                conds[dose] = r
            if not ok:
                continue
            for dlo, dhi in PAIRS:
                flo, nlo, xlo = conds[dlo]
                fhi, nhi, xhi = conds[dhi]
                for strat_key, strat_tag in (("mean_size", "size"), ("base_mol", "baseMol")):
                    pool = np.concatenate([flo[strat_key], fhi[strat_key]])
                    labels = tercile_labels(pool)
                    nL = len(flo[strat_key])
                    l_lo, l_hi = labels[:nL], labels[nL:]
                    chans = CHANNELS_MAIN + (CHANNELS_ROB if strat_key == "mean_size" else [])
                    for ch in chans:
                        d, ap, mt = delta_strat(flo[ch], fhi[ch], l_lo, l_hi)
                        ci_lo = ci_hi = np.nan
                        if dur == MAIN_DUR and strat_key == "size" and ch in CHANNELS_MAIN:
                            ci_lo, ci_hi = boot_delta(flo[ch], fhi[ch], l_lo, l_hi, rng)
                        rows.append({
                            "promoter": prom, "duration": dur,
                            "pair": f"{dlo}-{dhi}", "strat": strat_tag, "channel": ch,
                            "auc_pool": ap, "auc_terc": mt, "delta": d,
                            "ci_lo": ci_lo, "ci_hi": ci_hi,
                            "n_lo": nlo - xlo, "n_hi": nhi - xhi,
                            "excl_lo": xlo, "excl_hi": xhi,
                            "resp_lo": float(np.mean(flo["responded"])),
                            "resp_hi": float(np.mean(fhi["responded"])),
                        })
            log(f"[done] {prom} DM {dur}")

    # ---- v0.1.1 guardrail: zero data loaded → abort, do not print an empty ruling table (an empty run is not a ruling) ----
    if not rows:
        log("")
        log("!" * 72)
        log("Abort: no condition data were loaded. Check that ROOT points to the directory")
        log("that directly contains the promoter folders (ALD3 etc.). This run is an empty run, not a ruling.")
        log("!" * 72)
        _log_fh.close()
        return

    # ---------------- ruling (main analysis: 50min, size stratification, main decoders) ----------------
    log("")
    log("=" * 72)
    log(f"Ruling (main-analysis duration {MAIN_DUR}, size tercile stratification)")
    log("-" * 72)
    log(f"{'promoter':<8s} {'Δ_mol':>8s} {'Δ_AU':>8s} {'Δ_time':>8s} "
        f"{'P2-1':>6s} {'P2-2':>6s} {'P2-3':>6s}")
    verdict = {"seed": SEED, "n_boot": N_BOOT, "main_dur": MAIN_DUR,
               "margins": {"P2-1": MARGIN_P21, "P2-2": MARGIN_P22, "P2-3": MARGIN_P23},
               "promoters": {}, "clauses": {}}
    n1 = n2 = n3 = 0
    for prom in WT:
        dm = promoter_delta(rows, prom, "mol_peak")
        da = promoter_delta(rows, prom, "au_peak")
        dt3 = promoter_delta(rows, prom, "time_k3")
        p1 = bool(np.isfinite(dm) and dm >= MARGIN_P21)
        p2 = bool(np.isfinite(da) and np.isfinite(dm)
                  and da <= MARGIN_P22 and da < dm)
        p3 = bool(np.isfinite(dt3) and abs(dt3) <= MARGIN_P23)
        n1 += p1; n2 += p2; n3 += p3
        verdict["promoters"][prom] = {"d_mol": dm, "d_au": da, "d_time_k3": dt3,
                                      "P2-1": p1, "P2-2": p2, "P2-3": p3}
        f = lambda x: f"{x:+.4f}" if np.isfinite(x) else "   NaN"
        log(f"{prom:<8s} {f(dm):>8s} {f(da):>8s} {f(dt3):>8s} "
            f"{'v' if p1 else 'x':>6s} {'v' if p2 else 'x':>6s} {'v' if p3 else 'x':>6s}")
    c1, c2, c3 = n1 >= 6, n2 >= 6, n3 >= 6
    verdict["clauses"] = {"P2-1": {"pass_promoters": n1, "hit": bool(c1)},
                          "P2-2": {"pass_promoters": n2, "hit": bool(c2)},
                          "P2-3": {"pass_promoters": n3, "hit": bool(c3)}}
    overall = bool(c1 and c2 and c3)
    verdict["overall"] = overall
    log("-" * 72)
    log(f"P2-1 (amplitude paid)   :{n1}/7 promoters pass → {'hit' if c1 else '证伪'}")
    log(f"P2-2 (normalized redeem):{n2}/7 promoters pass → {'hit' if c2 else '证伪'}")
    log(f"P2-3 (event timing free):{n3}/7 promoters pass → {'hit' if c3 else '证伪'}")
    log(f"P2-F overall: {'prediction holds (all three clauses hit)' if overall else 'corresponding clause(s) 证伪, double-recorded in the ledger as-is'}")
    log("")

    # ---------------- details of the 21 adjudication units (main channels) ----------------
    log("Details of the 21 adjudication units (main-analysis duration, size stratification):")
    for ch in CHANNELS_MAIN:
        log(f"  [{ch}]")
        for r in rows:
            if (r["duration"] == MAIN_DUR and r["strat"] == "size"
                    and r["channel"] == ch and r["promoter"] in WT):
                ci = (f"[{r['ci_lo']:+.3f},{r['ci_hi']:+.3f}]"
                      if np.isfinite(r["ci_lo"]) else "[ point est ]")
                log(f"    {r['promoter']:<8s} {r['pair']:<9s} "
                    f"AUC_pool={r['auc_pool']:.3f} terc_mean={r['auc_terc']:.3f} "
                    f"Δ={r['delta']:+.4f} {ci} n={r['n_lo']}/{r['n_hi']}")
    log("")

    # ---------------- robust arm and alternative proxies (point-estimate overview) ----------------
    log("robust arm k∈{2,4} and other durations (promoter-level Δ, point estimates):")
    for ch in ["time_k2", "time_k4"]:
        ds = {p: promoter_delta(rows, p, ch) for p in WT}
        log(f"  [{ch}] " + " ".join(
            f"{p}:{ds[p]:+.3f}" if np.isfinite(ds[p]) else f"{p}:NaN" for p in WT))
    for dur in DURS:
        if dur == MAIN_DUR:
            continue
        dm = {p: promoter_delta(rows, p, "mol_peak", dur=dur) for p in WT}
        log(f"  [duration {dur} Δ_mol] " + " ".join(
            f"{p}:{dm[p]:+.3f}" if np.isfinite(dm[p]) else f"{p}:NaN" for p in WT))
    da_alt = {p: promoter_delta(rows, p, "mol_peak", strat="baseMol") for p in WT}
    log("  [alternative-proxy stratified Δ_mol] " + " ".join(
        f"{p}:{da_alt[p]:+.3f}" if np.isfinite(da_alt[p]) else f"{p}:NaN" for p in WT))
    log("")

    # ---------------- mutants (descriptive) ----------------
    log("Mutants (descriptive, not adjudicated):")
    for prom in MUT:
        dm = promoter_delta(rows, prom, "mol_peak")
        da = promoter_delta(rows, prom, "au_peak")
        dt3 = promoter_delta(rows, prom, "time_k3")
        f = lambda x: f"{x:+.4f}" if np.isfinite(x) else "NaN"
        log(f"  {prom}: Δ_mol={f(dm)} Δ_AU={f(da)} Δ_time={f(dt3)}")
    log("")

    # ---------------- secondary arm: FM/FM4 count decoding + CFP probe (descriptive) ----------------
    log("Secondary arm (descriptive, not adjudicated): FM/FM4 decoding")
    for prom in WT:
        fm_descriptive(prom)
    log("Secondary arm (descriptive): CFP gain probe")
    for prom in WT:
        cfp_probe(prom)

    # ---------------- save outputs ----------------
    import csv
    csv_path = os.path.join(outdir, "代码44_单元table.csv")
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    json_path = os.path.join(outdir, "代码44_ruling.json")
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(sanitize(verdict), fh, ensure_ascii=False, indent=2)
    log("")
    log(f"unit table → {csv_path}")
    log(f"ruling     → {json_path}")
    log(f"log        → {LOG_PATH}")
    log(f"total elapsed {_time.time() - t_start:.0f} s")
    _log_fh.close()


if __name__ == "__main__":
    main()
