# -*- coding: utf-8 -*-
"""
================================================================================
code 45 - P3 blind ruling: free-statistic sufficiency (Sung 2022, NF-kB spatial gradient)
version v0.1.1 | 2026-08-14 | run on the user machine
v0.1.1 (ledger entry 78): first-run guardrail interception - R measured as (5, n_chambers) (MATLAB 24x5, only the first 2 rows enter the analysis),
  get_R relaxed from "exactly 2 rows" to ">=2 rows", value-taking logic unchanged. An empty run constitutes no ruling; adjudication logic untouched.
================================================================================
pre-registration: 预注册_P3_免费统计量充分性_NFkB空间梯度_v01.md (ledger entry 75, criteria already pinned)
theory: theorem T' survivor class - free statistics (event times/ordinals/counts/durations) need no cross-cell calibration to suffice for decoding.
blind state: no adjudication statistic was computed before this script was written. Field semantics come entirely from the authors' public MATLAB script comments;
      （DeepAnalysisForFig4.m / HeatMap&Traces_Dose.m / HeatMap&Traces_Duration.m），
  trace values were never touched. The first formal run of this script is the ruling run, recorded in the ledger either way.

implementation decisions (all pinned before opening the data; basis = authors' script comments and the pre-registration, not results):
D1  TNF10 time alignment uses osD=[45,152]+osI_10. Rationale: the 10 ng section of HeatMap&Traces_Dose.m
    first declares [45,152]+osI_10; the later [39,148]+osI_30 override block, compared against DeepAnalysisForFig4.m
    (lines 150-165 use [45,152]+osI_10 with no override block), is judged a copy-paste leftover;
    and the script comments' wall-clock times independently corroborate (imaging 19:04:27, first feed 19:49:11 -> ~45 min;
    second feed 21:36:33 -> ~152 min). The alternative alignment Option B ([39,148]+osI_30) is reported only as a robust-arm point estimate.
D2  chamber inclusion mirrors the authors' scripts: 100/30 ng -> aa=1..12; 10 ng -> aa=[1,2,3,4,9,10,11,12];
    duration series -> aa 2-4=15 min, 5-8=30 min, 9-12=60 min (chamber aa=1 unused by the authors).
D3  analysis window (0,120] min: tau / upward-crossing count / above-threshold duration all defined within this window; never crossing = right-censored, recorded as +inf.
D4  invalid values = non-finite or <0; missing rate >10% (over all recorded frames of the cell) removed;
    baseline (x<=0) valid frames <2 removed; valid trace end x<120 removed; sigma_0=0 removed (std with ddof=1).
D5  no smoothing: threshold crossing is judged directly on the raw ratio trace (the authors' lowess smoothing is display-only and does not enter this decoder).
D6  P3-2 tertile boundaries = 33.33/66.67 percentiles of distances over all included TNF100 cells (fixed after merging;
    labels travel with the cells during bootstrap, boundaries not re-estimated).
D7  AUC = Mann-Whitney (average ranks; +inf participates in ranking). Directions pinned:
    P3-1: P(tau_low-dose > tau_high-dose); P3-2: P(tau_farther-band > tau_nearer-band);
    P3-3: P(above-threshold duration_long-stimulus > above-threshold duration_short-stimulus).
D8  clause verdict table (pinned, two-unit clauses):
    - any unit n<30 -> that unit "insufficient data", clause = insufficient data (neither hit nor falsified, recorded as-is);
    - both units AUC >= margin -> hit;
    - any unit AUC < 0.5 -> falsified (direction reversal);
    - both units < margin with no reversal -> falsified (all units below the line, P3-F literal);
    - otherwise (partially at line) -> miss-intermediate (per the pre-registration literal not a falsification, recorded as-is).
    margins (pre-registration frozen): P3-1=0.60, P3-2=0.55, P3-3=0.60.

guardrails (code-44 v0.1.1 precedent): missing files / failed structure assertions / any condition group with 0 included cells
-> abort the run, no ruling table produced; an empty run constitutes no ruling.

seed 20260814; bootstrap 2000 rounds (per-cell resampling with replacement within each group); 95% CI from the first round
(lesson from the code-44 CI triggering bug).
run environment: pure CPU rank statistics (thousands of cells x 2000 permutations, minute-scale), no matrix-parallel load,
          no GPU needed; depends on numpy/scipy/h5py (the .mat files are v7.3, pip install h5py first).

input (located recursively under ROOT by filename, tolerant of intermediate directory levels):
    FixedSource_TNF100_DC3.mat / FixedSource_TNF30_DC3.mat / FixedSource_TNF10_DC3.mat
    DifferentDuration_TNF100_DC2.mat
output (written to the script's own directory):
    代码45_判定日志.txt / 代码45_单元table.csv / 代码45_robust臂.csv /
    代码45_描述臂.csv / 代码45_ruling.json
================================================================================
"""
import os, sys, json, csv, time, datetime
import numpy as np
from scipy.stats import rankdata

try:
    import h5py
except ImportError:
    print("[ABORT] h5py not installed. Sung's .mat files are MATLAB v7.3; scipy cannot read them.")
    print("       run in the Spyder console:  pip install h5py   then rerun.")
    raise SystemExit(1)

# ----------------------------- config section --------------------------------------
ROOT = r"C:/Users/lihua/Desktop/NC/data/Sung2022"   # <- put the four .mat files here (subdirectories allowed)

SEED   = 20260814
N_BOOT = 2000
N_MIN  = 30                 # minimum cells per group per unit (guardrail D8)
K_MAIN = 3                  # main-analysis threshold fold number
K_ROB  = (2, 4)             # robust arms
WINDOW = (0.0, 120.0)       # analysis window (0,120] min
EXTENT_MIN = 120.0          # traces must extend at least to 120 min
MISS_MAX = 0.10             # missing-rate cap
MARGIN = {"P3-1": 0.60, "P3-2": 0.55, "P3-3": 0.60}
OVERLAP_BAND = (0.5, 1.5)   # P3-1 distance overlap-band robust arm (mm)

DOSE_FILES = {"100": "FixedSource_TNF100_DC3.mat",
              "30":  "FixedSource_TNF30_DC3.mat",
              "10":  "FixedSource_TNF10_DC3.mat"}
DUR_FILE = "DifferentDuration_TNF100_DC2.mat"

# time-alignment constants transcribed verbatim from the authors' scripts (DeepAnalysisForFig4.m / HeatMap&Traces_Dose.m)
OSD = {"100": [38.0, 148.0], "30": [39.0, 148.0], "10": [45.0, 152.0]}
OSI = {
 "100": [0,0,0,0,.7,.7,.7,.7,1.4,1.4,1.4,1.4,2.1,2.1,2.1,2.1,2.8,2.8,2.8,2.8,3.5,3.5,3.5,3.5],
 "30":  [0,0,0,0,.7,.7,.7,.7,1.4,1.4,1.4,1.4,3.5,3.5,3.5,3.5,2.8,2.8,2.8,2.8,2.1,2.1,2.1,2.1],
 "10":  [0,0,0,0,.7,.7,.7,.7,1.4,1.4,1.4,1.4,2.1,2.1,2.1,2.1,2.8,2.8,3.15,3.15,3.5,3.5,3.5,3.5],
}
OSD_10_ALT = [39.0, 148.0]          # Option B robust arm (the override block in the authors' script)
OSI_10_ALT = OSI["30"]

CHAMBERS = {"100": list(range(1, 13)), "30": list(range(1, 13)),
            "10": [1, 2, 3, 4, 9, 10, 11, 12]}
DUR_GROUPS = {"15min": [2, 3, 4], "30min": [5, 6, 7, 8], "60min": [9, 10, 11, 12]}
ADJ_DUR = -42.0                      # duration-series global offset (authors' script adj=-42)

try:
    OUT_DIR = os.path.dirname(os.path.abspath(__file__))
except NameError:
    OUT_DIR = os.getcwd()

_LOG_LINES = []
def log(msg=""):
    line = str(msg)
    print(line)
    _LOG_LINES.append(line)

def abort(msg):
    log("")
    log("[ABORT] " + msg)
    log("This run is an empty run; it constitutes no ruling and produces no ruling table. Please send this log back.")
    try:
        with open(os.path.join(OUT_DIR, "代码45_判定日志.txt"), "w", encoding="utf-8") as fh:
            fh.write("\n".join(_LOG_LINES))
    except Exception:
        pass
    raise SystemExit(2)

# ------------------------- v7.3 .mat load（h5py） ---------------------------
def find_file(name):
    hits = []
    for dp, _dn, fn in os.walk(ROOT):
        if name in fn:
            hits.append(os.path.join(dp, name))
    return hits

def read_label(f, ref):
    arr = np.asarray(f[ref][()]).flatten()
    return "".join(chr(int(c)) for c in arr if int(c) != 0)

def get_R(f, tag):
    if "R" not in f:
        abort(f"{tag}: no variable R in the file (top-level keys={list(f.keys())})")
    R = f["R"]
    if not isinstance(R, h5py.Dataset) or R.ndim != 2 or R.shape[0] < 2:
        abort(f"{tag}: R shape exception {getattr(R, 'shape', None)} (expected (>=2, n_chambers))")
    # v0.1.1: measured (5, n_chambers) - MATLAB R is (n_chambers,5); the authors' script uses only rows 1 and 2
    log(f"  [{tag}] R shape {R.shape} (row 0 = chamber labels, row 1 = per-cell blocks; other rows unused by the authors' script)")
    return R

def get_interval(f, tag):
    if "interval" in f:
        try:
            val = float(np.asarray(f["interval"][()]).flatten()[0])
            log(f"  [{tag}] interval = {val} min (read from the file)")
            return val
        except Exception:
            pass
    log(f"  [{tag}] no interval in the file; taking 6.0 min per the authors' script")
    return 6.0

def chamber_entry(f, R, aa):
    """Return (label, cells_ref); empty chambers return (None, None). aa is the 1-based chamber number."""
    lref = R[0, aa - 1]
    cref = R[1, aa - 1]
    if not lref:
        return None, None
    label = read_label(f, lref)
    if not cref:
        return label, None
    return label, cref

def iter_cells(f, cells_ref, tag, aa):
    """Per-cell yields (frames, values, dist_mm). Structure assertion failure -> abort (guardrail)."""
    S = f[cells_ref]
    if not isinstance(S, h5py.Dataset) or S.ndim != 2 or S.shape[0] < 3:
        abort(f"{tag} chamber {aa}: per-cell block shape exception {getattr(S,'shape',None)}"
              f" (expected a reference matrix of (>=3, n_cells))")
    ncell = S.shape[1]
    out = []
    n_skip = 0
    for j in range(ncell):
        tref = S[1, j]
        if not tref:
            n_skip += 1
            continue
        T = np.asarray(f[tref][()], dtype=float)
        if T.ndim != 2 or T.shape[0] < 2:
            abort(f"{tag} chamber {aa} cell {j+1}: trace block shape exception {T.shape}"
                  f" (expected (>=2, n_frames)): row 1 = frame number, row 2 = NF-kB")
        frames = T[0, :]
        values = T[1, :]
        fin_f = frames[np.isfinite(frames)]
        if fin_f.size == 0 or fin_f.min() < 1:
            abort(f"{tag} chamber {aa} cell {j+1}: frame-number exception (minimum value "
                  f"{fin_f.min() if fin_f.size else 'NaN'}, expected >=1) - structure inconsistent with decoding")
        dist = np.nan
        mref = S[2, j]
        if mref:
            M = np.asarray(f[mref][()], dtype=float)
            if M.ndim == 2 and M.shape[0] >= 4:
                drow = M[3, :]
                drow = drow[np.isfinite(drow)]
                if drow.size:
                    dist = float(drow.mean()) / 1000.0   # um -> mm (same as the authors' script)
        out.append((frames, values, dist))
    if n_skip:
        log(f"  [{tag}] chamber {aa}: {n_skip} empty cell slots skipped")
    return out

# ------------------------- include/remove + features --------------------------------
def featurize(frames, values, dist, interval, offset):
    """Return (status, feature dict). Status: 'ok' or a removal reason. All rules in header D3/D4/D5.
    x = (frame-1)*interval + offset (offset: dose series = chamber feed offset os; duration series = adj)."""
    x = (frames - 1.0) * interval + offset
    v = values.astype(float)
    valid = np.isfinite(v) & (v >= 0) & np.isfinite(x)   # missing frame numbers count toward the missing rate as well
    if v.size == 0:
        return "empty轨迹", None
    miss = 1.0 - float(valid.mean())
    if miss > MISS_MAX:
        return f"missing>{int(MISS_MAX*100)}%", None
    order = np.argsort(x, kind="stable")
    x, v, valid = x[order], v[order], valid[order]
    xv, vv = x[valid], v[valid]
    base = vv[xv <= 0.0]
    if base.size < 2:
        return "baseline帧<2", None
    if xv.size == 0 or xv[-1] < EXTENT_MIN:
        return "末端<120min", None
    mu0 = float(base.mean())
    sd0 = float(base.std(ddof=1))
    if sd0 == 0.0:
        return "sigma0=0", None
    win = (xv > WINDOW[0]) & (xv <= WINDOW[1])
    xw, vw = xv[win], vv[win]
    feat = {"dist": dist, "mu0": mu0, "sd0": sd0}
    for k in (K_MAIN,) + tuple(K_ROB):
        thr = mu0 + k * sd0
        hit = vw >= thr
        feat[f"tau{k}"] = float(xw[hit][0]) if hit.any() else np.inf
    thr3 = mu0 + K_MAIN * sd0
    hit3 = vw >= thr3
    feat["dur3"] = float(hit3.sum()) * interval          # total above-threshold duration (min)
    if vw.size >= 2:
        feat["count3"] = int(np.sum((vw[:-1] < thr3) & (vw[1:] >= thr3)))
    else:
        feat["count3"] = 0
    feat["peak"] = float(vw.max()) if vw.size else np.nan
    feat["responded"] = bool(np.isfinite(feat[f"tau{K_MAIN}"]))
    return "ok", feat

# ------------------------- AUC and bootstrap --------------------------------
def auc_mw(a, b):
    """P(a>b)+0.5*P(a=b), average ranks; NaN dropped, +inf participates in ranking (same convention as code 44)."""
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    n1, n2 = a.size, b.size
    if n1 == 0 or n2 == 0:
        return np.nan
    r = rankdata(np.concatenate([a, b]))
    u1 = r[:n1].sum() - n1 * (n1 + 1) / 2.0
    return float(u1 / (n1 * n2))

def boot_ci(a, b, rng, n_boot=N_BOOT):
    a = np.asarray(a, float); b = np.asarray(b, float)
    a = a[~np.isnan(a)]; b = b[~np.isnan(b)]
    n1, n2 = a.size, b.size
    stats = np.empty(n_boot)
    for i in range(n_boot):
        ia = rng.integers(0, n1, n1)
        ib = rng.integers(0, n2, n2)
        stats[i] = auc_mw(a[ia], b[ib])
    lo, hi = np.percentile(stats, [2.5, 97.5])
    return float(lo), float(hi)

# ------------------------- clause verdict table (D8, pinned) ---------------------------
def clause_verdict(clause, units):
    """units: list of dict(n1, n2, auc, margin). Returns (verdict, note)."""
    margin = MARGIN[clause]
    if any(u["n1"] < N_MIN or u["n2"] < N_MIN for u in units):
        return "数据不足", "至少一单元 n<%d，不构成hit也不构成证伪" % N_MIN
    aucs = [u["auc"] for u in units]
    if all(a >= margin for a in aucs):
        return "hit", f"全部单元 ≥ {margin}"
    if any(a < 0.5 for a in aucs):
        return "证伪", "existsdirection反转单元（AUC<0.5）"
    if all(a < margin for a in aucs):
        return "证伪", f"全单元不达线（均 < {margin}，P3-F 字面）"
    return "miss-中间态", "部分达线但未全达线；按pre-registration字面不记证伪，如实入账"

# ------------------------- assemble condition groups --------------------------------------
def collect_group(path, tag, chambers, offset_fn, log_head=True):
    """offset_fn(label, aa) -> float. Returns a list of feat dicts (included cells only)."""
    f = h5py.File(path, "r")
    interval = get_interval(f, tag)
    R = get_R(f, tag)
    feats, reasons = [], {}
    n_ch_used, n_raw = 0, 0
    for aa in chambers:
        if aa > R.shape[1]:
            log(f"  [{tag}] chamber {aa} beyond the R range ({R.shape[1]} chambers total), skipped")
            continue
        label, cref = chamber_entry(f, R, aa)
        if label is None or cref is None:
            continue
        n_ch_used += 1
        cells = iter_cells(f, cref, tag, aa)
        n_raw += len(cells)
        off = offset_fn(label, aa)
        for frames, values, dist in cells:
            st, ft = featurize(frames, values, dist, interval, off)
            if st == "ok":
                ft["chamber"] = aa
                feats.append(ft)
            else:
                reasons[st] = reasons.get(st, 0) + 1
    f.close()
    log(f"  [{tag}] chambers {n_ch_used}, raw cells {n_raw}, included {len(feats)}; "
        f"removeddetails {reasons if reasons else '{}'}")
    return feats

def dose_offset_fn(osd, osi):
    def fn(label, aa):
        sec = 0 if "1-" in label else 1     # authors' script: contains(label,'1-') -> osD(1)
        return -osd[sec] + osi[aa - 1]
    return fn

def dur_offset_fn(label, aa):
    return ADJ_DUR

# ------------------------- main pipeline -------------------------------------------
def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    log("=" * 74)
    log("code 45 - P3 blind ruling: free-statistic sufficiency (Sung 2022 NF-kB spatial gradient) v0.1.1")
    log("run time: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    log(f"seed={SEED}  bootstrap={N_BOOT}  adjudication margins={MARGIN}  analysis window={WINDOW} min")
    log(f"numpy={np.__version__}  scipy loaded  h5py={h5py.__version__}")
    log("pre-registration: 预注册_P3_免费统计量充分性_NFkB空间梯度_v01.md (entry 75)")
    log("=" * 74)

    # ---- guardrail 1: all files present ----
    paths = {}
    missing = []
    for key, fname in list(DOSE_FILES.items()) + [("DUR", DUR_FILE)]:
        hits = find_file(fname)
        if not hits:
            missing.append(fname)
        else:
            paths[key] = sorted(hits)[0]
    if missing:
        log("[missing files] the following files were not found under ROOT:")
        for m in missing:
            log("   - " + m)
        log(f"ROOT = {ROOT}")
        log("Please confirm the four .mat files (note: each is nested one extra zip layer inside the package, unzip once more) are in place.")
        abort("data files incomplete")
    for k, p in paths.items():
        log(f"  file[{k}] = {p}")

    # ---- dose series ----
    log("\n[load] dose series (FixedSource, per-chamber feed offset os=-osD+osI)")
    dose = {}
    for d in ("10", "30", "100"):
        dose[d] = collect_group(paths[d], f"TNF{d}",
                                CHAMBERS[d], dose_offset_fn(OSD[d], OSI[d]))
    # ---- duration series ----
    log("\n[load] duration series (DifferentDuration, adj=-42)")
    dur = {}
    for gname, chs in DUR_GROUPS.items():
        dur[gname] = collect_group(paths["DUR"], f"duration{gname}", chs, dur_offset_fn)

    # ---- guardrail 2: empty-run protection ----
    empties = ([f"TNF{d}" for d in ("10", "30", "100") if len(dose[d]) == 0]
               + [f"duration{g}" for g in DUR_GROUPS if len(dur[g]) == 0])
    if empties:
        log("[empty groups] the following condition groups have 0 included cells: " + ", ".join(empties))
        abort("empty condition groups exist (path/structure needs checking)")

    def col(group, key):
        return np.array([c[key] for c in group], float)

    # ---- P3-2 tertiles (D6) ----
    d100 = [c for c in dose["100"] if np.isfinite(c["dist"]) and c["dist"] > 0]
    log(f"\n[P3-2] TNF100 cells with valid distance {len(d100)}/{len(dose['100'])}")
    if len(d100) == 0:
        abort("TNF100 has no valid distances; P3-2 cannot be constructed")
    q1, q2 = np.percentile([c["dist"] for c in d100], [100/3, 200/3])
    band = {"near": [c for c in d100 if c["dist"] <= q1],
            "mid":  [c for c in d100 if q1 < c["dist"] <= q2],
            "far":  [c for c in d100 if c["dist"] > q2]}
    log(f"  tertile boundaries: {q1:.4f} / {q2:.4f} mm; "
        f"near/mid/far = {len(band['near'])}/{len(band['mid'])}/{len(band['far'])}")
    if any(len(band[b]) == 0 for b in band):
        abort("empty distance band exists")

    # ---- the 6 adjudication units ----
    units = [
        {"clause": "P3-1", "unit": "P3-1a", "contrast": "τ(3): TNF10 > TNF30",
         "a": col(dose["10"], f"tau{K_MAIN}"),  "b": col(dose["30"], f"tau{K_MAIN}")},
        {"clause": "P3-1", "unit": "P3-1b", "contrast": "τ(3): TNF30 > TNF100",
         "a": col(dose["30"], f"tau{K_MAIN}"),  "b": col(dose["100"], f"tau{K_MAIN}")},
        {"clause": "P3-2", "unit": "P3-2a", "contrast": "τ(3): 中带 > 近带",
         "a": col(band["mid"], f"tau{K_MAIN}"), "b": col(band["near"], f"tau{K_MAIN}")},
        {"clause": "P3-2", "unit": "P3-2b", "contrast": "τ(3): 远带 > 中带",
         "a": col(band["far"], f"tau{K_MAIN}"), "b": col(band["mid"], f"tau{K_MAIN}")},
        {"clause": "P3-3", "unit": "P3-3a", "contrast": "越threshold时程: 30min > 15min",
         "a": col(dur["30min"], "dur3"),        "b": col(dur["15min"], "dur3")},
        {"clause": "P3-3", "unit": "P3-3b", "contrast": "越threshold时程: 60min > 30min",
         "a": col(dur["60min"], "dur3"),        "b": col(dur["30min"], "dur3")},
    ]

    log("\n[adjudication] 6 units (AUC directions in header D7; CI = bootstrap 2000 percentile)")
    log(f"{'unit':<7}{'contrast':<26}{'n_A':>6}{'n_B':>6}{'AUC':>8}{'CI_lo':>8}{'CI_hi':>8}{'margin':>6}  at-line")
    rows = []
    for u in units:
        u["n1"], u["n2"] = int(u["a"].size), int(u["b"].size)
        if u["n1"] == 0 or u["n2"] == 0:
            abort(f"unit {u['unit']} empty group")
        u["auc"] = auc_mw(u["a"], u["b"])
        if u["n1"] >= N_MIN and u["n2"] >= N_MIN:
            u["ci"] = boot_ci(u["a"], u["b"], rng)
        else:
            u["ci"] = (np.nan, np.nan)
        ok = (u["auc"] >= MARGIN[u["clause"]]) if u["n1"] >= N_MIN and u["n2"] >= N_MIN else None
        u["pass"] = ok
        log(f"{u['unit']:<7}{u['contrast']:<26}{u['n1']:>6}{u['n2']:>6}"
            f"{u['auc']:>8.4f}{u['ci'][0]:>8.4f}{u['ci'][1]:>8.4f}"
            f"{MARGIN[u['clause']]:>6.2f}  "
            + ("v" if ok else ("x" if ok is False else "insufficient-data")))
        rows.append({"clause": u["clause"], "unit": u["unit"], "contrast": u["contrast"],
                     "n_a": u["n1"], "n_b": u["n2"], "auc": u["auc"],
                     "ci_lo": u["ci"][0], "ci_hi": u["ci"][1],
                     "margin": MARGIN[u["clause"]],
                     "pass": {True: "Y", False: "N", None: "INSUFFICIENT"}[ok]})

    # ---- clause verdicts (D8) ----
    log("\n[clause verdicts]")
    verdicts = {}
    for clause in ("P3-1", "P3-2", "P3-3"):
        us = [u for u in units if u["clause"] == clause]
        v, why = clause_verdict(clause, us)
        verdicts[clause] = {"verdict": v, "reason": why,
                            "units": [u["unit"] for u in us]}
        log(f"  {clause}: {v} —— {why}")
    any_fals = any(v["verdict"] == "证伪" for v in verdicts.values())
    all_hit = all(v["verdict"] == "hit" for v in verdicts.values())
    log("\n[P3-F overall falsification] " + ("triggered: at least one clause falsified" if any_fals else
                              ("not triggered; all three clauses hit" if all_hit else "not triggered; miss/insufficient-data clauses exist")))
    overall = ("P3 hit（三条款全达线）" if all_hit else
               "P3 证伪（exists证伪条款）" if any_fals else
               "P3 miss-中间态/数据不足（按判定table如实入账）")
    log("[overall ruling] " + overall)

    # ---- robust arms (point estimates, not adjudicated) ----
    log("\n[robust arms] point estimates (not adjudicated)")
    rob_rows = []
    def rob(tag, a, b):
        v = auc_mw(a, b)
        log(f"  {tag:<44} n={len(a)}/{len(b)}  AUC={v:.4f}")
        rob_rows.append({"arm": tag, "n_a": len(a), "n_b": len(b), "auc": v})
        return v
    for k in K_ROB:
        rob(f"k={k} P3-1a τ: TNF10>TNF30", col(dose["10"], f"tau{k}"), col(dose["30"], f"tau{k}"))
        rob(f"k={k} P3-1b τ: TNF30>TNF100", col(dose["30"], f"tau{k}"), col(dose["100"], f"tau{k}"))
        rob(f"k={k} P3-2a tau: mid>near", col(band["mid"], f"tau{k}"), col(band["near"], f"tau{k}"))
        rob(f"k={k} P3-2b tau: far>mid", col(band["far"], f"tau{k}"), col(band["mid"], f"tau{k}"))
    lo_b, hi_b = OVERLAP_BAND
    ov = {d: [c for c in dose[d]
              if np.isfinite(c["dist"]) and lo_b <= c["dist"] <= hi_b] for d in ("10", "30", "100")}
    log(f"  cells within the distance overlap band {OVERLAP_BAND} mm: "
        + ", ".join(f"TNF{d}={len(ov[d])}" for d in ("10", "30", "100")))
    if all(len(ov[d]) >= N_MIN for d in ("10", "30", "100")):
        rob(f"overlap-band P3-1a tau(3): TNF10>TNF30", col(ov["10"], f"tau{K_MAIN}"), col(ov["30"], f"tau{K_MAIN}"))
        rob(f"overlap-band P3-1b tau(3): TNF30>TNF100", col(ov["30"], f"tau{K_MAIN}"), col(ov["100"], f"tau{K_MAIN}"))
    else:
        log("  some group in the overlap band has n<30; overlap-band arm skipped (recorded as-is)")
        rob_rows.append({"arm": "重叠带 P3-1", "n_a": "", "n_b": "", "auc": "SKIPPED_n<30"})
    # Option B alignment robust arm (P3-1 only, point estimate)
    log("  [Option B] TNF10 re-aligned with the overridden block [39,148]+osI_30 (robust arm only)")
    dose10B = collect_group(paths["10"], "TNF10-OptB",
                            CHAMBERS["10"], dose_offset_fn(OSD_10_ALT, OSI_10_ALT))
    if len(dose10B) >= N_MIN:
        rob(f"OptB P3-1a τ(3): TNF10>TNF30", col(dose10B, f"tau{K_MAIN}"), col(dose["30"], f"tau{K_MAIN}"))
        rob(f"OptB P3-1b τ(3): TNF30>TNF100（control）", col(dose["30"], f"tau{K_MAIN}"), col(dose["100"], f"tau{K_MAIN}"))
    else:
        log("  Option B group n<30, skipped (recorded as-is)")
        rob_rows.append({"arm": "OptB P3-1a", "n_a": "", "n_b": "", "auc": "SKIPPED_n<30"})

    # ---- descriptive arm (not adjudicated) ----
    log("\n[descriptive arm] paid-type control (peak amplitude) + count decoder + response rate (none adjudicated)")
    desc_rows = []
    def desc(tag, a, b):
        v = auc_mw(a, b)
        log(f"  {tag:<44} n={len(a)}/{len(b)}  AUC={v:.4f}")
        desc_rows.append({"arm": tag, "n_a": len(a), "n_b": len(b), "auc": v})
    desc("peak P3-1a: TNF30>TNF10", col(dose["30"], "peak"), col(dose["10"], "peak"))
    desc("peak P3-1b: TNF100>TNF30", col(dose["100"], "peak"), col(dose["30"], "peak"))
    desc("peak P3-2a: near>mid", col(band["near"], "peak"), col(band["mid"], "peak"))
    desc("peak P3-2b: mid>far", col(band["mid"], "peak"), col(band["far"], "peak"))
    desc("peak P3-3a: 30min>15min", col(dur["30min"], "peak"), col(dur["15min"], "peak"))
    desc("peak P3-3b: 60min>30min", col(dur["60min"], "peak"), col(dur["30min"], "peak"))
    desc("count P3-1a: TNF30>TNF10", col(dose["30"], "count3"), col(dose["10"], "count3"))
    desc("count P3-1b: TNF100>TNF30", col(dose["100"], "count3"), col(dose["30"], "count3"))
    log("  response rate (uncensored fraction of tau(3)):")
    for name, g in [("TNF10", dose["10"]), ("TNF30", dose["30"]), ("TNF100", dose["100"]),
                    ("近带", band["near"]), ("中带", band["mid"]), ("远带", band["far"]),
                    ("15min", dur["15min"]), ("30min", dur["30min"]), ("60min", dur["60min"])]:
        rr = np.mean([c["responded"] for c in g]) if g else np.nan
        log(f"    {name:<8} n={len(g):<5} response rate={rr:.3f}")
        desc_rows.append({"arm": f"响应率 {name}", "n_a": len(g), "n_b": "", "auc": rr})

    # ---- write output ----
    def write_csv(path, rows_, fields):
        with open(path, "w", newline="", encoding="utf-8-sig") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            for r in rows_:
                w.writerow(r)

    p_unit = os.path.join(OUT_DIR, "代码45_单元table.csv")
    p_rob  = os.path.join(OUT_DIR, "代码45_robust臂.csv")
    p_desc = os.path.join(OUT_DIR, "代码45_描述臂.csv")
    p_json = os.path.join(OUT_DIR, "代码45_ruling.json")
    p_log  = os.path.join(OUT_DIR, "代码45_判定日志.txt")
    write_csv(p_unit, rows, ["clause", "unit", "contrast", "n_a", "n_b", "auc",
                             "ci_lo", "ci_hi", "margin", "pass"])
    write_csv(p_rob, rob_rows, ["arm", "n_a", "n_b", "auc"])
    write_csv(p_desc, desc_rows, ["arm", "n_a", "n_b", "auc"])
    with open(p_json, "w", encoding="utf-8") as fh:
        json.dump({"version": "v0.1.1", "seed": SEED, "n_boot": N_BOOT,
                   "margins": MARGIN, "units": rows, "verdicts": verdicts,
                   "overall": overall,
                   "any_clause_falsified": any_fals, "all_clauses_hit": all_hit},
                  fh, ensure_ascii=False, indent=2, default=str)
    log(f"\nelapsed {time.time()-t0:.1f}s; output:")
    for p in (p_unit, p_rob, p_desc, p_json, p_log):
        log("  " + p)
    with open(p_log, "w", encoding="utf-8") as fh:
        fh.write("\n".join(_LOG_LINES))

if __name__ == "__main__":
    main()
