# -*- coding: utf-8 -*-
"""
code 46 v0.1.5 - P4 blind ruling: ERK re-adjudication (GPCR-KTR raw-trace level)
=================================================================
v0.1.1 (2026-08-14): first empty-run fix - actual files have no Experiment column, cell key changed to
  Unique_Object (falls back to a Date+Slide+Position+Object composite key when missing); X,Y are per-frame
  plasma-membrane centroid coordinates, excluded from the key design; the batch-purity arm uses Date.
v0.1.2 (2026-08-14, never delivered, i.e. void, registered as-is): after the second empty run it was misjudged as
  "six wells with 0.5 min phase staggering, field-of-view-level anchor" - schema validation found it would anchor t0 wrongly to
  5.5 min，void。
v0.1.3 (2026-08-14, anchor fixed): schema validation confirmed the authors had linearly interpolated the raw acquisition (3.5 min/frame,
  six well positions staggered by 0.5 min delay) onto a **globally uniform 0.5 min grid** - all
  7,393 cells have full 162 frames, 2.5-83.0 min; t = Time_in_min - 23.5 (time point 7 =
  2.5+6x3.5 = 23.5 = the sample-addition point t=0, matching the authors' AUC window tp9-18 <-> 7-38.5 min post-stimulus
  verbatim); the tp7 (t=0) frame is discarded; baseline t<0 expected 42 frames. Criteria/window definitions unchanged.
v0.1.4 (2026-08-14): third empty-run fix - the actual dose level count is uniformly 6 (validated over all ligands x all
  inhibitor conditions; the "7 levels" frozen in v0.1.0 was a wrong hypothesis). Pair re-anchoring:
  low zone = pairs (1,2)(2,3), high zone = pairs (4,5)(5,6); the middle pair (3,4) = half-saturation transition band, excluded from adjudication
  and moved to the descriptive arm. Unit counts and count criteria recomputed for 6 levels (P4-1: 4 units; P4-2: 4+4;
  P4-3: 10 units, >=8/10 hit, >=3/10 falsify); adjudication thresholds 0.60/0.05/0.01/0.03
  untouched. Pre-registration bumped to v0.3.
The three empty runs (guardrail hits, none constituting a ruling) are all in the ledger.
v0.1.5 (2026-08-14): the formal ruling was produced (the verdict printed completely before the crash; the ruling holds,
  following the code-42 precedent); then the json-save bug was fixed - np.int64/np.float64 (auc_mw returns
  numpy types propagated via sum/abs) cannot be JSON-serialized directly; a default converter +
  explicit int() was added. Adjudication logic, seed, and rng call order unchanged -> main analysis reproduces bit-for-bit.
=================================================================
Pre-registration: 预注册_P4_ERK重裁_GPCR-KTR轨迹级_v01.md (ledger entry 81)
Dataset  : Chavez-Abiega & Goedhart (2022) J Cell Sci 135(6):jcs259685
          Zenodo doi:10.5281/zenodo.5836623
Run: user machine Spyder  %runfile this file   (CPU-light, no GPU needed)

Frozen decisions (pinned before opening the data, verbatim from the pre-registration section 2):
  D1 time benchmark (anchored in v0.1.3): the authors linearly interpolated the raw acquisition (3.5 min/frame, six well positions
     staggered by 0.5 min delay) onto a globally uniform 0.5 min grid (all cells full 162 frames,
     2.5–83.0 min）；t = Time_in_min − 23.5 min（time point 7 = 2.5+6×3.5 =
     23.5 = sample-addition point t=0; authors' AUC window tp9-18 <-> 7-38.5 min post-stimulus validated verbatim).
     baseline = t<0 frames (expected 42); tp7 (t=0, ambiguous ownership of the addition frame) discarded -
     neither into baseline nor into the response window. Guardrail: global grid start 2.5+/-0.1, step median
     0.5+/-0.05, end >=60 min; any violation -> empty run aborts.
  D2 baseline: all t<0 frames (expected 6); cells with <3 baseline frames removed; sigma_0=0 removed (ddof=1).
  D3 window: tau/count/duration window = (0, 38.5] min post-stim; censoring recorded as +inf, Mann-Whitney average ranks.
  D4 include/remove: frames with non-finite or <0 C/N are removed frame-wise (counted missing); cells with missing >10% removed;
     trace end must reach >=38.5 min post-stim. Main-analysis inhibitor = DMSO.
  D5 nonesmoothed。
  D6 meanArea three-stratum boundaries = 33.33/66.67 percentiles of the adjudication-unit merged cells; labels follow the bootstrap.
  D7 AUC direction frozen: tau-unit AUC = P(tau_low-dose > tau_high-dose) + 0.5*ties (earlier at higher dose is the predicted direction);
     P4-3 stratification penalty delta = AUC_pooled - mean(AUC_within-strata).
  D8 verdict table: hit / falsified / miss-intermediate / insufficient-data (any group n<30 in a unit), four states
     + P4-F overall falsification (any clause landing in the falsified cell -> overall ruling falsified); recorded in the ledger either way.

Clause criteria (frozen; unit counts re-anchored for the actual 6 levels in v0.1.4, thresholds untouched):
  P4-1: his+UK x low-zone 2 pairs = 4 units, tau(3) AUC >= 0.60 all reaching the line -> hit;
        any unit AUC<0.5 (direction reversal) or all units below the line -> falsified; in between -> intermediate.
  P4-2: high-zone 4 units, mean AUC_H < 0.60 and D = low-zone mean - AUC_H >= 0.05 -> hit;
        AUC_H >= 0.60 or D <= 0 -> falsified; in between -> intermediate. Logically independent of P4-1.
  P4-3: 10 units (2 ligands x 5 pairs), meanArea stratification penalty delta;
        >=8/10 units |delta|<=0.01 -> hit; >=3/10 units |delta|>=0.03 -> falsified; in between -> intermediate.

guardrail chain (code-44/45 precedent): cell key not unique / frame-order exception / time grid mismatch /
  dose levels != 7 / required columns missing -> [ABORT] empty run, no ruling constituted.
"""

import os, re, sys, json, datetime
import numpy as np
import pandas as pd

# ---------------- config (frozen) ----------------
ROOT      = r"C:/Users/lihua/Desktop/NC/data/ChavezAbiega2022"
SEED      = 20260815
N_BOOT    = 2000
N_MIN     = 30
K_MAIN    = 3
K_ROB     = (2, 4)

M41       = 0.60          # P4-1 unit criterion
M42_H     = 0.60          # P4-2 high-zone mean line
M42_D     = 0.05          # P4-2 gradient line
P43_LO, P43_HI = 0.01, 0.03
P43_HIT, P43_KILL = 8, 3   # hit/falsify count lines for 10 units (2 ligands x 5 pairs)

T_STIM    = 23.5          # min (time point 7 = sample-addition point t=0)
DT_GLOB   = 0.5           # min/frame (the authors' interpolated global uniform grid)
T_START   = 2.5           # grid start
T_END_MIN = 60.0          # grid end lower bound (expected 83.0)
WIN_LO, WIN_HI = 0.0, 38.5    # free-statistic window (0, 38.5]
PAID_LO, PAID_HI = 7.0, 38.5  # paid integral window [7, 38.5] (authors' AUC window)

LIG_MAIN  = ["Histamine", "UK"]     # main-analysis ligands
LIG_DESC  = "S1P"                    # descriptive-arm ligand
FILE_F1   = os.path.join("Figure_1", "Data", "All_Ex_{lig}_DMSO.csv")
FILE_F3   = os.path.join("Figure_3", "Data", "ALL_{lig3}.csv")   # inhibitor descriptive arm (may be absent)
LIG3_NAME = {"Histamine": "His", "UK": "UK", "S1P": "S1P"}

OUT_DIR   = os.path.dirname(os.path.abspath(__file__))
LOG_LINES = []

def log(s):
    print(s)
    LOG_LINES.append(str(s))

def abort(msg):
    log("[ABORT] " + msg + " - empty run, no ruling constituted.")
    flush_log()
    sys.exit(2)

def flush_log():
    with open(os.path.join(OUT_DIR, "代码46_判定日志.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(LOG_LINES) + "\n")

# ---------------- loading and guardrails ----------------
COLMAP = {
    "time_in_min": "Time", "time": "Time",
    "cn_erk": "ERK", "erk": "ERK",
    "cn_aktrb": "Akt", "cn_akt": "Akt", "akt": "Akt",
    "meanarea": "Area", "condition": "Dose", "experiment": "Exp",
    "x": "X", "y": "Y", "inhibitor": "Inh",
    "unique_object": "UID", "object": "Obj", "original_object": "OObj",
    "date": "Date", "slide": "Slide", "position": "Pos",
}

def parse_dose(v):
    m = re.search(r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?", str(v))
    return float(m.group()) if m else np.nan

def load_long(path, desc=False):
    """Read long table -> normalize columns -> guardrails -> per-cell dict list.
    desc=False: main analysis, guardrail failure = empty-run abort;
    desc=True : descriptive arm, guardrail failure = raise ValueError (caller warns and skips)."""
    def fail(msg):
        if desc:
            raise ValueError(msg)
        abort(msg)
    if not os.path.isfile(path):
        fail(f"file missing: {path}")
    df = pd.read_csv(path)
    ren = {}
    for c in df.columns:
        k = str(c).strip().lower()
        if k in COLMAP and COLMAP[k] not in ren.values():
            ren[c] = COLMAP[k]
    df = df.rename(columns=ren)
    need = ["Time", "ERK", "Area", "Dose", "Date"]
    miss = [c for c in need if c not in df.columns]
    if miss:
        fail(f"{os.path.basename(path)}: required columns missing {miss} (actual columns {list(df.columns)})")
    # cell key: prefer Unique_Object; otherwise the (Date,Slide,Position,Object) composite key.
    # X,Y are per-frame plasma-membrane centroid coordinates (cells move); strictly excluded from the key.
    if "UID" in df.columns and df["UID"].notna().all():
        key_cols = ["UID"]
    elif all(c in df.columns for c in ["Date", "Slide", "Pos", "Obj"]):
        key_cols = ["Date", "Slide", "Pos", "Obj"]
    else:
        fail(f"{os.path.basename(path)}: no usable cell key (need Unique_Object or Date/Slide/Position/Object; actual columns {list(df.columns)})")
    df["DoseV"] = df["Dose"].apply(parse_dose)
    if df["DoseV"].isna().any():
        fail(f"{os.path.basename(path)}: Condition column has unparseable dose values")
    # time-grid guardrail: the authors' interpolated global uniform 0.5 min grid
    tg = np.sort(df["Time"].unique())
    if len(tg) < 2 or abs(tg.min() - T_START) > 0.1 or tg.max() < T_END_MIN or \
       abs(np.median(np.diff(tg)) - DT_GLOB) > 0.05:
        fail(f"{os.path.basename(path)}: global grid exception ({tg.min() if len(tg) else float('nan'):.2f}-"
             f"{tg.max() if len(tg) else float('nan'):.2f} min, median step "
             f"{np.median(np.diff(tg)) if len(tg)>1 else float('nan'):.3f}; expected 2.5->=60 min, 0.5 step)")
    # dose-level guardrail (v0.1.4: actually uniformly 6 levels; schema validated over all ligands x all inhibitors)
    doses = np.sort(df["DoseV"].unique())
    if len(doses) != 6:
        fail(f"{os.path.basename(path)}: dose level count {len(doses)} != 6 (levels {doses.tolist()})")
    # cellkeyuniqueness guardrail
    if df.duplicated(subset=key_cols + ["Time"]).any():
        fail(f"{os.path.basename(path)}: cell key {key_cols}+Time not unique")
    cells = []
    n_drop_dose = 0
    for _k, g in df.groupby(key_cols, sort=False):
        if g["DoseV"].nunique() != 1:
            n_drop_dose += 1
            continue
        g = g.sort_values("Time")
        cells.append({
            "batch": str(g["Date"].iloc[0]),
            "dose": float(g["DoseV"].iloc[0]),
            "t": g["Time"].values.astype(float) - T_STIM,
            "erk": g["ERK"].values.astype(float),
            "akt": g["Akt"].values.astype(float) if "Akt" in g.columns else None,
            "area": g["Area"].values.astype(float),
        })
    if n_drop_dose:
        log(f"  [guardrail] {os.path.basename(path)}: {n_drop_dose} within-key duplicate doses removed")
    log(f"  load {os.path.basename(path)}: {len(cells)} cells, {len(doses)} dose levels {np.round(doses,4).tolist()}, frame grid {tg.min():.2f}-{tg.max():.2f} min ({len(tg)} frames)")
    return cells, doses

# ---------------- featurization (D2-D5) ----------------
def featurize(c, k=K_MAIN):
    """Return None=removed; otherwise a feature dict. All rules frozen.
    dt takes this cell's actual frame interval (main analysis 0.5 interpolated grid / inhibitor arm 3.5 acquisition grid)."""
    t, v = c["t"], c["erk"]
    if t.size < 2 or t.max() < WIN_HI - 0.1:
        return None
    dt = float(np.median(np.diff(t)))
    if not np.isfinite(dt) or dt <= 0:
        return None
    bad = (~np.isfinite(v)) | (v < 0) | (~np.isfinite(t))
    if bad.mean() > 0.10:
        return None
    vb = v[(t < 0) & ~bad]
    if vb.size < 3:
        return None
    mu0 = float(np.mean(vb)); sd0 = float(np.std(vb, ddof=1))
    if not np.isfinite(mu0) or not np.isfinite(sd0) or sd0 == 0.0:
        return None
    thr = mu0 + k * sd0
    win = (t > WIN_LO) & (t <= WIN_HI) & ~bad
    tw, vw = t[win], v[win]
    hit = vw >= thr
    tau = float(tw[hit][0]) if hit.any() else np.inf
    # count: number of upward crossings within the window (below -> >= transitions)
    ups = 0
    prev = False
    for val in vw:
        cur = val >= thr
        if cur and not prev:
            ups += 1
        prev = cur
    dur = float(hit.sum() * dt)
    pw = (t >= PAID_LO) & (t <= PAID_HI) & ~bad
    peak = float(np.max(v[pw])) if pw.any() else np.nan
    auci = float(np.sum(v[pw]) * dt) if pw.any() else np.nan
    aw = c["area"][np.isfinite(c["area"]) & (c["area"] > 0)]
    if aw.size == 0:
        return None
    return {"tau": tau, "cnt": ups, "dur": dur, "peak": peak, "auci": auci,
            "area": float(np.mean(aw)), "resp": bool(hit.any())}

# ---------------- AUC (Mann-Whitney, average ranks for ties/censoring) ----------------
def auc_mw(low, high):
    """AUC = P(tau_low > tau_high) + 0.5 P(tie). inf enters the ranking naturally."""
    a = np.asarray(low, float); b = np.asarray(high, float)
    allv = np.concatenate([a, b])
    order = allv.argsort(kind="mergesort")
    ranks = np.empty(len(allv), float)
    sv = allv[order]
    i = 0
    while i < len(sv):
        j = i
        while j + 1 < len(sv) and sv[j + 1] == sv[i]:
            j += 1
        ranks[order[i:j + 1]] = 0.5 * (i + j) + 1.0
        i = j + 1
    ra = ranks[:len(a)]
    U = ra.sum() - len(a) * (len(a) + 1) / 2.0
    return U / (len(a) * len(b))

def boot_auc(low, high, rng, n=N_BOOT):
    low = np.asarray(low, float); high = np.asarray(high, float)
    pt = auc_mw(low, high)
    bs = np.empty(n)
    for i in range(n):
        la = low[rng.integers(0, len(low), len(low))]
        hb = high[rng.integers(0, len(high), len(high))]
        bs[i] = auc_mw(la, hb)
    return pt, float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))

# ---------------- adjudication-unit assembly ----------------
def cells_by_dose(cells, feats):
    d = {}
    for c, f in zip(cells, feats):
        if f is not None:
            d.setdefault(c["dose"], []).append(f)
    return d

def pair_unit(by_d, dl, dh, key="tau"):
    low = [f[key] for f in by_d.get(dl, [])]
    high = [f[key] for f in by_d.get(dh, [])]
    return low, high

# ---------------- main pipeline ----------------
def main():
    log(f"code 46 v0.1.5 - P4 blind ruling: ERK re-adjudication (GPCR-KTR trace level)")
    log(f"run time {datetime.datetime.now().isoformat(timespec='seconds')} | SEED={SEED} | bootstrap={N_BOOT} (first round, i.e. execute)")
    log(f"ROOT={ROOT}")
    if not os.path.isdir(ROOT):
        abort(f"data root directory missing: {ROOT} (prepare data per pre-registration section 5 item 2 first)")
    rng = np.random.default_rng(SEED)

    unit_rows, rob_rows, desc_rows = [], [], []
    store = {}   # (lig, zone, pair_idx, side) -> tau array, for P4-2 aggregate bootstrap
    main_data = {}  # lig -> (cells, doses), reused by the Akt descriptive arm

    for lig in LIG_MAIN:
        path = os.path.join(ROOT, FILE_F1.format(lig=lig))
        log(f"\n=== main-analysis ligand {lig} (DMSO) ===")
        cells, doses = load_long(path)
        main_data[lig] = (cells, doses)
        feats = [featurize(c) for c in cells]
        n_ok = sum(f is not None for f in feats)
        log(f"  included {n_ok}/{len(cells)} cells (removed {len(cells)-n_ok}: baseline frames <3 / sigma_0=0 / end <38.5 / missing >10%)")
        by_d = cells_by_dose(cells, feats)
        for d in doses:
            log(f"    level {d:g}: n={len(by_d.get(d, []))}")

        # 6 levels, 5 pairs anchored: first 2 pairs low zone, last 2 pairs high zone; middle pair (3,4) = half-saturation transition band,
        # excluded from adjudication, moved to the descriptive arm (pre-registration v0.3)
        low_pairs  = [(doses[0], doses[1]), (doses[1], doses[2])]
        high_pairs = [(doses[3], doses[4]), (doses[4], doses[5])]
        mid_pair   = (doses[2], doses[3])

        # ---- P4-1 / P4-2 units ----
        for zone, pairs in (("低区", low_pairs), ("高区", high_pairs)):
            for pi, (dl, dh) in enumerate(pairs):
                low, high = pair_unit(by_d, dl, dh)
                suff = len(low) >= N_MIN and len(high) >= N_MIN
                if suff:
                    pt, lo, hi = boot_auc(low, high, rng)
                    hit = pt >= M41
                else:
                    pt = lo = hi = float("nan"); hit = False
                store[(lig, zone, pi, "low")] = np.asarray(low, float)
                store[(lig, zone, pi, "high")] = np.asarray(high, float)
                unit_rows.append({"clause": "P4-1" if zone == "低区" else "P4-2",
                                  "ligand": lig, "zone": zone,
                                  "pair": f"{dl:g}↔{dh:g}", "n_low": len(low), "n_high": len(high),
                                  "auc": pt, "ci_lo": lo, "ci_hi": hi,
                                  "line": M41, "pass": hit, "sufficient": suff})
                log(f"  [{zone}] tau({K_MAIN}) {dl:g}<->{dh:g}: AUC={pt:.4f} [{lo:.3f},{hi:.3f}] n={len(low)}/{len(high)} at-line={hit if suff else 'insufficient data'}")

        # ---- robust arm k=2/4 (point estimates) ----
        for kk in K_ROB:
            feats_k = [featurize(c, k=kk) for c in cells]
            by_k = cells_by_dose(cells, feats_k)
            for zone, pairs in (("低区", low_pairs), ("高区", high_pairs)):
                for dl, dh in pairs:
                    low, high = pair_unit(by_k, dl, dh)
                    if len(low) >= N_MIN and len(high) >= N_MIN:
                        rob_rows.append({"arm": f"k={kk}", "ligand": lig, "zone": zone,
                                         "pair": f"{dl:g}↔{dh:g}", "auc": auc_mw(low, high)})

        # ---- batch-purity arm (largest single Experiment) ----
        exps = pd.Series([c["batch"] for c in cells]).value_counts()
        big = exps.index[0]
        sub_cells = [c for c in cells if c["batch"] == big]
        sub = [featurize(c) for c in sub_cells]
        by_b = cells_by_dose(sub_cells, sub)
        for zone, pairs in (("低区", low_pairs), ("高区", high_pairs)):
            for dl, dh in pairs:
                low, high = pair_unit(by_b, dl, dh)
                if len(low) >= N_MIN and len(high) >= N_MIN:
                    rob_rows.append({"arm": f"batch纯净({big})", "ligand": lig, "zone": zone,
                                     "pair": f"{dl:g}↔{dh:g}", "auc": auc_mw(low, high)})

        # ---- P4-3: meanArea three-stratum penalty (all 5 pairs, including the middle transition-band pair) ----
        all_pairs = [(doses[i], doses[i + 1]) for i in range(5)]
        for dl, dh in all_pairs:
            pool = [(f["tau"], f["area"]) for d in (dl, dh) for f in by_d.get(d, [])]
            if len(pool) < 2 * N_MIN:
                unit_rows.append({"clause": "P4-3", "ligand": lig, "zone": "全",
                                  "pair": f"{dl:g}↔{dh:g}", "n_low": len(pool), "n_high": 0,
                                  "auc": float("nan"), "ci_lo": float("nan"), "ci_hi": float("nan"),
                                  "line": P43_HI, "pass": False, "sufficient": False})
                continue
            taus = np.array([p[0] for p in pool]); areas = np.array([p[1] for p in pool])
            q1, q2 = np.percentile(areas, [33.33, 66.67])
            lab = np.where(areas <= q1, 0, np.where(areas <= q2, 1, 2))
            # rebuild grouping: pool order = first dl then dh
            n_dl = len(by_d.get(dl, []))
            is_low = np.array([True] * n_dl + [False] * (len(pool) - n_dl))
            auc_pool = auc_mw(taus[is_low], taus[~is_low])
            layer_aucs, layer_ns = [], []
            for L in range(3):
                m = lab == L
                if m.sum() >= 2 * N_MIN and (m & is_low).sum() >= N_MIN // 2 and (m & ~is_low).sum() >= N_MIN // 2:
                    layer_aucs.append(auc_mw(taus[m & is_low], taus[m & ~is_low]))
                    layer_ns.append(int(m.sum()))
            delta = auc_pool - float(np.mean(layer_aucs)) if layer_aucs else float("nan")
            unit_rows.append({"clause": "P4-3", "ligand": lig, "zone": "全",
                              "pair": f"{dl:g}↔{dh:g}", "n_low": n_dl, "n_high": len(pool) - n_dl,
                              "auc": delta, "ci_lo": float("nan"), "ci_hi": float("nan"),
                              "line": P43_HI, "pass": abs(delta) <= P43_LO if np.isfinite(delta) else False,
                              "sufficient": bool(layer_aucs)})
            log(f"  [P4-3] {dl:g}<->{dh:g}: delta={delta:+.4f} (pooled {auc_pool:.4f} - within-strata mean {np.mean(layer_aucs) if layer_aucs else float('nan'):.4f}; strata n={layer_ns})")

        # ---- descriptive arm: paid peak/integral + count/duration + response rate ----
        for dl, dh in all_pairs:
            for key, nm in (("peak", "peakCN"), ("auci", "integralCN"), ("cnt", "count"), ("dur", "duration")):
                low, high = pair_unit(by_d, dl, dh, key)
                low = [x for x in low if np.isfinite(x)]; high = [x for x in high if np.isfinite(x)]
                if len(low) >= N_MIN and len(high) >= N_MIN:
                    # paid-type direction = larger at higher dose -> AUC = P(high > low)
                    desc_rows.append({"arm": nm, "ligand": lig, "pair": f"{dl:g}↔{dh:g}",
                                      "auc": auc_mw(high, low)})
        for d in doses:
            fs = by_d.get(d, [])
            if fs:
                desc_rows.append({"arm": "响应率", "ligand": lig, "pair": f"档{d:g}",
                                  "auc": float(np.mean([f["resp"] for f in fs]))})
        # middle-pair (half-saturation transition band) tau readings: excluded from adjudication, into the descriptive arm
        low, high = pair_unit(by_d, mid_pair[0], mid_pair[1])
        if len(low) >= N_MIN and len(high) >= N_MIN:
            desc_rows.append({"arm": "中间对τ(过渡带)", "ligand": lig,
                              "pair": f"{mid_pair[0]:g}↔{mid_pair[1]:g}", "auc": auc_mw(low, high)})

    # ---------------- P4-2 aggregate adjudication statistic ----------------
    log("\n=== P4-2 aggregation (saturation-zone disordering test) ===")
    aucL_list, aucH_list = [], []
    for lig in LIG_MAIN:
        for pi in range(2):
            low = store[(lig, "低区", pi, "low")]; high = store[(lig, "低区", pi, "high")]
            if len(low) >= N_MIN and len(high) >= N_MIN:
                aucL_list.append(auc_mw(low, high))
            low = store[(lig, "高区", pi, "low")]; high = store[(lig, "高区", pi, "high")]
            if len(low) >= N_MIN and len(high) >= N_MIN:
                aucH_list.append(auc_mw(low, high))
    AUC_L = float(np.mean(aucL_list)); AUC_H = float(np.mean(aucH_list))
    D_pt = AUC_L - AUC_H
    # bootstrap of D: resample all 8 units synchronously
    d_bs = np.empty(N_BOOT)
    for i in range(N_BOOT):
        ls, hs = [], []
        for lig in LIG_MAIN:
            for pi in range(2):
                for zone, acc in (("低区", ls), ("高区", hs)):
                    low = store[(lig, zone, pi, "low")]; high = store[(lig, zone, pi, "high")]
                    if len(low) >= N_MIN and len(high) >= N_MIN:
                        acc.append(auc_mw(low[rng.integers(0, len(low), len(low))],
                                          high[rng.integers(0, len(high), len(high))]))
        d_bs[i] = np.mean(ls) - np.mean(hs)
    D_lo, D_hi = float(np.percentile(d_bs, 2.5)), float(np.percentile(d_bs, 97.5))
    log(f"  AUC_L(low-zone mean)={AUC_L:.4f} | AUC_H(high-zone mean)={AUC_H:.4f} | D={D_pt:+.4f} [{D_lo:+.3f},{D_hi:+.3f}]")

    # ---------------- D8 adjudication (frozen table executed verbatim) ----------------
    log("\n=== D8 adjudication (frozen table executed verbatim) ===")
    verdicts = {}

    p41 = [r for r in unit_rows if r["clause"] == "P4-1"]
    if any(not r["sufficient"] for r in p41):
        verdicts["P4-1"] = "数据不足"
    elif all(r["pass"] for r in p41):
        verdicts["P4-1"] = "hit"
    elif any(np.isfinite(r["auc"]) and r["auc"] < 0.5 for r in p41) or not any(r["pass"] for r in p41):
        verdicts["P4-1"] = "证伪"
    else:
        verdicts["P4-1"] = "miss-中间态"

    if len(aucH_list) < 4 or len(aucL_list) < 4:
        verdicts["P4-2"] = "数据不足"
    elif AUC_H < M42_H and D_pt >= M42_D:
        verdicts["P4-2"] = "hit"
    elif AUC_H >= M42_H or D_pt <= 0:
        verdicts["P4-2"] = "证伪"
    else:
        verdicts["P4-2"] = "miss-中间态"

    p43 = [r for r in unit_rows if r["clause"] == "P4-3" and r["sufficient"]]
    n43_lo = sum(abs(r["auc"]) <= P43_LO for r in p43)
    n43_hi = sum(abs(r["auc"]) >= P43_HI for r in p43)
    if len(p43) < 10:
        verdicts["P4-3"] = "数据不足"
    elif n43_hi >= P43_KILL:
        verdicts["P4-3"] = "证伪"
    elif n43_lo >= P43_HIT:
        verdicts["P4-3"] = "hit"
    else:
        verdicts["P4-3"] = "miss-中间态"

    overall = "P4 证伪（exists证伪条款）" if any(v == "证伪" for v in verdicts.values()) else \
              ("P4 hit（none证伪条款且existshit条款）" if any(v == "hit" for v in verdicts.values())
               else "P4 miss-中间态")
    for k, v in verdicts.items():
        log(f"  {k}：{v}")
    log(f"  overall ruling: {overall}")
    log(f"  P4-3 count：|Δ|≤0.01 → {n43_lo}/10；|Δ|≥0.03 → {n43_hi}/10")

    # ---------------- main output (saved before the descriptive arm: the main ruling archive is unaffected by descriptive-arm exceptions) ----------------
    pd.DataFrame(unit_rows).to_csv(os.path.join(OUT_DIR, "代码46_单元table.csv"),
                                   index=False, encoding="utf-8-sig")
    pd.DataFrame(rob_rows).to_csv(os.path.join(OUT_DIR, "代码46_robust臂.csv"),
                                  index=False, encoding="utf-8-sig")
    with open(os.path.join(OUT_DIR, "代码46_ruling.json"), "w", encoding="utf-8") as f:
        json.dump({"script": "代码46 v0.1.5", "seed": int(SEED), "n_boot": int(N_BOOT),
                   "verdicts": {k: str(v) for k, v in verdicts.items()}, "overall": str(overall),
                   "P4-2": {"AUC_L": float(AUC_L), "AUC_H": float(AUC_H),
                            "D": float(D_pt), "D_ci": [float(D_lo), float(D_hi)]},
                   "P4-3": {"n_le_001": int(n43_lo), "n_ge_003": int(n43_hi)},
                   "timestamp": datetime.datetime.now().isoformat(timespec="seconds")},
                  f, ensure_ascii=False, indent=2,
                  default=lambda o: int(o) if isinstance(o, np.integer)
                         else float(o) if isinstance(o, np.floating) else str(o))
    log("\nmain output saved: 代码46_单元table.csv / 代码46_robust臂.csv / 代码46_ruling.json")
    flush_log()

    # ---------------- descriptive arm (not adjudicated; any exception only warns and skips) ----------------
    log("\n=== descriptive arm (not adjudicated) ===")
    try:
        path = os.path.join(ROOT, FILE_F1.format(lig=LIG_DESC))
        cells, doses = load_long(path, desc=True)
        feats = [featurize(c) for c in cells]
        by_d = cells_by_dose(cells, feats)
        for i in range(5):
            low, high = pair_unit(by_d, doses[i], doses[i + 1])
            if len(low) >= N_MIN and len(high) >= N_MIN:
                desc_rows.append({"arm": "S1P-τ", "ligand": "S1P", "pair": f"{doses[i]:g}↔{doses[i+1]:g}",
                                  "auc": auc_mw(low, high)})
        for d in doses:
            fs = by_d.get(d, [])
            if fs:
                desc_rows.append({"arm": "S1P-响应率", "ligand": "S1P", "pair": f"档{d:g}",
                                  "auc": float(np.mean([f["resp"] for f in fs]))})
        log(f"  S1P descriptive arm done ({len([r for r in desc_rows if r['ligand']=='S1P'])} rows)")
    except Exception as e:
        log(f"  [warning] S1P descriptive arm skipped: {e}")

    # Akt descriptive arm (reuses the main-analysis cache, switches to the CN_AktRB channel)
    for lig in LIG_MAIN:
        try:
            cells, doses = main_data[lig]
            akt_cells = [c for c in cells if c["akt"] is not None]
            if not akt_cells:
                raise ValueError("no Akt column or all empty")
            feats_a = []
            for c in akt_cells:
                cc = dict(c); cc["erk"] = c["akt"]
                feats_a.append(featurize(cc))
            by_a = cells_by_dose(akt_cells, feats_a)
            for i in range(5):
                low, high = pair_unit(by_a, doses[i], doses[i + 1])
                if len(low) >= N_MIN and len(high) >= N_MIN:
                    desc_rows.append({"arm": "Akt-τ", "ligand": lig, "pair": f"{doses[i]:g}↔{doses[i+1]:g}",
                                      "auc": auc_mw(low, high)})
            log(f"  Akt descriptive arm done ({lig})")
        except Exception as e:
            log(f"  [warning] Akt descriptive arm {lig} skipped: {e}")

    # inhibitor descriptive arm (Figure_3, may be absent)
    for lig in LIG_MAIN:
        p3 = os.path.join(ROOT, FILE_F3.format(lig3=LIG3_NAME[lig]))
        if not os.path.isfile(p3):
            log(f"  [warning] inhibitor arm file missing, skipped: {p3}")
            continue
        try:
            df = pd.read_csv(p3)
            ren = {}
            for c in df.columns:
                k = str(c).strip().lower()
                if k in COLMAP and COLMAP[k] not in ren.values():
                    ren[c] = COLMAP[k]
            df = df.rename(columns=ren)
            if "Inh" not in df.columns:
                log(f"  [warning] {os.path.basename(p3)} has no Inhibitor column, inhibitor arm skipped")
                continue
            # grid check: both the 0.5 min interpolated grid (Figure_1 form) and the 3.5 min acquisition grid
            # (Figure_3 form) are legal; t0=23.5 applies to both (tp7)
            tgi = np.sort(df["Time"].unique())
            sti = np.median(np.diff(tgi)) if len(tgi) > 1 else float("nan")
            ok_grid = len(tgi) >= 2 and abs(tgi.min() - T_START) <= 0.1 and tgi.max() >= T_END_MIN and \
                      (abs(sti - DT_GLOB) <= 0.05 or abs(sti - 3.5) <= 0.2)
            if not ok_grid:
                log(f"  [warning] {os.path.basename(p3)} gridexception（{tgi.min() if len(tgi) else float('nan'):.1f}–"
                    f"{tgi.max() if len(tgi) else float('nan'):.1f}, step median {sti:.3f}), inhibitor arm skipped")
                continue
            df["DoseV"] = df["Dose"].apply(parse_dose)
            for inh, g0 in df.groupby("Inh"):
                doses_i = np.sort(g0["DoseV"].dropna().unique())
                if len(doses_i) < 2:
                    continue
                key_i = ["UID"] if "UID" in g0.columns else \
                        [c for c in ["Date", "Slide", "Pos", "Obj"] if c in g0.columns]
                if not key_i:
                    log(f"  [warning] inhibitor arm {lig}/{inh}: no usable cell key, skipped")
                    continue
                cells_i = []
                for _k, g in g0.groupby(key_i, sort=False):
                    if g["DoseV"].nunique() != 1:
                        continue
                    g = g.sort_values("Time")
                    cells_i.append({"batch": "NA", "dose": float(g["DoseV"].iloc[0]),
                                    "t": g["Time"].values.astype(float) - T_STIM,
                                    "erk": g["ERK"].values.astype(float), "akt": None,
                                    "area": g["Area"].values.astype(float)})
                feats_i = [featurize(c) for c in cells_i]
                by_i = cells_by_dose(cells_i, feats_i)
                for j in range(len(doses_i) - 1):
                    low, high = pair_unit(by_i, doses_i[j], doses_i[j + 1])
                    if len(low) >= N_MIN and len(high) >= N_MIN:
                        desc_rows.append({"arm": f"inhibition剂-{inh}", "ligand": lig,
                                          "pair": f"{doses_i[j]:g}↔{doses_i[j+1]:g}",
                                          "auc": auc_mw(low, high)})
            log(f"  inhibitor arm done ({lig})")
        except Exception as e:
            log(f"  [warning] inhibitor arm {lig} exception, skipped: {e}")

    pd.DataFrame(desc_rows).to_csv(os.path.join(OUT_DIR, "代码46_描述臂.csv"),
                                   index=False, encoding="utf-8-sig")
    log("\noutput: 代码46_单元table.csv / 代码46_robust臂.csv / 代码46_描述臂.csv / 代码46_ruling.json / 代码46_判定日志.txt")
    log("- recorded either way: please return the five outputs; they enter the ledger as-is regardless of outcome. -")
    flush_log()

if __name__ == "__main__":
    main()
