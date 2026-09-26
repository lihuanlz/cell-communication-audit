#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code 48 · P6 blind ruling: FCD attribution — NF-κB sequential grid fold vs absolute dose
Pre-registration: 预注册_P6_FCD归属裁决_NFkB序贯网格_v01.md (registration no. 95, frozen)
Data: Wang AG, Son M, et al. 2022 Cell Rep 40:111159
     GitHub tay-lab/Sequential_NF-kB_stim
     scmat_sequentialstim.mat + featmat_sequentialstim.mat
v0.1.1 · 2026-08-15 (schema alignment fix, following the empty-run rerun precedent) · seed 20260815 · bootstrap 2000
Criterion (pre-registration §3, frozen): on the 36-unit S2 decay A_c
  ΔR2 = R2(A~log10 fold) - R2(A~log10 abs1)
  >= +0.15 → H1 hit (fold carrier); <= -0.15 → H-abs (absolute dose); in between → intermediate state
  Additional: the fold slope must be positive (higher preceding dose → stronger decay); if reversed, downgrade to intermediate state.
Guardrails (pre-registration §5): schema core check / cells per condition >= 50 / naive response rate <20% removed /
  <24 units power limitation / <18 abort / A<=0 removed.
Run: python3 代码48.py [DATA_DIR]
"""
import os, sys, json, hashlib, datetime
import numpy as np

SEED = 20260815
N_BOOT = 2000
DELTA_LINE = 0.15
MIN_N_COND = 50
MIN_NAIVE_RATE = 0.20
MIN_UNITS_PWR = 24
MIN_UNITS_ABORT = 18
INTERVAL_MIN = 6
STIM_MIN = [120, 240, 360, 480]          # §1 schema verification item to be checked

# dose table (ng/mL, paper methods; ligand order T,I,L,P × group H,M,L)
LIGANDS = ["TNF", "IL1", "LPS", "PAM"]
DOSE = {  # (ligand_idx, group) -> ng/mL ; group 1=H 2=M 3=L
    (0,1):90.0, (0,2):30.0, (0,3):3.0,
    (1,1):3.0,  (1,2):0.2,  (1,3):0.05,
    (2,1):400.0,(2,2):100.0,(2,3):12.5,
    (3,1):1.0,  (3,2):0.1,  (3,3):0.01,
}

LOG_PATH = None
def log(msg=""):
    print(msg)
    if LOG_PATH:
        with open(LOG_PATH, "a", encoding="utf-8") as f:
            f.write(str(msg) + "\n")

def linfit_r2(x, y):
    """Univariate linear regression R2 (with intercept)."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    X = np.column_stack([np.ones_like(x), x])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    ss_res = float(resid @ resid)
    ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return r2, float(beta[1])

def main():
    global LOG_PATH
    root = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
    outdir = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
    LOG_PATH = os.path.join(outdir, "代码48_判定日志.txt")
    if os.path.exists(LOG_PATH):
        os.remove(LOG_PATH)

    log("=" * 74)
    log("Code 48 · P6 blind ruling: FCD attribution — NF-κB sequential grid fold vs absolute dose  v0.1.1")
    log(f"Run time: {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
    log(f"seed={SEED}  bootstrap={N_BOOT}  criterion ΔR2=±{DELTA_LINE}")
    log(f"Pre-registration: 预注册_P6_FCD归属裁决_NFkB序贯网格_v01.md (registration no. 95, frozen)")
    log(f"Data root directory: {root}")
    log("=" * 74)

    import scipy.io as sio
    p_scmat = os.path.join(root, "scmat_sequentialstim.mat")
    p_feat  = os.path.join(root, "featmat_sequentialstim.mat")
    for p in (p_scmat, p_feat):
        if not os.path.isfile(p):
            log(f"[ABORT] missing file {p} — empty run, does not constitute a ruling"); sys.exit(2)
    scmat = sio.loadmat(p_scmat)["scmatcomb_norm"]
    feat  = sio.loadmat(p_feat)
    fkeys = [k for k in feat if not k.startswith("__")]
    log(f"scmat {scmat.shape}; featmat key {fkeys}")
    feat = feat[fkeys[0]]
    log(f"featmat {feat.shape}")

    # ---------- §5-1 schema core check ----------
    # v0.1.1: featmat rows are aligned to scmat after removing NaN rows (same filter as the authors' script)
    if scmat.shape[1] != 172 or feat.shape[1] < 28:
        log("[ABORT] schema mismatch (column count) — empty run"); sys.exit(2)
    nanrows = np.isnan(scmat).any(axis=1)   # authors' definition: a row is dropped if any of the 172 columns is NaN
    log(f"[schema] scmat contains {int(nanrows.sum())} NaN rows, removed per the authors' definition to align with featmat")
    scmat = scmat[~nanrows]
    if feat.shape[0] != scmat.shape[0]:
        log(f"[ABORT] row counts still mismatch scmat{scmat.shape[0]} vs feat{feat.shape[0]} — empty run"); sys.exit(2)
    log(f"[schema] cell count after alignment {scmat.shape[0]}")
    grp = scmat[:, 171].astype(int)
    codes = scmat[:, 167:171].astype(int)     # S1..S4 dose codes
    ngrp = sorted(set(grp.tolist()))
    if ngrp != [1, 2, 3]:
        log(f"[ABORT] dose groups {ngrp} are not [1,2,3] — empty run"); sys.exit(2)
    # dose code → (ligand, group): within-group codes 0-3/4-7/8-11 map to T/I/L/P
    def code2lg(cd, g):
        base = (g - 1) * 4
        if cd >= 12:
            return None
        return cd - base
    bad = 0
    for g in (1, 2, 3):
        cs = set(codes[grp == g].ravel().tolist())
        if not cs.issubset(set(range((g-1)*4, (g-1)*4+4)) | {12}):
            log(f"[ABORT] group {g} dose codes out of range {sorted(cs)} — empty run"); sys.exit(2)
        for cd in cs:
            if cd != 12 and code2lg(cd, g) not in (0, 1, 2, 3):
                bad += 1
    if bad:
        log("[ABORT] dose-code mapping failure — empty run"); sys.exit(2)
    log("[schema] dose-group and dose-code mapping core check passed (within-group 0-3/4-7/8-11 → T/I/L/P)")

    # main condition = (S1,S2,S3,S4,group); conditions containing 12 are excluded (same-ligand quadruple stimulation/control → descriptive arm)
    is_ctrl = (codes == 12).any(axis=1)
    log(f"Main-grid cells {int((~is_ctrl).sum())}; control/quadruple-stimulation cells {int(is_ctrl.sum())}")

    # ---------- Authors' features, main ruling ----------
    PEAK_COLS = [0, 7, 14, 21]               # featmat per-pulse max (1-based 1,8,15,22)
    peak = feat[:, PEAK_COLS]

    # noise floor: 95th percentile of control-condition S1-position peak (§5-3)
    floor = float(np.nanpercentile(peak[is_ctrl, 0], 95))
    log(f"Noise floor (control S1 peak 95th percentile) = {floor:.4f}")

    # naive reference: for each (j,g): all main-grid cells with S1 position = ligand j
    naive_med, naive_rate, naive_n = {}, {}, {}
    for j in range(4):
        for g in (1, 2, 3):
            m = (~is_ctrl) & (grp == g) & np.array([code2lg(c, g) == j for c, g2 in zip(codes[:, 0], grp)])
            v = peak[m, 0]
            v = v[~np.isnan(v)]
            naive_med[(j, g)] = float(np.median(v)) if len(v) else np.nan
            naive_rate[(j, g)] = float((v > floor).mean()) if len(v) else np.nan
            naive_n[(j, g)] = len(v)
    log("\n[naive reference] (ligand, group) median peak / response rate / n:")
    for j in range(4):
        for g in (1, 2, 3):
            log(f"  {LIGANDS[j]:4s} g{g}: med={naive_med[(j,g)]:.4f} rate={naive_rate[(j,g)]:.3f} n={naive_n[(j,g)]}")

    # §5-3 removal: (j,g) with naive response rate < 20%
    excl_jg = [k for k in naive_rate if not (naive_rate[k] >= MIN_NAIVE_RATE)]
    if excl_jg:
        log(f"[removed] naive response rate <{MIN_NAIVE_RATE}: {[(LIGANDS[j], g) for j, g in excl_jg]}")

    # ---------- 36 units ----------
    units = []   # dict per (i,j,g)
    drop_log = []
    for g in (1, 2, 3):
        for i in range(4):
            for j in range(4):
                if i == j:
                    continue
                if (j, g) in excl_jg:
                    continue
                # condition: S1=code i, S2=code j, group g (both S3/S4 arrangements merged)
                m = (~is_ctrl) & (grp == g) \
                    & (codes[:, 0] == (g-1)*4 + i) & (codes[:, 1] == (g-1)*4 + j)
                v = peak[m, 1]
                v = v[~np.isnan(v)]
                if len(v) < MIN_N_COND:
                    drop_log.append((i, j, g, f"n={len(v)}"))
                    continue
                if not (naive_med[(j, g)] > 0):   # §5-5: naive median = 0 → division by zero, removed
                    drop_log.append((i, j, g, f"naive_med={naive_med[(j,g)]:.4f}"))
                    continue
                A = float(np.median(v)) / naive_med[(j, g)]
                if not np.isfinite(A) or A <= 0:
                    drop_log.append((i, j, g, f"A={A}"))
                    continue
                d1 = DOSE[(i, g)]; d2 = DOSE[(j, g)]
                units.append(dict(i=i, j=j, g=g, n=len(v), A=A,
                                  d1=d1, d2=d2, fold=d2/d1))
    log(f"\n[units] valid {len(units)}/36; removed {len(drop_log)}: {drop_log}")
    if len(units) < MIN_UNITS_ABORT:
        log("[ABORT] valid units <18 — empty run"); sys.exit(2)
    power_note = len(units) < MIN_UNITS_PWR
    if power_note:
        log(f"[LIMITATION] valid units {len(units)} < {MIN_UNITS_PWR} — reduced-power limitation recorded")

    lgfold = np.array([np.log10(u["fold"]) for u in units])
    lgabs1 = np.array([np.log10(u["d1"]) for u in units])
    lgabs2 = np.array([np.log10(u["d2"]) for u in units])
    A      = np.array([u["A"] for u in units])

    r2_fold, b_fold = linfit_r2(lgfold, A)
    r2_abs1, b_abs1 = linfit_r2(lgabs1, A)
    r2_abs2, b_abs2 = linfit_r2(lgabs2, A)
    dR2 = r2_fold - r2_abs1

    log("\n" + "=" * 74)
    log("Ruling (main: authors' features, 36-unit S2 decay)")
    log("-" * 74)
    log(f"  R2(A~log fold) = {r2_fold:.4f} (slope {b_fold:+.4f})")
    log(f"  R2(A~log abs1) = {r2_abs1:.4f} (slope {b_abs1:+.4f})")
    log(f"  R2(A~log abs2) = {r2_abs2:.4f} (slope {b_abs2:+.4f}, control predictor)")
    log(f"  ΔR2 = {dR2:+.4f} (criterion ±{DELTA_LINE})")

    sign_ok = b_fold > 0
    if dR2 >= DELTA_LINE and sign_ok:
        verdict = "P6 hit（fold 载波）"
    elif dR2 >= DELTA_LINE and not sign_ok:
        verdict = "P6 中间态（ΔR2 达标但 fold slopereverse，按pre-registration降级）"
    elif dR2 <= -DELTA_LINE:
        verdict = "P6 证伪（绝对dose载波）"
    else:
        verdict = "P6 中间态"
    log(f"  Ruling: {verdict}")
    log("=" * 74)

    # ---------- bootstrap CI (descriptive) ----------
    rng = np.random.default_rng(SEED)
    boots = []
    nU = len(units)
    for _ in range(N_BOOT):
        idx = rng.integers(0, nU, nU)
        r2f, _ = linfit_r2(lgfold[idx], A[idx])
        r2a, _ = linfit_r2(lgabs1[idx], A[idx])
        boots.append(r2f - r2a)
    ci = np.percentile(boots, [2.5, 97.5])
    log(f"[CI] ΔR2 bootstrap 95% = [{ci[0]:+.4f}, {ci[1]:+.4f}] (descriptive, not part of the adjudication)")

    # ---------- §4 robust arm ----------
    rob = []
    # ① self-built features (local baseline)
    fr = [t // INTERVAL_MIN for t in STIM_MIN]   # onset frames 20/40/60/80
    if scmat.shape[1] < 83:
        log("[robust arm ①] insufficient trajectory columns — skip"); selfbuilt = None
    else:
        tr = scmat[:, :83]
        sb_peak = np.full((tr.shape[0], 4), np.nan)
        for p_, o in enumerate(fr):
            base = np.nanmedian(tr[:, max(0, o-5):o], axis=1)   # preceding 30 min (5 frames × 6 min)
            w = tr[:, o:min(83, o+20)]                           # 0..114 min window
            sb_peak[:, p_] = np.nanmax(w, axis=1) - base
        A_sb = []
        keep = []
        for k, u in enumerate(units):
            i, j, g = u["i"], u["j"], u["g"]
            m = (~is_ctrl) & (grp == g) & (codes[:, 0] == (g-1)*4+i) & (codes[:, 1] == (g-1)*4+j)
            v = sb_peak[m, 1]; v = v[~np.isnan(v)]
            mn = (~is_ctrl) & (grp == g) & np.array([code2lg(c, g) == j for c, g2 in zip(codes[:, 0], grp)])
            vn = sb_peak[mn, 0]; vn = vn[~np.isnan(vn)]
            if len(v) < MIN_N_COND or len(vn) == 0:
                keep.append(False); A_sb.append(np.nan); continue
            med_n = float(np.median(vn))
            a = float(np.median(v))/med_n if med_n > 0 else np.nan
            A_sb.append(a); keep.append(np.isfinite(a) and a > 0)
        keep = np.array(keep)
        if keep.sum() >= MIN_UNITS_ABORT:
            A_sb = np.array(A_sb)[keep]
            r2f_sb, bf_sb = linfit_r2(lgfold[keep], A_sb)
            r2a_sb, ba_sb = linfit_r2(lgabs1[keep], A_sb)
            d_sb = r2f_sb - r2a_sb
            rob.append(("自建特征(局部baseline)", int(keep.sum()), r2f_sb, r2a_sb, d_sb))
            log(f"[robust①] self-built features n={keep.sum()}: R2fold={r2f_sb:.4f} R2abs1={r2a_sb:.4f} ΔR2={d_sb:+.4f} (slope {bf_sb:+.4f})")
            selfbuilt = dict(r2f=r2f_sb, r2a=r2a_sb, d=d_sb, n=int(keep.sum()))
        else:
            log(f"[robust①] self-built features: too few valid units ({keep.sum()}) — skip"); selfbuilt = None
    # ② dual-predictor increment
    X2 = np.column_stack([np.ones(nU), lgfold, lgabs1])
    beta, *_ = np.linalg.lstsq(X2, A, rcond=None)
    resid = A - X2 @ beta
    r2_both = 1 - float(resid @ resid) / float(((A-A.mean())**2).sum())
    log(f"[robust②] dual predictors R2={r2_both:.4f} (fold increment {r2_both-r2_abs1:+.4f}; abs1 increment {r2_both-r2_fold:+.4f})")
    rob.append(("双prediction子", nU, r2_fold, r2_abs1, r2_both))
    # ③ cell-count weighting
    w = np.sqrt(np.array([u["n"] for u in units], float))
    def wlin_r2(x, y, w):
        X = np.column_stack([np.ones_like(x), x]) * w[:, None]
        beta, *_ = np.linalg.lstsq(X, y*w, rcond=None)
        r = y - np.column_stack([np.ones_like(x), x]) @ beta
        return 1 - float((w*r) @ (w*r)) / float(((w*(y - np.average(y, weights=w**2)))**2).sum())
    r2f_w = wlin_r2(lgfold, A, w); r2a_w = wlin_r2(lgabs1, A, w)
    log(f"[robust③] weighted: ΔR2={r2f_w-r2a_w:+.4f}")
    rob.append(("cell数加权", nU, r2f_w, r2a_w, r2f_w-r2a_w))
    # ④ Spearman
    from scipy import stats as sst
    rho_f = sst.spearmanr(lgfold, A).statistic
    rho_a = sst.spearmanr(lgabs1, A).statistic
    log(f"[robust④] Spearman: ρ_fold={rho_f:+.4f} ρ_abs1={rho_a:+.4f} (Δρ={rho_f-rho_a:+.4f})")
    rob.append(("Spearman", nU, rho_f**2, rho_a**2, rho_f-rho_a))

    # ---------- §6 descriptive arm ----------
    log("\n[descriptive arm] S3/S4 decay attribution (same-method extrapolation, not part of the adjudication):")
    for pos, pname in ((2, "S3"), (3, "S4")):
        u2 = []
        for g in (1, 2, 3):
            for i in range(4):
                for j in range(4):
                    if i == j or (j, g) in excl_jg:
                        continue
                    m = (~is_ctrl) & (grp == g) & (codes[:, pos-1] == (g-1)*4+i) & (codes[:, pos] == (g-1)*4+j)
                    v = peak[m, pos]; v = v[~np.isnan(v)]
                    # naive reference uses the same-group S1 position instead (simplified descriptive definition)
                    if len(v) < MIN_N_COND or naive_med[(j, g)] <= 0:
                        continue
                    u2.append((np.log10(DOSE[(j, g)]/DOSE[(i, g)]), np.log10(DOSE[(i, g)]),
                               float(np.median(v))/naive_med[(j, g)]))
        if len(u2) >= MIN_UNITS_ABORT:
            u2 = np.array(u2)
            r2f2, _ = linfit_r2(u2[:, 0], u2[:, 2]); r2a2, _ = linfit_r2(u2[:, 1], u2[:, 2])
            log(f"  {pname}: n={len(u2)} R2fold={r2f2:.4f} R2abs1={r2a2:.4f} ΔR2={r2f2-r2a2:+.4f}")
        else:
            log(f"  {pname}: insufficient valid units ({len(u2)})")

    log("\n[descriptive arm] per-pulse median peaks for same-ligand quadruple-stimulation/control conditions:")
    for g in (1, 2, 3):
        m = is_ctrl & (grp == g)
        if m.sum() >= 10:
            med = [float(np.nanmedian(peak[m, p_])) for p_ in range(4)]
            log(f"  group{g} n={int(m.sum())}: " + " ".join(f"{x:.3f}" for x in med))

    # ---------- output ----------
    import csv
    with open(os.path.join(outdir, "代码48_单元table.csv"), "w", newline="", encoding="utf-8") as f:
        wcsv = csv.writer(f)
        wcsv.writerow(["S1ligand","S2ligand","组","d1_ngml","d2_ngml","fold","log10fold","log10abs1","n","A_作者特征"])
        for u in units:
            wcsv.writerow([LIGANDS[u["i"]], LIGANDS[u["j"]], u["g"], u["d1"], u["d2"],
                           f"{u['fold']:.4f}", f"{np.log10(u['fold']):.4f}",
                           f"{np.log10(u['d1']):.4f}", u["n"], f"{u['A']:.6f}"])
    with open(os.path.join(outdir, "代码48_robust臂.csv"), "w", newline="", encoding="utf-8") as f:
        wcsv = csv.writer(f); wcsv.writerow(["臂","n","R2fold_or_rho2","R2abs1_or_rho2","delta"])
        for r in rob:
            wcsv.writerow([r[0], r[1], f"{r[2]:.4f}", f"{r[3]:.4f}", f"{r[4]:+.4f}"])
    verdict_json = {
        "code": "代码48 v0.1.0", "prereg": "P6 v0.1（注册·九十五）", "seed": SEED,
        "n_units": len(units), "power_note": bool(power_note),
        "excluded_jg": [[LIGANDS[j], g] for j, g in excl_jg],
        "r2_fold": round(r2_fold, 4), "r2_abs1": round(r2_abs1, 4),
        "slope_fold": round(b_fold, 4), "delta_R2": round(dR2, 4),
        "ci95_delta": [round(float(ci[0]), 4), round(float(ci[1]), 4)],
        "verdict": verdict,
        "selfbuilt_arm": selfbuilt,
    }
    with open(os.path.join(outdir, "代码48_ruling.json"), "w", encoding="utf-8") as f:
        json.dump(verdict_json, f, ensure_ascii=False, indent=2)
    log(f"\nOutput: {outdir}/代码48_{{判定日志.txt, 单元table.csv, robust臂.csv, ruling.json}}")

if __name__ == "__main__":
    main()
