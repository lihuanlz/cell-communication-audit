#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Code 49b · P7 erratum re-adjudication: FCD/Weber adjudication — chemotaxis single-cell FRET same-ligand background steps
pre-registration: 预注册_P7_FCD-Weber判决_趋化单细胞FRET_v01.md (registry No. 97, frozen)
Data: Moore et al. 2024 Cell Syst 15:628; Dryad doi:10.5061/dryad.nvx0k6dzz (CC0)
v0.2.0 · 2026-08-15 (erratum re-adjudication: B taken from metadata, within-row segmentation switched to positional criterion; spec = pre-registration §10 v0.2 frozen) · seed 20260815 · bootstrap 2000
Criterion (§3, frozen): over 30 condition units R(B,F)
  ΔR2 = R2(R~log10 r) - R2(R~log10 T)，r=(B+F)/B，T=B+F
  >=+0.15 and fold slope>0 → P7 hit; <= -0.15 → falsified; in between → intermediate
Guardrails (§5): separability |ρ(log r,log T)|>=0.95 abort; unit cell count <10 removed;
  <24 power-limited; <18 abort; response direction reversal >50% abort (empty run).
run: python3 代码49.py [DATA_DIR]
"""
import os, sys, json, datetime
import numpy as np

SEED = 20260815
N_BOOT = 2000
DELTA_LINE = 0.15
MIN_CELL_UNIT = 10
MIN_UNITS_PWR = 24
MIN_UNITS_ABORT = 18
SEP_ABORT = 0.95
MIN_LEVELS_CELL = 3
COLLAPSE_DLOGR = 0.15
COLLAPSE_BRATIO = 10.0
COLLAPSE_SPREAD = 0.10

FILES = {  # filename (minus .mat) -> BackgroundLIst background µM (§1; final per-row s measured recheck)
 "210802_FOV1":0,"210802_FOV2":0,"210805_FOV1":0,"210805_FOV2":0,"220106_FOV1":0,"230417_FOV1":0,
 "230815_FOV1":0.01,"230815_FOV2":0.01,"230816_FOV1":0.01,"230816_FOV2":0.01,
 "230830_FOV1":0.1,"230830_FOV2":0.1,"230831_FOV1":0.1,"230831_FOV2":0.1,  # FOV2 filed exception (s=100-180)
 "220615_FOV1":0.3,"230410_FOV1":0.3,
 "230428_FOV1":1.0,"230429_FOV1":1.0,
 "220302_FOV1":10.0,"220303_FOV1":10.0,
 "210816_FOV1":100.0,"210816_FOV2":100.0,"230717_FOV1":100.0,"230718_FOV1":100.0,
}

LOG_PATH = None
def log(msg=""):
    print(msg)
    if LOG_PATH:
        with open(LOG_PATH, "a", encoding="utf-8") as f:
            f.write(str(msg) + "\n")

def linfit(x, y):
    X = np.column_stack([np.ones_like(x), x])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    r = y - X @ beta
    sst = float(((y - y.mean()) ** 2).sum())
    r2 = 1 - float(r @ r) / sst if sst > 0 else np.nan
    return r2, float(beta[1])

def main():
    global LOG_PATH
    root = sys.argv[1] if len(sys.argv) > 1 else os.path.dirname(os.path.abspath(__file__))
    outdir = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
    LOG_PATH = os.path.join(outdir, "代码49b_判定日志.txt")
    if os.path.exists(LOG_PATH):
        os.remove(LOG_PATH)
    log("=" * 74)
    log("Code 49b · P7 erratum re-adjudication: FCD/Weber adjudication — chemotaxis single-cell FRET background steps  v0.1.0")
    log(f"run time: {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
    log(f"seed={SEED} bootstrap={N_BOOT} criterion ΔR2=±{DELTA_LINE}")
    log("pre-registration: 预注册_P7_FCD-Weber判决_趋化单细胞FRET_v01.md (registry No. 97, frozen)")
    log(f"data root directory: {root}")
    log("=" * 74)

    import scipy.io as sio
    # ---------- load + §5-3 schema guardrail + row-level s measured recheck ----------
    cell_rec = []   # (B_file_measured, F, Δa_i, ΔA_i, cell_key)
    file_log = []
    for name, bg_list in FILES.items():
        p = os.path.join(root, name + ".mat")
        if not os.path.isfile(p):
            log(f"[missing file] {name}.mat — removed and recorded"); file_log.append((name, "missing")); continue
        try:
            m = sio.loadmat(p)
            rd = m["reorgData"]["resp_data"][0, 0]
        except Exception as e:
            log(f"[guardrail] {name} schema exception ({str(e)[:50]}) — removed and recorded"); file_log.append((name, "schema")); continue
        ncells = rd.shape[1]
        # row-level s measured: baseline = first-column median; file background B taken as measured (§1 filing)
        # §10 v0.2: B taken from BackgroundLIst metadata; unique exception 230831_FOV2 → B=100
        B_file = 100.0 if name == "230831_FOV2" else float(bg_list)
        log(f"[loaded] {name}: cell {ncells}, background B={B_file:g} µM (metadata basis)")
        file_log.append((name, f"ok n={ncells} B={B_file:g}"))
        for ci in range(ncells):
            a = rd["a"][0, ci].astype(float)
            A = rd["A"][0, ci].astype(float)
            s = rd["s"][0, ci].astype(float)
            if a.shape != (35, 20) or s.shape != (35, 20):
                log(f"  [guardrail] {name} cell{ci} shape {a.shape}≠(35,20) — whole cell removed"); continue
            Bc = B_file
            lev_resp, lev_respA = {}, {}
            for row in range(35):
                sr = s[row]
                stim_val = float(sr.max())
                stim = np.where(sr >= stim_val - 1e-9)[0]
                if len(stim) == 0:
                    continue
                first_stim = stim[0]
                pre = np.arange(0, first_stim)   # positional criterion: all columns before the first stimulus column are baseline
                if len(pre) < 4 or len(stim) < 6:
                    continue
                pre_c, stim_c = pre[-4:], stim[-6:]
                da = float(np.median(a[row, pre_c])) - float(np.median(a[row, stim_c]))
                dA = float(np.median(A[row, pre_c])) - float(np.median(A[row, stim_c]))
                F = stim_val - Bc
                if F <= 0:
                    continue
                key = round(stim_val, 4)
                lev_resp.setdefault(key, []).append(da)
                lev_respA.setdefault(key, []).append(dA)
            good = {k for k, v in lev_resp.items() if len(v) >= 4}
            if len(good) < MIN_LEVELS_CELL:
                continue
            for k in sorted(good):
                da_i = float(np.median(lev_resp[k]))
                dA_i = float(np.median(lev_respA[k]))
                if not (np.isfinite(da_i) and np.isfinite(dA_i)):
                    continue
                cell_rec.append((Bc, k - Bc, da_i, dA_i, f"{name}#{ci}"))
    log(f"\n[total] cell×level records: {len(cell_rec)}")

    # ---------- condition units ----------
    units = {}
    for Bc, F, da, dA, ck in cell_rec:
        Bq = round(Bc, 4); Fq = round(F, 4)
        units.setdefault((Bq, Fq), []).append((da, dA))
    log(f"[units] raw (B,F) combinations: {len(units)} (incl. B=0 descriptive group)")

    rows = []
    for (B, F), lst in sorted(units.items()):
        das = np.array([x[0] for x in lst]); dAs = np.array([x[1] for x in lst])
        rows.append(dict(B=B, F=F, n=len(das), R=float(np.median(das)), RA=float(np.median(dAs))))
    main_units = [r for r in rows if r["B"] > 0 and r["n"] >= MIN_CELL_UNIT]
    dropped = [r for r in rows if r["B"] > 0 and r["n"] < MIN_CELL_UNIT]
    if dropped:
        log(f"[removed] units with cell count<{MIN_CELL_UNIT}: {len(dropped)}: " +
            str([(r["B"], r["F"], r["n"]) for r in dropped]))
    log(f"[main adjudication units] {len(main_units)}/30 (B>0)")
    if len(main_units) < MIN_UNITS_ABORT:
        log("[abort] valid units <18 — empty run"); sys.exit(2)
    power_note = len(main_units) < MIN_UNITS_PWR
    if power_note:
        log(f"[limited] valid units {len(main_units)} < {MIN_UNITS_PWR} — reduced-power limitation on record")

    neg = sum(1 for r in main_units if r["R"] < 0)
    if neg > len(main_units) / 2:
        log(f"[abort] response-direction-reversed units {neg}/{len(main_units)} >50% — check sign convention, empty run"); sys.exit(2)

    r_lg = np.array([np.log10((r["B"] + r["F"]) / r["B"]) for r in main_units])
    T_lg = np.array([np.log10(r["B"] + r["F"]) for r in main_units])
    R = np.array([r["R"] for r in main_units])
    RA = np.array([r["RA"] for r in main_units])

    # §5-1 separability
    rho_sep = float(np.corrcoef(r_lg, T_lg)[0, 1])
    log(f"[guardrail] separability Pearson(log r, log T) = {rho_sep:+.4f} (abort if |ρ|≥{SEP_ABORT})")
    if abs(rho_sep) >= SEP_ABORT:
        log("[abort] design cannot distinguish fold from absolute total — empty run"); sys.exit(2)

    r2_fold, b_fold = linfit(r_lg, R)
    r2_abs, b_abs = linfit(T_lg, R)
    dR2 = r2_fold - r2_abs

    log("\n" + "=" * 74)
    log("ruling (main adjudication: a feature, 30 condition units R(B,F))")
    log("-" * 74)
    for r in main_units:
        log(f"  B={r['B']:>7g} F={r['F']:>7g}  r={(r['B']+r['F'])/r['B']:>8.2f}  T={r['B']+r['F']:>8.3g}  n={r['n']:>4d}  R={r['R']:+.4f}")
    log(f"\n  R2(R~log fold) = {r2_fold:.4f} (slope {b_fold:+.4f})")
    log(f"  R2(R~log T)    = {r2_abs:.4f} (slope {b_abs:+.4f})")
    log(f"  ΔR2 = {dR2:+.4f} (criterion ±{DELTA_LINE})")
    sign_ok = b_fold > 0
    if dR2 >= DELTA_LINE and sign_ok:
        verdict = "P7 hit（fold/Weber 载波）"
    elif dR2 >= DELTA_LINE and not sign_ok:
        verdict = "P7 中间态（ΔR2 达标但 fold slopereverse，按pre-registration降级）"
    elif dR2 <= -DELTA_LINE:
        verdict = "P7 证伪（绝对总量载波）"
    else:
        verdict = "P7 中间态"
    log(f"  ruling：{verdict}")
    log("=" * 74)

    # bootstrap CI (descriptive)
    rng = np.random.default_rng(SEED)
    nU = len(main_units); boots = []
    for _ in range(N_BOOT):
        idx = rng.integers(0, nU, nU)
        r2f, _ = linfit(r_lg[idx], R[idx]); r2a, _ = linfit(T_lg[idx], R[idx])
        boots.append(r2f - r2a)
    ci = np.percentile(boots, [2.5, 97.5])
    log(f"[CI] ΔR2 bootstrap 95% = [{ci[0]:+.4f}, {ci[1]:+.4f}] (descriptive)")

    # ---------- §4 robust arm ----------
    r2f_A, bf_A = linfit(r_lg, RA); r2a_A, _ = linfit(T_lg, RA)
    log(f"[robust①] A feature: R2fold={r2f_A:.4f} R2abs={r2a_A:.4f} ΔR2={r2f_A-r2a_A:+.4f}")
    # ② cell-level + B dummies
    rec = [x for x in cell_rec if round(x[0],4) > 0]
    Bl = np.array([x[0] for x in rec]); Fl = np.array([x[1] for x in rec])
    yl = np.array([x[2] for x in rec])
    rl = np.log10((Bl + Fl) / Bl); Tl = np.log10(Bl + Fl)
    Bu = sorted(set(np.round(Bl, 4)))
    Dum = np.column_stack([(np.round(Bl, 4) == b).astype(float) for b in Bu[1:]])
    X1 = np.column_stack([np.ones_like(rl), rl, Dum])
    X2 = np.column_stack([np.ones_like(rl), Tl, Dum])
    def r2m(X, y):
        beta, *_ = np.linalg.lstsq(X, y, rcond=None)
        res = y - X @ beta
        return 1 - float(res @ res) / float(((y - y.mean()) ** 2).sum())
    r2f_c, r2a_c = r2m(X1, yl), r2m(X2, yl)
    log(f"[robust②] cell-level+B dummies: R2fold={r2f_c:.4f} R2abs={r2a_c:.4f} ΔR2={r2f_c-r2a_c:+.4f} (n={len(yl)})")
    # ③ matched-fold group collapse
    groups = []
    used = set()
    for i_, r1 in enumerate(main_units):
        if i_ in used: continue
        rr = (r1["B"] + r1["F"]) / r1["B"]
        mates = [j_ for j_, r2_ in enumerate(main_units)
                 if abs(np.log10((r2_["B"]+r2_["F"])/r2_["B"]) - np.log10(rr)) <= COLLAPSE_DLOGR]
        Bs = [main_units[j_]["B"] for j_ in mates]
        if len(set(Bs)) >= 2 and max(Bs)/min(Bs) >= COLLAPSE_BRATIO:
            groups.append(mates); used.update(mates)
    log(f"[robust③] matched-fold groups: {len(groups)} (r±{COLLAPSE_DLOGR} log10 and B spanning ≥{COLLAPSE_BRATIO:g}×):")
    ncol = 0
    for mates in groups:
        Rs = [main_units[j_]["R"] for j_ in mates]
        spread = max(Rs) - min(Rs)
        ok = spread <= COLLAPSE_SPREAD
        ncol += ok
        r0 = (main_units[mates[0]]["B"] + main_units[mates[0]]["F"]) / main_units[mates[0]]["B"]
        log(f"   r≈{r0:.2f}: Bs={[main_units[j_]['B'] for j_ in mates]} Rs={[round(x,3) for x in Rs]} spread={spread:.3f} {'collapsed' if ok else 'not collapsed'}")
    if groups:
        log(f"   collapse ratio {ncol}/{len(groups)} (≥2/3 supports fold, descriptive)")
    # ⑤ removed B=100
    keep = [k for k, r in enumerate(main_units) if r["B"] < 100]
    if len(keep) >= MIN_UNITS_ABORT:
        r2f5, _ = linfit(r_lg[keep], R[keep]); r2a5, _ = linfit(T_lg[keep], R[keep])
        log(f"[robust⑤] removed B=100：n={len(keep)} ΔR2={r2f5-r2a5:+.4f}")

    # ---------- §6 descriptive arm ----------
    log("\n[descriptive arm] B=0 group absolute dose-response:")
    for r in rows:
        if r["B"] == 0 and r["n"] >= MIN_CELL_UNIT:
            log(f"  F={r['F']:>5g}  n={r['n']:>4d}  R={r['R']:+.4f}")
    F_lg_only = np.array([np.log10(r["F"]) for r in main_units])
    r2_inc, b_inc = linfit(F_lg_only, R)
    log(f"[descriptive arm] third model R~log10(F) (pure foreground): R2={r2_inc:.4f} (slope {b_inc:+.4f})")

    # ---------- output ----------
    import csv
    with open(os.path.join(outdir, "代码49b_单元table.csv"), "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["B_uM", "F_uM", "fold_r", "T_uM", "log10r", "log10T", "n_cells", "R_a", "R_A"])
        for r in main_units:
            rr = (r["B"] + r["F"]) / r["B"]
            w.writerow([r["B"], r["F"], f"{rr:.4f}", r["B"] + r["F"],
                        f"{np.log10(rr):.4f}", f"{np.log10(r['B']+r['F']):.4f}",
                        r["n"], f"{r['R']:.6f}", f"{r['RA']:.6f}"])
    vjson = {
        "code": "代码49b v0.2.0", "prereg": "P7 v0.1（注册·九十七）", "seed": SEED,
        "n_units": len(main_units), "power_note": bool(power_note),
        "sep_rho": round(rho_sep, 4),
        "r2_fold": round(r2_fold, 4), "r2_abs": round(r2_abs, 4),
        "slope_fold": round(b_fold, 4), "delta_R2": round(dR2, 4),
        "ci95_delta": [round(float(ci[0]), 4), round(float(ci[1]), 4)],
        "verdict": verdict,
        "robust": {"A_feature_dR2": round(r2f_A - r2a_A, 4),
                   "celllevel_dR2": round(r2f_c - r2a_c, 4),
                   "collapse_groups": len(groups), "collapse_ok": ncol},
        "file_log": file_log,
    }
    with open(os.path.join(outdir, "代码49b_ruling.json"), "w", encoding="utf-8") as f:
        json.dump(vjson, f, ensure_ascii=False, indent=2)
    log(f"\noutput: {outdir}/代码49_{{判定日志.txt, 单元table.csv, ruling.json}}")

if __name__ == "__main__":
    main()
