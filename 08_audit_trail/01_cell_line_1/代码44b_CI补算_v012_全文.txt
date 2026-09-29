# -*- coding: utf-8 -*-
"""
代码44b · P2 附属补算：bootstrap CI（v0.1.2）
================================================
性质：附属补算脚本。**判决已冻结于注册·七十二（2026-08-14），本脚本不含、
不打印、不改变任何判决逻辑。** 唯一任务：补交付预注册 v0.1 承诺的
bootstrap 2000 次（逐细胞重抽样）置信区间——该交付物因代码44 的 CI 触发
条件笔误（strat_key/strat_tag 混用）未能产出。

范围（与预注册一致）：7 个野生型启动子 × 3 剂量对 × 3 主通道
（mol_peak / au_peak / time_k3），主分析时长 50min，尺寸三分层，共 63 行。
种子仍为 20260814（锁）。

复现一致性检查：本脚本在写 CI 前，先重算 63 行的点估计并与同目录下
代码44 的《代码44_单元表.csv》（正式裁决运行的落盘）逐格比对（容差 1e-12）。
不一致则大声报警——那意味着补算环境偏离了正式运行，CI 不可入账。
若同目录无该 CSV，跳过比对并注明（CI 仍算，但须人工核对点估计）。

输出：代码44b_日志.txt、代码44b_CI单元表.csv（写于本脚本所在目录）。
运行：同代码44（ROOT 指到直接含 ALD3 等启动子文件夹的那层）；纯 CPU 约 1–3 分钟。
"""
import os
import sys
import csv
import time as _time
import numpy as np
from scipy.io import loadmat
from scipy.stats import rankdata

# ---------------- 用户配置区 ----------------
ROOT = r"C:/Users/lihua/Desktop/NC/data/HansenZechner2021"
if len(sys.argv) > 1 and os.path.isdir(sys.argv[1]):
    ROOT = sys.argv[1]
# ---------------------------------------------

SEED = 20260814
N_BOOT = 2000
WT = ["ALD3", "DCS2", "DDR2", "HXK1", "RTN2", "SIP18", "TKL2"]
DOSE_TAG = {100: "100nM", 275: "275nM", 690: "690nM", 3000: "3uM"}
PAIRS = [(100, 275), (275, 690), (690, 3000)]
MAIN_DUR = "50min"
K_MAIN = 3
CHANNELS_MAIN = ["mol_peak", "au_peak", "time_k3"]

_log_fh = None


def log(msg=""):
    print(msg)
    if _log_fh is not None:
        _log_fh.write(str(msg) + "\n")
        _log_fh.flush()


# ---------------- 以下函数与代码44 完全一致（判决逻辑所在，逐字保留） ----------------
def cond_path(prom, series, dur, dose_tag):
    return os.path.join(ROOT, prom, f"{prom}_{series}_{dur}_{dose_tag}_size.mat")


def extract(fpath):
    m = loadmat(fpath)
    t = m["time"].ravel().astype(float)
    Y = m["YFP"].astype(float)
    M = m["YFP_molecules"].astype(float)
    S = m["cell_size_pixels"].astype(float)
    n, nF = Y.shape
    b_idx = np.where(t <= 0)[0]
    w_idx = np.where(t >= 0)[0]
    a_idx = np.where(t > 0)[0]
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
    mean_size = np.nanmean(S, axis=1)
    peak_mol = np.nanmax(M[:, w_idx], axis=1)
    peak_au = np.nanmax(Y[:, w_idx], axis=1)

    def cross(k):
        thr = (mu0 + k * sd0)[:, None]
        above = Y[:, a_idx] > thr
        has = above.any(1)
        first = above.argmax(1)
        return np.where(has, t[a_idx][first], np.inf)

    feat = {"mean_size": mean_size, "mol_peak": peak_mol, "au_peak": peak_au,
            "time_k3": cross(K_MAIN)}
    return {k: v[keep] for k, v in feat.items()}


def auc_mw(lo, hi):
    lo = lo[~np.isnan(lo)]
    hi = hi[~np.isnan(hi)]
    nL, nH = len(lo), len(hi)
    if nL < 2 or nH < 2:
        return np.nan
    r = rankdata(np.concatenate([lo, hi]))
    return float((r[nL:].sum() - nH * (nH + 1) / 2.0) / (nL * nH))


def tercile_labels(x):
    q1, q2 = np.quantile(x, [1.0 / 3.0, 2.0 / 3.0])
    return (x > q1).astype(int) + (x > q2).astype(int)


def delta_strat(v_lo, v_hi, l_lo, l_hi):
    ap = auc_mw(v_lo, v_hi)
    ts = []
    for g in range(3):
        a = auc_mw(v_lo[l_lo == g], v_hi[l_hi == g])
        if not np.isnan(a):
            ts.append(a)
    if np.isnan(ap) or len(ts) < 2:
        return np.nan, ap, np.nan
    return ap - float(np.mean(ts)), ap, float(np.mean(ts))


# ---------------- 唯一的修复：bootstrap 现在真正执行 ----------------
def boot_delta(v_lo, v_hi, l_lo, l_hi, rng):
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


# ---------------- 主入口 ----------------
def main():
    t_start = _time.time()
    global _log_fh
    outdir = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
    _log_fh = open(os.path.join(outdir, "代码44b_日志.txt"), "w", encoding="utf-8")

    rng = np.random.default_rng(SEED)
    log("=" * 72)
    log("代码44b · P2 附属补算：bootstrap CI（v0.1.2）")
    log("判决已冻结于注册·七十二；本脚本不打印、不改变任何判决。")
    log(f"种子 {SEED} ｜ bootstrap {N_BOOT} ｜ 范围 7启动子×3剂量对×3主通道 @ {MAIN_DUR} 尺寸分层")
    log(f"数据根目录：{ROOT}")
    log("=" * 72)

    rows = []
    for prom in WT:
        conds = {}
        ok = True
        for dose in (100, 275, 690, 3000):
            p = cond_path(prom, "DM", MAIN_DUR, DOSE_TAG[dose])
            if not os.path.isfile(p):
                log(f"[警告] 缺文件：{prom} {MAIN_DUR} {DOSE_TAG[dose]}")
                ok = False
                break
            conds[dose] = extract(p)
        if not ok:
            continue
        for dlo, dhi in PAIRS:
            flo, fhi = conds[dlo], conds[dhi]
            pool = np.concatenate([flo["mean_size"], fhi["mean_size"]])
            labels = tercile_labels(pool)
            nL = len(flo["mean_size"])
            l_lo, l_hi = labels[:nL], labels[nL:]
            for ch in CHANNELS_MAIN:
                d, ap, mt = delta_strat(flo[ch], fhi[ch], l_lo, l_hi)
                ci_lo, ci_hi = boot_delta(flo[ch], fhi[ch], l_lo, l_hi, rng)
                rows.append({"promoter": prom, "pair": f"{dlo}-{dhi}", "channel": ch,
                             "auc_pool": ap, "auc_terc": mt, "delta": d,
                             "ci_lo": ci_lo, "ci_hi": ci_hi,
                             "n_lo": nL, "n_hi": len(fhi[ch])})
        log(f"[完成] {prom}")

    if not rows:
        log("!" * 72)
        log("中止：未加载到任何数据，检查 ROOT。")
        log("!" * 72)
        _log_fh.close()
        return

    # ---------------- 复现一致性检查（对正式运行的单元表） ----------------
    ref_path = os.path.join(outdir, "代码44_单元表.csv")
    if os.path.isfile(ref_path):
        ref = {}
        with open(ref_path, encoding="utf-8-sig") as fh:
            for r in csv.DictReader(fh):
                if r["duration"] == MAIN_DUR and r["strat"] == "size" and r["channel"] in CHANNELS_MAIN:
                    ref[(r["promoter"], r["pair"], r["channel"])] = (
                        float(r["auc_pool"]), float(r["auc_terc"]), float(r["delta"]))
        maxdiff = 0.0
        nchk = 0
        for r in rows:
            key = (r["promoter"], r["pair"], r["channel"])
            if key in ref:
                nchk += 1
                for a, b in zip((r["auc_pool"], r["auc_terc"], r["delta"]), ref[key]):
                    maxdiff = max(maxdiff, abs(a - b))
        verdict = "通过（点估计与正式运行逐格一致）" if (nchk == 63 and maxdiff < 1e-12) \
            else f"**不一致**（比对 {nchk}/63 行，最大差 {maxdiff:.2e}）——CI 不可入账，需排查"
        log("")
        log(f"复现一致性检查：比对 {nchk}/63 行，最大偏差 {maxdiff:.2e} → {verdict}")
    else:
        log("")
        log("复现一致性检查：同目录未找到《代码44_单元表.csv》，跳过比对（请人工核对点估计）。")

    # ---------------- 落盘与打印 ----------------
    csv_path = os.path.join(outdir, "代码44b_CI单元表.csv")
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    log("")
    log("63 单元 CI 明细（Δ 的 bootstrap 2000 次 95% CI）：")
    for ch in CHANNELS_MAIN:
        log(f"  [{ch}]")
        for r in rows:
            if r["channel"] == ch:
                log(f"    {r['promoter']:<8s} {r['pair']:<9s} Δ={r['delta']:+.4f} "
                    f"CI=[{r['ci_lo']:+.4f}, {r['ci_hi']:+.4f}] n={r['n_lo']}/{r['n_hi']}")
    log("")
    log(f"单元表 → {csv_path}")
    log(f"总耗时 {_time.time() - t_start:.0f} s")
    _log_fh.close()


if __name__ == "__main__":
    main()
