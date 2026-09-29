# -*- coding: utf-8 -*-
"""
代码44 · P2 盲裁决：幅度分层衰减预言（Msn2，Hansen & Zechner 2021）
====================================================================
预注册书：05_主线纲领与设计/预注册_P2_幅度分层衰减预言_v01.md（2026-08-14 注册）
数据集  ：Hansen & Zechner 2021 (MSB)，Zenodo doi:10.5281/zenodo.2755026

判定线（注册时已锁死，禁止事后修改）：
  P2-1  Δ_mol  = AUC_混合 − mean(AUC_尺寸三分层内) ≥ +0.03，≥6/7 野生型启动子
  P2-2  Δ_AU   ≤ 0.01 且逐启动子 Δ_AU < Δ_mol，≥6/7 野生型启动子
  P2-3  |Δ_time| ≤ 0.01（k=3 主分析；k∈{2,4} 稳健臂），≥6/7 野生型启动子
  P2-F  任一方向反转 → 对应条款记证伪，死活双录入账

预注册未细化、由本脚本在开数据前写定的实现口径（同属盲态锁定，登记于总账）：
  · 启动子级 Δ = 三个相邻剂量对（100↔275、275↔690、690nM↔3µM）Δ 的算术平均；
    21 个判定单元（3 对 × 7 启动子）逐对另报于单元表
  · 三分层边界取自该单元观测合并样本（两剂量合并），bootstrap 重抽样时
    细胞携带其原始层标签（边界不随重抽样重算）
  · 事件时间解码器基于 YFP(AU) 通道：逐细胞基线 μ₀+kσ₀（t≤0 帧），
    首个 t>0 越阈时刻；未越阈 = 右删失（记 +inf），AUC 按 Mann–Whitney
    平均秩处理（删失细胞彼此并列、且大于一切已越阈者）
  · 缺失判定：YFP 或 cell_size 任一缺失 >10% 帧（64 帧中 ≥7 帧）→ 剔除
  · bootstrap（2000 次，逐细胞重抽样）CI 只给三条主判定通道
    （mol_peak / au_peak / time_k3）在主分析时长上的 Δ；
    备择臂（AUC 解码器、k∈{2,4}、备择代理分层、其余时长）报点估计
  · 主分析时长 = 50 min（全部启动子均有；确定性规则：50 在则取之，否则取最长）

运行方式（用户机，Windows + Spyder）：
  1) 解压 HansenZechner_RawData.zip
  2) 把下面 ROOT 改成解压后 2018_DataForPaper 目录的路径
  3) Spyder 中 %runfile 本文件；纯 CPU，约 5–15 分钟
依赖：numpy、scipy。随机种子 20260814（锁）。
输出：代码44_判定日志.txt、代码44_单元表.csv、代码44_裁决.json（写于本脚本所在目录）
"""
import os
import re
import sys
import json
import time as _time
import numpy as np
from scipy.io import loadmat
from scipy.stats import rankdata

# ---------------- 用户配置区 ----------------
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


# ---------------- 数据读取与逐细胞特征 ----------------
def cond_path(prom, series, dur, dose_tag):
    return os.path.join(ROOT, prom, f"{prom}_{series}_{dur}_{dose_tag}_size.mat")


def extract(fpath):
    """读取一个 _size.mat，返回逐细胞特征 dict + 剔除统计。全部确定性规则。"""
    m = loadmat(fpath)
    t = m["time"].ravel().astype(float)
    Y = m["YFP"].astype(float)            # AU：尺寸归一化浓度
    M = m["YFP_molecules"].astype(float)  # 绝对丰度（未归一化）
    S = m["cell_size_pixels"].astype(float)
    n, nF = Y.shape
    b_idx = np.where(t <= 0)[0]
    w_idx = np.where(t >= 0)[0]
    a_idx = np.where(t > 0)[0]
    dt = float(np.median(np.diff(t)))

    # --- 剔除规则（预注册 §2.2 + 本脚本锁定口径） ---
    miss_y = np.isnan(Y).sum(1)
    miss_s = np.isnan(S).sum(1)
    keep = (miss_y <= int(0.10 * nF)) & (miss_s <= int(0.10 * nF))  # 任一通道缺失 >10% 帧剔除
    dS = np.abs(np.diff(S, axis=1))
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = dS / np.abs(S[:, :-1])
    rel[~np.isfinite(rel)] = np.inf                 # 前帧为 0/NaN 视为跳变
    jump = rel.max(1)
    keep &= ~(jump > 0.40)                          # 相邻帧尺寸相对跳变 >40%

    mu0 = np.nanmean(Y[:, b_idx], axis=1)
    sd0 = np.nanstd(Y[:, b_idx], axis=1, ddof=1)
    keep &= np.isfinite(mu0) & np.isfinite(sd0) & (mu0 > 0)   # 校正后基线 AU ≤ 0 剔除

    # --- 逐细胞特征 ---
    mean_size = np.nanmean(S, axis=1)
    base_mol = np.nanmean(M[:, b_idx], axis=1)      # 备择分层代理
    peak_mol = np.nanmax(M[:, w_idx], axis=1)
    peak_au = np.nanmax(Y[:, w_idx], axis=1)
    auc_mol = np.nansum(M[:, w_idx], axis=1) * dt   # 备择幅度解码器
    auc_au = np.nansum(Y[:, w_idx], axis=1) * dt

    def cross(k):
        thr = (mu0 + k * sd0)[:, None]
        above = Y[:, a_idx] > thr                   # NaN 比较为 False，天然安全
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


# ---------------- 统计量：Mann–Whitney AUC 与分层罚 Δ ----------------
def auc_mw(lo, hi):
    """AUC = P(hi > lo) + 0.5·P(tie)，平均秩处理并列；NaN 丢弃，+inf（删失）保留。"""
    lo = lo[~np.isnan(lo)]
    hi = hi[~np.isnan(hi)]
    nL, nH = len(lo), len(hi)
    if nL < 2 or nH < 2:
        return np.nan
    r = rankdata(np.concatenate([lo, hi]))
    return float((r[nL:].sum() - nH * (nH + 1) / 2.0) / (nL * nH))


def tercile_labels(x):
    """按合并样本 1/3、2/3 分位数分三层，返回 0/1/2 标签。"""
    q1, q2 = np.quantile(x, [1.0 / 3.0, 2.0 / 3.0])
    return (x > q1).astype(int) + (x > q2).astype(int)


def delta_strat(v_lo, v_hi, l_lo, l_hi):
    """Δ = AUC_混合 − mean(AUC_层内)；层内某层任一组 <2 个有效细胞则跳过该层（<2 层可用 → NaN）。"""
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
    """逐细胞 bootstrap 2000 次；细胞携带原层标签（锁定口径）。"""
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


# ---------------- 主分析管线 ----------------
CHANNELS_MAIN = ["mol_peak", "au_peak", "time_k3"]
CHANNELS_ROB = ["time_k2", "time_k4", "mol_auc", "au_auc"]


def run_condition(prom, series, dur, dose):
    p = cond_path(prom, series, dur, DOSE_TAG[dose])
    if not os.path.isfile(p):
        return None
    return extract(p)


def promoter_delta(rows, prom, ch, strat="size", dur=MAIN_DUR):
    """启动子级 Δ = 三剂量对 Δ 的算术平均；任一对不可用 → NaN（按不通过计，保守、确定）。"""
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


# ---------------- 描述性臂（不进判定） ----------------
def fm_features(prom, series, dur):
    """FM/FM4 条件：越阈计数（μ₀+3σ₀，AU 通道）、首次越阈时间、峰值。"""
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
        ("FM", "8_5min", "FM4", "15minINT", "FM8(8脉冲) vs FM4-15min(4脉冲)"),
        ("FM", "4_5min", "FM4", "15minINT", "FM4(4脉冲) vs FM4-15min(4脉冲,同计数异间隔)"),
        ("FM", "4_5min", "FM4", "75minINT", "FM4(4脉冲) vs FM4-75min(4脉冲,同计数异间隔)"),
    ]
    for s1, d1, s2, d2, label in contrasts:
        f1 = fm_features(prom, s1, d1)
        f2 = fm_features(prom, s2, d2)
        if f1 is None or f2 is None:
            log(f"    [FM 描述] {prom} {label}：缺文件，跳过")
            continue
        a_cnt = auc_mw(f2["n_cross"], f1["n_cross"])
        a_tim = auc_mw(f1["first_time"], f2["first_time"])
        a_pk = auc_mw(f2["peak_au"], f1["peak_au"])
        log(f"    [FM 描述] {prom} {label}: "
            f"计数AUC={a_cnt:.3f} 首时AUC={a_tim:.3f} 峰值AUC={a_pk:.3f} "
            f"(n={len(f1['n_cross'])}/{len(f2['n_cross'])})")


def cfp_probe(prom):
    """CFP 通道作增益探针的可行性：基线 CFP 与时均尺寸的相关（描述性）。"""
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
        log(f"    [CFP 探针] {prom}: 有效细胞不足")
        return
    r = float(np.corrcoef(cb[ok], ms[ok])[0, 1])
    log(f"    [CFP 探针] {prom}: corr(基线CFP, 尺寸)={r:+.3f} (n={int(ok.sum())})")


# ---------------- 主入口 ----------------
def main():
    t_start = _time.time()
    global LOG_PATH, _log_fh
    outdir = os.path.dirname(os.path.abspath(__file__)) if "__file__" in globals() else os.getcwd()
    LOG_PATH = os.path.join(outdir, "代码44_判定日志.txt")
    _log_fh = open(LOG_PATH, "w", encoding="utf-8")

    rng = np.random.default_rng(SEED)
    log("=" * 72)
    log("代码44 · P2 盲裁决：幅度分层衰减预言")
    log(f"预注册：2026-08-14 v0.1 ｜ 种子 {SEED} ｜ bootstrap {N_BOOT}")
    log(f"数据根目录：{ROOT}")
    log(f"判定线：P2-1 Δ_mol≥+{MARGIN_P21}；P2-2 Δ_AU≤{MARGIN_P22} 且 <Δ_mol；"
        f"P2-3 |Δ_time|≤{MARGIN_P23}（k={K_MAIN}）；各 ≥6/7 启动子")
    log("=" * 72)

    rows = []
    for prom in WT + MUT:
        for dur in DURS:
            conds = {}
            ok = True
            for dose in (100, 275, 690, 3000):
                r = run_condition(prom, "DM", dur, dose)
                if r is None:
                    log(f"[警告] 缺文件：{prom} DM {dur} {DOSE_TAG[dose]}，跳过 {prom} {dur}")
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
            log(f"[完成] {prom} DM {dur}")

    # ---- v0.1.1 护栏：零数据加载 → 中止，不打印空裁决表（空跑不算裁决） ----
    if not rows:
        log("")
        log("!" * 72)
        log("中止：未加载到任何条件数据。请检查 ROOT 是否指向【直接包含 ALD3 等")
        log("启动子文件夹】的那一层目录。本次运行为空跑，不构成裁决。")
        log("!" * 72)
        _log_fh.close()
        return

    # ---------------- 裁决（主分析：50min、尺寸分层、主解码器） ----------------
    log("")
    log("=" * 72)
    log(f"裁决（主分析时长 {MAIN_DUR}，尺寸三分层）")
    log("-" * 72)
    log(f"{'启动子':<8s} {'Δ_mol':>8s} {'Δ_AU':>8s} {'Δ_time':>8s} "
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
            f"{'✓' if p1 else '✗':>6s} {'✓' if p2 else '✗':>6s} {'✓' if p3 else '✗':>6s}")
    c1, c2, c3 = n1 >= 6, n2 >= 6, n3 >= 6
    verdict["clauses"] = {"P2-1": {"pass_promoters": n1, "hit": bool(c1)},
                          "P2-2": {"pass_promoters": n2, "hit": bool(c2)},
                          "P2-3": {"pass_promoters": n3, "hit": bool(c3)}}
    overall = bool(c1 and c2 and c3)
    verdict["overall"] = overall
    log("-" * 72)
    log(f"P2-1（幅度付费）  ：{n1}/7 启动子通过 → {'命中' if c1 else '证伪'}")
    log(f"P2-2（归一化赎回）：{n2}/7 启动子通过 → {'命中' if c2 else '证伪'}")
    log(f"P2-3（事件时间免费）：{n3}/7 启动子通过 → {'命中' if c3 else '证伪'}")
    log(f"P2-F 总评：{'预言成立（三条款全命中）' if overall else '对应条款证伪，如实双录入账'}")
    log("")

    # ---------------- 21 判定单元明细（主通道） ----------------
    log("21 判定单元明细（主分析时长，尺寸分层）：")
    for ch in CHANNELS_MAIN:
        log(f"  [{ch}]")
        for r in rows:
            if (r["duration"] == MAIN_DUR and r["strat"] == "size"
                    and r["channel"] == ch and r["promoter"] in WT):
                ci = (f"[{r['ci_lo']:+.3f},{r['ci_hi']:+.3f}]"
                      if np.isfinite(r["ci_lo"]) else "[  点估计  ]")
                log(f"    {r['promoter']:<8s} {r['pair']:<9s} "
                    f"AUC混合={r['auc_pool']:.3f} 层均={r['auc_terc']:.3f} "
                    f"Δ={r['delta']:+.4f} {ci} n={r['n_lo']}/{r['n_hi']}")
    log("")

    # ---------------- 稳健臂与备择代理（点估计速览） ----------------
    log("稳健臂 k∈{2,4} 与其余时长（启动子级 Δ，点估计）：")
    for ch in ["time_k2", "time_k4"]:
        ds = {p: promoter_delta(rows, p, ch) for p in WT}
        log(f"  [{ch}] " + " ".join(
            f"{p}:{ds[p]:+.3f}" if np.isfinite(ds[p]) else f"{p}:NaN" for p in WT))
    for dur in DURS:
        if dur == MAIN_DUR:
            continue
        dm = {p: promoter_delta(rows, p, "mol_peak", dur=dur) for p in WT}
        log(f"  [时长 {dur} Δ_mol] " + " ".join(
            f"{p}:{dm[p]:+.3f}" if np.isfinite(dm[p]) else f"{p}:NaN" for p in WT))
    da_alt = {p: promoter_delta(rows, p, "mol_peak", strat="baseMol") for p in WT}
    log("  [备择代理分层 Δ_mol] " + " ".join(
        f"{p}:{da_alt[p]:+.3f}" if np.isfinite(da_alt[p]) else f"{p}:NaN" for p in WT))
    log("")

    # ---------------- 突变体（描述性） ----------------
    log("突变体（描述性，不进判定）：")
    for prom in MUT:
        dm = promoter_delta(rows, prom, "mol_peak")
        da = promoter_delta(rows, prom, "au_peak")
        dt3 = promoter_delta(rows, prom, "time_k3")
        f = lambda x: f"{x:+.4f}" if np.isfinite(x) else "NaN"
        log(f"  {prom}: Δ_mol={f(dm)} Δ_AU={f(da)} Δ_time={f(dt3)}")
    log("")

    # ---------------- 次要臂：FM/FM4 计数解码 + CFP 探针（描述性） ----------------
    log("次要臂（描述性，不进判定）：FM/FM4 解码")
    for prom in WT:
        fm_descriptive(prom)
    log("次要臂（描述性）：CFP 增益探针")
    for prom in WT:
        cfp_probe(prom)

    # ---------------- 落盘 ----------------
    import csv
    csv_path = os.path.join(outdir, "代码44_单元表.csv")
    with open(csv_path, "w", newline="", encoding="utf-8-sig") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        for r in rows:
            w.writerow(r)
    json_path = os.path.join(outdir, "代码44_裁决.json")
    with open(json_path, "w", encoding="utf-8") as fh:
        json.dump(sanitize(verdict), fh, ensure_ascii=False, indent=2)
    log("")
    log(f"单元表 → {csv_path}")
    log(f"裁决   → {json_path}")
    log(f"日志   → {LOG_PATH}")
    log(f"总耗时 {_time.time() - t_start:.0f} s")
    _log_fh.close()


if __name__ == "__main__":
    main()
