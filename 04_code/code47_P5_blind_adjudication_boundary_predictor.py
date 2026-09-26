#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
代码47_P5盲裁决_边界预测器_GPCR钙瞬态.py
版本 v0.2.0（2026-08-14）

【版本史】v0.1.0 首版 → v0.1.1：脉冲定时修订（预注册升 v0.2，零统计量窗口）——
原"中位轨迹直接检 35 峰"不可行（背景上漂淹没低剂量峰、低剂量档中位轨迹本无峰，
首次运行 guardrail 中止=空跑①）；改为去漂移残差检高剂量峰 + 锚定 120 s 均匀网格外推。
→ v0.2.0：预注册升 v0.3（空跑②"检出峰未全部落格"触发深诊）——均匀网格假设被证伪
（换档边界 85–240 s 不规则间隔、低剂量块群体不可见、末端 ionomycin 附加峰）；
§6 全替换为"逐细胞峰时刻聚类 → 35 位链组装"两级管线，含逐实验计时排除（纯元数据判据，
排除原因入日志；收录 <8 实验或有效细胞 <150 空跑中止；[150,200) 挂功效缩减限定语）。
判定线/单元/特征/条款/方向仍一字未动。

【案件】P5：边界预测器首次盲测——Keshelava 2018（Nat Commun 9:876）M3R-GPCR→Ca²⁺
瞬态尖峰，7 档 ACh 升序 × 每档 5 脉冲，within-cell 设计，27 实验 / 433 细胞（论文口径）。

【冻结依据】预注册_P5_边界预测器首测_GPCR钙瞬态_v01.md（注册·八十九，2026-08-14 冻结）。
本脚本所有判定线、单元定义、特征定义、guardrail 与该预注册逐字对应；任何修订只能发生在
零统计量窗口（guardrail 中止=空跑）并全文留痕（代码46 三次空跑先例）。

【执行地点】用户 2026-08-14 明示指令由我方沙箱执行（原"用户机执行"约束经用户变更，已入档）。

【冻结参数块 D1–D8】
D1 数据：27 个 experiment_XX.dot（行=时间帧 1 帧/s，列=trace1..N，无表头，ASCII 制表符）
D2 剂量映射：脉冲 1-5=D1(100nM) … 31-35=D7(10µM)；7 档升序；每档 5 脉冲
D3 脉冲定时（预注册 v0.3 §6）：逐细胞峰时刻(σ_diff口径,prominence≥4σ_diff,distance≥90)
   → 5s 分箱群体聚类(阈 max(3,⌈0.08N⌉),<60s 分裂簇加权合并) → T̂=规则间隔中位数∈[115,125]
   → 锚=升序首个过(尾窗/瞬态回落/前驱/ionomycin鉴别)的簇 → 回走 34 步(严格层[85,170]→中点净空
   →放宽层(170,250]，计分 2·count−0.1·|gap−T̂|) → 链校验(第36/37位±30s 无簇、D7块无空洞、
   空洞≤12 且连续≤3、基线 p1−8≥35s、匹配≥23)；onset=峰位−8s；计时排除逐实验登记
D4 特征（k=3）：基线=首 onset 前全部帧；μ₀/σ₀=基线均值/标准差；脉冲窗=[onset,onset+30s]；
   peak_p=窗内max−μ₀；responded_p=(peak_p≥kσ₀)；τ_p=窗内首穿时刻−onset（无穿越=inf）；
   每(细胞,剂量)：count=responded 数(0-5)；peak=5 脉冲 peak_p 中位数；τ=responded 脉冲 τ_p 中位数（全无=inf）
D5 单元：低区 (D1,D2)(D2,D3) + 高区 (D5,D6)(D6,D7) = 4 裁决单元；中区 (D3,D4)(D4,D5) 描述臂
D6 AUC 方向：τ-AUC=P(τ低>τ高)+½平；count/peak-AUC=P(高>低)+½平；MW 平均秩，inf 删失按平
D7 判定线：P5-1 命中=4/4 单元 τ-AUC∈[0.40,0.60)，证伪=≥2 单元出 [0.40,0.60]；
   P5-2 命中=4/4 count-AUC≥0.60，证伪=4/4<0.60；P5-3 校准=≥3/4 peak-AUC≥0.60（对照臂不进总裁决）；
   总裁决：获支持=P5-1命中∧P5-2命中；证伪=任一侧证伪；其余=中间态
D8 统计：种子 20260815；逐细胞配对 percentile bootstrap×2000；N_MIN=30；引擎交付前 scipy 逐位单测

【细胞有效性（独立于作者 config，冻结）】σ₀>0；基线帧≥30；全程有限；覆盖脉冲 35 窗。
【稳健臂】k=2 / k=4；脱敏校正（按块内脉冲位序扣全局中位漂移）；最大实验纯净臂。
【输出（主产物先于描述臂）】单元表.csv / 稳健臂.csv / 裁决.json / 判定日志.txt / 描述臂.csv
"""

import os
import sys
import json
import numpy as np

ROOT = "/mnt/agents/output/01_细胞线/公开数据/Keshelava2018/SD1/source_data_1"
OUTDIR = "/mnt/agents/output/01_细胞线/结果/Keshelava2018_P5"
os.makedirs(OUTDIR, exist_ok=True)
LOG_PATH = os.path.join(OUTDIR, "代码47_判定日志.txt")
_log = open(LOG_PATH, "w", encoding="utf-8")


def say(msg):
    print(msg)
    _log.write(str(msg) + "\n")
    _log.flush()


def fail(msg):
    say("【中止】" + msg)
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
MIN_BASELINE = 35       # 实验级：首 onset 前基线(s)，低于则计时排除（预注册 v0.3 §6-F）
MIN_EXP = 8             # 收录实验下限（预注册 v0.3 §7⑥）
MIN_CELLS = 150         # 有效细胞下限；[150,200) 挂功效缩减限定语
UNITS = [(0, 1), (1, 2), (4, 5), (5, 6)]      # 0-based；低区 (D1,D2)(D2,D3)，高区 (D5,D6)(D6,D7)
MID_UNITS = [(2, 3), (3, 4)]                   # 描述臂

say("代码47 v0.2.0 —— P5 盲裁决（边界预测器首测，Keshelava 2018 GPCR→Ca²⁺；预注册 v0.3）")
say(f"冻结参数：SEED={SEED} N_BOOT={N_BOOT} k={K_MAIN} 充分性线={M_AUC} 死区带=[{BAND_LO},{BAND_HI})")


# ============ AUC 引擎（代码44–46 同源实现：MW 平均秩，inf 删失按平局） ============
def auc_mw(x, y):
    """P(X>Y)+½P(X=Y)，平均秩。x=方向正向组。"""
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
    """guardrail⑦：对 scipy 逐位核对 + inf 删失方向 + 全删失平局"""
    from scipy.stats import mannwhitneyu
    rng = np.random.default_rng(0)
    for _ in range(5):
        a = rng.normal(0, 1, 40)
        b = rng.normal(0.3, 1, 35)
        u, _ = mannwhitneyu(a, b, alternative="two-sided")
        ref = u / (len(a) * len(b))
        if abs(auc_mw(a, b) - ref) > 1e-12:
            fail("AUC 引擎与 scipy 不一致（无删失）")
    a = [1.0, 2.0, np.inf, np.inf]
    b = [0.5, 1.5, 3.0, np.inf]
    # 手工：P(a>b)=4/16（a3,a4>b1,b2），平局 a1? 逐对： (1,.5)W (1,1.5)L (1,3)L (1,inf)L
    # (2,.5)W (2,1.5)W (2,3)L (2,inf)L (inf,.5)W (inf,1.5)W (inf,3)W (inf,inf)T ×2
    # W=5? 重数：a3=inf>b1,b2,b3=3W +T; a4=inf 同 3W+T → W=1+2+3+3=9? a1:1>.5=1W 余3L；a2:2>.5,1.5=2W，2<3，2<inf
    # → W=1+2+3+3=9，T=2，AUC=(9+1)/16=0.625
    if abs(auc_mw(a, b) - 0.625) > 1e-12:
        fail(f"AUC 引擎 inf 删失方向错误（得 {auc_mw(a, b)}）")
    if abs(auc_mw([np.inf, np.inf], [np.inf, np.inf]) - 0.5) > 1e-12:
        fail("AUC 引擎全删失平局错误")
    say("guardrail⑦ AUC 引擎 scipy 逐位单测通过（无删失 5/5 + inf 删失方向 + 全删失平局）")


# ============ 数据装载 + guardrail 链 ============
def load_all():
    files = sorted(f for f in os.listdir(ROOT) if f.endswith(".dot"))
    if len(files) != 27:
        fail(f"实验文件数 {len(files)} ≠ 27")
    exps = []
    tot_cols = 0
    for f in files:
        arr = np.loadtxt(os.path.join(ROOT, f), delimiter="\t")
        if arr.ndim != 2:
            fail(f"{f} 不是二维矩阵")
        t, c = arr.shape
        if not (4100 <= t <= 4800):
            fail(f"{f} 行数 {t} 出界 [4100,4800]")
        if not (5 <= c <= 40):
            fail(f"{f} 列数 {c} 出界 [5,40]")
        tot_cols += c
        exps.append((f, arr))
    if not (430 <= tot_cols <= 560):
        fail(f"总轨迹列 {tot_cols} 出界 [430,560]")
    say(f"guardrail①②③ 通过：27 实验，总轨迹列 {tot_cols}")
    return exps


def pulse_timeline(name, arr):
    """D3（v0.2.0，预注册 v0.3 §6）：两级管线——逐细胞峰时刻聚类 → 35 位链组装。
    纯定时元数据（峰时刻+群体聚合轨迹），不产出任何逐细胞×剂量读出统计量。
    返回 (onsets, info)；计时排除返回 (None, 原因)——排除为纯元数据判据，逐条入日志。"""
    from scipy.signal import find_peaks, medfilt
    nrows = arr.shape[0]
    # A. 逐细胞峰检测（仅取时刻）
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
        return None, "无可用峰"
    # B. 群体聚类（5s 分箱 → 强箱 → ≤10s 并簇 → <60s 分裂簇计数加权合并）
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

    # C. 周期估计
    d2 = np.diff(cs)
    reg = d2[(d2 >= 100) & (d2 <= 130)]
    if len(reg) < 5:
        return None, "规则间隔<5"
    That = float(np.median(reg))
    if not (115.0 <= That <= 125.0):
        return None, f"T̂={That:.1f}越界"

    # E. 回走（严格层 [85,170] → 中点净空检查 → 放宽层 (170,250]）
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

    # D. 锚定位（升序首个全过者）+ F. 链校验
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
    return None, "无有效锚/链"


# ============ 特征化 ============
def featurize(arr, onsets, k):
    """返回每细胞特征表：list[dict(count=[], peak=[], tau=[])]，按剂量索引。"""
    t_end = onsets[-1] + WIN
    base_end = onsets[0]
    cells = []
    n_excl = 0
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
        thr = mu0 + k * sd0
        cnt = np.zeros(N_DOSE)
        pk = np.zeros(N_DOSE)
        ta = np.full(N_DOSE, np.inf)
        pk_by_dose = [[] for _ in range(N_DOSE)]
        ta_by_dose = [[] for _ in range(N_DOSE)]
        for p, o in enumerate(onsets):
            w = v[o:o + WIN]
            d = p // PPD
            peak = w.max() - mu0
            pk_by_dose[d].append(peak)
            if peak >= k * sd0:
                cnt[d] += 1
                cross = np.where(w >= thr)[0]
                ta_by_dose[d].append(float(cross[0]) if len(cross) else np.inf)
        for d in range(N_DOSE):
            pk[d] = np.median(pk_by_dose[d])
            ta[d] = np.median(ta_by_dose[d]) if ta_by_dose[d] else np.inf
        cells.append({"cnt": cnt, "pk": pk, "ta": ta})
    return cells, n_excl


# ============ 单元 AUC + bootstrap ============
def unit_auc(cells, pair, stat, rng=None, boot=False):
    ia, ib = pair
    if stat == "ta":
        xa = np.array([c["ta"][ia] for c in cells])
        xb = np.array([c["ta"][ib] for c in cells])
        # 方向：P(τ低>τ高) → x=低档
        x, y = xa, xb
    else:
        key = "cnt" if stat == "cnt" else "pk"
        x = np.array([c[key][ib] for c in cells])   # 高档为正向
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


# ============ 主流程 ============
def main():
    engine_selftest()
    exps = load_all()

    rng = np.random.default_rng(SEED)
    # ---- 脉冲时刻表（纯定时元数据阶段；计时排除逐条登记） ----
    timelines = {}
    timing_excl = []
    for name, arr in exps:
        onsets, info = pulse_timeline(name, arr)
        if onsets is None:
            timing_excl.append({"experiment": name, "reason": info})
            say(f"计时排除：{name} —— {info}")
            continue
        timelines[name] = onsets
        say(f"计时收录：{name}（T̂={info['That']:.1f} 空洞{info['holes']} p1={info['p1']} p35={info['p35']}）")
    say(f"guardrail④⑤ 完成：收录 {len(timelines)}/27，计时排除 {len(timing_excl)}")
    if len(timelines) < MIN_EXP:
        fail(f"收录实验 {len(timelines)} < {MIN_EXP}")

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
        fail(f"有效细胞 {len(all_cells)} < {MIN_CELLS}")
    if len(all_cells) < 200:
        power_note = f"功效缩减限定：有效细胞 {len(all_cells)} ∈ [150,200)"
        say("【限定语】" + power_note)
    say(f"装载完成：有效细胞 {len(all_cells)}（无效排除 {tot_excl}；计时排除 {len(timing_excl)} 实验；论文口径 433，独立有效性规则出入如实报告）")

    # ---- 主分析（4 裁决单元 × 3 统计量，带 CI） ----
    main_res = run_units(all_cells, UNITS, rng, boot=True)

    say("\n===== 主分析：4 裁决单元（低区 D1↔D2、D2↔D3；高区 D5↔D6、D6↔D7） =====")
    header = f"{'统计量':<8}{'单元':<10}{'AUC':>8}{'CI95':>22}{'n':>6}"
    say(header)
    rows = []
    for stat, label in (("ta", "τ事件时"), ("cnt", "计数"), ("pk", "峰值(付费)")):
        for r in main_res[stat]:
            say(f"{label:<8}{r['pair']:<10}{r['auc']:>8.4f}   [{r['lo']:.4f},{r['hi']:.4f}]{r['n']:>6}")
            rows.append({"stat": stat, "pair": r["pair"], "auc": r["auc"],
                         "lo": r["lo"], "hi": r["hi"], "n": r["n"]})

    # ---- D7 判定（冻结表逐字执行） ----
    ta_aucs = [r["auc"] for r in main_res["ta"]]
    cnt_aucs = [r["auc"] for r in main_res["cnt"]]
    pk_aucs = [r["auc"] for r in main_res["pk"]]

    n_in = sum(1 for a in ta_aucs if BAND_LO <= a < BAND_HI)
    n_out = 4 - n_in
    if n_in == 4:
        v51 = "命中（死区确认：4/4 单元 τ-AUC∈[0.40,0.60)）"
    elif n_out >= 2:
        v51 = f"证伪（{n_out}/4 单元出 [0.40,0.60] 带）"
    else:
        v51 = f"中间态（{n_in}/4 单元在带内）"

    n_pass = sum(1 for a in cnt_aucs if a >= M_AUC)
    if n_pass == 4:
        v52 = "命中（活区确认：4/4 单元 count-AUC≥0.60）"
    elif n_pass == 0:
        v52 = "证伪（4/4 单元 count-AUC<0.60）"
    else:
        v52 = f"中间态（{n_pass}/4 单元达线）"

    n_cal = sum(1 for a in pk_aucs if a >= M_AUC)
    v53 = f"校准通过（{n_cal}/4 单元 peak-AUC≥0.60）" if n_cal >= 3 else f"数据集级弱解码限定（仅 {n_cal}/4 单元 peak-AUC≥0.60）"

    hit51 = n_in == 4
    fal51 = n_out >= 2
    hit52 = n_pass == 4
    fal52 = n_pass == 0
    if hit51 and hit52:
        overall = "P5 获支持（边界预测器：死区复核命中 ∧ 活区首测命中）"
    elif fal51 or fal52:
        overall = "P5 证伪（边界预测器在本形态上死亡，死亡侧见分款）"
    else:
        overall = "P5 中间态"

    say("\n===== 判定（冻结表 D7 逐字执行） =====")
    say(f"P5-1（R1 死区复核，τ）：{v51}　单元值 {['%.4f' % a for a in ta_aucs]}")
    say(f"P5-2（R2 活区首测，计数）：{v52}　单元值 {['%.4f' % a for a in cnt_aucs]}")
    say(f"P5-3（R3 校准，付费峰值，对照臂）：{v53}　单元值 {['%.4f' % a for a in pk_aucs]}")
    say(f"【总裁决】{overall}")

    # ---- 主产物落盘（先于描述臂/稳健臂结论性内容） ----
    import csv
    with open(os.path.join(OUTDIR, "代码47_单元表.csv"), "w", newline="", encoding="utf-8-sig") as f:
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
    with open(os.path.join(OUTDIR, "代码47_裁决.json"), "w", encoding="utf-8") as f:
        json.dump(verdict, f, ensure_ascii=False, indent=2,
                  default=lambda o: int(o) if isinstance(o, np.integer) else float(o) if isinstance(o, np.floating) else str(o))
    say("主产物落盘：单元表.csv / 裁决.json / 判定日志.txt")

    # ---- 稳健臂 ----
    say("\n===== 稳健臂 =====")
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
    # 脱敏校正臂：按块内位序扣全局中位 peak 漂移（论文同款），重算 cnt/pk
    # 收集所有 (cell, dose, j) peak
    say("脱敏校正臂：按块内脉冲位序 j∈{1..5} 扣群体中位漂移后重算")
    # 重新特征化时保留每脉冲 peak —— 简化：此处仅对 pk/cnt 近似校正（τ 不受幅度漂移影响）
    # 全局位序中位数
    # （严格实现：featurize 不返回逐脉冲值；此臂用 k=3 重特征化并在脉冲层校正）
    cells_dc = []
    for name, arr in exps:
        if name not in timelines:
            continue
        onsets = timelines[name]
        t_end = onsets[-1] + WIN
        base_end = onsets[0]
        # 第一遍：收集全部脉冲 peak 以估计位序漂移
        peaks_j = [[] for _ in range(PPD)]
        valid_cols = []
        for c in range(arr.shape[1]):
            v = arr[:, c]
            if not np.all(np.isfinite(v[:t_end])):
                continue
            b = v[:base_end]
            mu0 = b.mean()
            sd0 = b.std(ddof=1)
            if not (sd0 > 0) or len(b) < 30:
                continue
            pps = []
            for p, o in enumerate(onsets):
                pps.append(v[o:o + WIN].max() - mu0)
            for p, val in enumerate(pps):
                peaks_j[p % PPD].append(val)
            valid_cols.append((v, mu0, sd0, pps))
        med_j = np.array([np.median(x) for x in peaks_j])
        corr = med_j - med_j.mean()
        for v, mu0, sd0, pps in valid_cols:
            thr = mu0 + K_MAIN * sd0
            cnt = np.zeros(N_DOSE)
            pk = np.zeros(N_DOSE)
            ta = np.full(N_DOSE, np.inf)
            pk_d = [[] for _ in range(N_DOSE)]
            ta_d = [[] for _ in range(N_DOSE)]
            for p, o in enumerate(onsets):
                d = p // PPD
                val = pps[p] - corr[p % PPD]
                pk_d[d].append(val)
                if val >= K_MAIN * sd0:
                    cnt[d] += 1
                    w = v[o:o + WIN]
                    cross = np.where(w >= thr)[0]
                    ta_d[d].append(float(cross[0]) if len(cross) else np.inf)
            for d in range(N_DOSE):
                pk[d] = np.median(pk_d[d])
                ta[d] = np.median(ta_d[d]) if ta_d[d] else np.inf
            cells_dc.append({"cnt": cnt, "pk": pk, "ta": ta})
    res_dc = run_units(cells_dc, UNITS, rng, boot=False)
    for stat in ("ta", "cnt"):
        for r in res_dc[stat]:
            rob_rows.append({"arm": "脱敏校正", "stat": stat, "pair": r["pair"], "auc": r["auc"], "n": r["n"]})
            say(f"脱敏校正 {stat} {r['pair']} AUC={r['auc']:.4f}")
    # 最大实验纯净臂
    big = max(per_exp_cells, key=lambda x: len(x[1]))
    res_big = run_units(big[1], UNITS, rng, boot=False)
    for stat in ("ta", "cnt"):
        for r in res_big[stat]:
            rob_rows.append({"arm": f"单实验({big[0]},n={len(big[1])})", "stat": stat, "pair": r["pair"], "auc": r["auc"], "n": r["n"]})
            say(f"单实验纯净臂 {stat} {r['pair']} AUC={r['auc']:.4f}")
    with open(os.path.join(OUTDIR, "代码47_稳健臂.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["arm", "stat", "pair", "auc", "n"])
        w.writeheader()
        w.writerows(rob_rows)
    say("稳健臂落盘：稳健臂.csv")

    # ---- 描述臂（不进判定） ----
    say("\n===== 描述臂（不进判定） =====")
    desc_rows = []
    mid_res = run_units(all_cells, MID_UNITS, rng, boot=True)
    for stat in ("ta", "cnt", "pk"):
        for r in mid_res[stat]:
            desc_rows.append({"section": "中区对", "stat": stat, "key": r["pair"],
                              "auc": r["auc"], "lo": r["lo"], "hi": r["hi"], "n": r["n"]})
            say(f"中区对 {stat} {r['pair']} AUC={r['auc']:.4f} [{r['lo']:.4f},{r['hi']:.4f}]")
    # 逐剂量响应率（count≥1 的细胞比例）与平均计数
    for d in range(N_DOSE):
        rr = np.mean([1.0 if c["cnt"][d] >= 1 else 0.0 for c in all_cells])
        mc = np.mean([c["cnt"][d] for c in all_cells])
        desc_rows.append({"section": "响应率", "stat": "cnt", "key": f"D{d+1}({DOSES_NM[d]}nM)",
                          "auc": rr, "lo": mc, "hi": np.nan, "n": len(all_cells)})
        say(f"D{d+1}({DOSES_NM[d]:.0f}nM) 响应率={rr:.3f} 平均计数={mc:.2f}")
    with open(os.path.join(OUTDIR, "代码47_描述臂.csv"), "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=["section", "stat", "key", "auc", "lo", "hi", "n"])
        w.writeheader()
        w.writerows(desc_rows)
    say("描述臂落盘：描述臂.csv")

    say("\n===== 终稿判定复述（冻结表逐字） =====")
    say(f"P5-1：{v51}")
    say(f"P5-2：{v52}")
    say(f"P5-3：{v53}")
    say(f"【总裁决】{overall}")
    say("（死活双录：本日志与裁决.json 同步入总账，冻结裁决事后永不挪动）")
    _log.close()


if __name__ == "__main__":
    main()
