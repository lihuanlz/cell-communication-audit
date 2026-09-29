# -*- coding: utf-8 -*-
"""
代码46 v0.1.5 —— P4 盲裁决：ERK 重裁（GPCR-KTR 原始轨迹级）
=================================================================
v0.1.1（2026-08-14）：首次空跑修复——实际文件无 Experiment 列，细胞键改
  Unique_Object（缺失时退化 Date+Slide+Position+Object 复合键）；X,Y 为逐帧
  质心坐标退出键设计；批次纯净臂改用 Date。
v0.1.2（2026-08-14，未交付运行即废止，如实登记）：第二次空跑后误判为
  「六孔 0.5 min 相位交错、视野级锚定」——schema 验证发现会把 t0 错锚到
  5.5 min，作废。
v0.1.3（2026-08-14，定锚）：schema 验证确认作者已将原始采集（3.5 min/帧、
  六孔位 0.5 min 延迟交错）线性插值对齐到**全局统一 0.5 min 网格**——全部
  7,393 细胞满帧 162、2.5–83.0 min；t = Time_in_min − 23.5（time point 7 =
  2.5+6×3.5 = 23.5 = 加样点 t=0，与作者 AUC 窗 tp9–18↔刺激后 7–38.5 min
  逐字吻合）；tp7（t=0）帧弃用；基线 t<0 预期 42 帧。判定线/窗定义零改动。
v0.1.4（2026-08-14）：第三次空跑修复——实际剂量档数统一为 6（全配体×全
  抑制剂条件 schema 验证；v0.1.0 冻结的"7 档"为错误假设）。配对重锚：
  低区=对(1,2)(2,3)、高区=对(4,5)(5,6)，中间对(3,4)=半饱和过渡带剔出判定
  进描述臂。判定单元数与计数线按 6 档重算（P4-1: 4 单元；P4-2: 4+4；
  P4-3: 10 单元 ≥8/10 命中 ≥3/10 证伪）；判定阈值 0.60/0.05/0.01/0.03
  一字未动。预注册书升 v0.3。
  三次空跑（guardrail 命中，均不构成裁决）均已入总账。
v0.1.5（2026-08-14）：正式裁决已产生（判定先于崩溃完整打印，裁决成立，
  沿代码42 先例）后修复 json 落盘 bug——np.int64/np.float64（auc_mw 返回
  np 类型经 sum/abs 传播）不可直接 JSON 序列化，加 default 转换器 +
  显式 int()。判定逻辑、种子、rng 调用序列零改动 → 主分析逐位复现。
=================================================================
预注册书：预注册_P4_ERK重裁_GPCR-KTR轨迹级_v01.md（总账 注册·八十一）
数据集  ：Chavez-Abiega & Goedhart (2022) J Cell Sci 135(6):jcs259685
          Zenodo doi:10.5281/zenodo.5836623
运行方式：用户机 Spyder  %runfile 本文件   （CPU 轻量，无需 GPU）

冻结决策（开数据前写死，逐字沿预注册书 §2）：
  D1 时间基准（v0.1.3 定锚）：作者已将原始采集（3.5 min/帧、六孔位 0.5 min
     延迟交错）线性插值对齐到全局统一 0.5 min 网格（全细胞满帧 162，
     2.5–83.0 min）；t = Time_in_min − 23.5 min（time point 7 = 2.5+6×3.5 =
     23.5 = 加样点 t=0；作者 AUC 窗 tp9–18 ↔ 刺激后 7–38.5 min 逐字验证）。
     基线 = t<0 帧（预期 42 帧）；tp7（t=0，加样当帧归属模糊）弃用——
     不进基线、不进响应窗。Guardrail：全局网格起点 2.5±0.1、步长中位
     0.5±0.05、终点 ≥60 min，任一不符 → 空跑中止。
  D2 基线：t<0 全部帧（预期 6 帧）；基线帧 <3 剔细胞；σ₀=0 剔（ddof=1）。
  D3 窗：τ/计数/时程窗 = (0, 38.5] min post-stim；删失记 +inf，Mann–Whitney 平均秩。
  D4 纳入/剔除：非有限或 <0 的 C/N 帧剔除该帧（计缺失）；缺失 >10% 剔细胞；
     轨迹末端须 ≥38.5 min post-stim。主分析抑制剂 = DMSO。
  D5 无平滑。
  D6 meanArea 三分层边界 = 判定单元合并细胞的 33.33/66.67 百分位，标签随 bootstrap。
  D7 AUC 方向冻结：τ 单元 AUC = P(τ_低档 > τ_高档) + ½ties（高剂量更早为预言方向）；
     P4-3 分层罚 Δ = AUC_混合 − mean(AUC_层内)。
  D8 判定表：命中 / 证伪 / 未命中-中间态 / 数据不足（单元任一组 n<30）四态
     + P4-F 总证伪（任一条款落入证伪格 → 总裁决证伪）；死活双录入总账。

条款判定线（冻结；v0.1.4 按实际 6 档重锚单元数，阈值零改动）：
  P4-1：his+UK × 低区 2 对 = 4 单元，τ(3) AUC ≥ 0.60 全达线 → 命中；
        任一单元 AUC<0.5（方向反转）或全单元不达线 → 证伪；之间 → 中间态。
  P4-2：高区 4 单元均值 AUC_H < 0.60 且 D = 低区均值 − AUC_H ≥ 0.05 → 命中；
        AUC_H ≥ 0.60 或 D ≤ 0 → 证伪；之间 → 中间态。逻辑独立于 P4-1。
  P4-3：10 单元（2 配体 × 5 对）meanArea 分层罚 Δ；
        ≥8/10 单元 |Δ|≤0.01 → 命中；≥3/10 单元 |Δ|≥0.03 → 证伪；之间 → 中间态。

guardrail 链（代码44/45 先例）：细胞键不唯一 / 帧序异常 / 时间网格不符 /
  剂量档 ≠7 / 必需列缺失 → 【中止】空跑，不构成裁决。
"""

import os, re, sys, json, datetime
import numpy as np
import pandas as pd

# ---------------- 配置（冻结） ----------------
ROOT      = r"C:/Users/lihua/Desktop/NC/data/ChavezAbiega2022"
SEED      = 20260815
N_BOOT    = 2000
N_MIN     = 30
K_MAIN    = 3
K_ROB     = (2, 4)

M41       = 0.60          # P4-1 单元判定线
M42_H     = 0.60          # P4-2 高区均值线
M42_D     = 0.05          # P4-2 梯度线
P43_LO, P43_HI = 0.01, 0.03
P43_HIT, P43_KILL = 8, 3   # 10 单元（2 配体 × 5 对）的命中/证伪计数线

T_STIM    = 23.5          # min（time point 7 = 加样点 t=0）
DT_GLOB   = 0.5           # min/帧（作者插值后的全局统一网格）
T_START   = 2.5           # 网格起点
T_END_MIN = 60.0          # 网格终点下限（预期 83.0）
WIN_LO, WIN_HI = 0.0, 38.5    # 免费统计量窗 (0, 38.5]
PAID_LO, PAID_HI = 7.0, 38.5  # 付费积分窗 [7, 38.5]（作者 AUC 窗）

LIG_MAIN  = ["Histamine", "UK"]     # 主分析配体
LIG_DESC  = "S1P"                    # 描述臂配体
FILE_F1   = os.path.join("Figure_1", "Data", "All_Ex_{lig}_DMSO.csv")
FILE_F3   = os.path.join("Figure_3", "Data", "ALL_{lig3}.csv")   # 抑制剂描述臂（可缺）
LIG3_NAME = {"Histamine": "His", "UK": "UK", "S1P": "S1P"}

OUT_DIR   = os.path.dirname(os.path.abspath(__file__))
LOG_LINES = []

def log(s):
    print(s)
    LOG_LINES.append(str(s))

def abort(msg):
    log("【中止】" + msg + " —— 空跑，不构成裁决。")
    flush_log()
    sys.exit(2)

def flush_log():
    with open(os.path.join(OUT_DIR, "代码46_判定日志.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(LOG_LINES) + "\n")

# ---------------- 读取与 guardrail ----------------
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
    """读长表 → 规范化列 → guardrail → 逐细胞字典列表。
    desc=False：主分析，guardrail 失败 = 空跑中止；
    desc=True ：描述臂，guardrail 失败 = raise ValueError（调用方转警告跳过）。"""
    def fail(msg):
        if desc:
            raise ValueError(msg)
        abort(msg)
    if not os.path.isfile(path):
        fail(f"文件缺失：{path}")
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
        fail(f"{os.path.basename(path)}：必需列缺失 {miss}（实际列 {list(df.columns)}）")
    # 细胞键：优先 Unique_Object；否则 (Date,Slide,Position,Object) 复合键。
    # X,Y 为逐帧质心坐标（细胞移动），严禁进键。
    if "UID" in df.columns and df["UID"].notna().all():
        key_cols = ["UID"]
    elif all(c in df.columns for c in ["Date", "Slide", "Pos", "Obj"]):
        key_cols = ["Date", "Slide", "Pos", "Obj"]
    else:
        fail(f"{os.path.basename(path)}：无可用细胞键（需 Unique_Object 或 Date/Slide/Position/Object；实际列 {list(df.columns)}）")
    df["DoseV"] = df["Dose"].apply(parse_dose)
    if df["DoseV"].isna().any():
        fail(f"{os.path.basename(path)}：Condition 列存在无法解析的剂量值")
    # 时间网格 guardrail：作者插值后的全局统一 0.5 min 网格
    tg = np.sort(df["Time"].unique())
    if len(tg) < 2 or abs(tg.min() - T_START) > 0.1 or tg.max() < T_END_MIN or \
       abs(np.median(np.diff(tg)) - DT_GLOB) > 0.05:
        fail(f"{os.path.basename(path)}：全局网格异常（{tg.min() if len(tg) else float('nan'):.2f}–"
             f"{tg.max() if len(tg) else float('nan'):.2f} min，中位步长 "
             f"{np.median(np.diff(tg)) if len(tg)>1 else float('nan'):.3f}；预期 2.5–≥60 min、0.5 步长）")
    # 剂量档 guardrail（v0.1.4：实际统一 6 档，schema 已验证全配体×全抑制剂）
    doses = np.sort(df["DoseV"].unique())
    if len(doses) != 6:
        fail(f"{os.path.basename(path)}：剂量档数 {len(doses)} ≠ 6（档位 {doses.tolist()}）")
    # 细胞键唯一性 guardrail
    if df.duplicated(subset=key_cols + ["Time"]).any():
        fail(f"{os.path.basename(path)}：细胞键 {key_cols}+Time 不唯一")
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
        log(f"  [guardrail] {os.path.basename(path)}：{n_drop_dose} 个键内剂量不唯一，已剔除")
    log(f"  读取 {os.path.basename(path)}：{len(cells)} 细胞，{len(doses)} 档剂量 {np.round(doses,4).tolist()}，帧网格 {tg.min():.2f}–{tg.max():.2f} min（{len(tg)} 帧）")
    return cells, doses

# ---------------- 特征化（D2–D5） ----------------
def featurize(c, k=K_MAIN):
    """返回 None=剔除；否则特征字典。规则全部冻结。
    dt 取该细胞实际帧间隔（主分析 0.5 插值网格 / 抑制剂臂 3.5 采集网格）。"""
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
    # 计数：窗内上穿次数（below→>= 跳变）
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

# ---------------- AUC（Mann–Whitney，ties/删失平均秩） ----------------
def auc_mw(low, high):
    """AUC = P(τ_low > τ_high) + ½ P(tie)。inf 自然入秩。"""
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

# ---------------- 判定单元装配 ----------------
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

# ---------------- 主流程 ----------------
def main():
    log(f"代码46 v0.1.5 —— P4 盲裁决：ERK 重裁（GPCR-KTR 轨迹级）")
    log(f"运行时刻 {datetime.datetime.now().isoformat(timespec='seconds')} ｜ SEED={SEED} ｜ bootstrap={N_BOOT}（首轮即执行）")
    log(f"ROOT={ROOT}")
    if not os.path.isdir(ROOT):
        abort(f"数据根目录不存在：{ROOT}（请先按预注册书 §5② 备数据）")
    rng = np.random.default_rng(SEED)

    unit_rows, rob_rows, desc_rows = [], [], []
    store = {}   # (lig, zone, pair_idx, side) -> tau 数组，供 P4-2 聚合 bootstrap
    main_data = {}  # lig -> (cells, doses)，供 Akt 描述臂复用

    for lig in LIG_MAIN:
        path = os.path.join(ROOT, FILE_F1.format(lig=lig))
        log(f"\n=== 主分析配体 {lig}（DMSO） ===")
        cells, doses = load_long(path)
        main_data[lig] = (cells, doses)
        feats = [featurize(c) for c in cells]
        n_ok = sum(f is not None for f in feats)
        log(f"  纳入 {n_ok}/{len(cells)} 细胞（剔除 {len(cells)-n_ok}：基线帧<3 / σ₀=0 / 末端<38.5 / 缺失>10%）")
        by_d = cells_by_dose(cells, feats)
        for d in doses:
            log(f"    档 {d:g}：n={len(by_d.get(d, []))}")

        # 6 档 5 对分锚：低区前 2 对、高区后 2 对；中间对(3,4)=半饱和过渡带，
        # 剔出判定进描述臂（预注册 v0.3）
        low_pairs  = [(doses[0], doses[1]), (doses[1], doses[2])]
        high_pairs = [(doses[3], doses[4]), (doses[4], doses[5])]
        mid_pair   = (doses[2], doses[3])

        # ---- P4-1 / P4-2 单元 ----
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
                log(f"  [{zone}] τ({K_MAIN}) {dl:g}↔{dh:g}：AUC={pt:.4f} [{lo:.3f},{hi:.3f}] n={len(low)}/{len(high)} 达线={hit if suff else '数据不足'}")

        # ---- 稳健臂 k=2/4（点估计） ----
        for kk in K_ROB:
            feats_k = [featurize(c, k=kk) for c in cells]
            by_k = cells_by_dose(cells, feats_k)
            for zone, pairs in (("低区", low_pairs), ("高区", high_pairs)):
                for dl, dh in pairs:
                    low, high = pair_unit(by_k, dl, dh)
                    if len(low) >= N_MIN and len(high) >= N_MIN:
                        rob_rows.append({"arm": f"k={kk}", "ligand": lig, "zone": zone,
                                         "pair": f"{dl:g}↔{dh:g}", "auc": auc_mw(low, high)})

        # ---- 批次纯净臂（最大单一 Experiment） ----
        exps = pd.Series([c["batch"] for c in cells]).value_counts()
        big = exps.index[0]
        sub_cells = [c for c in cells if c["batch"] == big]
        sub = [featurize(c) for c in sub_cells]
        by_b = cells_by_dose(sub_cells, sub)
        for zone, pairs in (("低区", low_pairs), ("高区", high_pairs)):
            for dl, dh in pairs:
                low, high = pair_unit(by_b, dl, dh)
                if len(low) >= N_MIN and len(high) >= N_MIN:
                    rob_rows.append({"arm": f"批次纯净({big})", "ligand": lig, "zone": zone,
                                     "pair": f"{dl:g}↔{dh:g}", "auc": auc_mw(low, high)})

        # ---- P4-3：meanArea 三分层罚（全部 5 对，含中间过渡带对） ----
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
            # 重建分组：pool 顺序 = 先 dl 后 dh
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
            log(f"  [P4-3] {dl:g}↔{dh:g}：Δ={delta:+.4f}（混合 {auc_pool:.4f} − 层内均值 {np.mean(layer_aucs) if layer_aucs else float('nan'):.4f}；层 n={layer_ns}）")

        # ---- 描述臂：付费峰值/积分 + 计数/时程 + 响应率 ----
        for dl, dh in all_pairs:
            for key, nm in (("peak", "峰值CN"), ("auci", "积分CN"), ("cnt", "计数"), ("dur", "时程")):
                low, high = pair_unit(by_d, dl, dh, key)
                low = [x for x in low if np.isfinite(x)]; high = [x for x in high if np.isfinite(x)]
                if len(low) >= N_MIN and len(high) >= N_MIN:
                    # 付费类方向 = 高剂量更大 → AUC = P(high > low)
                    desc_rows.append({"arm": nm, "ligand": lig, "pair": f"{dl:g}↔{dh:g}",
                                      "auc": auc_mw(high, low)})
        for d in doses:
            fs = by_d.get(d, [])
            if fs:
                desc_rows.append({"arm": "响应率", "ligand": lig, "pair": f"档{d:g}",
                                  "auc": float(np.mean([f["resp"] for f in fs]))})
        # 中间对（半饱和过渡带）τ 读数：剔出判定，进描述臂
        low, high = pair_unit(by_d, mid_pair[0], mid_pair[1])
        if len(low) >= N_MIN and len(high) >= N_MIN:
            desc_rows.append({"arm": "中间对τ(过渡带)", "ligand": lig,
                              "pair": f"{mid_pair[0]:g}↔{mid_pair[1]:g}", "auc": auc_mw(low, high)})

    # ---------------- P4-2 聚合判定统计量 ----------------
    log("\n=== P4-2 聚合（饱和区失序检验） ===")
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
    # D 的 bootstrap：同步重抽全部 8 单元
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
    log(f"  AUC_L(低区均值)={AUC_L:.4f} ｜ AUC_H(高区均值)={AUC_H:.4f} ｜ D={D_pt:+.4f} [{D_lo:+.3f},{D_hi:+.3f}]")

    # ---------------- D8 判定（冻结表逐字） ----------------
    log("\n=== D8 判定（冻结表逐字执行） ===")
    verdicts = {}

    p41 = [r for r in unit_rows if r["clause"] == "P4-1"]
    if any(not r["sufficient"] for r in p41):
        verdicts["P4-1"] = "数据不足"
    elif all(r["pass"] for r in p41):
        verdicts["P4-1"] = "命中"
    elif any(np.isfinite(r["auc"]) and r["auc"] < 0.5 for r in p41) or not any(r["pass"] for r in p41):
        verdicts["P4-1"] = "证伪"
    else:
        verdicts["P4-1"] = "未命中-中间态"

    if len(aucH_list) < 4 or len(aucL_list) < 4:
        verdicts["P4-2"] = "数据不足"
    elif AUC_H < M42_H and D_pt >= M42_D:
        verdicts["P4-2"] = "命中"
    elif AUC_H >= M42_H or D_pt <= 0:
        verdicts["P4-2"] = "证伪"
    else:
        verdicts["P4-2"] = "未命中-中间态"

    p43 = [r for r in unit_rows if r["clause"] == "P4-3" and r["sufficient"]]
    n43_lo = sum(abs(r["auc"]) <= P43_LO for r in p43)
    n43_hi = sum(abs(r["auc"]) >= P43_HI for r in p43)
    if len(p43) < 10:
        verdicts["P4-3"] = "数据不足"
    elif n43_hi >= P43_KILL:
        verdicts["P4-3"] = "证伪"
    elif n43_lo >= P43_HIT:
        verdicts["P4-3"] = "命中"
    else:
        verdicts["P4-3"] = "未命中-中间态"

    overall = "P4 证伪（存在证伪条款）" if any(v == "证伪" for v in verdicts.values()) else \
              ("P4 命中（无证伪条款且存在命中条款）" if any(v == "命中" for v in verdicts.values())
               else "P4 未命中-中间态")
    for k, v in verdicts.items():
        log(f"  {k}：{v}")
    log(f"  总裁决：{overall}")
    log(f"  P4-3 计数：|Δ|≤0.01 → {n43_lo}/10；|Δ|≥0.03 → {n43_hi}/10")

    # ---------------- 主输出（先于描述臂落盘：主裁决存档不受描述臂异常影响） ----------------
    pd.DataFrame(unit_rows).to_csv(os.path.join(OUT_DIR, "代码46_单元表.csv"),
                                   index=False, encoding="utf-8-sig")
    pd.DataFrame(rob_rows).to_csv(os.path.join(OUT_DIR, "代码46_稳健臂.csv"),
                                  index=False, encoding="utf-8-sig")
    with open(os.path.join(OUT_DIR, "代码46_裁决.json"), "w", encoding="utf-8") as f:
        json.dump({"script": "代码46 v0.1.5", "seed": int(SEED), "n_boot": int(N_BOOT),
                   "verdicts": {k: str(v) for k, v in verdicts.items()}, "overall": str(overall),
                   "P4-2": {"AUC_L": float(AUC_L), "AUC_H": float(AUC_H),
                            "D": float(D_pt), "D_ci": [float(D_lo), float(D_hi)]},
                   "P4-3": {"n_le_001": int(n43_lo), "n_ge_003": int(n43_hi)},
                   "timestamp": datetime.datetime.now().isoformat(timespec="seconds")},
                  f, ensure_ascii=False, indent=2,
                  default=lambda o: int(o) if isinstance(o, np.integer)
                         else float(o) if isinstance(o, np.floating) else str(o))
    log("\n主输出已落盘：代码46_单元表.csv / 代码46_稳健臂.csv / 代码46_裁决.json")
    flush_log()

    # ---------------- 描述臂（不进判定；任何异常仅警告跳过） ----------------
    log("\n=== 描述臂（不进判定） ===")
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
        log(f"  S1P 描述臂完成（{len([r for r in desc_rows if r['ligand']=='S1P'])} 行）")
    except Exception as e:
        log(f"  [警告] S1P 描述臂跳过：{e}")

    # Akt 描述臂（复用主分析缓存，换 CN_AktRB 通道）
    for lig in LIG_MAIN:
        try:
            cells, doses = main_data[lig]
            akt_cells = [c for c in cells if c["akt"] is not None]
            if not akt_cells:
                raise ValueError("无 Akt 列或全为空")
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
            log(f"  Akt 描述臂完成（{lig}）")
        except Exception as e:
            log(f"  [警告] Akt 描述臂 {lig} 跳过：{e}")

    # 抑制剂描述臂（Figure_3，可缺）
    for lig in LIG_MAIN:
        p3 = os.path.join(ROOT, FILE_F3.format(lig3=LIG3_NAME[lig]))
        if not os.path.isfile(p3):
            log(f"  [警告] 抑制剂臂文件缺失，跳过：{p3}")
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
                log(f"  [警告] {os.path.basename(p3)} 无 Inhibitor 列，抑制剂臂跳过")
                continue
            # 网格检查：0.5 min 插值网格（Figure_1 形态）或 3.5 min 采集网格
            # （Figure_3 形态）均合法；t0=23.5 对两者通用（tp7）
            tgi = np.sort(df["Time"].unique())
            sti = np.median(np.diff(tgi)) if len(tgi) > 1 else float("nan")
            ok_grid = len(tgi) >= 2 and abs(tgi.min() - T_START) <= 0.1 and tgi.max() >= T_END_MIN and \
                      (abs(sti - DT_GLOB) <= 0.05 or abs(sti - 3.5) <= 0.2)
            if not ok_grid:
                log(f"  [警告] {os.path.basename(p3)} 网格异常（{tgi.min() if len(tgi) else float('nan'):.1f}–"
                    f"{tgi.max() if len(tgi) else float('nan'):.1f}，步长中位 {sti:.3f}），抑制剂臂跳过")
                continue
            df["DoseV"] = df["Dose"].apply(parse_dose)
            for inh, g0 in df.groupby("Inh"):
                doses_i = np.sort(g0["DoseV"].dropna().unique())
                if len(doses_i) < 2:
                    continue
                key_i = ["UID"] if "UID" in g0.columns else \
                        [c for c in ["Date", "Slide", "Pos", "Obj"] if c in g0.columns]
                if not key_i:
                    log(f"  [警告] 抑制剂臂 {lig}/{inh}：无可用细胞键，跳过")
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
                        desc_rows.append({"arm": f"抑制剂-{inh}", "ligand": lig,
                                          "pair": f"{doses_i[j]:g}↔{doses_i[j+1]:g}",
                                          "auc": auc_mw(low, high)})
            log(f"  抑制剂臂完成（{lig}）")
        except Exception as e:
            log(f"  [警告] 抑制剂臂 {lig} 异常跳过：{e}")

    pd.DataFrame(desc_rows).to_csv(os.path.join(OUT_DIR, "代码46_描述臂.csv"),
                                   index=False, encoding="utf-8-sig")
    log("\n输出：代码46_单元表.csv / 代码46_稳健臂.csv / 代码46_描述臂.csv / 代码46_裁决.json / 代码46_判定日志.txt")
    log("—— 死活双录：请将五份输出回传，无论结果如何如实入总账。——")
    flush_log()

if __name__ == "__main__":
    main()
