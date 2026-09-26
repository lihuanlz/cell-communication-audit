# -*- coding: utf-8 -*-
"""
================================================================================
代码45 —— P3 盲裁决：免费统计量充分性（Sung 2022, NF-κB 空间梯度）
版本 v0.1.1 ｜ 2026-08-14 ｜ 用户机运行
v0.1.1（注册·七十八）：首跑护栏拦截——R 实测 (5,腔室数)（MATLAB 24×5，仅前2列入分析），
  get_R 由"恰好2行"放宽为"≥2行"，取值逻辑不变。空跑不构成裁决，判定逻辑零改动。
================================================================================
预注册：预注册_P3_免费统计量充分性_NFkB空间梯度_v01.md（总账 注册·七十五，判定线已写死）
理论：定理 T′ 幸存类——免费统计量（事件时间/序数/计数/时程）无需跨细胞校准即足以解码。
盲态：本脚本写成前未计算任何判定统计量。字段语义全部来自作者公开 MATLAB 脚本文本
      （DeepAnalysisForFig4.m / HeatMap&Traces_Dose.m / HeatMap&Traces_Duration.m），
      未触碰轨迹数值。本脚本的第一次正式运行即为裁决运行，死活双录入账。

实现决定（全部在开数据前写死；依据=作者脚本文本与预注册，非结果）：
D1  TNF10 时间对齐采用 osD=[45,152]+osI_10。理由：HeatMap&Traces_Dose.m 的 10ng 段
    先声明 [45,152]+osI_10，其后出现的 [39,148]+osI_30 覆盖块经与 DeepAnalysisForFig4.m
    对照（该脚本第150-165行用 [45,152]+osI_10 且无覆盖块），判定为复制粘贴残留；
    且脚本注释挂钟时间独立印证（成像 19:04:27，首次馈液 19:49:11→≈45 min；
    二次馈液 21:36:33→≈152 min）。备选对齐 Option B（[39,148]+osI_30）仅作稳健臂报点估计。
D2  腔室纳入镜像作者脚本：100/30 ng→aa=1..12；10 ng→aa=[1,2,3,4,9,10,11,12]；
    时长系列→aa 2-4=15min、5-8=30min、9-12=60min（aa=1 号腔室作者未用）。
D3  分析窗 (0,120] min：τ/上穿计数/越阈时程均在此窗内定义；未越阈=右删失记 +inf。
D4  无效值=非有限或<0；缺失率>10%（按该细胞全部记录帧）剔除；
    基线（x≤0）有效帧<2 剔除；有效轨迹末端 x<120 剔除；σ₀=0 剔除（std 用 ddof=1）。
D5  无平滑：直接对原始比值轨迹判定越阈（作者 lowess 平滑仅用于展示，不进本解码器）。
D6  P3-2 三分带边界=TNF100 全体纳入细胞距离的 33.33/66.67 百分位（合并后固定，
    bootstrap 时标签随细胞携带，不重估边界）。
D7  AUC=Mann–Whitney（平均秩；+inf 参与排名）。方向写死：
    P3-1: P(τ_低剂量 > τ_高剂量)；P3-2: P(τ_较远带 > τ_较近带)；
    P3-3: P(越阈时程_长刺激 > 越阈时程_短刺激)。
D8  条款判定表（写死，两单元条款）：
    · 任一单元 n<30 → 该单元"数据不足"，条款=数据不足（不命中、不证伪，如实入账）；
    · 两单元 AUC 全 ≥ 边际 → 命中；
    · 任一单元 AUC < 0.5 → 证伪（方向反转）；
    · 两单元全 < 边际且无反转 → 证伪（全单元不达线，P3-F 字面）；
    · 其余（部分达线） → 未命中-中间态（按预注册字面不构成证伪，如实入账）。
    边际（预注册冻结）：P3-1=0.60，P3-2=0.55，P3-3=0.60。

护栏（沿用代码44 v0.1.1 先例）：缺文件 / 结构断言失败 / 任一条件组 0 个纳入细胞
→ 中止运行，不产出裁决表；空跑不构成裁决。

种子 20260814；bootstrap 2000 次（各组内逐细胞有放回重抽样），首轮即出 95% CI
（吸取代码44 CI 触发 bug 教训）。
运行环境：纯 CPU 秩统计（几千细胞 × 2000 次重排，分钟级），无矩阵并行负载，
          不需要 GPU；依赖 numpy/scipy/h5py（.mat 为 v7.3，务必先 pip install h5py）。

输入（在 ROOT 下递归按文件名定位，容忍中间目录层）：
    FixedSource_TNF100_DC3.mat / FixedSource_TNF30_DC3.mat / FixedSource_TNF10_DC3.mat
    DifferentDuration_TNF100_DC2.mat
输出（写入脚本所在目录）：
    代码45_判定日志.txt / 代码45_单元表.csv / 代码45_稳健臂.csv /
    代码45_描述臂.csv / 代码45_裁决.json
================================================================================
"""
import os, sys, json, csv, time, datetime
import numpy as np
from scipy.stats import rankdata

try:
    import h5py
except ImportError:
    print("[中止] 未安装 h5py。Sung 的 .mat 是 MATLAB v7.3 格式，scipy 读不了。")
    print("       请在 Spyder 控制台执行:  pip install h5py   然后重跑。")
    raise SystemExit(1)

# ----------------------------- 配置区 --------------------------------------
ROOT = r"C:/Users/lihua/Desktop/NC/data/Sung2022"   # ← 四个 .mat 放这里（可在子目录）

SEED   = 20260814
N_BOOT = 2000
N_MIN  = 30                 # 单元每组最少细胞数（护栏 D8）
K_MAIN = 3                  # 主分析阈值倍数
K_ROB  = (2, 4)             # 稳健臂
WINDOW = (0.0, 120.0)       # 分析窗 (0,120] min
EXTENT_MIN = 120.0          # 轨迹至少延伸到 120 min
MISS_MAX = 0.10             # 缺失率上限
MARGIN = {"P3-1": 0.60, "P3-2": 0.55, "P3-3": 0.60}
OVERLAP_BAND = (0.5, 1.5)   # P3-1 距离重叠带稳健臂（mm）

DOSE_FILES = {"100": "FixedSource_TNF100_DC3.mat",
              "30":  "FixedSource_TNF30_DC3.mat",
              "10":  "FixedSource_TNF10_DC3.mat"}
DUR_FILE = "DifferentDuration_TNF100_DC2.mat"

# 作者脚本逐字转录的时间对齐常数（DeepAnalysisForFig4.m / HeatMap&Traces_Dose.m）
OSD = {"100": [38.0, 148.0], "30": [39.0, 148.0], "10": [45.0, 152.0]}
OSI = {
 "100": [0,0,0,0,.7,.7,.7,.7,1.4,1.4,1.4,1.4,2.1,2.1,2.1,2.1,2.8,2.8,2.8,2.8,3.5,3.5,3.5,3.5],
 "30":  [0,0,0,0,.7,.7,.7,.7,1.4,1.4,1.4,1.4,3.5,3.5,3.5,3.5,2.8,2.8,2.8,2.8,2.1,2.1,2.1,2.1],
 "10":  [0,0,0,0,.7,.7,.7,.7,1.4,1.4,1.4,1.4,2.1,2.1,2.1,2.1,2.8,2.8,3.15,3.15,3.5,3.5,3.5,3.5],
}
OSD_10_ALT = [39.0, 148.0]          # Option B 稳健臂（作者脚本中的覆盖块）
OSI_10_ALT = OSI["30"]

CHAMBERS = {"100": list(range(1, 13)), "30": list(range(1, 13)),
            "10": [1, 2, 3, 4, 9, 10, 11, 12]}
DUR_GROUPS = {"15min": [2, 3, 4], "30min": [5, 6, 7, 8], "60min": [9, 10, 11, 12]}
ADJ_DUR = -42.0                      # 时长系列全局偏移（作者脚本 adj=-42）

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
    log("【中止】" + msg)
    log("本次运行为空跑，不构成裁决，不产出裁决表。请把本日志发回。")
    try:
        with open(os.path.join(OUT_DIR, "代码45_判定日志.txt"), "w", encoding="utf-8") as fh:
            fh.write("\n".join(_LOG_LINES))
    except Exception:
        pass
    raise SystemExit(2)

# ------------------------- v7.3 .mat 读取（h5py） ---------------------------
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
        abort(f"{tag}: 文件内没有变量 R（顶层键={list(f.keys())}）")
    R = f["R"]
    if not isinstance(R, h5py.Dataset) or R.ndim != 2 or R.shape[0] < 2:
        abort(f"{tag}: R 形状异常 {getattr(R, 'shape', None)}（预期 (≥2, 腔室数)）")
    # v0.1.1：实测 (5, 腔室数)——MATLAB R 为 (腔室数,5)，作者脚本仅用第1、2列
    log(f"  [{tag}] R 形状 {R.shape}（行0=腔室标签，行1=逐细胞块；其余行作者脚本未用）")
    return R

def get_interval(f, tag):
    if "interval" in f:
        try:
            val = float(np.asarray(f["interval"][()]).flatten()[0])
            log(f"  [{tag}] interval = {val} min（读自文件）")
            return val
        except Exception:
            pass
    log(f"  [{tag}] 文件内无 interval，按作者脚本取 6.0 min")
    return 6.0

def chamber_entry(f, R, aa):
    """返回 (label, cells_ref)；空腔室返回 (None, None)。aa 为 1 基腔室号。"""
    lref = R[0, aa - 1]
    cref = R[1, aa - 1]
    if not lref:
        return None, None
    label = read_label(f, lref)
    if not cref:
        return label, None
    return label, cref

def iter_cells(f, cells_ref, tag, aa):
    """逐细胞产出 (frames, values, dist_mm)。结构断言失败→中止（护栏）。"""
    S = f[cells_ref]
    if not isinstance(S, h5py.Dataset) or S.ndim != 2 or S.shape[0] < 3:
        abort(f"{tag} 腔室{aa}: 逐细胞块形状异常 {getattr(S,'shape',None)}"
              f"（预期 (≥3, 细胞数) 的引用矩阵）")
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
            abort(f"{tag} 腔室{aa} 细胞{j+1}: 轨迹块形状异常 {T.shape}"
                  f"（预期 (≥2, 帧数)）：第1行=帧号，第2行=NF-κB")
        frames = T[0, :]
        values = T[1, :]
        fin_f = frames[np.isfinite(frames)]
        if fin_f.size == 0 or fin_f.min() < 1:
            abort(f"{tag} 腔室{aa} 细胞{j+1}: 帧号异常（最小值 "
                  f"{fin_f.min() if fin_f.size else 'NaN'}，预期 ≥1）——结构与解码不符")
        dist = np.nan
        mref = S[2, j]
        if mref:
            M = np.asarray(f[mref][()], dtype=float)
            if M.ndim == 2 and M.shape[0] >= 4:
                drow = M[3, :]
                drow = drow[np.isfinite(drow)]
                if drow.size:
                    dist = float(drow.mean()) / 1000.0   # µm → mm（作者脚本同款）
        out.append((frames, values, dist))
    if n_skip:
        log(f"  [{tag}] 腔室{aa}: {n_skip} 个空细胞槽跳过")
    return out

# ------------------------- 纳入/剔除 + 特征 --------------------------------
def featurize(frames, values, dist, interval, offset):
    """返回 (状态, 特征dict)。状态: 'ok' 或剔除原因。全部规则见头部 D3/D4/D5。
    x = (帧号-1)*interval + offset（offset：剂量系列=腔室馈液偏移 os；时长系列=adj）。"""
    x = (frames - 1.0) * interval + offset
    v = values.astype(float)
    valid = np.isfinite(v) & (v >= 0) & np.isfinite(x)   # 帧号缺失同样计入缺失率
    if v.size == 0:
        return "空轨迹", None
    miss = 1.0 - float(valid.mean())
    if miss > MISS_MAX:
        return f"缺失>{int(MISS_MAX*100)}%", None
    order = np.argsort(x, kind="stable")
    x, v, valid = x[order], v[order], valid[order]
    xv, vv = x[valid], v[valid]
    base = vv[xv <= 0.0]
    if base.size < 2:
        return "基线帧<2", None
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
    feat["dur3"] = float(hit3.sum()) * interval          # 越阈总时长(min)
    if vw.size >= 2:
        feat["count3"] = int(np.sum((vw[:-1] < thr3) & (vw[1:] >= thr3)))
    else:
        feat["count3"] = 0
    feat["peak"] = float(vw.max()) if vw.size else np.nan
    feat["responded"] = bool(np.isfinite(feat[f"tau{K_MAIN}"]))
    return "ok", feat

# ------------------------- AUC 与 bootstrap --------------------------------
def auc_mw(a, b):
    """P(a>b)+0.5·P(a=b)，平均秩；NaN 丢弃，+inf 参与排名（代码44 同款口径）。"""
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

# ------------------------- 条款判定表（D8，写死） ---------------------------
def clause_verdict(clause, units):
    """units: list of dict(n1, n2, auc, margin)。返回 (判定, 说明)。"""
    margin = MARGIN[clause]
    if any(u["n1"] < N_MIN or u["n2"] < N_MIN for u in units):
        return "数据不足", "至少一单元 n<%d，不构成命中也不构成证伪" % N_MIN
    aucs = [u["auc"] for u in units]
    if all(a >= margin for a in aucs):
        return "命中", f"全部单元 ≥ {margin}"
    if any(a < 0.5 for a in aucs):
        return "证伪", "存在方向反转单元（AUC<0.5）"
    if all(a < margin for a in aucs):
        return "证伪", f"全单元不达线（均 < {margin}，P3-F 字面）"
    return "未命中-中间态", "部分达线但未全达线；按预注册字面不记证伪，如实入账"

# ------------------------- 组装条件组 --------------------------------------
def collect_group(path, tag, chambers, offset_fn, log_head=True):
    """offset_fn(label, aa) -> float。返回 list of feat dict（仅纳入细胞）。"""
    f = h5py.File(path, "r")
    interval = get_interval(f, tag)
    R = get_R(f, tag)
    feats, reasons = [], {}
    n_ch_used, n_raw = 0, 0
    for aa in chambers:
        if aa > R.shape[1]:
            log(f"  [{tag}] 腔室{aa} 超出 R 范围（共{R.shape[1]}腔室），跳过")
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
    log(f"  [{tag}] 腔室 {n_ch_used} 个，原始细胞 {n_raw}，纳入 {len(feats)}；"
        f"剔除明细 {reasons if reasons else '{}'}")
    return feats

def dose_offset_fn(osd, osi):
    def fn(label, aa):
        sec = 0 if "1-" in label else 1     # 作者脚本: contains(label,'1-')→osD(1)
        return -osd[sec] + osi[aa - 1]
    return fn

def dur_offset_fn(label, aa):
    return ADJ_DUR

# ------------------------- 主流程 -------------------------------------------
def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    log("=" * 74)
    log("代码45 — P3 盲裁决：免费统计量充分性（Sung 2022 NF-κB 空间梯度） v0.1.1")
    log("运行时间: " + datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    log(f"种子={SEED}  bootstrap={N_BOOT}  判定边际={MARGIN}  分析窗={WINDOW} min")
    log(f"numpy={np.__version__}  scipy 已载入  h5py={h5py.__version__}")
    log("预注册: 预注册_P3_免费统计量充分性_NFkB空间梯度_v01.md（注册·七十五）")
    log("=" * 74)

    # ---- 护栏1：文件齐备 ----
    paths = {}
    missing = []
    for key, fname in list(DOSE_FILES.items()) + [("DUR", DUR_FILE)]:
        hits = find_file(fname)
        if not hits:
            missing.append(fname)
        else:
            paths[key] = sorted(hits)[0]
    if missing:
        log("【缺文件】以下文件未在 ROOT 下找到：")
        for m in missing:
            log("   - " + m)
        log(f"ROOT = {ROOT}")
        log("请确认四个 .mat（注意压缩包内还各嵌套一层 zip，需再解一层）已放好。")
        abort("数据文件不齐")
    for k, p in paths.items():
        log(f"  文件[{k}] = {p}")

    # ---- 剂量系列 ----
    log("\n[装载] 剂量系列（FixedSource, 逐腔室馈液偏移 os=-osD+osI）")
    dose = {}
    for d in ("10", "30", "100"):
        dose[d] = collect_group(paths[d], f"TNF{d}",
                                CHAMBERS[d], dose_offset_fn(OSD[d], OSI[d]))
    # ---- 时长系列 ----
    log("\n[装载] 时长系列（DifferentDuration, adj=-42）")
    dur = {}
    for gname, chs in DUR_GROUPS.items():
        dur[gname] = collect_group(paths["DUR"], f"时长{gname}", chs, dur_offset_fn)

    # ---- 护栏2：空跑防护 ----
    empties = ([f"TNF{d}" for d in ("10", "30", "100") if len(dose[d]) == 0]
               + [f"时长{g}" for g in DUR_GROUPS if len(dur[g]) == 0])
    if empties:
        log("【空组】以下条件组纳入细胞为 0：" + ", ".join(empties))
        abort("存在空条件组（路径/结构需核查）")

    def col(group, key):
        return np.array([c[key] for c in group], float)

    # ---- P3-2 三分带（D6）----
    d100 = [c for c in dose["100"] if np.isfinite(c["dist"]) and c["dist"] > 0]
    log(f"\n[P3-2] TNF100 距离有效细胞 {len(d100)}/{len(dose['100'])}")
    if len(d100) == 0:
        abort("TNF100 无有效距离，P3-2 无法构造")
    q1, q2 = np.percentile([c["dist"] for c in d100], [100/3, 200/3])
    band = {"near": [c for c in d100 if c["dist"] <= q1],
            "mid":  [c for c in d100 if q1 < c["dist"] <= q2],
            "far":  [c for c in d100 if c["dist"] > q2]}
    log(f"  三分带边界: {q1:.4f} / {q2:.4f} mm；"
        f"近/中/远 = {len(band['near'])}/{len(band['mid'])}/{len(band['far'])}")
    if any(len(band[b]) == 0 for b in band):
        abort("存在空距离带")

    # ---- 6 个判定单元 ----
    units = [
        {"clause": "P3-1", "unit": "P3-1a", "contrast": "τ(3): TNF10 > TNF30",
         "a": col(dose["10"], f"tau{K_MAIN}"),  "b": col(dose["30"], f"tau{K_MAIN}")},
        {"clause": "P3-1", "unit": "P3-1b", "contrast": "τ(3): TNF30 > TNF100",
         "a": col(dose["30"], f"tau{K_MAIN}"),  "b": col(dose["100"], f"tau{K_MAIN}")},
        {"clause": "P3-2", "unit": "P3-2a", "contrast": "τ(3): 中带 > 近带",
         "a": col(band["mid"], f"tau{K_MAIN}"), "b": col(band["near"], f"tau{K_MAIN}")},
        {"clause": "P3-2", "unit": "P3-2b", "contrast": "τ(3): 远带 > 中带",
         "a": col(band["far"], f"tau{K_MAIN}"), "b": col(band["mid"], f"tau{K_MAIN}")},
        {"clause": "P3-3", "unit": "P3-3a", "contrast": "越阈时程: 30min > 15min",
         "a": col(dur["30min"], "dur3"),        "b": col(dur["15min"], "dur3")},
        {"clause": "P3-3", "unit": "P3-3b", "contrast": "越阈时程: 60min > 30min",
         "a": col(dur["60min"], "dur3"),        "b": col(dur["30min"], "dur3")},
    ]

    log("\n[判定] 6 单元（AUC 方向见头部 D7；CI=bootstrap 2000 百分位）")
    log(f"{'单元':<7}{'对比':<26}{'n_A':>6}{'n_B':>6}{'AUC':>8}{'CI低':>8}{'CI高':>8}{'边际':>6}  达线")
    rows = []
    for u in units:
        u["n1"], u["n2"] = int(u["a"].size), int(u["b"].size)
        if u["n1"] == 0 or u["n2"] == 0:
            abort(f"单元 {u['unit']} 空组")
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
            + ("✓" if ok else ("✗" if ok is False else "数据不足")))
        rows.append({"clause": u["clause"], "unit": u["unit"], "contrast": u["contrast"],
                     "n_a": u["n1"], "n_b": u["n2"], "auc": u["auc"],
                     "ci_lo": u["ci"][0], "ci_hi": u["ci"][1],
                     "margin": MARGIN[u["clause"]],
                     "pass": {True: "Y", False: "N", None: "INSUFFICIENT"}[ok]})

    # ---- 条款判定（D8）----
    log("\n[条款判定]")
    verdicts = {}
    for clause in ("P3-1", "P3-2", "P3-3"):
        us = [u for u in units if u["clause"] == clause]
        v, why = clause_verdict(clause, us)
        verdicts[clause] = {"verdict": v, "reason": why,
                            "units": [u["unit"] for u in us]}
        log(f"  {clause}: {v} —— {why}")
    any_fals = any(v["verdict"] == "证伪" for v in verdicts.values())
    all_hit = all(v["verdict"] == "命中" for v in verdicts.values())
    log("\n[P3-F 总证伪] " + ("触发：至少一条款证伪" if any_fals else
                              ("未触发；三条款全命中" if all_hit else "未触发；存在未命中/数据不足条款")))
    overall = ("P3 命中（三条款全达线）" if all_hit else
               "P3 证伪（存在证伪条款）" if any_fals else
               "P3 未命中-中间态/数据不足（按判定表如实入账）")
    log("[总裁决] " + overall)

    # ---- 稳健臂（点估计，不进判定）----
    log("\n[稳健臂] 点估计（不进判定）")
    rob_rows = []
    def rob(tag, a, b):
        v = auc_mw(a, b)
        log(f"  {tag:<44} n={len(a)}/{len(b)}  AUC={v:.4f}")
        rob_rows.append({"arm": tag, "n_a": len(a), "n_b": len(b), "auc": v})
        return v
    for k in K_ROB:
        rob(f"k={k} P3-1a τ: TNF10>TNF30", col(dose["10"], f"tau{k}"), col(dose["30"], f"tau{k}"))
        rob(f"k={k} P3-1b τ: TNF30>TNF100", col(dose["30"], f"tau{k}"), col(dose["100"], f"tau{k}"))
        rob(f"k={k} P3-2a τ: 中>近", col(band["mid"], f"tau{k}"), col(band["near"], f"tau{k}"))
        rob(f"k={k} P3-2b τ: 远>中", col(band["far"], f"tau{k}"), col(band["mid"], f"tau{k}"))
    lo_b, hi_b = OVERLAP_BAND
    ov = {d: [c for c in dose[d]
              if np.isfinite(c["dist"]) and lo_b <= c["dist"] <= hi_b] for d in ("10", "30", "100")}
    log(f"  距离重叠带 {OVERLAP_BAND} mm 内细胞数: "
        + ", ".join(f"TNF{d}={len(ov[d])}" for d in ("10", "30", "100")))
    if all(len(ov[d]) >= N_MIN for d in ("10", "30", "100")):
        rob(f"重叠带 P3-1a τ(3): TNF10>TNF30", col(ov["10"], f"tau{K_MAIN}"), col(ov["30"], f"tau{K_MAIN}"))
        rob(f"重叠带 P3-1b τ(3): TNF30>TNF100", col(ov["30"], f"tau{K_MAIN}"), col(ov["100"], f"tau{K_MAIN}"))
    else:
        log("  重叠带内某组 n<30，重叠带臂跳过（如实记录）")
        rob_rows.append({"arm": "重叠带 P3-1", "n_a": "", "n_b": "", "auc": "SKIPPED_n<30"})
    # Option B 对齐稳健臂（仅 P3-1，点估计）
    log("  [Option B] TNF10 改用被覆盖块 [39,148]+osI_30 重对齐（仅稳健臂）")
    dose10B = collect_group(paths["10"], "TNF10-OptB",
                            CHAMBERS["10"], dose_offset_fn(OSD_10_ALT, OSI_10_ALT))
    if len(dose10B) >= N_MIN:
        rob(f"OptB P3-1a τ(3): TNF10>TNF30", col(dose10B, f"tau{K_MAIN}"), col(dose["30"], f"tau{K_MAIN}"))
        rob(f"OptB P3-1b τ(3): TNF30>TNF100（对照）", col(dose["30"], f"tau{K_MAIN}"), col(dose["100"], f"tau{K_MAIN}"))
    else:
        log("  Option B 组 n<30，跳过（如实记录）")
        rob_rows.append({"arm": "OptB P3-1a", "n_a": "", "n_b": "", "auc": "SKIPPED_n<30"})

    # ---- 描述臂（不进判定）----
    log("\n[描述臂] 付费类对照（峰值幅度）+ 计数解码器 + 响应率（均不进判定）")
    desc_rows = []
    def desc(tag, a, b):
        v = auc_mw(a, b)
        log(f"  {tag:<44} n={len(a)}/{len(b)}  AUC={v:.4f}")
        desc_rows.append({"arm": tag, "n_a": len(a), "n_b": len(b), "auc": v})
    desc("峰值 P3-1a: TNF30>TNF10", col(dose["30"], "peak"), col(dose["10"], "peak"))
    desc("峰值 P3-1b: TNF100>TNF30", col(dose["100"], "peak"), col(dose["30"], "peak"))
    desc("峰值 P3-2a: 近>中", col(band["near"], "peak"), col(band["mid"], "peak"))
    desc("峰值 P3-2b: 中>远", col(band["mid"], "peak"), col(band["far"], "peak"))
    desc("峰值 P3-3a: 30min>15min", col(dur["30min"], "peak"), col(dur["15min"], "peak"))
    desc("峰值 P3-3b: 60min>30min", col(dur["60min"], "peak"), col(dur["30min"], "peak"))
    desc("计数 P3-1a: TNF30>TNF10", col(dose["30"], "count3"), col(dose["10"], "count3"))
    desc("计数 P3-1b: TNF100>TNF30", col(dose["100"], "count3"), col(dose["30"], "count3"))
    log("  响应率（τ(3) 未删失比例）：")
    for name, g in [("TNF10", dose["10"]), ("TNF30", dose["30"]), ("TNF100", dose["100"]),
                    ("近带", band["near"]), ("中带", band["mid"]), ("远带", band["far"]),
                    ("15min", dur["15min"]), ("30min", dur["30min"]), ("60min", dur["60min"])]:
        rr = np.mean([c["responded"] for c in g]) if g else np.nan
        log(f"    {name:<8} n={len(g):<5} 响应率={rr:.3f}")
        desc_rows.append({"arm": f"响应率 {name}", "n_a": len(g), "n_b": "", "auc": rr})

    # ---- 写输出 ----
    def write_csv(path, rows_, fields):
        with open(path, "w", newline="", encoding="utf-8-sig") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            for r in rows_:
                w.writerow(r)

    p_unit = os.path.join(OUT_DIR, "代码45_单元表.csv")
    p_rob  = os.path.join(OUT_DIR, "代码45_稳健臂.csv")
    p_desc = os.path.join(OUT_DIR, "代码45_描述臂.csv")
    p_json = os.path.join(OUT_DIR, "代码45_裁决.json")
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
    log(f"\n耗时 {time.time()-t0:.1f}s；输出：")
    for p in (p_unit, p_rob, p_desc, p_json, p_log):
        log("  " + p)
    with open(p_log, "w", encoding="utf-8") as fh:
        fh.write("\n".join(_LOG_LINES))

if __name__ == "__main__":
    main()
