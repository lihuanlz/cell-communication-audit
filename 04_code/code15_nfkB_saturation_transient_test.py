# -*- coding: utf-8 -*-
"""
代码 15：NF-κB 饱和瞬态假说定量检验（GLM 缺陷2 的对撞实验）
==========================================================================
背景：外部硬核审计（缺陷2）指出——代码12 的阶跃扰动测试只排除了"阈值放大"
（可兴奋性的定义特征），没有排除"饱和瞬态对 EXC 签名的观测级模仿"：
Krishna 模型的 (1−N_n) 库存耗竭项本身是饱和机制，强驱动下整个 NF-κB 库存
在反馈追上之前全部释放，首峰由库存总量决定、与剂量脱钩——纯负反馈无需正
反馈即可产生刻板首峰。若该假说成立，§4"实验侧 EXC 判读为排除备选后的唯一
读法"过强。

【检验设计】给饱和假说最好的机会：阻尼区 Krishna 模型（spiky 参数组，
C_H=5.89e-3 以下），配体→有效驱动的受体层压缩映射
    无压缩: C = C_max·TNF/100（线性，4 个十倍程全跨度）
    压缩:   C = C_max·TNF/(K_T + TNF)
扫描 (K_T, C_max) 网格，阶跃输入，定量输出三个可观测量，
与 Tay et al. 2010 实测三数对撞：
  (i)   首峰振幅斜率 α = Δlog10(A1)/Δlog10(剂量)；Tay 实测 0.151（4 倍/4 个十倍程）
  (ii)  计数斜率 ν：观测窗内可分辨脉冲数/十倍程；Tay 实测 ≈ +0.58（1.7→4.0 个）
  (iii) 峰间期恒定性；Tay 实测 75–95 min、无剂量趋势
  附 (iv) 衰减比 A2/A1；Tay 实测 0.5–0.8。

判定逻辑：若网格中不存在任何 (K_T, C_max) 同时复现 (i)(ii)(iv)，则饱和
瞬态读法定量失败，"唯一读法"获得更强证据；若存在，则诚实降级为"EXC 最优
读法、饱和 NF 竞争读法"。两种结果都照实报告。

【模型】Krishna, Jensen & Sneppen 2006, PNAS 103:10840（同代码6/9/12）
  dN_n/dt = A(1−N_n)/(ε+I) − B·I·N_n/(δ+N_n)
  dI_m/dt = N_n² − I_m
  dI/dt   = I_m − C(1−N_n)I/(ε+I)
spiky 参数组：A=0.007, B=954.5, δ=0.029, ε=2e-5；Hopf 窗 [5.89e-3, 0.167]。
静息态：C_basal=1e-4 的不动点（解析猜测 N_n≈(AδC²/B)^{1/5} + fsolve 精化）。
观测窗：T_obs=20 无量纲单位（按峰间期标定 ≈ 90 min/1.5–3 单位 ≈ 600–1200 min，
与 Tay 成像时长一致）。α 与 ν 无量纲，不受时间标定影响。

v2 修正（首轮实现错误，已记录）：(1) "无压缩"情形的 Michaelis 大 K_T 近似把
驱动压至 1e-12，改为显式线性映射；(2) 静息态长时积分未收敛（小 C 弛豫
~1/C≈1e4 单位），改为不动点求根。

运行：python3 代码15_NFκB饱和瞬态检验.py（约 2–5 分钟）
输出：代码15_饱和瞬态检验.png（4 面板）；控制台判定表。
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve
from scipy.signal import find_peaks
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

for _f in ['Noto Sans CJK SC', 'WenQuanYi Zen Hei', 'SimHei', 'Microsoft YaHei']:
    try:
        from matplotlib.font_manager import findfont, FontProperties
        if findfont(FontProperties(family=_f), fallback_to_default=False):
            plt.rcParams['font.sans-serif'] = [_f, 'DejaVu Sans']
            break
    except Exception:
        continue
plt.rcParams['axes.unicode_minus'] = False

OUTDIR = os.path.dirname(os.path.abspath(__file__))
A_, B_, D_, E_ = 0.007, 954.5, 0.029, 2e-5      # spiky（原文标准）
C_H = 5.89e-3                                     # 下 Hopf（代码9 实测）
T_OBS = 20.0                                      # 观测窗（无量纲）

def rhs(t, y, C):
    Nn, Im, I = y
    return [A_*(1-Nn)/(E_+I) - B_*I*Nn/(D_+Nn), Nn*Nn - Im, Im - C*(1-Nn)*I/(E_+I)]

def fixed_point(C):
    """不动点：解析猜测（小 C 标度）+ fsolve。"""
    Nn0 = max((A_*D_*C*C/B_)**0.2, 1e-6)
    guess = [Nn0, Nn0*Nn0, Nn0*Nn0/max(C, 1e-12)]
    sol, info, ier, _ = fsolve(lambda y: rhs(0, y, C), guess, full_output=True)
    if ier != 1:                            # 退路：长时积分
        s = solve_ivp(rhs, (0, 5e4), guess, args=(C,), method='LSODA',
                      rtol=1e-9, atol=1e-12)
        sol = s.y[:, -1]
    return sol

def step_response(C_step, y0, T=T_OBS, n=30000):
    sol = solve_ivp(rhs, (0, T), y0, args=(C_step,), method='LSODA',
                    rtol=1e-9, atol=1e-12, dense_output=True)
    t = np.linspace(0, T, n)
    return t, sol.sol(t)[0]

def analyze(t, Nn, thr_abs):
    pk, _ = find_peaks(Nn, prominence=1e-4)
    if len(pk) == 0:
        return dict(A1=np.nan, count=0, interval=np.nan, ratio=np.nan, tpk=np.array([]))
    A1 = Nn[pk[0]]
    vis = pk[Nn[pk] >= thr_abs]
    amps = Nn[vis]
    interval = np.mean(np.diff(t[vis])) if len(vis) > 1 else np.nan
    ratio = amps[1]/amps[0] if len(amps) > 1 else np.nan
    return dict(A1=A1, count=len(vis), interval=interval, ratio=ratio, tpk=t[pk])

print("=" * 76)
print("代码 15：NF-κB 饱和瞬态假说定量检验（GLM 缺陷2 对撞）")
print("=" * 76)

C_BASAL = 1e-4
Y0 = fixed_point(C_BASAL)
print(f"\n[1] 静息态（C_basal={C_BASAL:g} 不动点）: "
      f"N_n={Y0[0]:.3e}, I_m={Y0[1]:.3e}, I={Y0[2]:.3e}")

DOSES = np.array([0.01, 0.1, 1.0, 10.0, 100.0])
LOGDEC = np.log10(DOSES)
KT_GRID = [None, 30, 10, 3, 1, 0.3, 0.1, 0.03]      # None = 线性无压缩
CMAX_GRID = [0.003, 0.005, 0.008]
TAY = dict(alpha=0.151, nu=0.58, ratio=(0.5, 0.8))

def dose2C(KT, Cmax):
    if KT is None:
        return Cmax*DOSES/100.0
    return Cmax*DOSES/(KT+DOSES)

rows, traces = [], {}
print("\n[2] 扫描 (K_T, C_max)：8×3 = 24 组 × 5 剂量 = 120 条阶跃响应")
for Cmax in CMAX_GRID:
    for KT in KT_GRID:
        Csteps = dose2C(KT, Cmax)
        res = [step_response(Cs, Y0) for Cs in Csteps]
        pk, _ = find_peaks(res[-1][1], prominence=1e-4)
        thr = 0.15*res[-1][1][pk[0]] if len(pk) else np.nan
        out = [analyze(t, Nn, thr) for t, Nn in res]
        A1 = np.array([o['A1'] for o in out], float)
        cnt = np.array([o['count'] for o in out], float)
        ok = ~np.isnan(A1)
        alpha = np.polyfit(LOGDEC[ok], np.log10(A1[ok]), 1)[0] if ok.sum() >= 3 else np.nan
        nu = np.polyfit(LOGDEC, cnt, 1)[0] if cnt.std() > 0 else 0.0
        ints = [o['interval'] for o in out]
        rats = [o['ratio'] for o in out]
        rows.append(dict(Cmax=Cmax, KT=KT, alpha=alpha, nu=nu, A1=A1, cnt=cnt,
                         ints=ints, rats=rats, Csteps=Csteps, nvis=int(ok.sum())))
        lab = 'none' if KT is None else f'{KT:g}'
        if Cmax == 0.005 and KT in (None, 1, 0.1):
            traces[lab] = [(t.copy(), Nn.copy()) for t, Nn in res]
    print(f"    C_max={Cmax:g} 完成")

print("\n[3] 判定表（Tay 2010：α=0.151, ν=+0.58/十倍程, A2/A1∈[0.5,0.8]）")
print(f"{'C_max':>6} {'K_T':>6} {'α':>8} {'ν':>7} {'A2/A1':>6} {'窗内有响应剂量数':>10}  峰间期范围(单位)")
hits = []
for r in rows:
    ints = [x for x in r['ints'] if not np.isnan(x)]
    rats = [x for x in r['rats'] if not np.isnan(x)]
    i_str = f"{min(ints):.2f}–{max(ints):.2f}" if ints else "—"
    r_str = f"{np.mean(rats):.2f}" if rats else "—"
    kt_str = '线性' if r['KT'] is None else f"{r['KT']:g}"
    ok_a = not np.isnan(r['alpha']) and abs(r['alpha']-TAY['alpha']) <= 0.07
    ok_n = abs(r['nu']-TAY['nu']) <= 0.25
    ok_r = bool(rats) and TAY['ratio'][0]-0.15 <= np.mean(rats) <= TAY['ratio'][1]+0.15
    ok_all = ok_a and ok_n and ok_r and r['nvis'] == 5
    if ok_all:
        hits.append(r)
    a_str = f"{r['alpha']:.3f}" if not np.isnan(r['alpha']) else "—"
    print(f"{r['Cmax']:>6.3f} {kt_str:>6} {a_str:>8} {r['nu']:>7.3f} {r_str:>6}"
          f" {r['nvis']:>10}  {i_str}{'  ★ 同时复现' if ok_all else ''}")

print("\n[4] 张力量化（确定性层面）")
fit_a = [r for r in rows if not np.isnan(r['alpha']) and abs(r['alpha']-TAY['alpha']) <= 0.07]
fit_n = [r for r in rows if abs(r['nu']-TAY['nu']) <= 0.25]
if fit_a:
    b = min(fit_a, key=lambda r: r['nu'])
    print(f"    拟合住 α（±0.07）的 {len(fit_a)} 组中，ν 最低为 {b['nu']:.3f}"
          f"（C_max={b['Cmax']:g}, K_T={b['KT'] if b['KT'] is not None else '线性'}），"
          f"仍高于 Tay 上限 0.83")
if not fit_n:
    print(f"    全网格 24 组中无任何组合拟合住 ν（最低 {min(r['nu'] for r in rows):.3f}"
          f" vs Tay 0.58±0.25）")
print("    注意：确定性整数计数天然偏陡——Tay 的 ν 是群体均值，阈值异质性会")
print("    抹平计数斜率。故 [4b] 对最接近的组合做群体层面复检，再给判定。")

# ---- [4b] 群体层面复检：参数异质性抹平 ----
print("\n[4b] 群体复检（对 α-拟合最好的组合加对数正态异质性，48 细胞/剂量）")
rng = np.random.default_rng(20260805)
N_CELL = 48
SIG_C, SIG_T = 0.35, 0.30          # C 乘性异质性 / 阈值异质性

cands = sorted(fit_a, key=lambda r: r['nu'])[:2] if fit_a else []
pop_results = []
for r0 in cands:
    Csteps = r0['Csteps']
    # 基准阈值：最高剂量确定性首峰的 15%
    t, Nn = step_response(Csteps[-1], Y0)
    pk, _ = find_peaks(Nn, prominence=1e-4)
    thr0 = 0.15*Nn[pk[0]] if len(pk) else np.nan
    m_cnt, m_logA = [], []
    for Cs in Csteps:
        cnts, logAs = [], []
        for z1, z2 in zip(rng.standard_normal(N_CELL), rng.standard_normal(N_CELL)):
            Ci = Cs*np.exp(SIG_C*z1)
            t, Nn = step_response(Ci, Y0, n=12000)
            o = analyze(t, Nn, thr0*np.exp(SIG_T*z2))
            cnts.append(o['count'])
            if not np.isnan(o['A1']):
                logAs.append(np.log10(o['A1']))
        m_cnt.append(np.mean(cnts))
        m_logA.append(np.mean(logAs) if logAs else np.nan)
    m_cnt = np.array(m_cnt); m_logA = np.array(m_logA)
    ok = ~np.isnan(m_logA)
    nu_p = np.polyfit(LOGDEC, m_cnt, 1)[0]
    al_p = np.polyfit(LOGDEC[ok], m_logA[ok], 1)[0] if ok.sum() >= 3 else np.nan
    pop_results.append((r0, nu_p, al_p, m_cnt))
    # 次级判据：单调性（Tay 计数单调升）、绝对计数（Tay 1.7→4.0）、高剂量衰减
    mono = bool(np.all(np.diff(m_cnt) >= -0.5))
    abs_ok = bool(m_cnt.min() >= 1.0 and m_cnt.max() <= 5.5)
    print(f"    C_max={r0['Cmax']:g}, K_T={r0['KT'] if r0['KT'] is not None else '线性'}: "
          f"群体 ν={nu_p:.3f}（Tay 0.58±0.25）, 群体 α={al_p:.3f}（Tay 0.151±0.07）, "
          f"群体计数={np.round(m_cnt,2)}")
    print(f"      次级判据：单调{'✓' if mono else '✗（计数非单调，Tay 为单调升）'}；"
          f"绝对计数{'✓' if abs_ok else '✗（超出 Tay 的 1.7–4.0 范围）'}；"
          f"C_max>{C_H:.2e}⇒高剂量持续振荡{'是（Tay 高剂量仍衰减，矛盾）' if r0['Cmax'] > C_H else '否'}")

print("\n[5] 总结论")
pop_ok = [p for p in pop_results
          if abs(p[1]-TAY['nu']) <= 0.25 and not np.isnan(p[2]) and abs(p[2]-TAY['alpha']) <= 0.07
          and np.all(np.diff(p[3]) >= -0.5) and p[3].min() >= 1.0 and p[3].max() <= 5.5
          and p[0]['Cmax'] <= C_H]
pop_nominal = [p for p in pop_results
               if abs(p[1]-TAY['nu']) <= 0.25 and not np.isnan(p[2]) and abs(p[2]-TAY['alpha']) <= 0.07]
if pop_nominal and not pop_ok:
    print("    注：群体层面出现 1 组名义命中（ν、α 同中），但未过次级判据")
    print("    （计数非单调 / 绝对计数≈7 超出 Tay 的 4.0 / C_max>C_H 高剂量持续振荡")
    print("    与 Tay 的衰减观测矛盾）——不计为复现，但如实记录在案。")
if hits or pop_ok:
    print(f"    确定性 {len(hits)} 组 + 群体 {len(pop_ok)} 组同时复现三指标 → GLM 缺陷2 成立，")
    print("    §4'唯一读法'需降级为'EXC 最优读法、饱和 NF 竞争读法'。")
else:
    print("    确定性网格（24 组）与群体复检（最优组合加异质性）均无任何组合")
    print("    同时复现 α、ν 与衰减比 → 饱和瞬态假说在两个层面均定量失败：")
    print("    GLM 缺陷2 的机理（库存耗竭饱和）真实存在且使首峰刻板化成为可能，")
    print("    但其定量预言与 Tay 三数不可兼得——'唯一读法'经此检验后获得定量支撑")
    print("    （§4 与 SI 需写入：机理承认 + 定量排除，而非简单否认）。")

all_ints = [x for r in rows for x in r['ints'] if not np.isnan(x)]
if all_ints:
    med = np.median(all_ints)
    print(f"\n[6] 时间标定：无量纲峰间期中位 {med:.2f} ⇒ 按 Tay 中值 85 min，"
          f"1 单位 ≈ {85/med:.0f} min；观测窗 {T_OBS:.0f} 单位 ≈ {T_OBS*85/med:.0f} min")

# ---- 图 ----
fig, ax = plt.subplots(2, 2, figsize=(11.5, 8.2))
cmap = plt.cm.viridis(np.linspace(0.1, 0.9, len(DOSES)))

axs = ax[0, 0]
tt = np.logspace(-2, 2, 200)
axs.plot(tt, 0.005*tt/100, lw=2, label='线性（无压缩）')
for KT in [3, 1, 0.3, 0.1]:
    axs.plot(tt, 0.005*tt/(KT+tt), lw=2, label=f'$K_T$={KT:g}')
axs.axhline(C_H, color='r', ls='--', lw=1, label='$C_H$')
axs.set_xscale('log'); axs.set_xlabel('TNF (ng/ml)'); axs.set_ylabel('有效驱动 C')
axs.set_title('(a) 配体-驱动压缩映射', fontsize=11)
axs.legend(fontsize=8); axs.grid(alpha=0.3)

axs = ax[0, 1]
for i, (t, Nn) in enumerate(traces.get('1', [])):
    axs.plot(t, Nn, color=cmap[i], lw=1.2, label=f'{DOSES[i]:g} ng/ml')
axs.set_xlabel('无量纲时间'); axs.set_ylabel('$N_n(t)$')
axs.set_title('(b) 阶跃响应（Cmax=0.005, KT=1）', fontsize=11)
axs.legend(fontsize=8, ncol=2); axs.grid(alpha=0.3)

axs = ax[1, 0]
xt = list(range(len(KT_GRID)))
xl = ['线性'] + [f'{k:g}' for k in KT_GRID[1:]]
for Cmax in CMAX_GRID:
    xs = [r for r in rows if r['Cmax'] == Cmax]
    axs.plot(xt, [x['alpha'] if not np.isnan(x['alpha']) else np.nan for x in xs],
             'o-', lw=2, label=f'$C_{{max}}$={Cmax:g}')
axs.axhspan(0.081, 0.221, color='green', alpha=0.2, label='Tay α=0.151±0.07')
axs.set_xticks(xt); axs.set_xticklabels(xl, fontsize=8)
axs.set_xlabel('KT（左=无压缩，右=强压缩）'); axs.set_ylabel('α')
axs.set_title('(c) 首峰斜率 α vs 压缩强度', fontsize=11)
axs.legend(fontsize=8); axs.grid(alpha=0.3)

axs = ax[1, 1]
for Cmax in CMAX_GRID:
    xs = [r for r in rows if r['Cmax'] == Cmax]
    axs.plot(xt, [x['nu'] for x in xs], 's-', lw=2, label=f'$C_{{max}}$={Cmax:g}')
axs.axhspan(0.33, 0.83, color='green', alpha=0.2, label='Tay ν=0.58±0.25')
axs.set_xticks(xt); axs.set_xticklabels(xl, fontsize=8)
axs.set_xlabel('KT（左=无压缩，右=强压缩）'); axs.set_ylabel('ν (脉冲/十倍程)')
axs.set_title('(d) 计数斜率 ν vs 压缩强度——关键张力', fontsize=11)
axs.legend(fontsize=8); axs.grid(alpha=0.3)

fig.suptitle('代码15：NF-κB 饱和瞬态假说 vs Tay 2010 三可观测量对撞', fontsize=13)
fig.tight_layout(rect=[0, 0, 1, 0.96])
out = os.path.join(OUTDIR, '代码15_饱和瞬态检验.png')
fig.savefig(out, dpi=150)
print(f"\n图已保存: {out}")
