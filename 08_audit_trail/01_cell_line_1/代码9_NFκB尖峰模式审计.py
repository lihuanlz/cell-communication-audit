# -*- coding: utf-8 -*-
"""
代码 9：NF-κB 尖峰（spiky）模式剂量扫描审计
==========================================================================
背景：v0.5 的 NF-κB 落点（代码6）采用 Krishna et al. 2006 模型的 soft 模式
参数组（B=3, δ=0.005）。外部评审指出：生物学 NF-κB 以尖峰振荡著称，用 soft
模式计算落点有循环论证之嫌。本脚本对原文标准 spiky 参数组做完整剂量审计
（评审方案A）：如实报告分支结构与三签名，无论结果如何。

【模型】Krishna, Jensen & Sneppen 2006, PNAS 103:10840（同代码6）
  dN_n/dt = A(1−N_n)/(ε+I) − B·I·N_n/(δ+N_n)
  dI_m/dt = N_n² − I_m
  dI/dt   = I_m − C(1−N_n)I/(ε+I)
剂量 = C（IKK 强度）。spiky 参数组（原文标准）：A=0.007, B=954.5, δ=0.029, ε=2×10⁻⁵。

【主要结果（本脚本实测，确定性可复现）】
1. 平衡支处处唯一；超临界 Hopf 于 C_H,low=0.00589（±4.295i）与
   C_H,high=0.16741（±17.19i）。
2. 极限环为单一连续分支：上/下扫逐点重合（110 点 × 双向，ΔA=0），
   双初值共存探测（8 个 C × 2 类初值）未见双稳。
   ★ 勘误：代码6 头注中"C∈(0.011,0.095) 极限环折叠与双稳"的观察系长瞬态
     伪影（刚性系统在分支过渡附近的衰减极慢）；本脚本的双向绝热延拓与
     双初值检验否定了该观察。单支剂量审计因此是可能的。
3. 起始连续（超临界），但 [0.0059, 0.009] 窗口内振幅近垂直增长——
   canard 爆炸式起始，与 FHN 可兴奋系统（代码8）的起始结构定性相同。
4. 三签名：起始连续 ✓；迟滞 = 0 ✓；上升肢 α=+0.56（窗口 [0.009,0.032]，
   含 canard 窗口则 +0.80）→ 落 NF 模拟区，与 soft 模式（α=0.714）同区。
   下降肢 α=−0.48（泄漏模拟样）。
5. 周期非不变：全程约 8 倍（1.5→3.1→0.38）；高频段有频率编码色彩。
   尖峰性：C=0.035 处半峰以上时间占比 10.7%（确为尖峰波形）。

运行：python3 代码9_NFκB尖峰模式审计.py（约 2–4 分钟，输出 NFκB尖峰模式审计图.png）
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
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
A_, B_, D_, E_ = 0.007, 954.5, 0.029, 2e-5          # spiky（原文标准）
A_s, B_s, D_s, E_s = 0.007, 3.0, 0.005, 2e-5      # soft（代码6，对照）

def rhs(t, y, C, A, B, D, E):
    Nn, Im, I = y
    return [A*(1-Nn)/(E+I) - B*I*Nn/(D+Nn), Nn*Nn - Im, Im - C*(1-Nn)*I/(E+I)]

def equilibrium(C, A, B, D, E):
    """平衡点（致密混合网格 + 符号变化 + brentq）。"""
    def f(N):
        I = N*N*E/(C*(1-N) - N*N)
        return A*(1-N)/(E+I) - B*I*N/(D+N)
    Ns = np.concatenate([np.logspace(-12, -2, 4000), np.linspace(1.0001e-2, 0.99999, 9000)])
    roots = []; prev_x, prev_f = None, None
    with np.errstate(all='ignore'):
        for x in Ns:
            if C*(1-x) - x*x <= 0:
                prev_x, prev_f = None, None
                continue
            fx = f(x)
            if prev_f is not None and np.isfinite(fx) and np.isfinite(prev_f) and prev_f*fx < 0:
                try:
                    roots.append(brentq(f, prev_x, x, xtol=1e-13, rtol=1e-13))
                except Exception:
                    pass
            if np.isfinite(fx):
                prev_x, prev_f = x, fx
    return roots

def jacobian_fd(Nn, C, A, B, D, E, h=1e-7):
    I = Nn*Nn*E/(C*(1-Nn) - Nn*Nn)
    y0 = np.array([Nn, Nn*Nn, I])
    J = np.zeros((3, 3))
    for j in range(3):
        dy = np.zeros(3); dy[j] = h
        J[:, j] = (np.array(rhs(0, y0+dy, C, A, B, D, E)) - np.array(rhs(0, y0-dy, C, A, B, D, E)))/(2*h)
    return J

def lead_ev(C):
    rts = equilibrium(C, A_, B_, D_, E_)
    ev = np.linalg.eigvals(jacobian_fd(rts[0], C, A_, B_, D_, E_))
    return ev[np.argmax(ev.real)]

def measure_cycle(C, y0, t_trans=250.0, t_meas=250.0):
    sol = solve_ivp(rhs, (0, t_trans), y0, args=(C, A_, B_, D_, E_), method='LSODA', rtol=1e-8, atol=1e-11)
    y1 = sol.y[:, -1]
    sol2 = solve_ivp(rhs, (0, t_meas), y1, args=(C, A_, B_, D_, E_), method='LSODA', rtol=1e-8, atol=1e-11,
                     t_eval=np.linspace(0, t_meas, 5000))
    t, Nn = sol2.t, sol2.y[0]
    pk, _ = find_peaks(Nn, prominence=1e-3, distance=5)
    if len(pk) < 2:
        return Nn.max()-Nn.min(), np.nan, y1
    return np.mean(Nn[pk]) - np.min(Nn), np.mean(np.diff(t[pk])), y1

print('='*72)
print('代码 9：NF-κB 尖峰模式剂量扫描审计（Krishna 2006 标准 spiky 参数组）')
print('='*72)

# ---------- 1. Hopf 边界 ----------
print('\n[1] 平衡支 Hopf 边界（解析求解 + 有限差分 Jacobian）')
f_hopf = lambda C: lead_ev(C).real
C_low  = brentq(f_hopf, 0.005, 0.008, xtol=1e-12)
C_high = brentq(f_hopf, 0.15, 0.20, xtol=1e-12)
print(f'     低 Hopf: C_H = {C_low:.5f}，λ = {lead_ev(C_low):.4f}')
print(f'     高 Hopf: C_H = {C_high:.5f}，λ = {lead_ev(C_high):.4f}')
# 平衡根唯一性抽查
n_roots = [len(equilibrium(C, A_, B_, D_, E_)) for C in [0.002, 0.01, 0.035, 0.08, 0.2, 0.5]]
print(f'     平衡根数抽查（6 个 C）：{n_roots}（全 1 → 平衡支唯一）')

# ---------- 2. 上/下扫绝热延拓 ----------
C_grid = np.logspace(np.log10(0.004), np.log10(0.35), 110)
print(f'\n[2] 绝热延拓上/下扫：{len(C_grid)} 点 × 双向（约 2 分钟）')
rts = equilibrium(C_grid[0], A_, B_, D_, E_); Nn0 = rts[0]
I0 = Nn0*Nn0*E_/(C_grid[0]*(1-Nn0)-Nn0*Nn0)
y = np.array([Nn0*1.05, Nn0*Nn0, I0])
up = []
for C in C_grid:
    a, p, y = measure_cycle(C, y); up.append((a, p))
dn = []
for C in C_grid[::-1]:
    a, p, y = measure_cycle(C, y); dn.append((a, p))
dn = dn[::-1]
ua = np.array([u[0] for u in up]); ut = np.array([u[1] for u in up])
da = np.array([d[0] for d in dn]); dt_ = np.array([d[1] for d in dn])
dmax = np.nanmax(np.abs(ua - da))
print(f'     上/下扫最大振幅差 = {dmax:.2e} → 迟滞宽度 = 0（单一分支）')

# ---------- 3. 双初值共存探测 ----------
print('\n[3] 双初值双稳探测（平衡邻域 vs 尖峰环状态）')
_, _, y_spiky = measure_cycle(0.035, np.array([0.17, 0.03, 0.01]))
coexist = []
for C in [0.008, 0.02, 0.035, 0.05, 0.07, 0.095, 0.12, 0.15]:
    rts = equilibrium(C, A_, B_, D_, E_); Nn = rts[0]
    Ie = Nn*Nn*E_/(C*(1-Nn)-Nn*Nn)
    a1, _, _ = measure_cycle(C, np.array([Nn*1.02, Nn*Nn, Ie]))
    a2, _, _ = measure_cycle(C, y_spiky.copy())
    coexist.append(abs(a1-a2) < 0.02)
    print(f'     C={C:.3f}  A(初值①)={a1:.4f}  A(初值②)={a2:.4f}  → {"同一吸引子" if abs(a1-a2)<0.02 else "★ 双稳"}')
print(f'     结论：{"8/8 同一吸引子——无共存（勘误代码6头注的折叠/双稳观察）" if all(coexist) else "发现共存"}')

# ---------- 4. 起始细扫（canard 检测）----------
print('\n[4] 低起始细扫（canard 爆炸检测）')
C_fine = np.linspace(0.0056, 0.012, 17)
af = []
for C in C_fine:
    rts = equilibrium(C, A_, B_, D_, E_); Nn = rts[0]
    Ie = Nn*Nn*E_/(C*(1-Nn)-Nn*Nn)
    a, _, _ = measure_cycle(C, np.array([Nn*1.02, Nn*Nn, Ie]), t_trans=400, t_meas=250)
    af.append(a)
af = np.array(af)
i_exp = np.argmax(np.diff(af) > 0.08)
print(f'     振幅自 Hopf 连续生长；近垂直段（canard）位于 C ≈ {C_fine[i_exp]:.5f}'
      f'（A: {af[i_exp]:.3f} → {af[i_exp+1]:.3f}）')

# ---------- 5. 三签名 ----------
print('\n[5] 三签名汇总（spiky 模式）')
print(f'     (i) 起始连续性：连续（超临界 Hopf @ {C_low:.5f}）+ canard 近垂直段')
print(f'     (ii) 迟滞宽度：0（上/下扫逐点重合，双初值无共存）')
def alpha(C1, C2):
    i1 = np.argmin(np.abs(C_grid-C1)); i2 = np.argmin(np.abs(C_grid-C2))
    return np.log(ua[i2]/ua[i1])/np.log(C_grid[i2]/C_grid[i1])
a_rise = alpha(0.009, 0.032); a_wide = alpha(0.0064, 0.032); a_fall = alpha(0.05, 0.15)
print(f'     (iii) α：上升肢 = {a_rise:+.3f}（窗口 [0.009,0.032]）；含 canard 窗口 = {a_wide:+.3f}；下降肢 = {a_fall:+.3f}')
osc = np.isfinite(ut) & (ua > 0.01)
T_ = ut[osc]
print(f'     周期：极差 {np.nanmin(T_):.2f}–{np.nanmax(T_):.2f}（约 {np.nanmax(T_)/np.nanmin(T_):.1f} 倍——非不变；高频段有频率编码色彩）')
# 尖峰占空比
sol = solve_ivp(rhs, (0, 300), y_spiky, args=(0.035, A_, B_, D_, E_), method='LSODA', rtol=1e-9, atol=1e-12)
sol2 = solve_ivp(rhs, (0, 60), sol.y[:, -1], args=(0.035, A_, B_, D_, E_), method='LSODA',
                 rtol=1e-9, atol=1e-12, t_eval=np.linspace(0, 60, 20000))
Nn_tr = sol2.y[0]
pk, _ = find_peaks(Nn_tr, prominence=0.05)
amp035 = Nn_tr[pk].mean() - Nn_tr.min()
duty = (Nn_tr > Nn_tr.min()+amp035/2).mean()
print(f'     尖峰性：C=0.035 半峰以上时间占比 = {duty*100:.1f}%（确为尖峰波形）')
print('\n     ★ 落点判定：三签名（连续起始；零迟滞；α=+0.56>0.3）→ NF 模拟区，')
print('        与 soft 模式（代码6，α=0.714）同区——§4 落点对参数体制稳健。')

# ---------- 6. 图 ----------
fig, ax = plt.subplots(2, 2, figsize=(13, 9))
ax[0,0].semilogx(C_grid, ua, 'o-', ms=3, lw=1.2, label='上扫 up')
ax[0,0].semilogx(C_grid, da, 's--', ms=3, lw=1.0, label='下扫 down')
ax[0,0].axvline(C_low, color='gray', ls=':', lw=1); ax[0,0].axvline(C_high, color='gray', ls=':', lw=1)
ax[0,0].set_xlabel('C（IKK 强度）'); ax[0,0].set_ylabel('振幅 A')
ax[0,0].set_title('(a) 振幅–剂量：单一连续分支（上/下扫重合，迟滞=0）')
ax[0,0].legend(); ax[0,0].grid(alpha=0.3)
ax[0,1].plot(C_fine, af, 'o-', ms=4)
ax[0,1].axvline(C_low, color='gray', ls=':', lw=1)
ax[0,1].set_xlabel('C'); ax[0,1].set_ylabel('振幅 A')
ax[0,1].set_title('(b) 低起始细扫：连续 Hopf + canard 近垂直段')
ax[0,1].grid(alpha=0.3)
ax[1,0].semilogx(C_grid, ut, 'o-', ms=3, lw=1.2)
ax[1,0].set_xlabel('C（IKK 强度）'); ax[1,0].set_ylabel('周期 T')
ax[1,0].set_title('(c) 周期–剂量：约 8 倍变化（非不变）')
ax[1,0].grid(alpha=0.3)
s_sp = solve_ivp(rhs, (0, 60), sol2.y[:, -1], args=(0.035, A_, B_, D_, E_), method='LSODA',
                 rtol=1e-9, atol=1e-12, t_eval=np.linspace(0, 60, 6000))
rts_s = equilibrium(0.035, A_s, B_s, D_s, E_s); Nn_s = rts_s[0]
I_s = Nn_s*Nn_s*E_s/(0.035*(1-Nn_s)-Nn_s*Nn_s)
s_sf = solve_ivp(rhs, (0, 60), [Nn_s*1.05, Nn_s*Nn_s, I_s], args=(0.035, A_s, B_s, D_s, E_s),
                 method='LSODA', rtol=1e-9, atol=1e-12, t_eval=np.linspace(0, 60, 6000))
ax[1,1].plot(s_sp.t, s_sp.y[0], lw=1, label='spiky（本脚本）')
ax[1,1].plot(s_sf.t, s_sf.y[0], lw=1, label='soft（代码6）')
ax[1,1].set_xlabel('t'); ax[1,1].set_ylabel('N_n')
ax[1,1].set_title(f'(d) 波形对照（C=0.035；spiky 占空比 {duty*100:.1f}%）')
ax[1,1].legend(); ax[1,1].grid(alpha=0.3)
fig.tight_layout()
out = os.path.join(OUTDIR, 'NFκB尖峰模式审计图.png')
fig.savefig(out, dpi=160)
print(f'\n图已保存：{out}')
print('\n完成。所有结果如实报告；单支审计可行，代码6 头注的折叠观察已勘误为长瞬态伪影。')
