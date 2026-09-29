# -*- coding: utf-8 -*-
"""
代码 10：非对称米氏常数（K_m,1 ≠ K_m,2）下的 PdPC 审计
==========================================================================
背景：定理 1–3 在对称 K_m 假设下陈述（手稿 §3 范围注记："非对称见 SI"）。
外部评审断言"非对称时简并群降至 2 维"——本脚本先验证再落笔，结论：
  ★ 评审断言不成立。非对称情形参数空间 6→7 维（Km1, Km2 分列），
    可辨识组合 2→3 维（ξ, κ1, κ2），简并群维数保持 4（7−3）。
  ★ n_H 存在精确闭式推广（非近似）：
        n_H(κ1, κ2) = 4 / (4 − 1/(κ1+½) − 1/(κ2+½))
    对称退化：κ1=κ2=κ → 4/(4 − 2/(κ+½)) = 1 + 1/(2κ) ✓
  ★ Fisher 壁垒（形状可辨识 vs 绝对标度不可辨识）在全部不对称度下保持。

【模型】非对称 Goldbeter–Koshland 稳态（κ1 = K_m,1/S_T 挂 (1−u) 支，
κ2 = K_m,2/S_T 挂 u 支）：
        ξ = u/(1−u) · (κ1+1−u)/(κ2+u)
r=1（对称）时复现代码 4 全部基线（n_H 精确、群作用机器零、χ 分层）。

【结构】
P1 简并群维数：7 参数 (k1,E1T,k2,E2T,Km1,Km2,ST) → 3 可辨识组合 (ξ,κ1,κ2)
   → 4 条独立精确群作用（数值验证 max|Δu| = 0）+ 1 条反例（Km1 独变）
P2 n_H 闭式推广：数值 Hill 系数 vs 闭式，跨不对称度 r = κ1/κ2
P3 Fisher 壁垒：形状块 (ln c, ln κ1, ln κ2) vs 绝对块 (ln c, ln Km1, ln Km2, ln ST)
输出：非对称Km分析图.png（四联）+ SI 就绪表格
运行：python3 代码10_非对称Km简并群分析.py（约 1 分钟）
"""
import numpy as np
from scipy.optimize import brentq
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
SIG = 0.03

# ============================================================
# 主方程与解析导数（非对称）
# ============================================================
def solve_u(xi, k1_, k2_):
    xi = np.atleast_1d(np.asarray(xi, float))
    f = lambda u, x: x*(1-u)*(k2_+u) - u*(k1_+1-u)
    out = np.empty_like(xi)
    for i, x in enumerate(xi.flat):
        out.flat[i] = brentq(f, 1e-16, 1-1e-16, args=(x,), xtol=1e-14, rtol=1e-14)
    return out

def dxi_du(u, k1_, k2_):
    xi = u/(1-u)*(k1_+1-u)/(k2_+u)
    return xi*(1/u - 1/(k1_+1-u) + 1/(1-u) - 1/(k2_+u))

def du_dki(u, xi, k1_, k2_, which):
    """du/dκ_i（隐函数求导）。"""
    G = lambda u, a, b: xi*(1-u)*(b+u) - u*(a+1-u)
    Gu = xi*(1-k2_-2*u) - (k1_+1-2*u)
    if which == 1:
        Gi = -u
    else:
        Gi = xi*(1-u)
    return -Gi/Gu

def nH_closed(k1_, k2_):
    """n_H 精确闭式推广：u=1/2 处 logit 斜率。"""
    return 4.0/(4.0 - 1.0/(k1_+0.5) - 1.0/(k2_+0.5))

def nH_numeric(k1_, k2_):
    """数值 Hill 系数：u=1/2 处 d logit(u)/d ln ξ = 4ξ·du/dξ。"""
    xi50 = 0.5/0.5*(k1_+0.5)/(k2_+0.5)
    return 4.0*xi50/dxi_du(0.5, k1_, k2_)

print('='*72)
print('代码 10：非对称 K_m（κ1 ≠ κ2）PdPC 审计')
print('='*72)

# ============================================================
# P1 简并群维数：4 条独立精确群作用 + 1 反例
# ============================================================
print('\n[P1] 简并群维数（参数空间 7 维 → 可辨识组合 3 维 (ξ, κ1, κ2)）')
# 基准（非对称）：k1,E1T,k2,E2T,Km1,Km2,ST
base = dict(k1=1.0, E1T=1.0, k2=1.0, E2T=1.0, Km1=0.03, Km2=0.1, ST=1.0)
x_grid = np.logspace(-2, 2, 200)

def titration(p):
    xi = p['k1']*p['E1T']*x_grid/(p['k2']*p['E2T'])
    return solve_u(xi, p['Km1']/p['ST'], p['Km2']/p['ST'])

u_base = titration(base)
actions = {
    '(a) k1,k2 ×2':        {**base, 'k1': 2.0, 'k2': 2.0},           # ξ 不变
    "(a') E1T,E2T ×3":     {**base, 'E1T': 3.0, 'E2T': 3.0},         # ξ 不变（独立方向）
    '(b) Km1,Km2,ST ×5':   {**base, 'Km1': 0.15, 'Km2': 0.5, 'ST': 5.0},  # κ1,κ2 不变
    '(c) k1,E2T ×4':       {**base, 'k1': 4.0, 'E2T': 4.0},          # ξ 不变
}
print(f'     基准 (Km1,Km2,ST)=(0.03,0.1,1)，200 点对数滴定：')
for name, g in actions.items():
    d = np.max(np.abs(titration(g) - u_base))
    print(f'     群作用 {name:<22} max|Δu| = {d:.3e}  {"（精确 ✓）" if d < 1e-10 else "（破缺 ✗）"}')
# 反例：Km1 独变（κ1 改变 → 曲线应变）
g_bad = {**base, 'Km1': 0.06}
d_bad = np.max(np.abs(titration(g_bad) - u_base))
print(f'     反例   Km1 ×2 独变           max|Δu| = {d_bad:.3e}  （非群作用，曲线改变 ✓ 预期）')
print('     结论：4 条独立精确作用 → 简并群保持 4 维（7−3）；评审"降至 2 维"断言不成立。')

# ============================================================
# P2 n_H 闭式推广验证
# ============================================================
print('\n[P2] n_H 闭式推广：n_H = 4/(4 − 1/(κ1+½) − 1/(κ2+½))')
k2_fix = 0.01
ratios = [0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0]
print(f'     κ2 = {k2_fix} 固定；r = κ1/κ2 扫描')
print(f'     {"r":<8}{"n_H 数值":<12}{"n_H 闭式":<12}{"相对偏差":<12}{"对称式 1+1/(2κ2)":<18}')
rows = []
for r in ratios:
    k1_ = r*k2_fix
    nh_n, nh_c = nH_numeric(k1_, k2_fix), nH_closed(k1_, k2_fix)
    dev = abs(nh_n-nh_c)/nh_c
    nh_sym = 1 + 1/(2*k2_fix)
    rows.append((r, nh_n, nh_c, dev))
    print(f'     {r:<8}{nh_n:<12.4f}{nh_c:<12.4f}{dev:<12.2e}{nh_sym:<18.2f}')
print(f'     最大相对偏差 = {max(x[3] for x in rows):.2e}（机器精度 → 闭式精确成立）')
print(f'     注：r≠1 时对称式 1+1/(2κ2) 系统性失效（如 r=30：闭式 {rows[-1][2]:.1f} vs 对称式 51）')
# 对称基线复核
nh1_n, nh1_c = nH_numeric(k2_fix, k2_fix), 1+1/(2*k2_fix)
print(f'     对称基线复核（r=1）：数值 {nh1_n:.6f} vs 1+1/(2κ) = {nh1_c:.6f} ✓')

# ============================================================
# P3 Fisher 壁垒持续性
# ============================================================
print('\n[P3] Fisher 壁垒：形状块 (ln c, ln κ1, ln κ2) vs 绝对块 (ln c, ln Km1, ln Km2, ln ST)')
LNXI = np.linspace(-1.5*np.log(10), 1.5*np.log(10), 20)
XI = np.exp(LNXI)

def jac_shape(k1_, k2_):
    u = solve_u(XI, k1_, k2_)
    a = XI/dxi_du(u, k1_, k2_)
    b1 = k1_*du_dki(u, XI, k1_, k2_, 1)
    b2 = k2_*du_dki(u, XI, k1_, k2_, 2)
    return np.column_stack([a, b1, b2])

def jac_absolute(k1_, k2_):
    """绝对参数 (ln c, ln Km1, ln Km2, ln ST)：u 经 κ1=Km1/ST、κ2=Km2/ST 依赖。"""
    u = solve_u(XI, k1_, k2_)
    a = XI/dxi_du(u, k1_, k2_)
    b1 = k1_*du_dki(u, XI, k1_, k2_, 1)
    b2 = k2_*du_dki(u, XI, k1_, k2_, 2)
    return np.column_stack([a, b1, b2, -(b1+b2)])

print(f'     {"r":<8}{"χ(形状 3参数)":<16}{"χ(绝对 4参数)":<16}{"λ_min(绝对)":<14}')
fisher_rows = []
for r in ratios:
    k1_ = r*k2_fix
    Fs = jac_shape(k1_, k2_fix).T @ jac_shape(k1_, k2_fix)/SIG**2
    Fa = jac_absolute(k1_, k2_fix).T @ jac_absolute(k1_, k2_fix)/SIG**2
    es = np.linalg.eigvalsh(Fs); ea = np.linalg.eigvalsh(Fa)
    chi_s = np.sqrt(es.max()/es.min())
    chi_a = np.sqrt(ea.max()/max(ea.min(), 1e-300))
    fisher_rows.append((r, chi_s, chi_a, ea.min()))
    print(f'     {r:<8}{chi_s:<16.2f}{chi_a:<16.2e}{ea.min():<14.2e}')
print('     结论：形状 χ 保持 O(10)（可辨识），绝对块始终含数值零特征值（壁垒持续）✓')

# ============================================================
# 图：四联
# ============================================================
fig, ax = plt.subplots(2, 2, figsize=(13, 9))
xi_c = np.logspace(-1.5, 1.5, 400)
for r, col in zip([0.1, 0.3, 1.0, 3.0, 10.0], plt.cm.viridis(np.linspace(0, 0.9, 5))):
    ax[0,0].semilogx(xi_c, solve_u(xi_c, r*k2_fix, k2_fix), color=col, lw=1.5, label=f'r={r}')
ax[0,0].set_xlabel('ξ（激酶/磷酸酶活性比）'); ax[0,0].set_ylabel('u（修饰分数）')
ax[0,0].set_title(f'(a) 非对称滴定曲线（κ2={k2_fix}）'); ax[0,0].legend(); ax[0,0].grid(alpha=0.3)
rs = np.array([x[0] for x in rows])
ax[0,1].loglog(rs, [x[1] for x in rows], 'o', ms=7, label='数值')
ax[0,1].loglog(rs, [x[2] for x in rows], '-', lw=1.5, label='闭式 4/(4−1/(κ1+½)−1/(κ2+½))')
ax[0,1].axhline(1+1/(2*k2_fix), color='r', ls='--', lw=1, label='对称式 1+1/(2κ2)（失效）')
ax[0,1].set_xlabel('r = κ1/κ2'); ax[0,1].set_ylabel('n_H')
ax[0,1].set_title('(b) n_H 闭式推广 vs 数值'); ax[0,1].legend(fontsize=8); ax[0,1].grid(alpha=0.3)
fr = np.array(fisher_rows)
ax[1,0].semilogx(fr[:,0], fr[:,1], 'o-', label='χ 形状块 (lnc, lnκ1, lnκ2)')
ax[1,0].semilogx(fr[:,0], fr[:,2], 's--', label='χ 绝对块 (lnc, lnKm1, lnKm2, lnST)')
ax[1,0].set_yscale('log')
ax[1,0].set_xlabel('r = κ1/κ2'); ax[1,0].set_ylabel('条件数 χ')
ax[1,0].set_title('(c) Fisher 壁垒跨不对称度持续'); ax[1,0].legend(); ax[1,0].grid(alpha=0.3)
names = list(actions.keys()) + ['Km1 ×2 独变（反例）']
vals = [np.max(np.abs(titration(g)-u_base)) for g in actions.values()] + [d_bad]
colors = ['tab:blue']*4 + ['tab:red']
ax[1,1].bar(range(5), np.maximum(vals, 1e-17), color=colors)
ax[1,1].set_yscale('log'); ax[1,1].set_ylim(1e-17, 1)
ax[1,1].set_xticks(range(5)); ax[1,1].set_xticklabels(['(a)', "(a')", '(b)', '(c)', '反例'], fontsize=9)
ax[1,1].set_ylabel('max|Δu|')
ax[1,1].set_title('(d) 群作用精确性（4 精确 + 1 破缺）→ 群维数 4')
ax[1,1].grid(alpha=0.3, axis='y')
fig.tight_layout()
out = os.path.join(OUTDIR, '非对称Km分析图.png')
fig.savefig(out, dpi=160)
print(f'\n图已保存：{out}')

# SI 就绪表
print('\n[SI 就绪] 表 S1：非对称 K_m 审计汇总')
print('r\tκ1\tn_H 数值\tn_H 闭式\tχ形状\tχ绝对\tλmin绝对')
for (r, nh_n, nh_c, dev), (r2, cs, ca, lm) in zip(rows, fisher_rows):
    print(f'{r}\t{r*k2_fix:.4f}\t{nh_n:.3f}\t{nh_c:.3f}\t{cs:.1f}\t{ca:.2e}\t{lm:.1e}')
print('\n完成。')
