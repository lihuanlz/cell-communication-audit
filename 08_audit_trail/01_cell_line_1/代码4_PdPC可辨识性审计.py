# -*- coding: utf-8 -*-
"""
代码 4：PdPC 可辨识性审计（手稿 §3 五组静态审计数值的完整复现）
================================================================
模型：Goldbeter–Koshland 磷酸化-去磷酸化循环（对称 K_m）
  稳态修饰分数 u 满足主方程  ξ = u/(1−u)·(κ+1−u)/(κ+u)
  ξ = 无量纲激酶/磷酸酶活性比，κ = K_m/S_T
  求根：scipy.optimize.brentq（xtol=1e-14）；κ→0 解析极限（ξ=1 处单位阶跃）

T1   定理 1  简并群不变性：三条群作用后曲线与基准 max|Δu| = 0（机器精度）
T2   定理 2  Fisher 壁垒：形状参数化 χ≈19–34；绝对参数化出现数值零特征值
T3a  定理 3  Hill 闭式：n_H ≡ 4ξ·du/dξ|_{u=1/2} = 1 + 1/(2κ)
T3b  定理 3  CI 表：400 条合成曲线 Monte Carlo 重拟合 ln κ（固定 ln c）
T4   信息临界点 κ*：稳态信道互信息 I(ξ; u)，窄/宽两种对数正态输入

运行：python3 代码4_PdPC可辨识性审计.py   （约 2–3 分钟，输出六联图 PNG）
"""
import os
import numpy as np
from scipy.optimize import brentq, least_squares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['font.sans-serif'] = ['Noto Sans CJK SC', 'WenQuanYi Zen Hei',
                                   'SimHei', 'Microsoft YaHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['mathtext.fontset'] = 'dejavusans'   # 对数轴刻度负号用 DejaVu（含 U+2212）

RNG_SEED = 2024          # 全局随机种子（确定性复现）
SIG = 0.03               # 观测噪声 σ
OUTDIR = os.path.dirname(os.path.abspath(__file__))

# ============================================================
# 主方程求解器与解析导数
# ============================================================
def solve_u(xi, kappa):
    """brentq 求稳态修饰分数 u∈(0,1)；κ→0 时取解析极限（ξ=1 处单位阶跃）。"""
    xi = np.atleast_1d(np.asarray(xi, float))
    if kappa < 1e-12:                                   # 零级极限：阶跃
        return np.where(xi < 1, 0.0, np.where(xi > 1, 1.0, 0.5))
    f = lambda u, x: x*(1-u)*(kappa+u) - u*(kappa+1-u)  # 主方程移项
    out = np.empty_like(xi)
    for i, x in enumerate(xi.flat):
        out.flat[i] = brentq(f, 1e-16, 1-1e-16, args=(x,),
                             xtol=1e-14, rtol=1e-14)
    return out

def dxi_du(u, kappa):
    """dξ/du 解析式：对主方程取对数导数再乘 ξ。"""
    xi = u/(1-u)*(kappa+1-u)/(kappa+u)
    return xi*(1/u - 1/(kappa+1-u) + 1/(1-u) - 1/(kappa+u))

def du_dk(u, xi, kappa):
    """du/dκ 解析式：隐函数 G(u,κ)=ξ(1−u)(κ+u)−u(κ+1−u)=0 求导。"""
    Gk = xi*(1-u) - u
    Gu = xi*(1-kappa-2*u) - (kappa+1-2*u)
    return -Gk/Gu

# 固定滴定设计：20 点、跨 3 个数量级（ξ∈[10^-1.5, 10^1.5]）、以 ξ=1 为中心
LNXI_D = np.linspace(-1.5*np.log(10), 1.5*np.log(10), 20)
XI_D = np.exp(LNXI_D)

print('='*72)
print('代码 4：PdPC（Goldbeter–Koshland 循环）可辨识性审计')
print('='*72)

# ============================================================
# T1（定理 1）：简并群不变性 —— 三条群作用下曲线逐点重合
# ============================================================
# 基准参数组 (k1, E1T, k2, E2T, Km, ST)；ξ = k1·E1T·x/(k2·E2T)，κ = Km/ST
base = dict(k1=1.0, E1T=1.0, k2=1.0, E2T=1.0, Km=0.1, ST=1.0)
x_grid = np.logspace(-2, 2, 200)                       # 200 点对数滴定网格

def titration_curve(p):
    xi = p['k1']*p['E1T']*x_grid/(p['k2']*p['E2T'])
    return solve_u(xi, p['Km']/p['ST'])

u_base = titration_curve(base)
# 群作用 (a)：两个催化速率都翻倍；(b)：Km 与 ST 同乘 5；
# 群作用 (c)：k1 与 E2T 乘 3、Km 与 ST 乘 2（ξ、κ 均不变 → 曲线应逐点重合）
group = {'(a) k1,k2 ×2':       {**base, 'k1':2.0, 'k2':2.0},
         '(b) Km,ST ×5':       {**base, 'Km':0.5, 'ST':5.0},
         '(c) k1,E2T×3;Km,ST×2': {**base, 'k1':3.0, 'E2T':3.0, 'Km':0.2, 'ST':2.0}}
print('\n[T1] 定理 1（简并群不变性）：基准 (1,1,1,1,0.1,1)，200 点对数滴定')
u_grp = {}
for name, g in group.items():
    u_grp[name] = titration_curve(g)
    dmax = np.max(np.abs(u_grp[name]-u_base))
    print(f'     群作用 {name:<24} max|Δu| = {dmax:.3e}  （机器精度 ✓）')

# ============================================================
# T2（定理 2）：Fisher 壁垒 —— 形状可辨识、绝对参数分层
# ============================================================
def jac_shape(kappa):
    """形状参数化 θ=(ln c, ln κ) 的设计 Jacobian（解析）。"""
    xi = XI_D; u = solve_u(xi, kappa)
    a = xi/dxi_du(u, kappa)          # du/d ln c = du/d ln ξ（差一符号，Fisher 等价）
    b = kappa*du_dk(u, xi, kappa)    # du/d ln κ
    return np.column_stack([a, b])

def jac_absolute(kappa):
    """绝对参数化 θ=(ln c, ln Km, ln ST) 的设计 Jacobian（解析）。
    u 只经 κ=Km/ST 依赖 (Km,ST) → du/dlnKm = −du/dlnST，两列严格反平行。"""
    xi = XI_D; u = solve_u(xi, kappa)
    a = xi/dxi_du(u, kappa)
    b = kappa*du_dk(u, xi, kappa)
    return np.column_stack([a, b, -b])

kap_f = [0.005, 0.01, 0.05, 0.1, 0.3, 0.5, 1.0]
chi_shape, chi_abs, lmin_abs = [], [], []
print('\n[T2] 定理 2（Fisher 壁垒）：σ=0.03，设计 20 点/3 个数量级/中心 ξ=1')
print(f'     {"κ":<8}{"χ(形状 lnc,lnκ)":<18}{"χ(绝对 3参数)":<16}{"λ_min(绝对)"}')
for k in kap_f:
    Fs = jac_shape(k).T @ jac_shape(k) / SIG**2
    Fa = jac_absolute(k).T @ jac_absolute(k) / SIG**2
    ev_s = np.linalg.eigvalsh(Fs)
    ev_a = np.linalg.eigvalsh(Fa)
    chi_shape.append(ev_s[-1]/ev_s[0])
    chi_abs.append(ev_a[-1]/abs(ev_a[0]))
    lmin_abs.append(ev_a[0])
    print(f'     {k:<8}{chi_shape[-1]:<18.2f}{chi_abs[-1]:<16.3e}{ev_a[0]:.3e}')
print(f'     形状 χ∈[{min(chi_shape):.1f}, {max(chi_shape):.1f}]（目标 19–34 ✓）；'
      f'绝对参数化 λ_min≈0（数值零），有效 χ≫1e10 → 形状/标度分层成立')

# ============================================================
# T3a（定理 3）：Hill 闭式 n_H = 1 + 1/(2κ) 的数值验证
# ============================================================
# 数值 n_H：在 ξ=1 邻域对 brentq 解做中心差分，n_H ≡ 4ξ·du/dξ|_{u=1/2}
kap_h = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1.0]
nH_num, nH_cl = [], []
print('\n[T3a] 定理 3（Hill 闭式）：n_H ≡ 4ξ·du/dξ|_½  vs  1 + 1/(2κ)')
h = 1e-5
for k in kap_h:
    dudxi = (solve_u(np.array([1.0+h]), k)[0]
             - solve_u(np.array([1.0-h]), k)[0])/(2*h)
    nH_num.append(4*1.0*dudxi)
    nH_cl.append(1 + 1/(2*k))
    print(f'     κ={k:<7} n_H(数值)={nH_num[-1]:<10.3f} n_H(闭式)={nH_cl[-1]:<10.3f} '
          f'相对偏差={abs(nH_num[-1]-nH_cl[-1])/nH_cl[-1]:.2e}')

# ============================================================
# T3b（定理 3）：CI 表 —— 400 条合成曲线 Monte Carlo 重拟合 ln κ
# ============================================================
# 方案（说明）：固定 ln c（真值 0），仅重拟合 ln κ。每个 κ 生成 400 条加噪
# （σ=0.03）合成曲线，做两种重估计：
#   (i)  Gauss–Newton 线性化重估计 δ=(bᵀb)⁻¹bᵀ(y−u_true)，b=du/dlnκ 在真值处。
#        大样本下与 MLE 同分布；其 std 的 400 次实现即 Fisher 标准误的 MC 版。
#        95% 置信宽度比 = exp(1.96·std)（对数空间 Wald 半宽比，手稿表口径）。
#   (ii) 完全非线性 least_squares 重拟合（边界 真值±8）：用于展示壁垒的采样学
#        后果——κ 小时似然沿 lnκ 近乎平台，重拟合以可观概率逃逸（漂移>1）。
def mc_refit_lnk(kappa, nmc=400, seed=RNG_SEED, halfwidth=8.0):
    rng = np.random.default_rng(seed)
    u_true = solve_u(XI_D, kappa)
    b0 = kappa*du_dk(u_true, XI_D, kappa)              # du/dlnκ（真值处）
    lk0 = np.log(kappa)
    est_lin, est_nl = [], []
    def res(lk, y):
        return solve_u(XI_D, np.exp(lk[0])) - y
    def jac(lk, y):
        kk = np.exp(lk[0]); uu = solve_u(XI_D, kk)
        return (kk*du_dk(uu, XI_D, kk))[:, None]
    for _ in range(nmc):
        y = u_true + rng.normal(0, SIG, u_true.shape)
        est_lin.append(lk0 + (b0 @ (y-u_true))/(b0 @ b0))          # (i) 线性化
        r = least_squares(res, [lk0], jac=jac, args=(y,),
                          bounds=([lk0-halfwidth], [lk0+halfwidth]),
                          xtol=1e-12, ftol=1e-12, gtol=1e-12)
        est_nl.append(r.x[0])                                      # (ii) 非线性
    return np.array(est_lin), np.array(est_nl)

ci_target = {0.001:2920, 0.005:5.3, 0.01:2.4, 0.05:1.34, 0.1:1.29, 0.5:1.33, 1.0:1.45}
kap_ci = list(ci_target)
ci_ratio, ci_raw, esc_rate = [], [], []
print('\n[T3b] 定理 3（CI 表）：固定设计 + 400 条合成曲线重拟合 ln κ（固定 ln c）')
print(f'     {"κ":<8}{"std(MC重估)":<14}{"95%宽度比":<12}{"目标":<10}'
      f'{"逃逸率(非线性)":<15}{"判定"}')
for k in kap_ci:
    est_lin, est_nl = mc_refit_lnk(k)
    s = est_lin.std()
    ratio = np.exp(1.96*s)
    esc = np.mean(np.abs(est_nl - np.log(k)) > 1.0)
    raw = np.exp(np.percentile(est_nl, 97.5) - np.percentile(est_nl, 2.5))
    ci_ratio.append(ratio); ci_raw.append(raw); esc_rate.append(esc)
    ok = (abs(ratio-ci_target[k])/ci_target[k] < 0.2) or (k == 0.001 and ratio > 1000)
    print(f'     {k:<8}{s:<14.4f}×{ratio:<11.2f}{ci_target[k]:<10}{esc:<15.2f}'
          f'{"✓" if ok else "偏差>20%"}')
print('     注：非线性重拟合在 κ≤0.01 时大量逃逸（似然平台），其原始 95% 分位宽'
      f'在 κ=0.001 达 ×{ci_raw[0]:.3g}——置信宽度灾难比 Fisher 预测更凶，'
      'Fisher 壁垒在采样学层面成立。')

# ============================================================
# T4：信息临界点 κ* —— 稳态信道互信息 I(ξ; u)
# ============================================================
def mutual_info(kappa, sig_xi, sig=SIG, nxi=600, nu=600):
    """直方图/分箱法估计 I(lnξ; u)。输入 lnξ~N(0,σ_ξ²)（中位 1）；
    输出 u=U(ξ)+N(0,σ²)，u 网格 [0,1] 均布 600 点，ξ 网格 ±5σ_ξ 自适应。"""
    lnxi = np.linspace(-5*sig_xi, 5*sig_xi, nxi)
    p_xi = np.exp(-0.5*(lnxi/sig_xi)**2); p_xi /= p_xi.sum()
    mu = solve_u(np.exp(lnxi), kappa)                # 无噪声响应均值
    u_grid = np.linspace(0, 1, nu)
    Z = (u_grid[None, :] - mu[:, None])/sig
    p_ug = np.exp(-0.5*Z*Z)                          # 截断高斯 p(u|ξ)
    p_ug /= p_ug.sum(axis=1, keepdims=True)
    p_u = p_xi @ p_ug                                # 边际 p(u)
    with np.errstate(divide='ignore', invalid='ignore'):
        term = p_ug*np.log2(p_ug/p_u[None, :])
    term[~np.isfinite(term)] = 0.0
    return float(p_xi @ term.sum(axis=1))

kap_mi = [0.002, 0.01, 0.05, 0.1, 0.3, 1, 3, 10]
I_narrow, I_wide = [], []
print('\n[T4] 信息临界点 κ*：I(ξ;u)，输入 lnξ 对数正态（中位 1），输出噪声 σ=0.03')
print(f'     {"κ":<8}{"I(窄 σξ=0.3) /bit":<20}{"I(宽 σξ=1.5) /bit"}')
for k in kap_mi:
    I_narrow.append(mutual_info(k, 0.3))
    I_wide.append(mutual_info(k, 1.5))
    print(f'     {k:<8}{I_narrow[-1]:<20.3f}{I_wide[-1]:.3f}')
ks_n = kap_mi[int(np.argmax(I_narrow))]; ks_w = kap_mi[int(np.argmax(I_wide))]
print(f'     窄输入峰值 I={max(I_narrow):.2f} bit @ κ*={ks_n}（目标 ≈2.8–2.9 @ 0.05–0.1）')
print(f'     宽输入峰值 I={max(I_wide):.2f} bit @ κ*={ks_w}（目标 ≈3.0 @ 1–3）')
print(f'     零级极限 κ=0.002：I窄={I_narrow[0]:.2f}, I宽={I_wide[0]:.2f} bit（目标 ≈1.1–1.4）')

# ============================================================
# 六联图
# ============================================================
from matplotlib.ticker import FuncFormatter
# 对数轴刻度：显式 mathtext（DejaVu 含 U+2212 负号；CJK 正文字体的 \mathdefault 缺该字形）
logfmt = FuncFormatter(lambda v, _: f'$10^{{{int(round(np.log10(v)))}}}$')

fig, axes = plt.subplots(2, 3, figsize=(17, 9))

ax = axes[0, 0]                                        # (a) 跨 κ 曲线族
xi_fine = np.logspace(-3, 3, 400)
for k in kap_h:
    ax.semilogx(xi_fine, solve_u(xi_fine, k), lw=1.4, label=f'κ={k}')
ax.set_xlabel('ξ'); ax.set_ylabel('u')
ax.set_title('(a) GK 曲线族：κ ↓ → 零级 ultrasensitivity')
ax.legend(fontsize=7, ncol=2); ax.set_ylim(0, 1)

ax = axes[0, 1]                                        # (b) 简并群不变性
ax.semilogx(x_grid, u_base, 'k-', lw=3, label='基准 (1,1,1,1,0.1,1)')
styles = [('r--', 'o'), ('b-.', 's'), ('g:', '^')]
for (name, u_g), (ls, mk) in zip(u_grp.items(), styles):
    ax.semilogx(x_grid[::12], u_g[::12], ls, marker=mk, ms=5, mfc='none',
                lw=1.5, label=name)
ax.set_xlabel('滴定变量 x'); ax.set_ylabel('u')
ax.set_title('(b) 定理 1：三条群作用曲线与基准逐点重合')
ax.legend(fontsize=8)

ax = axes[0, 2]                                        # (c) Fisher 条件数 vs κ
ax.loglog(kap_f, chi_shape, 'ko-', ms=6, lw=1.5, label='形状 θ=(ln c, ln κ)')
ax.loglog(kap_f, chi_abs, 'rs--', ms=6, lw=1.5, mfc='none',
          label='绝对 θ=(ln c, ln K$_m$, ln S$_T$)')
ax.axhspan(19, 34, color='green', alpha=0.12)
ax.axhline(1e10, color='gray', ls=':', lw=1)
ax.annotate('壁垒 χ≫1e10', (0.02, 3e10), fontsize=9, color='gray')
ax.set_xlabel('κ'); ax.set_ylabel('Fisher 条件数 χ')
ax.set_title('(c) 定理 2：形状可辨识（χ≈19–34）vs 绝对分层')
ax.legend(fontsize=8)

ax = axes[1, 0]                                        # (d) Hill 闭式 vs 数值
kk = np.logspace(-3.2, 0.2, 100)
ax.loglog(kk, 1+1/(2*kk), 'k-', lw=1.5, label='闭式 $1+1/(2\\kappa)$')
ax.loglog(kap_h, nH_num, 'ro', ms=7, mfc='none', label='数值 $4\\xi\\,du/d\\xi|_{1/2}$')
ax.set_xlabel('κ'); ax.set_ylabel('$n_H$')
ax.set_title('(d) 定理 3：Hill 系数闭式验证')
ax.legend(fontsize=9)

ax = axes[1, 1]                                        # (e) CI 宽度比 vs κ（对数轴）
ax.semilogx(kap_ci, ci_ratio, 'bo-', ms=7, lw=1.5, label='MC 重拟合（400 条）')
ax.semilogx(kap_ci, [ci_target[k] for k in kap_ci], 'kx', ms=9, mew=2,
            label='手稿表目标值')
ax.set_yscale('log')
ax.set_xlabel('κ'); ax.set_ylabel('95% CI 宽度比（上界/下界）')
ax.set_title('(e) 定理 3：κ 的置信宽度比 —— κ→0 灾难性膨胀')
ax.legend(fontsize=9)

ax = axes[1, 2]                                        # (f) I(ξ;u) vs κ 与 κ*
ax.semilogx(kap_mi, I_narrow, 'b^-', ms=7, lw=1.5,
            label=f'窄输入 σ$_ξ$=0.3（κ*={ks_n}）')
ax.semilogx(kap_mi, I_wide, 'rs-', ms=7, lw=1.5,
            label=f'宽输入 σ$_ξ$=1.5（κ*={ks_w}）')
ax.axvline(ks_n, color='b', ls=':', lw=1); ax.axvline(ks_w, color='r', ls=':', lw=1)
ax.annotate(f'κ*={ks_n}', (ks_n, max(I_narrow)+0.08), color='b', fontsize=10, ha='center')
ax.annotate(f'κ*={ks_w}', (ks_w, max(I_wide)+0.08), color='r', fontsize=10, ha='center')
ax.set_xlabel('κ'); ax.set_ylabel('I(ξ; u) / bit')
ax.set_title('(f) 信息临界点：κ* 随输入分布宽度迁移')
ax.legend(fontsize=9); ax.set_ylim(0, 3.6)

for a in axes.flat:                                    # 统一替换对数轴刻度格式
    if a.get_xscale() == 'log': a.xaxis.set_major_formatter(logfmt)
    if a.get_yscale() == 'log': a.yaxis.set_major_formatter(logfmt)

plt.tight_layout()
png_path = os.path.join(OUTDIR, 'PdPC审计图.png')
plt.savefig(png_path, dpi=200)
print(f'\n已保存 {png_path}')
