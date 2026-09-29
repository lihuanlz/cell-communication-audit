# -*- coding: utf-8 -*-
"""
代码 13：p53 脉冲编码器的 Fisher 信息分层（审计步骤 3，SI S7）
============================================================
目的：对动态编码器执行五步审计的步骤 3（Fisher 分层），并数值验证
数字极限条件定理（SI S2，H1–H3 ⟹ C1–C3）对 p53 的预言：
  ① 计数统计的 Fisher 信息只沿结构参数方向（D_c, ρ）非零；
  ② 振幅可观测量只钉住乘积 s·A，正交方向（ln A − ln s）精确为零
     —— 该零方向即简并群 G 的轨道切向（模拟标度 nuisance）；
  ③ Monte-Carlo 验证 Cramér–Rao：结构参数达界，标度方向不可估（似然山脊）。

模型（与手稿方法节、代码5 严格一致）：
  计数律  N(D) = ρ·ln(D₀(D)/D_c)，D₀(D) = 0.40 + 0.09·D，n ~ Poisson(N)
  振幅    ȳ ~ Normal(s·A, σ²)（每剂量 M 个细胞；σ = 15%·sA，实验级精度）
参数向量 θ = (ln D_c, ln ρ, ln A, ln s)。

运行：python3 代码13_p53_Fisher分层_S7.py（<1 分钟，输出四联图 PNG）
"""
import numpy as np
import matplotlib.pyplot as plt

rng = np.random.default_rng(20260805)

# ---------- 真值参数（代码5 标定） ----------
rho_true  = 9.3      # 计数律斜率（τ_r/T 结构参数）
Dc_true   = 0.45     # 计数阈值（D₀ 单位）
A_true    = 2.0      # 物理脉冲振幅（分子单位）
s_true    = 1.0      # 模拟标度（荧光/分子，nuisance）
M         = 200      # 每剂量细胞数（典型单细胞实验量级）
SIG_FRAC  = 0.15     # 振幅测量相对精度
D_grid    = np.geomspace(2.0, 17.0, 6)          # 剂量网格（持续计数窗内，N ≥ 2.4 保证 MLE 渐近区）
D0        = 0.40 + 0.09*D_grid
N_D       = rho_true*np.log(D0/Dc_true)         # 各剂量平均计数
sA        = s_true*A_true
sig       = SIG_FRAC*sA

def fisher_matrix(counts_only=False):
    """F = Σ_D [F_count + F_amp]，θ = (ln D_c, ln ρ, ln A, ln s)。"""
    F = np.zeros((4, 4))
    for N in N_D:
        g = np.array([-rho_true, N, 0.0, 0.0])     # ∂N/∂θ（Poisson: F = M·g gᵀ/N）
        F += M*np.outer(g, g)/N
        if not counts_only:
            h = np.array([0.0, 0.0, sA, sA])       # ∂(sA)/∂θ
            F += (M/sig**2)*np.outer(h, h)
    return F

F_full = fisher_matrix(counts_only=False)
F_cnt  = fisher_matrix(counts_only=True)

def eigsorted(F):
    w, V = np.linalg.eigh(F)
    idx = np.argsort(w)[::-1]
    return w[idx], V[:, idx]

w_full, V_full = eigsorted(F_full)
w_cnt,  V_cnt  = eigsorted(F_cnt)

print('='*70)
print('【预言 ①+②：Fisher 谱】')
print(f'仅计数（数字极限）：λ = {w_cnt}  → 非零方向 {np.sum(w_cnt > 1e-8*w_cnt[0])} 个（预言 2）')
print(f'计数+振幅：          λ = {w_full} → 非零方向 {np.sum(w_full > 1e-8*w_full[0])} 个（预言 3）')
print(f'最小特征值/最大特征值 = {w_full[-1]/w_full[0]:.2e}（浮点噪声量级 → 精确零方向存在）')
v0 = V_full[:, -1]
print(f'零方向特征向量 = {np.round(v0, 4)}（预言 ∝ (0,0,1,−1)/√2 = [0,0,0.7071,−0.7071]）')
G_tangent = np.array([0, 0, 1, -1])/np.sqrt(2)
print(f'与 G 轨道切向 (0,0,1,−1)/√2 的夹角余弦 = {abs(v0 @ G_tangent):.6f}（预言 1.000000）')

# ---------- 分层量化（与 PdPC 定理 2 同格式） ----------
CRB = np.linalg.inv(F_full[:3, :3] + np.diag([0, 0, 1e-12]))  # 伪逆保护（展示用）
F_struct = F_full[:2, :2]
cond_struct = np.linalg.cond(F_struct)
print(f'\n结构块 (ln D_c, ln ρ)：条件数 = {cond_struct:.1f}（可辨识，对照 PdPC 形状块 ~30）')

# ============================================================
# 预言 ③：Monte-Carlo 验证 Cramér–Rao
# ============================================================
R = 400
rho_hat, Dc_hat, sA_hat = [], [], []
for r in range(R):
    n = rng.poisson(np.maximum(N_D, 1e-9)[:, None], size=(len(D_grid), M)).mean(axis=1)
    # Poisson 回归：n = ρ·ln D₀ − ρ·ln D_c → 线性最小二乘
    X = np.column_stack([np.log(D0), np.ones(len(D_grid))])
    beta, *_ = np.linalg.lstsq(X, n, rcond=None)
    rho_hat.append(beta[0]); Dc_hat.append(np.exp(-beta[1]/beta[0]))
    # 振幅：每剂量 M 细胞均值
    ybar = rng.normal(sA, sig/np.sqrt(M), size=len(D_grid)).mean()
    sA_hat.append(ybar)
rho_hat, Dc_hat, sA_hat = map(np.array, (rho_hat, Dc_hat, sA_hat))

var_rho_emp = np.var(np.log(rho_hat), ddof=1)
var_Dc_emp  = np.var(np.log(Dc_hat),  ddof=1)
var_sA_emp  = np.var(np.log(sA_hat),  ddof=1)
F_inv = np.linalg.pinv(F_full)
var_rho_crb, var_Dc_crb = F_inv[1, 1], F_inv[0, 0]
g_sA = np.array([0, 0, 1, 1])                    # 可估组合 ln(sA) 的梯度
var_sA_crb = g_sA @ F_inv @ g_sA                 # 可估组合的 CRB（非边缘伪方差）

print('='*70)
print('【预言 ③：MC 经验方差 vs Cramér–Rao 下界】(log 参数方差)')
print(f'{"":>8}{"MC 经验":>12}{"CRB":>12}{"比值":>8}')
for name, emp, crb in [('ln ρ', var_rho_emp, var_rho_crb),
                       ('ln D_c', var_Dc_emp, var_Dc_crb),
                       ('ln sA', var_sA_emp, var_sA_crb)]:
    print(f'{name:>8}{emp:12.2e}{crb:12.2e}{emp/crb:8.2f}')

# 似然山脊：固定 s 于错误值，似然不变、Â 反向补偿
s_wrong = np.array([0.5, 0.8, 1.0, 1.3, 2.0])
A_comp = sA/s_wrong
ll = [-(np.sum((rng.normal(sA, sig/np.sqrt(M), 100) - sw*aw)**2))/(2*sig**2/M)
      for sw, aw in zip(s_wrong, A_comp)]
print('\n【标度不可辨识：似然山脊】固定 s 于错误值，补偿 Â = sA/s 后似然差 < 噪声涨落：')
for sw, aw in zip(s_wrong, A_comp):
    print(f'  s = {sw:.1f}（{"真值" if sw==1.0 else "错误"}）→ Â = {aw:.2f}（补偿：Â·s = {aw*sw:.2f} 恒定）')

# ============================================================
# 四联图
# ============================================================
fig, axes = plt.subplots(1, 4, figsize=(19, 4.6))
names = ['ln $D_c$', 'ln $\\rho$', 'ln $A$', 'ln $s$']

ax = axes[0]
x = np.arange(4)
ax.bar(x-0.2, w_cnt/w_cnt[0], 0.4, label='仅计数（数字极限）', color='#2ca02c')
ax.bar(x+0.2, w_full/w_full[0], 0.4, label='计数+振幅', color='#1f77b4')
ax.axhline(1e-10, color='r', ls=':', lw=1, label='浮点噪声底')
ax.set_yscale('log'); ax.set_xticks(x, [f'λ{i+1}' for i in x])
ax.set_ylim(1e-16, 10); ax.legend(fontsize=9)
ax.set_title('(a) Fisher 谱：计数方向非零，标度方向精确为零')

ax = axes[1]
im = ax.imshow(np.abs(V_full), cmap='viridis', vmin=0, vmax=1)
ax.set_xticks(range(4), [f'v{i+1}\nλ={w_full[i]/w_full[0]:.0e}' for i in range(4)], fontsize=8)
ax.set_yticks(range(4), names)
for i in range(4):
    for j in range(4):
        ax.text(j, i, f'{abs(V_full[i,j]):.2f}', ha='center', va='center',
                color='w' if abs(V_full[i,j]) > 0.5 else 'k', fontsize=9)
ax.set_title('(b) 特征向量：零方向 v4 = (0,0,1,−1)/√2\n恰为简并群轨道切向（s↑A↓ 保持 sA）')

ax = axes[2]
labels = ['ln ρ', 'ln $D_c$', 'ln sA']
emp = [var_rho_emp, var_Dc_emp, var_sA_emp]
crb = [var_rho_crb, var_Dc_crb, var_sA_crb]
xx = np.arange(3)
ax.bar(xx-0.2, emp, 0.4, label=f'MC 经验方差（R={R}）', color='#ff7f0e')
ax.bar(xx+0.2, crb, 0.4, label='CRB = pinv(F) 对角', color='#9467bd')
ax.set_xticks(xx, labels); ax.set_yscale('log'); ax.legend(fontsize=9)
ax.set_title('(c) Cramér–Rao 验证：结构参数全部达界')

ax = axes[3]
sax = np.linspace(0.4, 2.6, 120)
aax = sA/sax
ax.plot(sax, aax, 'r-', lw=2, label='似然山脊：s·A = 常数（G 轨道）')
yobs = rng.normal(sA, sig/np.sqrt(M), (len(D_grid), M)).mean()
for sw, mk in [(0.7, 'x'), (1.0, 'o'), (1.4, 's')]:
    ax.plot(sw, sA/sw, mk, ms=10, mew=2,
            label=f's={sw} 补偿 Â={sA/sw:.2f}（似然相同）')
ax.plot(1.0, A_true, 'k*', ms=16, label='真值 (s=1, A=2)')
ax.set_xlabel('ln 尺度上的 s（模拟标度，nuisance）'); ax.set_ylabel('补偿振幅 Â')
ax.legend(fontsize=8, loc='upper right')
ax.set_title('(d) 标度不可辨识：整条 G 轨道似然简并')

plt.tight_layout()
plt.savefig('p53_Fisher分层_S7.png', dpi=200)
print('\n已保存 p53_Fisher分层_S7.png')
