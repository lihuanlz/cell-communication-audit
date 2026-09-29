# -*- coding: utf-8 -*-
"""
Created on Fri Apr 17 10:18:51 2026

@author: lihua
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Apr 13 08:43:36 2026

@author: lihua
"""

# -*- coding: utf-8 -*-
"""
Ω-κ Framework Simulation Code (精简版)
仅生成论文核心三图：
Figure 1: β-μ 空间响应曲面与 45° 斜坡
Figure 2: ξ 标度不变性验证
Figure 3: 定量准确性与 CV 范围
"""

import numpy as np
import matplotlib.pyplot as plt
from nature_style import apply_nature_style
apply_nature_style()  # 统一字体样式 + SVG/PNG 双输出（nature_style.py）
from mpl_toolkits.mplot3d import Axes3D
import os
from datetime import datetime
from scipy.stats import binom  # 新增：用于 LoB/LoD 计算
import warnings

import seaborn as sns


sns.set_style("white")   # 或 "ticks"，无网格背景





warnings.filterwarnings('ignore')

# 输出目录
output_dir = "omega_kappa_figures"
os.makedirs(output_dir, exist_ok=True)

# 绘图设置
plt.rcParams.update({
    # 'font.family': 'serif',
    
    
    'font.family': 'Times new roman',
    
    
    'font.size': 10,
    'savefig.dpi': 300,
    'axes.linewidth': 0.8,
})

# =============================================================================
# 核心函数（保留原计算函数，无需改动）
# =============================================================================
def calculate_p(mu, kappa, beta):
    A = kappa * beta + beta + mu
    discriminant = max(A**2 - 4 * beta * mu, 0.0)
    sqrt_disc = np.sqrt(discriminant)
    if mu < beta:
        p = (A - sqrt_disc) / (2 * beta)
    else:
        p = 2 * mu / (A + sqrt_disc)
    return np.clip(p, 0.0, 1.0)

def calculate_P_specific(p, beta):
    if p >= 1 - 1e-12:
        return 1.0
    if p <= 1e-12:
        return beta * p
    return -np.expm1(beta * np.log1p(-p))

def calculate_P_pos(mu, kappa, beta, b=0.01):
    p = calculate_p(mu, kappa, beta)
    P_specific = calculate_P_specific(p, beta)
    P_pos = b + (1 - b) * P_specific
    return P_pos, p, P_specific

def calculate_M_hat(P_pos_obs, kappa, beta, N, b=0.01, model='R1'):
    if P_pos_obs <= b:
        P_specific_hat = 0.0
    else:
        P_specific_hat = (P_pos_obs - b) / (1 - b)
    P_specific_hat = max(0.0, min(1.0 - 1e-12, P_specific_hat))
    if model == 'R1':
        term1 = kappa * ((1 - P_specific_hat) ** (-1 / beta) - 1)
        term2 = 1 - (1 - P_specific_hat) ** (1 / beta)
        M_hat = N * beta * (term1 + term2)
    elif model == 'R2':
        M_hat = -N * (kappa + 1) * np.log1p(-P_specific_hat)
    elif model == 'R3':
        M_hat = N * (kappa + 1) * P_specific_hat
    else:
        raise ValueError("Model must be 'R1', 'R2', or 'R3'")
    return max(0.0, M_hat)



def calculate_CV(M, kappa, beta, N, n, b=0.01, model='R1'):
    """
    计算 TCS 模型的变异系数 (CV)
    参数：
        M: 绝对分子数
        kappa: 平台特征常数
        beta: 每分区结合位点数（R1 需要；R2 和 R3 不依赖 beta 具体值）
        N: 总分区数
        n: 成像分区数
        b: 背景概率
        model: 模型选择 ('R1', 'R2', 'R3')
    返回：
        CV (无量纲，非百分比)
    """
    if model == 'R1' or model == 'R2':
        # 计算 P_pos, p, P_specific 用原函数
        P_pos, p, P_specific = calculate_P_pos(M/N, kappa, beta, b)
    
    if model == 'R1':
        term = kappa * (1 - P_specific) ** (-1/beta - 1) + (1 - P_specific) ** (1/beta - 1)
        CV = (N / (M * (1 - b))) * np.abs(term) * np.sqrt(P_pos * (1 - P_pos) / n)
    elif model == 'R2':
        CV = (N * (kappa + 1) / (M * (1 - b) * (1 - P_specific))) * np.sqrt(P_pos * (1 - P_pos) / n)
    # elif model == 'R3':
    #     # 按照论文 S4.25 公式
    #     gamma = 1.0 / (kappa + 1)          # γ = 1/(1+κ)
    #     mu = M / N                          # μ = M/N
    #     P_specific_approx = gamma * mu     # 线性近似下的 P_specific
    #     # 公式 CV ≈ 1/(γ(1-b)μ) * sqrt( (b + (1-b)γμ) / n )
    #     CV = 1.0 / (gamma * (1 - b) * mu) * np.sqrt((b + (1 - b) * P_specific_approx) / n)
    elif model == 'R3':
        mu = M / N                          # μ = M/N
        P_specific = mu / (kappa + 1)       # γ·μ = μ/(1+κ)
        P_pos = b + (1 - b) * P_specific    # 正向模型
        # 精确 S2f.25: CV = 1/[(1-b)·P_specific] · sqrt[P_pos·(1-P_pos)/n]
        CV = 1.0 / ((1 - b) * P_specific) * np.sqrt(P_pos * (1 - P_pos) / n)
    else:
        raise ValueError("model must be 'R1', 'R2', or 'R3'")
    return CV









# =============================================================================
# 新增：LoB 和 LoD 计算函数（严格按 S7 公式）
# =============================================================================
def calculate_LoB(n, b, method='exact'):
    """
    严格按 S7.1 公式计算 Limit of Blank (LoB)
    :param n: 观测分区数
    :param b: 背景阳性概率
    :param method: 'exact' (精确二项分位数, S7.0) 或 'normal' (正态近似, S7.1)
    :return: LoB (阳性率阈值)
    """
    if method == 'exact':
        # S7.0 精确公式: LoB = (1/n) * F_Bin^{-1}(0.95; n, b)
        k_95 = binom.ppf(0.95, n, b)
        return k_95 / n
    else:
        # S7.1 正态近似: LoB ≈ b + 1.645 * sqrt(b(1-b)/n)
        return b + 1.645 * np.sqrt(b * (1 - b) / n)

def calculate_LoD_exact_R1(kappa, beta, N, n, b, LoB, 
                           M_low_guess=1e-3, M_high_guess=1e5, tol=1e-2):
    """
    严格按 S7.2.1 精确 R1 方法计算 Limit of Detection (LoD)
    用二分法求解: Pr(P_pos_hat > LoB | M=LoD) = 0.95
    """
    def detection_probability(M):
        """计算给定 M 时的检测概率 (S7.2.1 求和公式)"""
        # 1. 从 M 正向计算 P_pos
        mu = M / N
        P_pos, _, P_specific = calculate_P_pos(mu, kappa, beta, b)
        
        # 2. 计算 Pr(P_pos_hat > LoB) = sum_{k=floor(n*LoB)+1}^n binom(n,k) P_pos^k (1-P_pos)^(n-k)
        k_threshold = int(np.floor(n * LoB)) + 1
        return binom.sf(k_threshold - 1, n, P_pos)  # sf(k-1) = Pr(X >= k)
    
    # 二分法求解 (利用 S7.2.1 提到的严格单调性)
    while M_high_guess - M_low_guess > tol:
        M_mid = (M_low_guess + M_high_guess) / 2
        prob = detection_probability(M_mid)
        if prob < 0.95:
            M_low_guess = M_mid  # 需要更大的 M
        else:
            M_high_guess = M_mid  # 可以尝试更小的 M
    
    return (M_low_guess + M_high_guess) / 2

# =============================================================================
# Figure 1: β-μ 空间响应（精简版，合并原 plot_contour_and_projection）
# =============================================================================



    
    
def figure1_Omega_M_space(
    kappa_values=[0.001, 0.1, 10],         # 覆盖强/中/弱耗竭
    Omega_range=(1e4, 1e14, 100),          # 总位点数从 1e6 到 1e14
    M_range=(1, 1e10, 100),                # 靶标分子数从 1 到 1e10
    N=1,
    b=0.01
):
    
    
    
    
    """
    Figure 1: Ω-M 空间响应。
    正方形子图，坐标轴刻度间距一致，共享色条置于右侧。
    黑色细等高线展示等 p 值分布，斜率恒为 1。
    """
    Omega_vals = np.logspace(np.log10(Omega_range[0]), np.log10(Omega_range[1]), Omega_range[2])
    M_vals = np.logspace(np.log10(M_range[0]), np.log10(M_range[1]), M_range[2])
    Omega_grid, M_grid = np.meshgrid(Omega_vals, M_vals)

    n_cols = len(kappa_values)
    fig_width = 5 * n_cols + 0.8
    fig, axes = plt.subplots(1, n_cols, figsize=(fig_width, 5))
    if n_cols == 1:
        axes = [axes]

    cmap = plt.cm.plasma
    levels = np.linspace(0, 1, 21)
    contour_levels = np.linspace(0.00001, 0.99999, 6)

    log_Omega_min, log_Omega_max = np.log10(Omega_range[0]), np.log10(Omega_range[1])
    log_M_min, log_M_max = np.log10(M_range[0]), np.log10(M_range[1])
    span_Omega = log_Omega_max - log_Omega_min
    span_M = log_M_max - log_M_min
    max_span = max(span_Omega, span_M)

    cf = None
    for ax, kappa in zip(axes, kappa_values):
        p_surface = np.zeros_like(Omega_grid)
        for i in range(len(M_vals)):
            for j in range(len(Omega_vals)):
                Omega = Omega_grid[i, j]
                M = M_grid[i, j]
                beta = Omega / N          # 由 Ω 和分区数反算 β
                _, p, _ = calculate_P_pos(M, kappa, beta, b)
                p_surface[i, j] = p

        cf = ax.contourf(np.log10(Omega_grid), np.log10(M_grid), p_surface,
                         levels=levels, cmap=cmap)
        ax.contour(np.log10(Omega_grid), np.log10(M_grid), p_surface,
                   levels=contour_levels, colors='black', linewidths=0.8, alpha=0.6)
        ax.contour(np.log10(Omega_grid), np.log10(M_grid), p_surface,
                   levels=[0.5], colors='red', linewidths=2, linestyles='--')

        ax.set_xlabel(r'$\log_{10}(\Omega)$ [total sites]', fontsize=14)
        ax.set_ylabel(r'$\log_{10}(M)$ [total targets]', fontsize=14)
        ax.set_title(f'κ = {kappa}', fontsize=14)
        ax.grid(False, alpha=0.3)

        Omega_center = (log_Omega_min + log_Omega_max) / 2
        M_center = (log_M_min + log_M_max) / 2
        ax.set_xlim(Omega_center - max_span/2, Omega_center + max_span/2)
        ax.set_ylim(M_center - max_span/2, M_center + max_span/2)
        ax.set_aspect('equal', adjustable='box')

        # 斜率 1 箭头
        mid_x = np.mean(ax.get_xlim())
        mid_y = np.mean(ax.get_ylim())
        ax.annotate('', xy=(mid_x+1.2, mid_y+1.2), xytext=(mid_x-1.2, mid_y-1.2),
                    arrowprops=dict(arrowstyle='->', color='white', lw=2))
        ax.text(mid_x+0.8, mid_y+1.5, 'slope = 1', color='black', fontsize=11,
                bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='gray', alpha=0.8))

    fig.subplots_adjust(right=0.88)
    cbar_ax = fig.add_axes([0.9, 0.15, 0.02, 0.7])
    cbar = fig.colorbar(cf, cax=cbar_ax, orientation='vertical')
    cbar.set_label('Site Occupancy Probability $p$', fontsize=14)
    cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])

    plt.suptitle('Fig. 1: Response in Ω-M Space — Scale Invariance', 
                 fontsize=16, fontweight='bold', y=0.98)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = os.path.join(output_dir, f"Figure1_Omega_M_space_{timestamp}.svg")
    plt.savefig(filename, format='svg', bbox_inches='tight')
    print(f"Saved Fig. 1: {filename}")
    plt.show()
    return fig



# ================== 主绘图函数（3×2布局，已删除图6，统一实验参数 + 统一字体） ==================
def figure2_li_number_collapse(K_nM=1.0, V_uL=100.0,
                               omega_pairs_physical=None, highlight_M=None):
    """
    3×2 布局（右下角空白）：
    (1) 左上：固定 κ 下不同 N, β 的数据坍缩（标度不变性）
    (2) 右上：不同 κ 对主曲线形态的影响（ξ 空间）
    (3) 左中：化学调控 —— 固定 Ω，变 K（亲和力）时的 p‑M 曲线
    (4) 右中：物理调控 —— 固定 K,V，变 Ω（包被量）时的 p‑M 曲线
    (5) 左下：体积缩放效应
    (6) 右下：空白（已删除）
    所有字体、字号、标题格式与原图完全一致。
    """
    # ---------- 全局物理常数 ----------
    N_A = 6.02214076e23
    K_M = K_nM * 1e-9
    V_L = V_uL * 1e-6
    KVNA = K_M * V_L * N_A               # ≈ 6.02e10 for K=1nM, V=100μL

    # ---------- 公用参数 ----------
    M_common = np.logspace(6, 14, 300)   # 分子数范围 1e6 ~ 1e12，对应 ~0.01nM~10nM
    colors_std = ['#1f77b4', '#ff7f0e', '#2ca02c']

    fig, axes = plt.subplots(3, 2, figsize=(12, 18))
    (ax1, ax2), (ax3, ax4), (ax5, ax6) = axes
    ax6.set_visible(False)               # 删除图(6)
    for ax in [ax1, ax2, ax3, ax4, ax5]:
        ax.grid(False)

    # ==================== (1) 标度坍缩：固定 κ = 1.0 ====================
    kappa_fixed = 1.0
    N_vals = [1e4, 1e5, 1e6]
    beta_vals = [6e6, 6e5, 6e4]          
    xi_theory = np.logspace(-2, 2, 200)
    p_theory = np.zeros_like(xi_theory)
    for i, xi in enumerate(xi_theory):
        p_guess = xi / (1 + xi)
        for _ in range(10):
            f = p_guess/(1 - p_guess) + p_guess/kappa_fixed - xi
            df = 1/(1 - p_guess)**2 + 1/kappa_fixed
            p_guess = np.clip(p_guess - f/df, 1e-12, 1-1e-12)
        p_theory[i] = p_guess
    ax1.plot(xi_theory, p_theory, 'k-', lw=1.5, label=f'Exact (κ={kappa_fixed})')

    markers = ['o', 's', '^', 'D']
    for idx, (N, beta) in enumerate(zip(N_vals, beta_vals)):
        Omega = N * beta
        xi = M_common / (kappa_fixed * Omega)
        mu = M_common / N
        p = np.array([calculate_p(mu_i, kappa_fixed, beta) for mu_i in mu])
        ax1.plot(xi, p, marker=markers[idx], linestyle='None', markersize=4,
                 alpha=0.7, markevery=20,
                 label=f'N={N:.0e}, β={beta:.0e}')
    xi_half = 1 + 0.5 / kappa_fixed
    ax1.scatter(xi_half, 0.5, color='red', s=60, zorder=10,
                label=r'$\xi = 1+0.5/\kappa$ (half‑sat)')
    ax1.set_xscale('log')
    ax1.set_xlabel(r'$\log_{10} \xi$' + '\n' + r'$\xi = M/(\kappa\Omega)$', fontsize=14)
    ax1.set_ylabel(r'$p$', fontsize=14)
    ax1.set_title('(1) Scale Collapse for Fixed κ', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=14)
    ax1.set_ylim(0, 1.05)
    ax1.set_xlim(1e-2, 1e2)
    ax1.axvline(x=1, color='gray', ls=':', lw=1)
    ax1.text(1, 0.92, r'$\xi=1$', transform=ax1.get_xaxis_transform(), fontsize=9, color='gray')
    ax1.set_box_aspect(1)   # 正方形框，使斜率看起来是 45°

    # ==================== (2) 不同 κ 对主曲线（ξ 空间） ====================
    xi_range = np.logspace(-2, 2, 200)
    kappa_compare = [0.1, 1.0, 10.0]
    Omega_ref = KVNA / np.array(kappa_compare)
    N_ref = 1e5
    beta_ref = Omega_ref / N_ref

    ax2.plot(xi_range, xi_range/(1+xi_range), 'k:', lw=2.5,
             label=r'$p = \xi/(\xi+1)$ (κ→∞)')
    for i, kappa in enumerate(kappa_compare):
        p_th = np.zeros_like(xi_range)
        for j, xi in enumerate(xi_range):
            p_guess = xi/(1+xi)
            for _ in range(10):
                f = p_guess/(1-p_guess) + p_guess/kappa - xi
                df = 1/(1-p_guess)**2 + 1/kappa
                p_guess = np.clip(p_guess - f/df, 1e-12, 1-1e-12)
            p_th[j] = p_guess
        ax2.plot(xi_range, p_th, color=colors_std[i], lw=2, label=f'κ = {kappa:.2f}')

        mu_sim = M_common / N_ref
        p_sim = np.array([calculate_p(mu_i, kappa, beta_ref[i]) for mu_i in mu_sim])
        xi_sim = M_common / (kappa * Omega_ref[i])
        ax2.plot(xi_sim, p_sim, 'o', color=colors_std[i], markersize=3, alpha=0.3, markevery=25)

        xi_half = 1 + 0.5 / kappa
        ax2.scatter(xi_half, 0.5, color=colors_std[i], marker='o', s=60, zorder=10,
                    edgecolors='black', linewidth=0.5)
    ax2.set_xscale('log')
    ax2.set_xlabel(r'$\log_{10} \xi$' + '\n' + r'$\xi = M/(\kappa\Omega)$', fontsize=14)
    ax2.set_ylabel(r'$p$', fontsize=14)
    ax2.set_title(r'(2) Effect of κ on the Master Curve ($\xi$ space)', fontsize=14, fontweight='bold')
    ax2.legend(fontsize=14)
    ax2.set_ylim(0, 1.05)
    ax2.set_xlim(1e-2, 1e2)
    ax2.axvline(x=1, color='gray', ls=':', lw=1)
    ax2.text(1, 0.92, r'$\xi=1$', transform=ax2.get_xaxis_transform(), fontsize=9, color='gray')
    ax2.set_box_aspect(1)   # 正方形框，使斜率看起来是 45°

    # ==================== (3) 化学调控：固定 Ω，变 K ====================
    Omega_chem = 1e8
    N_chem = 1e5
    beta_chem = Omega_chem / N_chem       # 1000
    K_vals = [0.1, 1.0, 10.0]            # nM
    kappa_chem = [(k * 1e-9 * V_L * N_A) / Omega_chem for k in K_vals]

    for k_val, kappa, col in zip(K_vals, kappa_chem, colors_std):
        mu = M_common / N_chem
        p = np.array([calculate_p(mu_i, kappa, beta_chem) for mu_i in mu])
        label = f'K={k_val:.2f} nM, κ={kappa:.2f}'
        ax3.semilogx(M_common, p, color=col, lw=2, label=label)
        M_half = (1 + 0.5/kappa) * kappa * Omega_chem
        if M_common[0] <= M_half <= M_common[-1]:
            ax3.scatter(M_half, 0.5, color=col, edgecolors='black', s=60, zorder=10)
    ax3.set_xscale('log')
    ax3.set_xlabel(r'$\log_{10} M$', fontsize=14)
    ax3.set_ylabel('p', fontsize=14)
    ax3.set_title('(3) Tuning κ via affinity (K): a chemical control knob', fontsize=14, fontweight='bold')
    ax3.legend(fontsize=14)
    # ax3.set_xlim(1e6, 1e13)
    ax3.set_ylim(0, 1.05)
    ax3.set_box_aspect(1)   


    
    
        # ==================== (4) 物理调控：固定 K,V，变 Ω ====================
    if omega_pairs_physical is None:
        omega_pairs_physical = [
            (1e11, 1e2),      # Ω=1e13,  κ→0   (强耗竭)
            (1e10, 1e2),      # Ω=1e12
            (5e9, 1e2),       # Ω=5e11
            (1e9, 1e2),       # Ω=1e3
            (1e4, 1e-2),      # Ω=1e2,  κ→∞   (零耗竭，AAI)
        ]
    if highlight_M is None:
        highlight_M = [1e8, 1e10, 1e12]
    
    colors_phys = plt.cm.tab10(np.linspace(0, 1, len(omega_pairs_physical)))
    markers_phys = ['o', 's', '^', 'D', 'v']
    
    # 找出极限情况的索引
    Omega_all = [N * beta for (N, beta) in omega_pairs_physical]
    idx_max_Omega = np.argmax(Omega_all)   # κ→0
    idx_min_Omega = np.argmin(Omega_all)   # κ→∞ (AAI)
    
    for idx, (N, beta) in enumerate(omega_pairs_physical):
        Omega = N * beta
        kappa = KVNA / Omega
        mu = M_common / N
        p = np.array([calculate_p(mu_i, kappa, beta) for mu_i in mu])
    
        # 高亮两种极限曲线
        if idx == idx_max_Omega:
            lw = 3
            label = f'Ω={Omega:.1e}, κ={kappa:.1f} (→0)'
        elif idx == idx_min_Omega:
            lw = 3
            label = f'Ω={Omega:.1e}, κ={kappa:.1f} (AAI, κ→∞)'
        else:
            lw = 2
            label = f'Ω={Omega:.1e}, κ={kappa:.1f}'
    
        ax4.loglog(M_common, p, color=colors_phys[idx], lw=lw, label=label)
    
        # 标记特定 M 值对应的点
        for M_val in highlight_M:
            if M_common[0] <= M_val <= M_common[-1]:
                p_exact = calculate_p(M_val / N, kappa, beta)
                ax4.scatter(M_val, p_exact, color=colors_phys[idx],
                            marker=markers_phys[idx % len(markers_phys)],
                            s=60, edgecolors='black', linewidth=0.8, zorder=10)
                ax4.axvline(x=M_val, color='gray', ls=':', lw=0.8, alpha=0.5)
                ax4.text(M_val*1.1, 2.96, f'M={M_val:.0e}', fontsize=8, color='gray')
    
    # # 添加 κ→0 理想参考线（以最大的 Ω 为例）
    # Omega_max = Omega_all[idx_max_Omega]
    # p_ref = M_common / Omega_max
    # ax4.loglog(M_common, p_ref, '--', color='darkred', lw=2,
    #            label=r'$\kappa \to 0$: $p = M/\Omega$ (slope 1)')
    # 不再单独画理想极限线——κ最小的那条实线本身就是近似极限展示

    
    # 双对数坐标已由 loglog 设置，无需再 set_xscale/set_yscale
    ax4.set_xlabel(r'$\log_{10} M$', fontsize=14)
    ax4.set_ylabel(r'$\log_{10} p$', fontsize=14)
    ax4.set_title('(4) Tuning κ via capture capacity (Ω): a physical control knob', fontsize=14, fontweight='bold')
    ax4.legend(fontsize=10)
    
    # 强制两轴数值跨度（decades）相等，保证斜率视觉为 1
    x_min, x_max = ax4.get_xlim()
    y_min, y_max = ax4.get_ylim()
    
    # 放在你调整完 xlim/ylim 之后
    y_upper = ax4.get_ylim()[1]
    if y_upper > 1:
        ax4.axhspan(1, y_upper, facecolor='gray', alpha=0.08, zorder=0)
        ax4.text(0.01, 0.98, 'No data above p=1\n(expanded to preserve\nslope visual)',
                 transform=ax4.transAxes, ha='left', va='top',
                 fontsize=9, color='gray', style='italic')
    
    
    
    
    
    x_decades = np.log10(x_max) - np.log10(x_min)
    y_decades = np.log10(y_max) - np.log10(y_min)
    max_decades = max(x_decades, y_decades)
    
    x_center = 0.5 * (np.log10(x_min) + np.log10(x_max))
    y_center = 0.5 * (np.log10(y_min) + np.log10(y_max))
    
    ax4.set_xlim(10**(x_center - max_decades/2), 10**(x_center + max_decades/2))
    ax4.set_ylim(10**(y_center - max_decades/2), 10**(y_center + max_decades/2))
    
    ax4.set_box_aspect(1)   # 正方形框，保证斜率 1 的视觉为 45°
   

    # ==================== (5) 体积缩放效应 ====================
    kappa_base = 1.0
    v_factors = [0.5, 1.0, 2.0]
    v_colors = ['#d62728', '#2ca02c', '#1f77b4']
    line_styles = ['--', '-', '-.']
    for v, col, ls in zip(v_factors, v_colors, line_styles):
        kappa_scaled = kappa_base * v
        p_curve = np.zeros_like(xi_range)
        for i, xi in enumerate(xi_range):
            p_guess = xi/(1+xi)
            for _ in range(15):
                f = p_guess/(1-p_guess) + p_guess/kappa_scaled - xi
                df = 1/(1-p_guess)**2 + 1/kappa_scaled
                p_guess = np.clip(p_guess - f/df, 1e-12, 1-1e-12)
            p_curve[i] = p_guess
        ax5.plot(xi_range, p_curve, color=col, lw=2.5, ls=ls,
                 label=f'$V/V_0 = {v}$')
        xi_half = 1 + 0.5 / kappa_scaled
        ax5.scatter(xi_half, 0.5, color=col, edgecolors='black', s=60)
    ax5.text(0.05, 0.9,
             r'$V \uparrow \;\Rightarrow\; \kappa \uparrow,\; \xi \downarrow$' + '\n' +
             r'$V \downarrow \;\Rightarrow\; \kappa \downarrow,\; \xi \uparrow$',
             transform=ax5.transAxes, fontsize=11, bbox=dict(facecolor='white', alpha=0.7))
    ax5.set_xscale('log')
    ax5.set_xlabel(r'$\log_{10} \xi$' + '\n' + r'$\xi = M/(\kappa\Omega)$', fontsize=14)
    ax5.set_ylabel(r'$p$', fontsize=14)
    ax5.set_title('(5) Volume Scaling Effect', fontsize=14, fontweight='bold')
    ax5.legend(fontsize=14)
    ax5.set_ylim(0, 1.05)
    ax5.set_xlim(1e-2, 1e2)
    ax5.axvline(x=1, color='gray', ls=':', lw=1)
    ax5.text(1, 0.92, r'$\xi=1$', transform=ax5.get_xaxis_transform(), fontsize=9, color='gray')
    ax5.set_box_aspect(1)   # 正方形框，使斜率看起来是 45°



    plt.suptitle('Fig. 1d: Universality of the TCS master equation: scale collapse, κ-dependence, and K/Ω tuning',
                 fontsize=18, fontweight='bold')
    plt.tight_layout(rect=[0, 0, 1, 0.97])

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = os.path.join(output_dir, f"Figure2_chem_phys_tuning_{timestamp}.svg")
    plt.savefig(filename, format='svg', bbox_inches='tight')
    print(f"Saved Figure 1d: {filename}")
    plt.show()
    return fig






def figure3_quantification_performance(kappa=1, beta=1, N=1e6, n=1e4, b=0.01):
    """
    生成 Figure 3：定量准确性对比与 CV 范围（保留a-e，删除f，前4个2×2，e居中）。
    """
    M_low, M_high = 1e1, 1e9
    M_test = np.logspace(np.log10(M_low), np.log10(M_high), 500)
    
    # 基准参数下的计算
    M_hat_R1 = np.zeros_like(M_test)
    M_hat_R2 = np.zeros_like(M_test)
    M_hat_R3 = np.zeros_like(M_test)
    CV_R1 = np.zeros_like(M_test)
    CV_R2 = np.zeros_like(M_test)
    CV_R3 = np.zeros_like(M_test)
    
    for i, M in enumerate(M_test):
        P_pos, _, P_specific = calculate_P_pos(M/N, kappa, beta, b)
        M_hat_R1[i] = calculate_M_hat(P_pos, kappa, beta, N, b, 'R1')
        M_hat_R2[i] = calculate_M_hat(P_pos, kappa, beta, N, b, 'R2')
        M_hat_R3[i] = calculate_M_hat(P_pos, kappa, beta, N, b, 'R3')
        
        if P_specific >= 1 - 1e-12:
            CV_R1[i] = 1e6
            CV_R2[i] = 1e6
            CV_R3[i] = 1e6
        else:
            CV_R1[i] = calculate_CV(M, kappa, beta, N, n, b, 'R1')
            CV_R2[i] = calculate_CV(M, kappa, beta, N, n, b, 'R2')
            CV_R3[i] = calculate_CV(M, kappa, beta, N, n, b, 'R3')
    
    # -------------------------- 计算 LoB 和 LoD --------------------------
    LoB_p_pos = calculate_LoB(n, b, method='exact')
    M_hat_LoB = calculate_M_hat(LoB_p_pos, kappa, beta, N, b, 'R1')
    LoD_M = calculate_LoD_exact_R1(kappa, beta, N, n, b, LoB_p_pos)
    # -------------------------------------------------------------------------
    
    # 寻找 LoQ 区间（基于 R1 的 CV ≤ 20%）
    valid = np.isfinite(CV_R1)
    M_valid = M_test[valid]
    CV_valid = CV_R1[valid]
    below = CV_valid <= 0.2
    if np.any(below):
        loq_low = M_valid[np.argmax(below)]
        loq_high = M_valid[len(below) - 1 - np.argmax(below[::-1])]
    else:
        loq_low = loq_high = M_test[-1]
    
    # 【核心修改】布局调整：先用2×2画前4个子图(a-d)，再单独画第5个子图(e)居中
    fig = plt.figure(figsize=(14, 18))  # 调整高度以容纳5个子图
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1])  # 3行2列网格
    
    
    # ========== 子图 1：准确性 (a) ==========
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.loglog(M_test, M_hat_R1, 'k-', lw=3, label='R1: General TCS theory')
    ax1.loglog(M_test, M_hat_R2, 'b--', lw=2, label='R2: Large-β')
    ax1.loglog(M_test, M_hat_R3, 'r:', lw=2, label='R3: low-concentration')
    ax1.loglog(M_test, M_test, 'k:', lw=1, alpha=0.5, label='y = x')
    
    ax1.axhline(M_hat_LoB, color='darkorange', ls='-.', lw=2.5, 
                label=f'LoB (Blank Limit): {M_hat_LoB:.1e}')
    ax1.axvline(LoD_M, color='crimson', ls='-.', lw=2.5, 
                label=f'LoD (Detection Limit): {LoD_M:.1e}')
    ax1.axvline(loq_low, color='gray', ls='--', lw=2, label=f'LoQ Lower: {loq_low:.1e}')
    ax1.axvline(loq_high, color='gray', ls='--', lw=2, label=f'Linear upper bound: {loq_high:.1e}')
    
    ax1.set_xlabel('True Molecule Count M', fontsize=14)
    ax1.set_ylabel('Estimated M̂', fontsize=14)
    ax1.set_title('(1) Quantification Accuracy', fontsize=14, fontweight='bold')
    ax1.legend(fontsize=9, loc='upper left')
    
    # ========== 子图 2：CV 与偏差（基准 n, b）(b) ==========
    ax2 = fig.add_subplot(gs[0, 1])
    bias_R2 = (M_hat_R2 - M_test) / M_test * 100
    bias_R3 = (M_hat_R3 - M_test) / M_test * 100
    
    ax2.semilogx(M_test, CV_R1 * 100, 'k-', lw=3, label='R1 CV (%)')
    ax2.semilogx(M_test, CV_R2 * 100, 'b--', lw=2, label='R2 CV (%)')
    ax2.semilogx(M_test, CV_R3 * 100, 'r--', lw=2, label='R3 CV (%)')
    ax2.semilogx(M_test, bias_R2, 'b:', lw=2, alpha=0.6, label='R2 Bias (%)')
    ax2.semilogx(M_test, bias_R3, 'r:', lw=2, alpha=0.6, label='R3 Bias (%)')
    
    ax2.axhline(20, color='dimgray', ls='--', lw=2, label='20% CV Threshold (LoQ Criterion)')
    ax2.axhline(10, color='darkgray', ls=':', lw=1.2, label='±10% Bias Acceptable Range')
    ax2.axhline(-10, color='darkgray', ls=':', lw=1.2)
    ax2.axvline(LoD_M, color='crimson', ls='-.', lw=2, alpha=0.8, label=f'LoD (Detection Limit): {LoD_M:.1e}')
    ax2.axvline(loq_low, color='gray', ls='--', lw=1.5, alpha=0.7, label=f'LoQ Lower: {loq_low:.1e}')
    ax2.axvline(loq_high, color='gray', ls='--', lw=1.5, alpha=0.7, label=f'Linear Upper Bound: {loq_high:.1e}')
    
    ax2.set_xlabel('True Molecule Count M', fontsize=14)
    ax2.set_ylabel('Bias / CV (%)', fontsize=14)
    ax2.set_title(f'(2) Precision and Bias (n={n:.0e}, b={b})', fontsize=14, fontweight='bold')
    ax2.set_ylim(-50, 200)
    ax2.legend(fontsize=8, ncol=2, loc='upper right')
    

    
    
    # ========== 子图 3：改变 n 对 R1 模型 CV 的影响 ==========
    ax3 = fig.add_subplot(gs[1, 0])
    n_values = [N/1000, N/100, N/10, N]
    n_colors = plt.cm.viridis(np.linspace(0.2, 0.9, len(n_values)))
    model = 'R1'
    
    for ni, color in zip(n_values, n_colors):
        CV_temp = np.zeros_like(M_test)
        for i, M in enumerate(M_test):
            P_pos, _, P_specific = calculate_P_pos(M/N, kappa, beta, b)
            if P_specific >= 1 - 1e-12:
                CV_temp[i] = 1e6
            else:
                CV_temp[i] = calculate_CV(M, kappa, beta, N, ni, b, model)
        
        # 修改：label 中同时显示 n 和 n/N
        ratio = ni / N
        if np.isclose(ni, N):
            label = f'R1, n=N (n/N=1)'
        else:
            label = f'R1, n={ni:.0e} (n/N={ratio:.3f})'
        
        ax3.semilogx(M_test, CV_temp * 100, color=color, lw=2.5, label=label)
    
    ax3.axhline(20, color='dimgray', ls='--', lw=2, label='20% CV Threshold')
    ax3.set_xlabel('True Molecule Count M', fontsize=14)
    ax3.set_ylabel('CV (%)', fontsize=14)
    ax3.set_title('(3) Effect of Replicate Number n on R1 CV', fontsize=14, fontweight='bold')
    ax3.legend(fontsize=9)
    ax3.set_ylim(0, 200)
        

    ax4 = fig.add_subplot(gs[1, 1])
    b_values = [0.001, 0.01, 0.05, 0.1]
    b_colors = plt.cm.plasma(np.linspace(0.2, 0.9, len(b_values)))
    model = 'R1'  # 【核心修改】固定为精确模型R1，和子图c对称
    
    for bi, color in zip(b_values, b_colors):
        CV_temp = np.zeros_like(M_test)
        for i, M in enumerate(M_test):
            P_pos, _, P_specific = calculate_P_pos(M/N, kappa, beta, bi)
            if P_specific >= 1 - 1e-12:
                CV_temp[i] = 1e6
            else:
                CV_temp[i] = calculate_CV(M, kappa, beta, N, n, bi, model)
        label = f'R1, b={bi}'
        ax4.semilogx(M_test, CV_temp * 100, color=color, lw=2.5, label=label)
    
    ax4.axhline(20, color='dimgray', ls='--', lw=2, label='20% CV Threshold')
    ax4.set_xlabel('True Molecule Count M', fontsize=14)
    ax4.set_ylabel('CV (%)', fontsize=14)
    ax4.set_title('(4) Effect of Background b on R1 CV', fontsize=14, fontweight='bold')
    ax4.legend(fontsize=9)
    ax4.set_ylim(0, 200)  # 【优化】去掉负的y轴范围，CV不会为负
    
    # ========== 子图 5：b的估计不确定性对比（精确公式+线性y轴）(e) —— 居中 ==========
    ax5 = fig.add_subplot(gs[2, 0])  # 【核心修改】占满整行，居中
    
    b_range = np.logspace(-4, -1, 200)
    n_blank_list = [N/1000, N/100, N/10, N]
    n_blank_colors = plt.cm.cividis(np.linspace(0.2, 0.9, len(n_blank_list)))
    
    for n_b, color in zip(n_blank_list, n_blank_colors):
        RSE_exact = np.sqrt( (1 - b_range) / (n_b * b_range) ) * 100
        ax5.semilogx(b_range, RSE_exact, color=color, lw=2.5, label=f'n={n_b:.0e}')
    
    ax5.axhline(20, color='dimgray', ls='--', lw=2, label='20% RSE')
    ax5.set_xlabel('True Background b', fontsize=14)
    ax5.set_ylabel('Relative Standard Error of b̂ (%)', fontsize=14)
    ax5.set_title('(5) Uncertainty of b', fontsize=14, fontweight='bold')
    ax5.legend(fontsize=9)
    ax5.set_ylim(0, 150)
 

    # plt.suptitle('Figure 3: Model Validation — Quantification Performance & Background Estimation', 
    #          fontsize=18, y=0.99, fontweight='bold')


    plt.suptitle('Extended Data Fig. 1b: Model Validation — Quantification Performance and Parameter Sensitivity.', 
             fontsize=16, y=0.99, fontweight='bold')
    plt.tight_layout()
    
    
    
    
    # 打印性能指标
    print(f"\n===== Performance Metrics (Strict S7 Formulas) =====")
    print(f"LoB (Exact Binomial): {LoB_p_pos:.4f} (pos. rate) → M̂_LoB = {M_hat_LoB:.1e}")
    print(f"LoD (Exact R1 Bisection): {LoD_M:.1e}")
    print(f"Quantification Range (R1 CV ≤ 20%): {loq_low:.1e} ≤ M ≤ {loq_high:.1e}")
    print(f"\n===== Direct Blank Measurement RSE (n_blank={n:.0e}) =====")
    for bi in b_values:
        rse = np.sqrt( (1 - bi) / (n * bi) ) * 100
        print(f"b={bi}: RSE={rse:.1f}%")
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = os.path.join(output_dir, f"Figure3_quantification_final_{timestamp}.svg")
    plt.savefig(filename, format='svg', bbox_inches='tight')
    print(f"\nSaved Figure 3 (final): {filename}")
    plt.show()
    return fig




def figure_3d_xi_relation_log_axes():
    """
    额外 3D 图：显示 ξ = p/(1-p) + p/κ 的关系曲面。
    x 轴：log10(ξ)，y 轴：log10(κ)，z 轴：p (线性刻度)。
    新增：ξ=1 关键曲线（白色粗线，强制浮在最上层）。
    """
    plt.rcParams.update({
    'font.family': 'serif',
    'font.serif': ['Times New Roman'],
    'mathtext.fontset': 'stix',  # 数学公式字体，和Times New Roman视觉完全统一
})
    # 对数采样 p (避免 p=0, p=1 导致发散)
    p_vals = np.logspace(-4, np.log10(0.9999), 1000)
    kappa_vals = np.logspace(-2, 3, 1000)
    P, K = np.meshgrid(p_vals, kappa_vals)
    Xi = P/(1-P) + P/K

    logXi = np.log10(Xi)
    logK = np.log10(K)
    z_P = P

    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111, projection='3d')
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis._axinfo["grid"]['color'] = (0.9, 0.9, 0.9, 0.5)   # 半透明浅灰
    ax.yaxis._axinfo["grid"]['color'] = (0.9, 0.9, 0.9, 0.5)
    ax.zaxis._axinfo["grid"]['color'] = (0.9, 0.9, 0.9, 0.5)

    # 绘制曲面（降低透明度，强制在最下层）
    surf = ax.plot_surface(logXi, logK, z_P, 
                            cmap='coolwarm', 
                            edgecolor='none', 
                            alpha=0.9,  # 降低透明度
                            zorder=1)   # 强制曲面在最底层

    # ---- 绘制 ξ=1 关键曲线 (白色粗线，强制浮在最上层) ----
    kappa_curve = np.logspace(-2, 3, 500)
    p_curve = (2*kappa_curve + 1 - np.sqrt(4*kappa_curve**2 + 1)) / 2
    x_curve = np.zeros_like(kappa_curve)
    y_curve = np.log10(kappa_curve)
    
    # 强制抬升z值，避免和曲面重合
    p_curve_elevated = p_curve + 0.008
    
    ax.plot(x_curve, y_curve, p_curve_elevated, 
            color='red',  # 白色在plasma色背景上对比度最高
            lw=5,           # 加粗
            zorder=15,      # 强制最高层级
            label=r'$\xi = 1$')
    
    # 文字标注（同样抬升z值）
    idx_text = np.argmin(np.abs(kappa_curve - 1.0))
    ax.text(0, 0, p_curve_elevated[idx_text]+0.03, 
            r'$\xi = 1$', 
            color='red', 
            fontsize=14, 
            weight='bold',
            zorder=20)

    # ---- 坐标轴设置 ----
    ax.set_xlabel(r'$\log_{10}\xi$', fontsize=14)
    xi_ticks = [0.01, 0.1, 1, 10, 100, 1000]
    xi_tick_log = np.log10(xi_ticks)
    ax.set_xticks(xi_tick_log)
    ax.set_xticklabels([f'{val:g}' for val in xi_ticks])

    ax.set_ylabel(r'$\log_{10}\kappa$', fontsize=14)
    kappa_ticks = [0.001, 0.01, 0.1, 1, 10, 100, 1000]
    kappa_tick_log = np.log10(kappa_ticks)
    ax.set_yticks(kappa_tick_log)
    ax.set_yticklabels([f'{val:g}' for val in kappa_ticks])

    p_ticks = np.arange(0, 1.01, 0.2)
    ax.set_zticks(p_ticks)
    ax.set_zticklabels([f'{t:.1f}' for t in p_ticks])
    ax.set_zlim(0, 1)

    # ax.set_title(r'Figure 1c. 3D View: $\xi = \frac{p}{1-p} + \frac{p}{\kappa}$', 
    #              fontsize=16, fontweight='bold', pad=0)
    ax.set_title(r'$\mathbf{Figure\ 1c.}\ \xi = \frac{p}{1-p} + \frac{p}{\kappa}$', 
             fontsize=16, 
             fontweight='normal',  # 【注意】这里设为normal，因为我们用\mathbf{}单独加粗了
             pad=0,
             fontfamily='Times New Roman')
    ax.view_init(elev=30, azim=-110)
    
    # 添加图例
    ax.legend(loc='upper left', fontsize=14)

    cbar = fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10, pad=0.01)
    cbar.set_label(r'$p$', fontsize=14)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    filename = os.path.join(output_dir, f"Figure_extra_3d_xi_relation_logaxes_{timestamp}.svg")
    plt.savefig(filename, format='svg', bbox_inches='tight')
    print(f"Saved Extra 3D Figure with CLEAR ξ=1 curve: {filename}")
    plt.show()
    return fig



# =============================================================================
# 主程序
# =============================================================================
if __name__ == "__main__":
    print("="*60)
    print("Generating three core figures for Ω-κ framework")
    print("="*60)
    
    print("\n[1/3] Figure 1: β-μ response surfaces...")
    fig1 = figure1_Omega_M_space(kappa_values=[0.1, 100])
    
    print("\n[2/3] Figure 2: ξ scale invariance...")
    fig2 = figure2_li_number_collapse()
    
    print("\n[3/3] Figure 3: Quantification performance...")
    fig3 = figure3_quantification_performance()
    
    print("\n" + "="*60)
    print("All figures saved in:", output_dir)
    print("="*60)
    print("\n[Extra] Generating 3D ξ relation surface with log p and log ξ axes...")
    fig_extra_log = figure_3d_xi_relation_log_axes()
    # 新增：真实物理情形 —— 固定 K·V·N_A，Ω 变化导致 κ 变化
    # figure_real_omega_effect(KVNA=1e6)   # 您可以根据自己的体系修改 KVNA 的值
    # figure_real_omega_effect()