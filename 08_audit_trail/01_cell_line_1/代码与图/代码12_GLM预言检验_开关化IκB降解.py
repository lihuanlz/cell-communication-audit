# -*- coding: utf-8 -*-
"""
代码 12（v2）：GLM 预言的计算检验——开关化 IκBα 降解是否把 Krishna 模型从 NF 移到 EXC
================================================================================
背景：v0.9 记录 NF-κB 模型-实验落差（Krishna α=0.714 vs Tay 2010 α=0.151）。
GLM 诊断：模型与实验处在不同泛类（NF 自持极限环 vs EXC 阻尼可兴奋），
预言：把 IκBα 的 IKK 驱动降解项从 MM（Hill 1）换为开关型（Hill 4），
模型应从 NF 类移到 EXC 类——α 降至 <0.3、出现 fold/canard 结构。

【须先声明的事实】Krishna 模型 ε=2×10⁻⁵：I/(ε+I) 在数值上已近阶跃。
GLM 修改（n=1→4, K=ε）的边际效果存疑——本代码不预判，如实输出。

方程（Krishna 2006；GLM 修改仅作用于 dI/dt 降解项 degr(I)）：
  dN_n/dt = A(1−N_n)/(ε+I) − B·I·N_n/(δ+N_n)
  dI_m/dt = N_n² − I_m
  dI/dt   = I_m − C(1−N_n)·degr(I)
  degr: n=1 → I/(ε+I)；n=4 → I⁴/(ε⁴+I⁴)

稳态化简：I_m=N²；I*(N) = ε·(N²/(C(1−N)−N²))^{1/n}（需 C(1−N)>N²）；
代回 dN_n=0 一维求根（N 稠密网格 2×10⁵，向量化）。

审计协议（同代码6/9）：平衡支多值（fold）扫描；Hopf 边界（解析 Jacobian）；
窗口内极限环振幅→α；起始 A²∝(C−C_H) 检验；可兴奋性探针（窗口外阶梯扰动，
刻版全或无 vs 渐变阻尼）。
运行：python3 代码12_GLM预言检验_开关化IκB降解.py（约 3–6 分钟）
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
from scipy.signal import find_peaks

A_, EPS = 0.007, 2e-5

def steady_roots(B, d, C, n_hill):
    """向量化一维求根，返回 [(N, Im, I), ...]"""
    Nmax = (-C + np.sqrt(C*C + 4*C))/2
    Ns = np.linspace(1e-9, Nmax*(1-1e-9), 200000)
    gap = Ns**2/(C*(1-Ns) - Ns**2)          # = target/(1-target)
    Iv = EPS*gap if n_hill == 1 else EPS*gap**(1.0/n_hill)
    f = A_*(1-Ns)/(EPS+Iv) - B*Iv*Ns/(d+Ns)
    sc = np.where(np.sign(f[:-1])*np.sign(f[1:]) < 0)[0]
    out = []
    for i in sc:
        g = lambda n: A_*(1-n)/(EPS + (EPS*(n*n/(C*(1-n)-n*n)) if n_hill==1
                     else EPS*(n*n/(C*(1-n)-n*n))**(1.0/n_hill))) \
            - B*(EPS*(n*n/(C*(1-n)-n*n)) if n_hill==1
                 else EPS*(n*n/(C*(1-n)-n*n))**(1.0/n_hill))*n/(d+n)
        try:
            N = brentq(g, Ns[i], Ns[i+1], xtol=1e-13)
            I = EPS*(N*N/(C*(1-N)-N*N)) if n_hill==1 else EPS*(N*N/(C*(1-N)-N*N))**(1.0/n_hill)
            out.append((N, N*N, I))
        except Exception:
            pass
    return out

def eigs(B, d, C, N, Im, I, n_hill):
    if n_hill == 1:
        degr = I/(EPS+I); ddegr = EPS/(EPS+I)**2
    else:
        n = n_hill; In = I**n
        degr = In/(EPS**n+In); ddegr = n*EPS**n*I**(n-1)/(EPS**n+In)**2
    J = np.array([
        [-A_/(EPS+I) - B*I*d/(d+N)**2, 0, -A_*(1-N)/(EPS+I)**2 - B*N/(d+N)],
        [2*N, -1, 0],
        [C*degr, 1, -C*(1-N)*ddegr]])
    return np.linalg.eigvals(J)

def make_rhs(B, d, C, n_hill):
    def rhs(t, y):
        N, Im, I = y
        if n_hill == 1:
            degr = I/(EPS+I)
        else:
            I4 = max(I, 0.0)**n_hill
            degr = I4/(EPS**n_hill + I4)
        return [A_*(1-N)/(EPS+I) - B*I*N/(d+N),
                N*N - Im,
                Im - C*(1-N)*degr]
    return rhs

def hopf_window(B, d, n_hill):
    """扫 C∈[1e-4, 1] 失稳区间，返回 [(Cin, Cout), ...] 与多值点数"""
    Cs = np.logspace(-4, 0, 150)
    flags, multi = [], 0
    for C in Cs:
        eqs = steady_roots(B, d, C, n_hill)
        if len(eqs) > 1: multi += 1
        unst = any(np.any(np.real(eigs(B, d, C, *eq, n_hill)) > 1e-9) for eq in eqs) if eqs else None
        flags.append(unst)
    wins, i = [], 0
    while i < len(Cs):
        if flags[i]:
            j = i
            while j < len(Cs) and flags[j]: j += 1
            wins.append((Cs[i-1] if i > 0 else Cs[0], Cs[j-1]))
            i = j
        else:
            i += 1
    return wins, multi

def limit_cycle_amp(B, d, C, n_hill, y0=None, T=2500.0):
    rhs = make_rhs(B, d, C, n_hill)
    if y0 is None:
        eqs = steady_roots(B, d, C, n_hill)
        if not eqs: return None, None, None
        N, Im, I = eqs[0]
        y0 = [N*1.03, Im, I*1.03]
    sol = solve_ivp(rhs, [0, T], y0, method='LSODA', rtol=1e-8, atol=1e-11, dense_output=True)
    t = np.linspace(T*0.6, T, 12000)
    y = sol.sol(t)
    N = y[0]
    pk, _ = find_peaks(N, prominence=max(N.max()*1e-4, 1e-9))
    if len(pk) < 3: return None, None, y0
    return N[pk[-3:]].mean(), np.diff(t[pk[-3:]]).mean(), [N[pk[-1]], y[1][pk[-1]], y[2][pk[-1]]]

def excitability_probe(B, d, C, n_hill):
    eqs = steady_roots(B, d, C, n_hill)
    if not eqs: return None
    N0, Im0, I0 = eqs[0]
    Ds = np.array([1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 0.1, 0.3])
    peaks = []
    rhs = make_rhs(B, d, C, n_hill)
    for D in Ds:
        sol = solve_ivp(rhs, [0, 3000], [N0+D, Im0, I0], method='LSODA',
                        rtol=1e-8, atol=1e-11)
        peaks.append(sol.y[0].max() - N0)
    return Ds, np.array(peaks), N0

def audit(B, d, n_hill, label):
    print(f"\n{'='*72}\n【{label}】B={B}, δ={d}, Hill n={n_hill}\n{'='*72}")
    wins, multi = hopf_window(B, d, n_hill)
    print(f"平衡支多值 C 网格点：{multi}/150（>0 → fold/双稳存在）")
    print(f"失稳窗口：{[(f'{a:.3e}', f'{b:.3e}') for a, b in wins]}")
    if not wins:
        print("无失稳窗口——该参数组无自持振荡区。")
    else:
        Cin, Cout = wins[0]
        hi = min(Cin*60, Cout*0.98)
        C_scan = np.logspace(np.log10(Cin*1.03), np.log10(hi), 12)
        amps, Cs_ok, pers, y0 = [], [], [], None
        for C in C_scan:
            amp, per, y0 = limit_cycle_amp(B, d, C, n_hill, y0=y0)
            if amp is not None:
                amps.append(amp); Cs_ok.append(C); pers.append(per)
        if len(amps) >= 4:
            la, lc = np.log(amps), np.log(Cs_ok)
            sl, ic = np.polyfit(lc, la, 1)
            r2 = 1 - np.sum((la-(sl*lc+ic))**2)/np.sum((la-la.mean())**2)
            print(f"极限环振幅 {len(amps)} 点：{min(amps):.3g}–{max(amps):.3g}"
                  f"（{max(amps)/min(amps):.2f} 倍跨度）")
            print(f"α = ∂lnA/∂lnC = {sl:+.3f}（R²={r2:.3f}）")
            print(f"周期跨度：{max(pers)/min(pers):.2f} 倍")
            k = min(6, len(amps))
            dc = np.array(Cs_ok[:k]) - Cin
            a2 = np.array(amps[:k])**2
            if np.all(dc > 0) and a2.std() > 0:
                s2, i2 = np.polyfit(dc, a2, 1)
                r2s = 1 - np.sum((a2-(s2*dc+i2))**2)/np.sum((a2-a2.mean())**2)
                print(f"起始检验（近起始 {k} 点 A²∝(C−C_in)）：R²={r2s:.3f}"
                      f"（→1 超临界 √ 起始；差→canard/不连续）")
        else:
            print(f"极限环点不足（{len(amps)}）。")
        # 可兴奋性探针：窗口下边界以下
        C_sub = Cin/3
        pr = excitability_probe(B, d, C_sub, n_hill)
        if pr:
            Ds, pk, N0 = pr
            print(f"可兴奋性探针（C={C_sub:.2e}=窗口下界/3，稳态 N*={N0:.4f}）：")
            print(f"  Δ     = {Ds}")
            print(f"  响应峰= {np.round(pk, 5)}")
            if pk[0] > 0:
                print(f"  最强/最弱响应比 = {pk[-1]/pk[0]:.1f}，扰动跨度比 = {Ds[-1]/Ds[0]:.0f}"
                      f"（两比接近→渐变阻尼螺线=NF 静息侧；响应比≪扰动比→阈值刻版=EXC 候选）")

print("GLM 预言检验：开关化 IκBα 降解（n=1→4）是否把 Krishna 模型移出 NF 类")
print("="*72)
audit(3, 0.005, 1, "原版 soft")
audit(3, 0.005, 4, "GLM 修改版 soft（n=4）")
audit(954.5, 0.029, 1, "原版 spiky（标准参数）")
audit(954.5, 0.029, 4, "GLM 修改版 spiky（n=4）")
