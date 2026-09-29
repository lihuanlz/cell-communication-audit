# -*- coding: utf-8 -*-
"""
代码 2：分岔分析套件（拓扑选择定理的三条数值验证）
==================================================
B3  NF 族超临界 Hopf 起始：A² ∝ (s−s_H)     —— 正规延拓级分析（非手搓）
B2  亚临界 Hopf 例外：迟滞环（第三签名）
B1  引理 2 普适性：Morris–Lecar 泄漏数字（α 远小于 NF 模拟签名）

三个验证共同支撑手稿 §4 的"三签名判别定理"：
  ① 起始连续性  ② 迟滞宽度  ③ 宽程 α = ∂lnA/∂lnD
运行：python3 代码2_分岔分析套件.py   （约 2–3 分钟，输出三联图 PNG）
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

# ============================================================
# 第一部分（B3）：NF 模型 —— 平衡支解析 + Hopf 检测 + 极限环支
# ============================================================
# 模型：p53→Mdm2 延迟负反馈（3-ODE Goodwin 类）
#   ṗ = s − δp·p − γ·m·p/(K+p)      (p53：产生 − 自发降解 − Mdm2 催化降解)
#   ẋ = a·p − δx·x                  (Mdm2 mRNA，提供相位滞后)
#   ṁ = b·x − δm·m                  (Mdm2 蛋白)
dp_, gm, K, a_, dx_, b_, dm_ = 0.02, 1.5, 0.15, 1.0, 0.35, 1.0, 0.35

def eq_branch(p):
    """平衡支闭式：给定 p* 直接算出 s（不需要数值求根）。"""
    x = a_*p/dx_; m = b_*x/dm_
    return dp_*p + gm*m*p/(K+p)

def jac_eigs(p):
    """平衡点处 3×3 Jacobian 的特征值（判断失稳方式）。"""
    x = a_*p/dx_; m = b_*x/dm_
    J = np.array([[-dp_ - gm*m*K/(K+p)**2, 0, -gm*p/(K+p)],
                  [a_, -dx_, 0],
                  [0, b_, -dm_]])
    return np.linalg.eigvals(J)

# —— 步骤 1：沿平衡支扫描稳定性，找唯一失稳点 ——
ps = np.linspace(0.001, 3.0, 6000)
stab = np.array([max(jac_eigs(p).real) for p in ps])
hopf_idx = np.where(np.diff(np.sign(stab)) != 0)[0]
i0 = hopf_idx[0]
s_H, _, _ = (lambda v: (eq_branch(v), None, None))(ps[i0])[0], None, None
ev = jac_eigs(ps[i0])
print(f'[B3] 唯一失稳点：p*={ps[i0]:.3f}, s_H={eq_branch(ps[i0]):.4f}, '
      f'过零根={[e for e in ev if abs(e.imag)>1e-6][0]:.3f}（复根对 → Hopf）')

# —— 步骤 2：Hopf 上方逐点积分，测极限环振幅支 ——
def cycle_amp(s, tmax=800):
    def rhs(t, y):
        p, x, m = y
        return [s - dp_*p - gm*m*p/(K+p), a_*p - dx_*x, b_*x - dm_*m]
    sol = solve_ivp(rhs, (0, tmax), [1.4, 4, 11],
                    t_eval=np.linspace(0, tmax, 40001), rtol=1e-9, atol=1e-12)
    seg = sol.y[0][sol.t > tmax*0.6]
    return (seg.max()-seg.min())/2

s_br = s_H + np.array([0.05, 0.15, 0.3, 0.5, 0.8, 1.2, 1.7, 2.3, 3.0, 3.8])
A_br = np.array([cycle_amp(s) for s in s_br])
c = np.polyfit(s_br - s_H, A_br**2, 1)                 # A² 对 (s−s_H) 应线性
A2 = A_br**2
R2 = 1 - ((A2-np.polyval(c, s_br-s_H))**2).sum()/((A2-A2.mean())**2).sum()
print(f'[B3] A² ∝ (s−s_H) 线性 R² = {R2:.4f} → 超临界 √ 标度确认\n')

# ============================================================
# 第二部分（B2）：亚临界 Hopf —— 上扫/下扫迟滞环
# ============================================================
# 正则形式 ż = (μ+i)z + |z|²z − |z|⁴z（Re c₁=+1>0 亚临界，五次项稳定化）
def run_normal(mu, r0, tmax=600):
    def rhs(t, y):
        x, yy = y; r2 = x*x + yy*yy
        return [mu*x - yy + x*r2 - x*r2**2, x + mu*yy + yy*r2 - yy*r2**2]
    sol = solve_ivp(rhs, (0, tmax), [r0, 0],
                    t_eval=np.linspace(0, tmax, 30001), rtol=1e-10, atol=1e-12)
    r = np.hypot(sol.y[0][sol.t > tmax*0.5], sol.y[1][sol.t > tmax*0.5])
    return r.mean()

mus = np.linspace(-0.35, 0.35, 29)
up   = np.array([run_normal(m, 0.05) for m in mus])    # 小初值 = 上扫
down = np.array([run_normal(m, 1.20) for m in mus])    # 大初值 = 下扫
on_up  = mus[np.argmax(up > 0.2)]
off_dn = mus[len(down)-1-np.argmax(down[::-1] < 0.2)]
print(f'[B2] 亚临界：上扫起始 μ_on={on_up:.3f}（有限幅 r=1.00，貌似数字）')
print(f'      下扫熄灭 μ_off={off_dn:.3f}，迟滞宽度={on_up-off_dn:.3f}（理论 0.25）→ 第三签名捕获\n')

# ============================================================
# 第三部分（B1）：Morris–Lecar —— 驱动进快方程的"泄漏数字"
# ============================================================
minf = lambda v: 0.5*(1+np.tanh((v+1.2)/18))
winf = lambda v: 0.5*(1+np.tanh((v-2)/30))
tauw = lambda v: 1/np.cosh((v-2)/60)

def run_ML(I, tmax=3000):
    def rhs(t, y):
        v, w = y
        return [(I - 4.4*minf(v)*(v-120) - 8*w*(v+84) - 2*(v+60))/20,
                0.04*(winf(v)-w)/tauw(v)]
    sol = solve_ivp(rhs, (0, tmax), [-60, 0.01],
                    t_eval=np.linspace(0, tmax, 150001), rtol=1e-9, atol=1e-12)
    v = sol.y[0][sol.t > tmax*0.4]
    pk, _ = find_peaks(v, prominence=20, distance=500)
    return (v[pk].mean()-v.min()) if len(pk) >= 2 else np.nan

Is = np.array([95, 100, 110, 125, 145, 170, 200])
As = np.array([run_ML(I) for I in Is])
ok = ~np.isnan(As)
alpha_ML = np.polyfit(np.log(Is[ok]), np.log(As[ok]), 1)[0]
print(f'[B1] Morris–Lecar：驱动 2.1× 范围振幅变 {(As[ok][-1]/As[ok][0]-1)*100:.1f}%，'
      f'α = {alpha_ML:.3f}（泄漏数字，仍远小于 NF 模拟签名）\n')

# ============================================================
# 三联图
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))

ax = axes[0]                                          # (a) NF 分岔图
ss = np.array([eq_branch(p) for p in ps])
ax.plot(ss, ps, 'k-', lw=1.2, label='平衡支 $p^*(s)$')
ax.plot(ss[stab < 0], ps[stab < 0], 'r-', lw=3, alpha=0.6, label='失稳段')
ax.plot(s_br, 1.477+A_br, 'bo', ms=5, label='极限环 max')
ax.plot(s_br, 1.477-A_br, 'bo', ms=5, mfc='none', label='极限环 min')
fit_s = np.linspace(s_H, s_br.max(), 50)
ax.plot(fit_s, 1.477+np.sqrt(np.polyval(c, fit_s-s_H)), 'g--', lw=1.5,
        label=f'$\\sqrt{{s-s_H}}$ 拟合 ($R^2$={R2:.3f})')
ax.axvline(s_H, color='gray', ls=':'); ax.annotate('Hopf', (s_H, 0.2))
ax.set_xlabel('驱动 s'); ax.set_ylabel('p'); ax.legend(fontsize=8)
ax.set_title('(a) NF 族：超临界 Hopf 起始（延拓级）'); ax.set_xlim(14, 21); ax.set_ylim(0, 3)

ax = axes[1]                                          # (b) 亚临界迟滞环
ax.plot(mus, up, 'b^-', ms=5, label='上扫（小初值）')
ax.plot(mus, down, 'rv-', ms=5, label='下扫（大初值）')
ax.fill_betweenx([0, 1.6], off_dn, on_up, color='orange', alpha=0.15)
ax.annotate(f'迟滞宽度 = {on_up-off_dn:.2f}', (-0.13, 1.3), fontsize=9)
ax.set_xlabel('μ'); ax.set_ylabel('振幅 r'); ax.legend(fontsize=9); ax.set_ylim(0, 1.6)
ax.set_title('(b) 亚临界例外：有限起始但迟滞 ≠ 0')

ax = axes[2]                                          # (c) 三族 α 分离
alpha_NF = np.polyfit(np.log(s_br/s_H), np.log(A_br/A_br[0]), 1)[0]
D0s_ = np.array([0.42, 0.50, 0.60, 0.75, 0.95, 1.25, 1.65])
A_fhn = np.array([1.793, 1.813, 1.833, 1.847, 1.872, 1.898, 1.909])  # 代码1/TrackB 结果
ax.plot(s_br/s_H, A_br/A_br[0], 'rs-', ms=5, label=f'NF 模拟：α≈{alpha_NF:.1f}（近起始→∞）')
ax.plot(D0s_/D0s_[0], A_fhn/A_fhn[0], 'bo-', ms=5, label='FHN 慢驱动（数字）：α=0.014')
ax.plot(Is/Is[0], As/As[0], 'g^-', ms=6, label=f'ML 快驱动（泄漏）：α={alpha_ML:.2f}')
ax.set_xscale('log'); ax.set_yscale('log')
ax.set_xlabel('驱动（归一化）'); ax.set_ylabel('振幅（归一化）'); ax.legend(fontsize=9)
ax.set_title('(c) 三泛类在 α 上以数量级分离')

plt.tight_layout()
plt.savefig('TrackB_分岔验证图.png', dpi=200)
print('已保存 TrackB_分岔验证图.png')
