# -*- coding: utf-8 -*-
"""
图 3 | 动态编码器（手稿图 3 实现，四联）
  (a) NF 分岔图与 √ 起始（Krishna soft 模式振幅支 + √(C−C_H) 拟合）
  (b) 亚临界迟滞环——第三签名（正则形式解析支）
  (c) 三类 α 分离：NF / FHN / Morris–Lecar（模型侧签名汇总）
  (d) 计数律 N ∝ ln $D_0$ 与周期自洽（FHN 仿真，代码5 协议精简版）
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

# ---------- (a) Krishna soft 模式振幅支 ----------
A_, EPS, B, DEL = 0.007, 2e-5, 3.0, 0.005
C_H = 2.62e-3
def krishna(t, y, C):
    Nn, Im, I = y
    return [A_*(1-Nn)/(EPS+I) - B*I*Nn/(DEL+Nn),
            Nn**2 - Im,
            Im - C*(1-Nn)*I/(EPS+I)]
Cs = np.array([3.0, 3.5, 4.5, 6, 8, 12, 20, 40])*1e-3
amps = []
for C in Cs:
    sol = solve_ivp(krishna, (0, 4000), [0.5, 0.1, 0.01], args=(C,),
                    dense_output=True, rtol=1e-8, atol=1e-10)
    t = np.linspace(3000, 4000, 20000)
    v = sol.sol(t)[0]
    pk, _ = find_peaks(v, prominence=1e-4)
    tr, _ = find_peaks(-v, prominence=1e-4)
    amps.append((v[pk].mean()-v[tr].mean()) if len(pk) and len(tr) else 0.0)
amps = np.array(amps)
ok = amps > 0
cf = np.polyfit(np.sqrt(Cs[ok]-C_H), amps[ok], 1)
R2 = 1 - ((amps[ok]-np.polyval(cf, np.sqrt(Cs[ok]-C_H)))**2).sum()/((amps[ok]-amps[ok].mean())**2).sum()

# ---------- (b) 亚临界正则形式迟滞环 ----------
c1 = 2.0
mu = np.linspace(-1.2, 1.5, 400)
disc = c1**2 + 4*mu
A_stable = np.sqrt((c1 + np.sqrt(np.maximum(disc, 0)))/2)
A_unstab = np.sqrt(np.maximum((c1 - np.sqrt(np.maximum(disc, 0)))/2, 0))
mu_fold = -c1**2/4

# ---------- (c) 三类 α 汇总 ----------
sys_names = ['NF-κB soft', 'NF-κB spiky\n(上升肢)', 'p53-FHN', 'Morris–Lecar']
alphas = [0.714, 0.549, 0.014, -0.258]
cols = ['#d62728', '#d62728', '#1f77b4', '#2ca02c']

# ---------- (d) 计数律（代码5 协议精简：12 剂量） ----------
a, b_, eps, sf, tau_r = 0.7, 0.8, 0.08, 20.0, 20.0
def fhn_rhs(t, y, D):
    v, w = y
    I = (0.40+0.09*D)*(t/(t+0.5))*np.exp(-t/tau_r)
    return [sf*(v - v**3/3 - w + I), sf*eps*(v + a - b_*w)]
y0 = solve_ivp(lambda t, y: [sf*(y[0]-y[0]**3/3-y[1]), sf*eps*(y[0]+a-b_*y[1])],
               (0, 10), [-1.0, 0.0], rtol=1e-10, atol=1e-12).y[:, -1]
Ds = np.geomspace(0.45, 17.0, 24)
Ns = []
for D in Ds:
    sol = solve_ivp(fhn_rhs, (0, 100), y0, args=(D,),
                    t_eval=np.linspace(0, 100, 100001), rtol=1e-9, atol=1e-12)
    pk, _ = find_peaks(sol.y[0], prominence=1.0, distance=1000)
    Ns.append(len(pk))
Ns = np.array(Ns)
lnD0 = np.log(0.40+0.09*Ds)
win = Ns >= 2
cN = np.polyfit(lnD0[win], Ns[win], 1)
R2N = 1 - ((Ns[win]-np.polyval(cN, lnD0[win]))**2).sum()/((Ns[win]-Ns[win].mean())**2).sum()

# ---------- 绘图 ----------
fig, axes = plt.subplots(2, 2, figsize=(12.5, 9))

ax = axes[0, 0]
eq_C = np.geomspace(5e-4, 0.05, 200)
ax.axvline(C_H, color='gray', ls=':', lw=1)
ax.plot(Cs[ok], amps[ok], 'o', color='#1f77b4', ms=7, label='极限环振幅（数值）')
xf = np.linspace(C_H*1.001, Cs[ok].max(), 100)
ax.plot(xf, np.polyval(cf, np.sqrt(xf-C_H)), 'k--', lw=1.4,
        label=f'√(C−C_H) 拟合，R²={R2:.3f}')
ax.set_xscale('log'); ax.set_xlabel('剂量 C（∝ IKK 强度）'); ax.set_ylabel('振幅 A')
ax.annotate(f'Hopf\nC_H={C_H:.4f}', xy=(C_H, 0), xytext=(1.6*C_H, amps.max()*0.55), fontsize=9,
            arrowprops=dict(arrowstyle='->', color='gray'))
ax.legend(fontsize=9); ax.set_title('(a) NF 超临界 Hopf 起始：A ∝ √(C−C_H)')

ax = axes[0, 1]
ax.plot(mu, A_stable**2, 'b-', lw=2, label='稳定极限环 A²')
m2 = mu >= mu_fold
ax.plot(mu[m2], A_unstab[m2]**2, 'r--', lw=2, label='不稳定环（共存）')
ax.plot(mu[mu<=0], np.zeros((mu<=0).sum()), 'g-', lw=2, label='稳定静息态')
ax.plot(mu[(mu>0)&(mu<1.5)], np.zeros(((mu>0)&(mu<1.5)).sum()), 'g--', lw=2)
ax.axvspan(mu_fold, 0, color='orange', alpha=0.15)
ax.annotate('迟滞窗\n（双稳共存）', xy=((mu_fold)/2, 0.5), ha='center', fontsize=10, color='#b36b00')
ax.axvline(mu_fold, color='r', ls=':', lw=1); ax.axvline(0, color='g', ls=':', lw=1)
ax.set_xlabel('驱动 μ'); ax.set_ylabel('A²'); ax.legend(fontsize=9, loc='upper left')
ax.set_title(f'(b) 亚临界起始：有限振幅 + 迟滞（宽度 {-mu_fold:.2f}）——第三签名')

ax = axes[1, 0]
bars = ax.bar(sys_names, alphas, color=cols, alpha=0.85)
ax.axhline(0.5, color='r', ls='--', lw=1.2); ax.text(3.45, 0.51, 'NF 签名 α≳0.5', color='r', fontsize=9, ha='right')
ax.axhline(0.3, color='k', ls=':', lw=1.2); ax.text(3.45, 0.305, '判别边界 0.3', fontsize=9, ha='right')
ax.axhline(0, color='k', lw=0.8)
for b, v in zip(bars, alphas):
    ax.text(b.get_x()+b.get_width()/2, v+0.02 if v>0 else v-0.05, f'{v:+.3f}',
            ha='center', fontsize=10, fontweight='bold')
ax.set_ylabel('α = ∂lnA/∂lnD'); ax.set_ylim(-0.45, 0.9)
ax.set_title('(c) 三类 α 分离：NF（红）≫ EXC（蓝/绿）≈ 0')

ax = axes[1, 1]
ax.plot(lnD0, Ns, 's', color='#2ca02c', ms=7, label='FHN 仿真脉冲数')
xf2 = np.linspace(lnD0[win].min(), lnD0[win].max(), 50)
ax.plot(xf2, np.polyval(cN, xf2), 'k--', lw=1.4,
        label=f'线性拟合：斜率={cN[0]:.1f}，R²={R2N:.3f}')
ax.set_xlabel('ln $D_0$'); ax.set_ylabel('脉冲数 N'); ax.legend(fontsize=9)
ax.set_title(f'(d) 计数律：N ∝ ln $D_0$（ν = {cN[0]/np.log(10):.1f} 次/十倍）')

plt.tight_layout()
plt.savefig('图3_动态编码器.png', dpi=200)
print(f'(a) √(C−C_H) 拟合 R²={R2:.3f}；(d) 计数律斜率={cN[0]:.1f}，R²={R2N:.3f}，ν={cN[0]/np.log(10):.2f}/十倍')
print('已保存 图3_动态编码器.png')
