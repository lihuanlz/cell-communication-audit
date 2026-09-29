# -*- coding: utf-8 -*-
"""
代码 5：FHN 剂量扫描与计数律（手稿动态侧缺失分析）
==================================================
目的：在确定性 FitzHugh–Nagumo 可兴奋约化上补齐三条数值验证
  ① 振幅不变性（引理 2）：宽剂量扫描，振幅缓升，α = ∂lnA/∂lnD ≈ 0.014
     —— 比 NF 模拟签名（α ≳ 0.3）低一个数量级以上
  ② 计数律（引理 3）：N ∝ ln D₀，斜率 ≈ τ_r/T ≈ 9.3，R² > 0.9
  ③ 周期自洽：T_推 = τ_r/斜率 与直接实测的中位脉冲间隔 T_测 对比（相对偏差 <10%）

模型（手稿方法节）：
  dv/dt = sf·(v − v³/3 − w + I(D,t))      ← 快变量（读作 p53）
  dw/dt = sf·ε·(v + a − b·w)              ← 慢恢复变量
  I(D,t) = (0.40+0.09D)·t/(t+0.5)·exp(−t/τ_r)，τ_r = 20（升起后指数修复）
参数：a=0.7, b=0.8, ε=0.08, sf=20；确定性积分（本脚本不加噪声）。
注：剂量只经"初始损伤幅值" D₀ = 0.40+0.09·D 进入模型，计数律自变量取 ln D₀。

运行：python3 代码5_FHN剂量扫描与计数律.py   （约 1–2 分钟，输出三联图 PNG）
"""
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

# ---------- 模型参数（与手稿方法节严格一致） ----------
a, b, eps = 0.7, 0.8, 0.08            # FHN 标准参数（可兴奋区）
sf   = 20.0                           # 时间尺度加速因子
tau_r = 20.0                          # 损伤修复时间常数
T_END, DT = 100.0, 0.001              # 积分窗（足够长使 I 衰减回阈值以下）与输出步长
D0_of = lambda D: 0.40 + 0.09*D       # 初始损伤幅值 D₀（计数律自变量）

def I_inj(D, t):
    """损伤输入：升起后指数修复。"""
    return D0_of(D) * (t/(t+0.5)) * np.exp(-t/tau_r)

# —— 初值平衡化：无输入下从 (v,w)=(−1,0) 弛豫到真静息态 ——
def rhs_rest(t, y):
    v, w = y
    return [sf*(v - v**3/3 - w), sf*eps*(v + a - b*w)]
y_rest = solve_ivp(rhs_rest, (0, 10), [-1.0, 0.0],
                   rtol=1e-11, atol=1e-13).y[:, -1]
print(f'平衡化静息态：v*={y_rest[0]:.4f}, w*={y_rest[1]:.4f}\n')

def simulate(D):
    """确定性积分一条剂量-D 轨迹并检测脉冲（prominence≥1.0，最小间隔 1.0 au）。"""
    def rhs(t, y):
        v, w = y
        return [sf*(v - v**3/3 - w + I_inj(D, t)),
                sf*eps*(v + a - b*w)]
    nt = int(T_END/DT) + 1
    sol = solve_ivp(rhs, (0, T_END), y_rest,
                    t_eval=np.linspace(0, T_END, nt), rtol=1e-9, atol=1e-12)
    tt, v = sol.t, sol.y[0]
    pk, _ = find_peaks(v, prominence=1.0, distance=int(1.0/DT))
    return tt, v, pk

# ============================================================
# 第一部分（引理 2）：振幅不变性 —— 宽剂量扫描
# ============================================================
# 振幅定义：手稿读数为脉冲峰高 A = 峰均值（v 以 0 为基线的峰读数；
# 峰–谷全幅 ≈ 3.77–3.90 同样缓升，仅作对照打印，α 结论不变）。
D_amp = np.geomspace(0.42, 55.0, 16)          # 跨度 ≈ 130 倍（≥3.5 倍要求）
amp, nsp, famp = [], [], []
for D in D_amp:
    tt, v, pk = simulate(D)
    nsp.append(len(pk))
    amp.append(v[pk].mean() if len(pk) else np.nan)        # 峰均值（峰高）
    famp.append(v[pk].mean()-v.min() if len(pk) else np.nan)  # 峰–谷全幅（对照）
amp, nsp, famp = np.array(amp), np.array(nsp), np.array(famp)
ok = ~np.isnan(amp)
cal = np.polyfit(np.log(D_amp[ok]), np.log(amp[ok]), 1)   # ln A = α·ln D + const
res = np.log(amp[ok]) - np.polyval(cal, np.log(D_amp[ok]))
R2_a = 1 - (res**2).sum()/(((np.log(amp[ok])-np.log(amp[ok]).mean())**2).sum())
alpha = cal[0]

print('【引理 2：振幅不变性】')
print(f'{"剂量 D":>8} {"脉冲数":>6} {"振幅A(峰高)":>12} {"全幅(对照)":>11}')
for D, n, A, F in zip(D_amp, nsp, amp, famp):
    print(f'{D:8.2f} {n:6d} {A:12.3f} {F:11.3f}')
print(f'振幅范围 {np.nanmin(amp):.2f} → {np.nanmax(amp):.2f}，'
      f'全范围变化 {(np.nanmax(amp)/np.nanmin(amp)-1)*100:.1f}%（目标 ≈6–7%）')
print(f'α = ∂lnA/∂lnD = {alpha:.4f}（ln–ln 拟合 R²={R2_a:.3f}；'
      f'目标 ≈0.014，|α|<0.05 通过，须比 NF 签名低一个数量级）')
print(f'判定：{"通过" if abs(alpha) < 0.05 else "未通过"}\n')

# ============================================================
# 第二部分（引理 3）：计数律 N ∝ ln D₀
# ============================================================
D_cnt = np.geomspace(0.45, 17.0, 24)          # 从刚好超阈值到 N≈15–17（饱和前）
ncnt, ipis = [], {}
for D in D_cnt:
    tt, v, pk = simulate(D)
    ncnt.append(len(pk))
    ipis[D] = np.diff(tt[pk])
ncnt = np.array(ncnt)
lnD0 = np.log(D0_of(D_cnt))
win = ncnt >= 2                               # 持续计数窗
cN = np.polyfit(lnD0[win], ncnt[win], 1)      # N = 斜率·ln D₀ + 截距
res = ncnt[win] - np.polyval(cN, lnD0[win])
R2_N = 1 - (res**2).sum()/(((ncnt[win]-ncnt[win].mean())**2).sum())
slope = cN[0]
nu = slope/np.log(10)                         # 每十倍剂量的脉冲数（手稿定义）

print('【引理 3：计数律 N ∝ ln D₀】')
print(f'{"剂量 D":>8} {"D₀":>7} {"ln D₀":>8} {"N":>4}')
for D, n in zip(D_cnt, ncnt):
    print(f'{D:8.2f} {D0_of(D):7.3f} {np.log(D0_of(D)):8.3f} {n:4d}')
print(f'持续计数窗（N≥2，{win.sum()} 点，N 范围 {ncnt[win].min()}–{ncnt[win].max()}）：')
print(f'斜率 = {slope:.2f}（目标 ≈9.3±1.5），R² = {R2_N:.4f}（目标 >0.9），'
      f'ν = 斜率/ln10 = {nu:.2f} 脉冲/十倍剂量（目标 ≈4.05）')
print(f'判定：{"通过" if abs(slope-9.3) <= 1.5 and R2_N > 0.9 else "未通过"}\n')

# ============================================================
# 第三部分：周期自洽（方法节"计数律"段）
# ============================================================
T_pred = tau_r/slope                          # 由计数律反推的周期
ipi_pool = np.concatenate([v for v in ipis.values() if len(v)])
T_meas = np.median(ipi_pool)                  # 直接实测的中位脉冲间隔
dev = abs(T_pred-T_meas)/T_meas

print('【周期自洽】')
print(f'T_推 = τ_r/斜率 = {tau_r:.0f}/{slope:.2f} = {T_pred:.2f} au（参考值 2.14）')
print(f'T_测 = 全部脉冲间隔中位数 = {T_meas:.2f} au（参考值 2.30）')
print(f'相对偏差 = {dev*100:.1f}%（<10% 通过）')
print(f'判定：{"通过" if dev < 0.10 else "未通过"}\n')

# ============================================================
# 三联图
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))

ax = axes[0]                                  # (a) 代表剂量轨迹
for D, c, lab in [(0.6, '#1f77b4', '低剂量 D=0.6（少脉冲）'),
                  (8.0, '#d62728', '高剂量 D=8.0（多脉冲）')]:
    tt, v, pk = simulate(D)
    ax.plot(tt, v, lw=0.9, color=c, label=f'{lab}，N={len(pk)}')
ax.set_xlabel('时间 (au)'); ax.set_ylabel('v（p53 读数）')
ax.legend(fontsize=9); ax.set_xlim(0, 40)
ax.set_title('(a) 确定性轨迹：低剂量 vs 高剂量')

ax = axes[1]                                  # (b) 振幅 vs 剂量（ln–ln 拟合）
ax.plot(D_amp[ok], amp[ok], 'o', color='#ff7f0e', ms=6, label='仿真振幅（峰高）')
Df = np.geomspace(D_amp[ok].min(), D_amp[ok].max(), 100)
ax.plot(Df, np.exp(cal[1])*Df**cal[0], 'k--', lw=1.2,
        label=f'幂律拟合：α = {alpha:.3f}')
ax.set_xscale('log'); ax.set_xlabel('损伤剂量 D（对数轴）'); ax.set_ylabel('振幅 A')
ax.legend(fontsize=9); ax.set_title(f'(b) 振幅不变性：α = {alpha:.3f} ≪ NF 签名')

ax = axes[2]                                  # (c) N vs ln D₀ + 线性拟合
ax.plot(lnD0, ncnt, 's', color='#2ca02c', ms=6, label='仿真脉冲数')
xf = np.linspace(lnD0[win].min(), lnD0[win].max(), 50)
ax.plot(xf, np.polyval(cN, xf), 'k--', lw=1.2,
        label=f'线性拟合：斜率 = {slope:.2f}，R² = {R2_N:.3f}')
ax.text(0.03, 0.97, f'$\\nu$ = {nu:.2f} 脉冲/十倍剂量\n'
                    f'T推 = {T_pred:.2f} au\nT测 = {T_meas:.2f} au（偏差 {dev*100:.1f}%）',
        transform=ax.transAxes, va='top', fontsize=9,
        bbox=dict(boxstyle='round', fc='#f0f0f0', ec='#999'))
ax.set_xlabel('ln $D_0$（$D_0$ = 0.40+0.09·D）'); ax.set_ylabel('脉冲数 N')
ax.legend(fontsize=9, loc='lower right'); ax.set_title('(c) 计数律：$N \\propto \\ln D_0$')

plt.tight_layout()
plt.savefig('FHN计数律验证图.png', dpi=200)
print('已保存 FHN计数律验证图.png')
