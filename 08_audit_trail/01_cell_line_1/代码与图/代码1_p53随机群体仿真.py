# -*- coding: utf-8 -*-
"""
代码 1：p53 随机可兴奋群体仿真
================================
目的：复现 Lahav 2004 / Batchelor 2008 的四条实验特征
  ① 脉冲数随剂量增加（数字编码）
  ② 振幅与剂量无关（全或无）
  ③ 脉冲间隔（IPI）精确 ≈ 5.5 h
  ④ 首脉冲时刻弥散（相位由噪声触发）

模型：FitzHugh–Nagumo 可兴奋约化（Mönke 2017 类）
  dv/dt = (v − v³/3 − w) × sf + 噪声        ← 快变量（读作 p53）
  dw/dt = ε(v + a − b·w) × sf                ← 慢恢复变量（读作 Mdm2 等）
损伤输入：I(D,t) = (0.40+0.09D)·t/(t+0.5)·exp(−t/20)   ← 升起后经修复衰减

噪声拆两份（论文要点）：
  · 动力学噪声 σ 调低 → 保节拍（IPI 精确）
  · 外源噪声（每细胞增益/起效时间不同）→ 造相位弥散（首脉冲时刻可变）
运行：python3 代码1_p53随机群体仿真.py   （约 1–2 分钟，输出四联图 PNG）
"""
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt

rng = np.random.default_rng(7)          # 固定种子，结果可复现

# ---------- 模型参数 ----------
a, b, eps = 0.7, 0.8, 0.08              # FHN 标准参数（可兴奋区）
sf   = 20.0                             # 时间尺度加速因子 → 周期 ≈ 2 任意单位
sig  = 0.10                             # 动力系统噪声（压低，保节拍）

T, dt = 24.0, 0.001                     # 仿真窗与步长（任意时间单位）
tt = np.arange(0, T, dt)
NC, doses = 200, [0.1, 0.3, 1.0, 3.0]   # 每剂量 200 个细胞，4 个剂量

results = {}
for D in doses:
    # 外源（细胞间）变异：输入增益 ~ 对数正态，起效时间 ~ 正态
    gain = rng.lognormal(0, 0.22, NC)
    ton  = np.clip(rng.normal(0.5, 0.18, NC), 0.1, None)
    v = np.full(NC, -1.0); w = np.full(NC, 0.0)   # 静息态出发
    traj = np.empty((NC, len(tt)), dtype=np.float32)
    for i, t in enumerate(tt):
        traj[:, i] = v
        I = gain * (0.40 + 0.09*D) * (t/(t+ton)) * np.exp(-t/20.0)
        v += sf*(v - v**3/3 - w + I)*dt + sig*np.sqrt(dt)*rng.standard_normal(NC)
        w += sf*eps*(v + a - b*w)*dt
    # 逐细胞脉冲检测与统计
    stats = []
    for c in range(NC):
        pk, prop = find_peaks(traj[c], prominence=1.0, distance=int(1.0/dt), width=0.1/dt)
        tpk = tt[pk]
        stats.append(dict(n=len(pk),
                          amp=traj[c][pk].mean() if len(pk) else np.nan,
                          width=prop['widths'].mean()*dt if len(pk) else np.nan,
                          first=tpk[0] if len(pk) else np.nan,
                          ipi=np.diff(tpk) if len(pk) > 1 else np.array([])))
    results[D] = (traj, stats)

# ---------- 时间标定：令中位 IPI = 5.5 h ----------
all_ipi = np.concatenate([s['ipi'] for _, ss in results.values() for s in ss if len(s['ipi'])])
hour = 5.5 / np.median(all_ipi)
print(f'中位 IPI = {np.median(all_ipi):.2f} 任意单位 → 1 单位 = {hour:.2f} h\n')
print(f'{"剂量":>5} {"脉冲数":>11} {"振幅(CV)":>13} {"脉宽(h)":>11} {"首脉冲(h)":>13} {"IPI(h)":>13}')
for D in doses:
    ss = results[D][1]
    ns   = np.array([s['n'] for s in ss]); amps = np.array([s['amp'] for s in ss])
    wid  = np.array([s['width'] for s in ss])*hour
    fst  = np.array([s['first'] for s in ss])*hour
    ipi  = np.concatenate([s['ipi'] for s in ss if len(s['ipi'])])*hour
    print(f'{D:5.1f} {ns.mean():4.1f}±{ns.std():<4.1f} {np.nanmean(amps):5.2f} '
          f'({np.nanstd(amps)/np.nanmean(amps)*100:3.1f}%) {np.nanmean(wid):5.1f}±{np.nanstd(wid):<4.1f} '
          f'{np.nanmean(fst):5.1f}±{np.nanstd(fst):<4.1f} {np.nanmean(ipi):5.2f}±{np.nanstd(ipi):<4.2f}')

# ---------- 四联图 ----------
fig, axes = plt.subplots(2, 2, figsize=(13, 8))
t_h = tt*hour
ax = axes[0, 0]                                   # (a) 示例轨迹
for j, D in enumerate([0.3, 3.0]):
    for c in range(3):
        ax.plot(t_h, results[D][0][c]+j*3.5, lw=0.7, alpha=0.85,
                color=['#888', '#1f77b4'][j], label=f'D={D}' if c == 0 else None)
ax.set_xlabel('时间 (h)'); ax.set_ylabel('p53 水平（错位显示）'); ax.legend()
ax.set_title('(a) 单细胞脉冲轨迹：低剂量 vs 高剂量'); ax.set_xlim(0, 24)

ax = axes[0, 1]                                   # (b) 计数 ∝ 剂量
bp = ax.boxplot([[s['n'] for s in results[D][1]] for D in doses],
                positions=doses, widths=0.15, patch_artist=True)
for p in bp['boxes']: p.set_facecolor('#aec7e8')
ax.set_xlabel('损伤剂量'); ax.set_ylabel('脉冲数 N'); ax.set_title('(b) 数字编码：N ∝ 剂量')

ax = axes[1, 0]                                   # (c) 振幅与剂量无关
bp = ax.boxplot([[s['amp'] for s in results[D][1] if not np.isnan(s['amp'])] for D in doses],
                positions=doses, widths=0.15, patch_artist=True)
for p in bp['boxes']: p.set_facecolor('#ffbb78')
ax.set_xlabel('损伤剂量'); ax.set_ylabel('脉冲振幅'); ax.set_ylim(1.5, 2.1)
ax.set_title('(c) 振幅与剂量无关（全或无）')

ax = axes[1, 1]                                   # (d) 相位弥散 vs 节拍精确
ss = results[3.0][1]
fst = np.array([s['first'] for s in ss if not np.isnan(s['first'])])*hour
ipi = np.concatenate([s['ipi'] for s in ss if len(s['ipi'])])*hour
ax.hist(fst, bins=30, alpha=0.6, color='#d62728', label=f'首脉冲时刻 (SD={fst.std():.1f} h)')
ax.hist(ipi, bins=30, alpha=0.6, color='#2ca02c', label=f'脉冲间隔 (SD={ipi.std():.2f} h)')
ax.set_xlabel('时间 (h)'); ax.set_ylabel('细胞数'); ax.legend()
ax.set_title('(d) D=3.0：相位弥散 vs 节拍精确')

plt.tight_layout()
plt.savefig('p53脉动仿真.png', dpi=200)
print('\n已保存 p53脉动仿真.png')
