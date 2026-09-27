# -*- coding: utf-8 -*-
"""
数值核查 S2.7 · Gamma 间隔更新过程族（工作笔记 §8 的代数防错，非注册性质）
=====================================================================
目的：钉死两个解析结论
  V3a 参数 Gamma(k,lam) 更新族：似然比在 (N, ΣX, ΠX) 纤维上严格恒定
       => (N,ΣX,ΠX) 是 R_y 的真粗化且充分 => R_y 不极小（极小性条款证伪）
  V3b 自由形状更新族：似然比对间隔置换不变
       => 间隔多重集充分（经典：全族非参数 iid 下次序统计量极小充分）
       => 极小充分 = 多重集 ≠ 有序序列（极小性条款再次证伪）
两行阴性 => H4 退役，定理按 T' 口径重述（见工作笔记 §8）
运行：python 数值核查_S27_Gamma族_V3.py   （秒级，纯 CPU）
"""
import numpy as np
from scipy.special import gammaln
from scipy.optimize import brentq
from scipy.stats import lognorm, weibull_min

rng = np.random.default_rng(20260814)
print("=" * 64)
print("V3a 参数 Gamma(k,lam) 更新族：LR 是否只依赖 (N, sumX, prodX)")
print("=" * 64)

def logL_gamma(X, k, lam):
    N = len(X)
    return N * k * np.log(lam) - N * gammaln(k) + (k - 1) * np.sum(np.log(X)) - lam * np.sum(X)

th1, th2 = (2.0, 1.0), (3.5, 0.7)
N = 8
X = rng.gamma(2.0, 1.0, N)

# 构造 X'：与 X 的 (sum, sum log) 全同，但内部模式/排序不同
a, b, c = X[0], X[1], X[2]
u = 0.3
g = lambda v: a * np.exp(u) + b * np.exp(v) + c * np.exp(-u - v) - (a + b + c)
vs = np.linspace(-2, 2, 4001)
gs = [g(v) for v in vs]
v_star = None
for i in range(1, len(vs)):
    if gs[i - 1] * gs[i] < 0:
        r = brentq(g, vs[i] - 0.005, vs[i] + 0.005)
        if abs(r - u) > 1e-3 and abs(r + u) > 1e-3 and abs(r) > 1e-3:
            v_star = r
            break
Xp = X.copy()
Xp[0], Xp[1], Xp[2] = a * np.exp(u), b * np.exp(v_star), c * np.exp(-u - v_star)
Xp = Xp[rng.permutation(N)]

print("sum 差      :", X.sum() - Xp.sum())
print("sum log 差  :", np.log(X).sum() - np.log(Xp).sum())
print("内部模式变了:", not np.allclose(np.sort(X), np.sort(Xp)))
LR_X = np.exp(logL_gamma(X, *th1) - logL_gamma(X, *th2))
LR_Xp = np.exp(logL_gamma(Xp, *th1) - logL_gamma(Xp, *th2))
print("LR(X)  =", LR_X)
print("LR(X') =", LR_Xp, "  差:", abs(LR_X - LR_Xp))
Xq = X.copy(); Xq[3] *= 1.01
LR_Xq = np.exp(logL_gamma(Xq, *th1) - logL_gamma(Xq, *th2))
print("对照：sum 动 1% => LR 相对变化:", abs(LR_Xq - LR_X) / LR_X)
ok_a = (abs(LR_X - LR_Xp) < 1e-9) and (abs(LR_Xq - LR_X) / LR_X > 1e-4)

print()
print("=" * 64)
print("V3b 自由形状更新族：LR 对间隔置换是否不变")
print("=" * 64)
f1 = lambda x: 0.6 * lognorm.pdf(x, 0.8, scale=1.0) + 0.4 * weibull_min.pdf(x, 1.7, scale=2.0)
f2 = lambda x: 0.3 * lognorm.pdf(x, 1.1, scale=1.5) + 0.7 * weibull_min.pdf(x, 2.5, scale=1.2)
X2 = rng.gamma(2.0, 1.0, 12)
LR_full = np.prod(f1(X2) / f2(X2))
perm = rng.permutation(12)
LR_perm = np.prod(f1(X2[perm]) / f2(X2[perm]))
print("LR 原序 :", LR_full)
print("LR 置换 :", LR_perm, "  差:", abs(LR_full - LR_perm))
X3 = X2.copy(); X3[5] *= 1.001
print("对照：动一个间隔 0.1% => LR 相对变化:",
      abs(np.prod(f1(X3) / f2(X3)) - LR_full) / LR_full)
ok_b = abs(LR_full - LR_perm) < 1e-12

print()
print("=" * 64)
print("裁决：V3a", "通过（极小性条款证伪确认）" if ok_a else "未通过",
      "| V3b", "通过（置换不变确认）" if ok_b else "未通过")
print("=" * 64)
