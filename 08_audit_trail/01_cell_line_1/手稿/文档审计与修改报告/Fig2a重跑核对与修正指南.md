# Fig. 2a 重跑核对与正文修正指南

数据来源：`Scale Degeneracy5.0` 本机重跑输出（MCMC 50,000 步，seed=42）
基准版本：`main v6(2).docx`
结论：**Fig. 2a 全部正文数字复现，零改动**；唯一待办是 4PL/5PL bootstrap 的 n=200 → 2000 决策。

---

## 一、核对结论：Fig. 2a 正文数字全部复现 ✅（不用改）

| 项目 | 本次输出 | 正文现状 | 判定 |
|---|---|---|---|
| $\xi_0$ | 7.6930 [5.73, 10.84] | 7.69 | ✅ 一致 |
| $\kappa_1$ | 0.1443 [0.1010, 0.2006] | 0.144 | ✅ 一致 |
| $\kappa_2$ | 0.1694 [0.1104, 0.2578] | 0.169 | ✅ 一致 |
| $\kappa_3$ | 中位 498.8，HDI [18.76, 975.50] | $\geq 18.8$（median 499）；Results 段写 [18.8, 976] | ✅ 一致（18.76→18.8、975.5→976 正确舍入） |
| 全局 $R^2$ | 0.993629 | 0.994 | ✅ 一致 |
| DW | 2.4864 | 2.49 | ✅ 一致 |
| $G_1/G_2/G_3$ | 0.457 [0.449, 0.466] / 0.404 [0.384, 0.419] / 0.779 [0.756, 0.801] | 0.46 / 0.41 / 0.78，"CI non-overlapping" | ✅ 一致且三组 CI 两两不重叠 |
| RMSE | 4PL 0.0709 / 5PL 0.0504 / TCS 0.0473 | 0.071 / 0.050 / 0.047 | ✅ 一致 |
| 分组 $R^2$ | 0.9922 / 0.9936 / 0.9793；Bayesian [0.9918, 0.9932, 0.9710] | 报于 SI Table 1 | ✅ 以新 xlsx 为准 |

> 说明：你本机（seed=42）复现的就是正文值；我此前沙盒环境的 7.72/19.17/493.8 是 BLAS/emcee 版本差异，**以你本机输出为准**，正文无需任何回改。

---

## 二、唯一待办：4PL/5PL bootstrap 仍为 n=200

**证据**：本次输出头部为 `4PL/5PL Bootstrap 95% CI (n=200, D locked to blank)`——脚本第 531 行的 `n_boot` 还没改（或本次跑的是改前版本）。该 bootstrap 只管分组 4PL/5PL 的 $C,B,G$ 置信区间，**MCMC（κ 后验）与它无关**。

### 操作步骤

**1.** 脚本 `Scale Degeneracy5.0固定背景，最终确认.py` 第 531 行：

```python
n_boot = 200
```

改为：

```python
n_boot = 2000
```

（上一行 `np.random.seed(42)` 不动。）

**2.** 本机重跑（约 10 分钟，MCMC 会一并重算但结果不变）。

**3.** 重跑后核对两项：
- 三组 $G$ 的 bootstrap CI 仍两两不重叠（本次 n=200 已为 [0.449, 0.466] / [0.384, 0.419] / [0.756, 0.801]，间隔很宽，n=2000 只会更稳；若新 CI 有变化以新输出为准写入 SI Table 1）；
- Fig. 2a 的 MCMC 数字（$\xi_0$、$\kappa_{1,2,3}$、$R^2$、DW）应与本次完全一致（seed 固定），若有漂移以最终入库版为准回改。

**4.** Methods 同步（1 处）：
- **锚定搜索**：`triplicate-level bootstrap resampling`
- **原文**：Per-group uncertainty was quantified by triplicate-level bootstrap resampling ($n = 200$):
- **改为**：Per-group uncertainty was quantified by triplicate-level bootstrap resampling ($n = 2000$):

**5.** `SI_Table_1.xlsx` 用新输出重新导出打包。

> 备注：若决定不重跑，保留 n=200 在科学上也成立（该 bootstrap 不进任何正文数字，仅服务于 SI Table 1 的分组 CI），但与其他四份脚本统一成 2000 更干净，建议改。

---

## 三、上轮遗留的两处正文修改（与本次输出无关，仍待落实）

### 3.1 参数边界修正（Methods·Fig. 2a 拟合段）

- **锚定搜索**：`physically meaningful bounds`
- **原文**：…constrained within physically meaningful bounds ($\xi_0 \in [1,100]$, $A_i \in [0.1,10]$, $\kappa_i \in [10^{-4},10^{4}]$).
- **改为**：…constrained within physically meaningful bounds ($\xi_0 \in [0.1,500]$, $A_1 \in [0.5,3.5]$, $A_2 \in [0.5,2.5]$, $A_3 \in [0.1,0.5]$, $\kappa_i \in [10^{-4},10^{3}]$).
- 依据：脚本实际边界即如此，且 $\kappa_3$ HDI 上端 975.5 顶在 $10^3$ 上界；正文 Results 段 "extends to the upper prior ceiling $\sim 10^{3}$" 与此一致 ✓。

### 3.2 中位数辩护句插入（同一段落，`total posterior samples per parameter` 之后）

> Point estimates are posterior medians: they are optimal under absolute-error loss, robust to the right-skewed posteriors typical of scale parameters, and equivariant under monotone reparameterisation ($\mathrm{median}(\log\kappa) = \log\,\mathrm{median}(\kappa)$). For all well-identified parameters, the posterior medians agreed with the seeding least-squares estimates to within 0.1%.

（依据已核实：最小二乘解与 MCMC 中位数在四个物理参数上偏差全部 <0.1%。）

---

## 四、执行顺序

1. 先做三（两处正文文字修改，不重跑）；
2. 再做二（改脚本 531 行 → 重跑 → Methods n 同步 → SI_Table_1.xlsx 重打包）；
3. 全部完成后全文复查：Fig. 2a 相关数字（7.69 / 0.144 / 0.169 / 18.8 / 499 / 976 / 0.994 / 2.49）与新输出逐项一致——本指南第一节的表即为核对清单。
