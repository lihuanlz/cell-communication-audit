# 先案检索：已有工作 vs 本框架（第一轮）

- 检索日期：2026-08-18
- 工具：Google Scholar（scholar 插件），10 条查询，原始命中 CSV 存于 `先案检索原始记录/`
- 局限声明（诚实双录）：仅 Google Scholar 单源；英文偏向；未交叉 Web of Science / Scopus；未查专利与工业软件文档。本表结论=**第一轮、非穷尽**。二轮欠账见文末。

---

## 一、三条线检索结果

### 线 A：数据协调 / 显著误差检测（data reconciliation & gross error detection, DR/GED）

成熟工程领域，1980 年代起于化工：

| 文献 | 年份 | 要点 |
|---|---|---|
| Crowe, AIChE J | 1988 | 线性数据协调中递归识别显著误差（经典） |
| Zhang & Chen, Comput. Chem. Eng. | 2014 | 动态系统同时协调+显著误差检测 |
| Yuan et al., AIChE J | 2015 | 贝叶斯同时检测与协调 |
| Fuente et al., IFAC | 2015 | 显著误差管理综述 |

**已进入生物/医药的痕迹（重要，必须引用划界）：**
- Medeiros et al., Crit. Rev. Anal. Chem. 2023：数据协调进分析化学（综述）
- Chaffee, AIChE 年会 2008：连续制药制造中的 DR/GED
- **代谢工程**：Leighty & Antoniewicz, Metab. Eng. 2011（DMFA，明确使用"data reconciliation and detection of gross errors"）；Jordà et al., Metabolites 2014（"data consistency verified using standard data reconciliation"）；Kim & Lee 2006（同位素平衡数据协调）

**与本框架的关系：哲学孪生。** 同为"不拟合全模型、用冗余+约束检验数据内在一致性"。
差异：约束来源不同——他们靠守恒律（物料/同位素平衡），本框架靠**模型结构本身（跨协议共享参数的恒等式）**；审计对象不同——他们协调过程测量值，本框架审计**已发表的参数集与原始数据**。
**威胁度：中。** "没人做过"的声称必须限定在"离子通道/药理单细胞数据"内。

### 线 B：model invalidation（模型作废）

控制论经典 + 系统生物学移植：

| 文献 | 年份 | 要点 |
|---|---|---|
| Smith & Doyle, IEEE TAC | 1992 | "只能作废、不能证实"的名言出处（505 引） |
| Rosa et al., CDC | 2010 | set-valued observers，LTV 故障检测 |
| Harirchi & Ozay, IFAC | 2015 | switched affine 系统的模型作废用于故障检测 |
| Anderson & Papachristodoulou, BMC Bioinf. | 2009 | 生物模型的 validation/invalidation（94 引） |
| Rumschinski et al., BMC Syst. Biol. | 2010 | **集合参数估计+模型作废，生化反应网络**（89 引） |
| Streif et al., Bioinformatics | 2012 | ADMIT 工具箱：保证性模型作废 |
| Bates & Cosentino, IET Syst. Biol. | 2011 | 鲁棒性分析做模型 validation/invalidation |

**与本框架的关系：最近的方法学邻居。** 本框架的"兼容带成员资格 / FIM 椭球判据 / 跨度演示"与 set-membership invalidation（Rumschinski 一脉）在精神上同构：都是"参数集是否落在与数据相容的集合内"。
差异：①用途——他们在**建模过程中**筛选/作废候选模型，本框架作为**外部 referee 审计他人已发表的数据与参数**；②制度层——本框架有冻结判线、合成对照电池（校准误报/漏报）、盲法预注册、公开自审计记录，邻居均无；③对象——他们判模型，本框架判数据+参数集。
**威胁度：高。** 引言必须显式引用 Rumschinski/ADMIT 一脉并划界，否则审稿人会替我们做这件事。

### 线 C：离子通道 / 药理的可辨识性、复现性审计、数据 QC

**可辨识性（建模工具，非审计）：**
- Ball & Sansom, Proc. R. Soc. 1989；Ball & Rice, Math. Biosci. 1992：单通道门控模型辨识（160/168 引）
- Gutenkunst et al., PLoS Comput. Biol. 2007：sloppy 普适性（1679 引）——本框架"简并"概念的直接前辈，必引
- Janzén et al., Front. Physiol. 2016：药效模型参数可辨识性（含离子通道阻滞剂案例）
- Mangold et al., PLoS Comput. Biol. 2021：通道动力学**结构**系统辨识（枚举 Markov 结构）
- Epstein et al., Biophys. J. 2016：通道贝叶斯推断+漏检事件校正
- （公认经典 Raue et al. 2009 profile-likelihood 本轮未单独检索，二轮补确认）

**hERG 拟合与变异性（拟合，不审计）：**
- Lei, Clerx et al., Biophys. J. 2019（I、II 两篇）：15 秒协议快速表征 hERG 动力学+参数变异性——**本框架 hERG 数据的来源方**，拟合变异性研究，不审计他人数据
- Li et al., J. Pharmacokinet. Pharmacodyn. 2016；Ridder et al., Toxicol. Appl. Pharmacol. 2020

**复现性审计（审计模型，非数据一致性）：**
- **Tiwari et al., Mol. Syst. Biol. 2021（181 引）：审计已发表系统生物学模型的可复现性——本轮发现的"审计他人已发表工作"的最近邻居。** 但对象=模型文件（SBML），手段=重实现重跑，**不是**零拟合的数据内在一致性检验；且依赖作者公开的模型代码。
- Papin et al., PLoS Comput. Biol. 2020；Karr et al., Front. 2022：复现性社区努力

**电生理数据 QC（原始记录级，与参数一致性无关）：**
- Harris et al., Nat. Neurosci. 2016；Nowotny et al., PLoS ONE 2013；Seibertz et al., Commun. Biol. 2022：封接阻抗/漏电/噪声/漂移级质控

**线 C 结论：对"零拟合、第三方、审计已发表药理/电生理**数据**内在一致性"这一精确生态位，本轮检索零命中。** gap 在第一轮范围内确认。

---

## 二、总对比表

| 工作线 | 代表文献 | 零拟合 | 约束/判据来源 | 用途 | 校准电池/盲测 | 自审计记录 | 与本框架关系 |
|---|---|---|---|---|---|---|---|
| DR/GED（化工） | Crowe 1988; Yuan 2015 | ✅ | 守恒律 | 协调自家过程测量 | 有仿真检验传统 | 无 | 哲学孪生（零件） |
| DR→代谢工程 | Leighty 2011; Jordà 2014 | ✅ | 同位素/物料平衡 | 自家通量估计 | χ² 一致性 | 无 | 跨域先例，必须引用 |
| model invalidation | Smith & Doyle 1992; Rumschinski 2010 | 部分（集合法免点估计） | 模型+误差界 | 建模过程中筛模型 | 无 | 无 | **最近邻居，必须显式划界** |
| sloppy/可辨识性 | Gutenkunst 2007; Janzén 2016 | — | 灵敏度谱 | 建模者自查 | 无 | 无 | 零件（简并概念前辈） |
| 通道结构辨识 | Mangold 2021 | ❌（贝叶斯拟合） | Markov 结构枚举 | 选模型结构 | 合成基准 | 无 | 相邻（判模型不判数据） |
| 复现性审计 | Tiwari 2021 | — | 模型文件重跑 | 审计已发表**模型** | 人工评级 | 无 | 相邻（审计对象与手段均不同） |
| 电生理 QC | Harris 2016; Seibertz 2022 | — | 记录质量指标 | 原始数据准入 | 有 | 无 | 不同层（记录级 vs 参数级） |
| **本框架** | — | ✅（审计步零自由参数） | 跨协议共享参数恒等式+绝对锚 | **第三方审计已发表数据+参数** | **合成电池+盲法+预注册** | **有（F4 全档）** | — |

---

## 三、创新定位（修订版，替换"没人系统造过"）

**旧声称（作废）：** "在药理学数据审计领域，之前没人这么系统地造过。"——太宽，DR/GED 与 model invalidation 两脉会被审稿人拍。

**新声称（可防守）：**
> 零拟合的数据一致性审计思想 1980 年代起于化工数据协调，2010 年代经 13C 代谢通量分析进入生物学；model invalidation 在控制论与系统生物学中成熟用于建模过程。本框架首次将二者与可辨识性分析融合并移植到**离子通道/药理学已发表数据的第三方审计**，其约束来自模型结构本身（跨协议共享参数恒等式）而非守恒律，并配套五个制度化构件：入场券、五步程序、校准电池（误报/漏报实测）、盲法预注册、公开自审计执行记录。

**三个独占差异点（逐条可防守）：**
1. **约束来源**：无守恒律领域的数据协调——约束由模型结构（共享参数）代数生成。
2. **用途转向**：从"建模者自查"到"外部 referee 审计他人已发表数据"，配冻结判线与注册制度。
3. **公开自审计记录**：F4 触发→尸检→v1.1 补丁→首考通过全档留存（§31–34）。检索范围内，DR/GED 与 model invalidation 两脉均无此制度证据。

---

## 四、A3 引言可用段落（中文草稿）

> 检验已发表定量数据的内在一致性，通常依赖重新拟合或人工判读，两者都引入审计者自身的主观自由度。化工领域的数据协调与显著误差检测（data reconciliation & gross error detection）自上世纪 80 年代起提供了一条免拟合路径：用冗余测量与守恒约束直接对账；该思想近年经 13C 代谢通量分析进入生物学，控制论与系统生物学中的模型作废（model invalidation）则以集合成员资格判定模型与数据的相容性。然而，这两条路径依赖守恒律或服务于建模者自身的模型筛选；对离子通道与药理学这类无守恒律、以模型拟合参数为主要数据载体的领域，尚不存在系统化的第三方数据审计工具。
>
> 本文提出的审计框架将上述思想移植到药理-电生理数据：审计约束不由守恒律提供，而由模型结构本身生成——同一组有效参数被多条电压协议共享，任何可重复计算的观测量必须满足相应恒等式；配合绝对定标锚与公开原始数据（三张入场券），审计在零自由参数下执行。框架附带五项制度化设计：冗余面/可辨识面的显式二分、五步操作程序、误报与漏报率实测校准的合成对照电池、盲法预注册、以及公开的自审计执行记录。我们在 hERG 通道公开数据上完成全流程演示，并预注册了一个预期静默的对照体系。

---

## 五、二轮检索欠账（排雷未完）

1. Web of Science / Scopus 交叉验证（Google Scholar 覆盖偏向）。
2. Rumschinski 式集合方法是否已被应用于离子通道/电生理数据（本轮未见，需定向复查）。
3. 工业自动膜片钳平台（Sophion、Nanion 等）内置 QC 是否含模型一致性检查（专利/白皮书层面）。
4. Raue 2009 profile-likelihood 等公认文献的正式确认引用。
5. 中文文献与学位论文库。


---

## 六、二轮检索（2026-08-20）

- 工具：scholar 插件（Google Scholar）4 组定向查询 + 公开网络源交叉（PMC/出版社页）；原始命中 CSV 存于 `先案检索原始记录/二轮_*.csv`（4 份）。
- 局限声明（双录）：WoS/Scopus 仍无法直连，欠账 1 只以公开网络源部分覆盖，**投稿前仍需机构数据库正式跑一轮——继续挂账**。

### 欠账 2：集合方法是否已进离子通道/电生理 —— 两轮定向，零命中

- "model invalidation + ion channel" 组命中大量基因敲除语义噪声，无集合法；定向复查 "set-based parameter estimation guaranteed ion channel" 组：最近邻为 Rauh & Kretzberg（神经元放电模型的集合仿真，前向保证包围，非审计）、Dang et al. 2019（生物建模集合分析综述）、Podlaski et al. 2017（eLife，2378 模型集合归类）。**Rumschinski 式集合作废应用于离子通道已发表数据：检索范围内仍为零。gap 维持。**
- **附带新发现（必引）**：
  - Pathmanathan & Gray, Front. Physiol. 2018（77 引）：心脏电生理多尺度模型的验证与可信度框架——本框架审计制度的**制度侧最近邻居**，讨论节必引。
  - **Sigg & Carnevale, Biophys. J. 2025："Markov models and long-term memory in ion channels: A contradiction in terms?"——与我们 §51 史依赖签名（固定速率模型类算术上不可调和）直接对话的 2025 年文献，通缉令方向的 citation 锚。**
  - Clerx, Beattie, Gavaghan & Mirams, Biophys. J. 2019（67 引）："Four ways to fit an ion channel model"——数据血脉方方法论自述。

### 欠账 3：工业自动膜片钳 QC —— 记录级，无模型一致性检查

- Danker & Möller 2014、Obergrussberger 2016/2018、Kutchinsky 2003（QPatch）、Dubin 2005（PatchXpress）、Guo & Guthrie 2005、Polonchuk 2012（PatchLiner 温控）、Jones 2009、Mathes 2006、Houtmann 2017、Bell & Fermini 2021：QC 全部是**封接阻抗/串联阻抗/电流稳定性**等记录级指标。
- Brinkwirth et al. 2020：基准与校准标准（best practices benchmarking）；**Baron et al. 2025：ICH S7B Q&A 2.1 最佳实践下手工膜片钳 hERG 数据变异性（Qube 384 对照）——监管级数据质量基准的现行文献，必引**。
- **无一例含模型级/参数级一致性检查。生态位空置确认。**

### 欠账 4：Raue 2009 正式确认 ✅

Raue, Kreutz, Maiwald, Bachmann et al., **Bioinformatics 25(15):1923–1929, 2009**（检索时点 1,882 引）。谱系一并确认：Raue et al. 2014（方法对比，266 引）、Kreutz et al. 2013（FEBS J 综述）、Wieland et al. 2021（Curr. Opin. Syst. Biol.，471 引）、Simpson & Maclaren 2023（profile-wise workflow，PLoS Comput. Biol.）。

### 欠账 5：中文文献 —— 零命中

参数辨识/似然剖面可辨识性在电力系统等工程领域有中文文献（陈润泽等 2015 等）；离子通道中文文献仅见单通道信号重构（韩晓东 2001；乔晓艳 2007）与门控动力学（方积乾等 1999）。**第三方数据一致性审计：中文库零命中。gap 维持。**

### 二轮总判

- 五条欠账销四条（2/3/4/5），欠账 1（WoS/Scopus 直连）继续挂账。
- 一轮"线 C 零命中"结论经二轮定向加固；新增两篇必引（Pathmanathan & Gray 2018 制度邻居；Sigg & Carnevale 2025 史依赖对话）与一篇监管数据基准（Baron 2025）。
- 创新定位（一轮第三节修订版）无需变动；A3 讨论节增加史依赖↔长时记忆文献的对话接口。
