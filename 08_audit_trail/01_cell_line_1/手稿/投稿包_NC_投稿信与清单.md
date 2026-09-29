# NC 投稿包：投稿信 + 投审清单（v1.0，2026-08-11）

## 一、投稿信（Cover Letter，英文——投审系统粘贴用）

Dear Editors,

We are pleased to submit our manuscript, **"Why cells count: an identifiability audit of cellular information circuits,"** for consideration at *Nature Communications*.

Living cells sense and encode information with biochemical circuits, yet what can in principle be inferred about a circuit from measurable outputs has lacked systematic characterization. Our manuscript transfers an identifiability audit—originally developed for bioanalytical measurement, where it proved that digital counting is calibration-free—to cellular information circuits, and asks whether the same rules hold, and whether the same digital limit reappears.

Three contributions result. **First**, for two canonical circuits—the phosphorylation–dephosphorylation cycle and the p53 pulse generator—the audit is fully constructive: we derive the exact degeneracy groups leaving all observables invariant, quantify the Fisher-information barrier separating identifiable circuit structure from unidentifiable molecular scale, and establish one-sided reporting protocols for parameters at the identifiability edge. **Second**, for the p53 pulser we prove a topology-selection theorem: uniform-amplitude, dose-counting pulses require an excitable architecture containing positive feedback, excluding pure negative-feedback models as sufficient explanations of the p53 data. **Third**, across both cases the calibration-free regime is the same mathematical limit—when analog scale factors cancel, information survives only in counting statistics—and experimental-side discrimination on three published single-cell datasets (p53, NF-κB, ERK) finds the calibration-free observables to be counting-type without exception.

Two features of the manuscript may be of particular interest to your editors and referees. Its claims are advanced under an explicit pre-registered falsification protocol: an earlier version's analogue-region prediction for ERK was falsified by public data and corrected openly, a falsification–correction cycle we document in the manuscript itself as evidence of the framework's adjudicative force. And every numerical claim is deterministically reproducible from analysis scripts that will be deposited, with SHA-256 fingerprints, in a public repository upon submission.

We believe the work will interest readers across systems biology, biological information theory, and synthetic biology—for whom the digital limit doubles as a design principle: synthetic circuits that must transmit dose under unresolvable parameter uncertainty should be built to count.

This manuscript is not under consideration elsewhere. The authors declare no competing interests. Data- and code-availability statements are included in the manuscript; a companion preprint (the bioanalytical-measurement audit that forms the cited foundation) is being deposited concurrently and its DOI will be provided.

Sincerely,
[作者姓名与通讯地址——投稿时填入]

---

## 二、投审清单（状态截至 2026-08-11）

| # | 事项 | 状态 | 责任人 |
|---|---|---|---|
| 1 | 标题定案《Why cells count: an identifiability audit of cellular information circuits》 | ✅ 拍板（注册·二十七） | — |
| 2 | 正稿 §1–6 + 反向预言小节（双稿 v0.9.6） | ✅ 完成 | — |
| 3 | 摘要压缩 ~150 词（EN） | ✅ 完成 | — |
| 4 | 文献表 39 条卷期全部核验（含 Micali Sci Rep 9:16898 终核通过） | ✅ 完成（注册·三十二） | — |
| 5 | 数据/代码可用性声明（含代码14 预留编号如实声明） | ✅ 完成 | — |
| 6 | SI 八节 + 引用终同步 | ✅ 完成 | — |
| 7 | 投稿信 | ✅ 本文件 | — |
| 8 | **公开仓平台拍板** | ✅ Zenodo（注册·三十五，用户「无偏好」→按推荐；GitHub 活仓可选，以 Zenodo 快照 DOI 为准） | — |
| 9 | **代码打包上传 Zenodo + DOI 获取** | ⬜ | 你（我可协助整理打包清单） |
| 10 | **TCS 姊妹篇预印本投 bioRxiv → DOI 补入文献 27** | ⬜ | 你 |
| 11 | **本稿同步投 bioRxiv**（NC 政策允许，建议同日） | ⬜ | 你 |
| 12 | **EndNote/文献管理器导入定稿**（39 条已核验，直接导入） | ⬜ | 你 |
| 13 | **图件终导出**（5 主图 PNG→TIFF/300dpi，NC 图件规范；现 PNG 均在位） | ⬜ | 你（我可给规格清单） |
| 14 | **md→投稿格式转换**（docx/PDF；NC 初投接受单 PDF） | ⬜ | 我（说一声即做） |
| 15 | 作者名单/单位/贡献声明/ORCID | ⬜ | 你 |
| 16 | 建议/排除审稿人（可选；建议避开 Tkačik/ten Wolde 信息论直系以免立场审） | ⬜ | 你 |

## 三、投稿信使用说明

- 方括号占位（作者/地址）在投稿系统填写时替换；
- 「companion preprint ... deposited concurrently」一句以第 10 项完成为前提——若 TCS 预印本暂缓，删该句并把文献 27 表述改为「manuscript in preparation」；
- 审稿人建议名单若需要，我可按文献表反查候选人（引我们的人+我们引的中立方法学家）。
