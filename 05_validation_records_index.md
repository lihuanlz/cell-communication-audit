# Validation & Audit Records Index

Paper: *Scale-free statistics are the currency of cellular communication: an audit across four signalling systems* (manuscript v06, 2026-09-28; v05 archived).

This index lists every validation/audit artefact cited in the manuscript and SI. Paths are relative to the workspace root (`D:\Kimi_Agent_细胞仿真工具包扩展以及具身智能20260911\`). All artefacts remain in their original archive locations; this package copies only the manuscript, figures, data tables and code. Original filenames are quoted verbatim (Chinese names included) for traceability.

## 1. Pre-registration texts (frozen before unblinding)

| Pre-registration | Archive path (relative to workspace root) |
|---|---|
| P2 amplitude-stratification prediction (frozen 2026-08-14) | `05_主线纲领与设计\预注册_P2_幅度分层衰减预言_v01.md` |
| P3 free-statistic sufficiency, NF-κB spatial gradient | `05_主线纲领与设计\预注册_P3_免费统计量充分性_NFkB空间梯度_v01.md` |
| P4 ERK re-adjudication, GPCR-KTR trajectory level (v0.3, re-anchored to 6 doses) | `05_主线纲领与设计\预注册_P4_ERK重裁_GPCR-KTR轨迹级_v01.md` |
| P5 boundary-predictor first test, GPCR Ca²⁺ transients (v0.4 with erratum appendix) | `05_主线纲领与设计\预注册_P5_边界预测器首测_GPCR钙瞬态_v01.md` |
| P6 fold-change-attribution adjudication, NF-κB sequential grid | `05_主线纲领与设计\预注册_P6_FCD归属裁决_NFkB序贯网格_v01.md` |
| P7 FCD/Weber verdict, chemotaxis single-cell FRET (frozen 2026-08-15) | `05_主线纲领与设计\预注册_P7_FCD-Weber判决_趋化单细胞FRET_v01.md` |
| P8 prospective cohort registration | `03_细胞线3\结果\预注册_P8前瞻队列登记.md` |
| Mixed-ligand validation arm P1/P2 (frozen 2026-08-16; duplicate also at `04_细胞线4\_归档\结果\预注册_混合配体验证臂_P1P2.md`) | `03_细胞线3\结果\预注册_混合配体验证臂_P1P2.md` |

## 2. Verdict cards / closure cards (frozen at first execution, never edited)

| Card | Archive path |
|---|---|
| p53 Fisher audit verdict card | `03_细胞线3\结果\判词卡_P53-Fisher审计_2026-09-18.md` |
| TrackC2 discrimination-map verdict card | `03_细胞线3\结果\判词卡_TrackC2_判别图_2026-09-18.md` |
| p53 signalling-line closure card | `03_细胞线3\结果\封卷卡_p53信号线_2026-09-18.md` |
| p53 Fisher audit pathfinder card | `03_细胞线3\结果\探路卡_P53-Fisher审计_2026-09-18.md` |
| Cell-line-3 detailed report (chemotaxis receptor-layer rupture ledger) | `03_细胞线3\结果\报告_细胞线3_尺度简并框架与趋化受体层断裂_详细报告.md` |
| Cell-line-4 GPCR bias fracture-audit detailed report | `04_细胞线4\_归档\结果\报告_细胞线4_GPCR偏向性断裂审计_详细报告.md` |
| Cell-line-4 GPCR analysis note | `04_细胞线4\_归档\结果\分析_GPCR偏向性_断裂审计_v01.md` |

The P2–P7 adjudication verdicts themselves are the `裁决.json` files listed in §4 below; each is accompanied by a judgement log (`判定日志.txt`) in the same folder.

## 3. Errata (appended after freezing under the never-move rule)

| Erratum | Scope | Archive path |
|---|---|---|
| 44b | P2 bootstrap-CI recomputation; frozen verdict untouched; 63/63-row reproduction check bit-exact | `01_细胞线\结果\P2P3P4_盲裁决_本机原件\代码44b_CI单元表.csv`, `代码44b_CI补算_v012_全文.txt`, `代码44b_日志.txt` (same folder) |
| 47b | P5 feature-implementation defect; original code-47 verdict voided-by-erratum and double-recorded; 47b re-adjudication authoritative (195 valid cells) | `01_细胞线\结果\Keshelava2018_P5_47b勘误重裁\` (full folder); voided original at `01_细胞线\结果\Keshelava2018_P5\` |
| 49b | P7 unit-definition/background-assignment correction; both passes intermediate | `02_细胞线2\结果\P7_Moore2024_FCD-Weber\代码49b_裁决.json` (+ unit table at `03_细胞线3\结果\P7_Moore2024_FCD-Weber\代码49b_单元表.csv`) |

## 4. Adjudication verdict records (P2–P7, blind, pre-registered)

| Adjudication | Verdict JSON | Unit table | Folder |
|---|---|---|---|
| P2 (code 44) | `代码44_裁决.json` | `代码44_单元表.csv` | `01_细胞线\结果\P2P3P4_盲裁决_本机原件\` |
| P3 (code 45) | `代码45_裁决.json` | `代码45_单元表.csv` (+ descriptive/robustness arms) | same |
| P4 (code 46) | `代码46_裁决.json` | `代码46_单元表.csv` (+ arms) | same |
| P5 (code 47b, authoritative) | `代码47b_裁决.json` | `代码47b_单元表.csv` | `01_细胞线\结果\Keshelava2018_P5_47b勘误重裁\` |
| P6 (code 48) | `代码48_裁决.json` | `代码48_单元表.csv` | `02_细胞线2\结果\P6_Wang2022_FCD归属\` |
| P7 (codes 49/49b) | `代码49_裁决.json`, `代码49b_裁决.json` | `代码49_单元表.csv`, `代码49b_单元表.csv` | `02_细胞线2\结果\P7_Moore2024_FCD-Weber\` |

## 5. Cloud re-execution records (independent re-runs, 2026-08-15)

| Scope | Archive path | Contents |
|---|---|---|
| P2 + P3 cloud re-execution | `01_细胞线\结果\P2P3_云端复核_2026-08-15\` | codes 44/45 with unit tables, descriptive/robustness arms, judgement logs, verdict JSONs; P2 promoter Δ values bit-identical to local |
| P4 cloud re-execution | `01_细胞线\结果\P4_云端复核_2026-08-15\` | code 46 unit table, judgement log, verdict JSON (AUC_L = 0.5668941796791128 identical to local card) |

No cloud–local discrepancy exists in the ledger (SI S7.7, S10.5).

## 6. Prior-art search (formal round, 2026-09-23)

`05_主线纲领与设计\论文_细胞通讯审计_2026-09-23\先案检索_正式轮_2026-09-23\`

- Search log: `检索记录_正式轮_v01.md`
- Per-query result CSVs (top-20 records each, individually checked): `Q01_scalefree_cellular_signaling.csv`, `Q02_foldchange_scaleinvariance.csv`, `Q03_identifiability_sloppy.csv`, `Q04_Fisher_p53_pulse.csv`, `Q05_blind_adjudication.csv`, `Q06_GPCR_audit.csv`, `Q07_chemotaxis_weber.csv`, `Q08_counting_digital_encoding.csv`, `Q09_DR_GED_invalidation.csv`, `Q10_cascade_information_loss.csv`, `Q11_pulse_counting_biology.csv`
- First-round map: `A3论文素材\先案检索_已有工作vs本框架.md`
- Verdict (SI S11): no collision on any of the four novelty pillars; the WoS/Scopus residual pass is flagged and the eleven query strings are archived for verbatim re-run.

## 7. Script-number reconciliation table (closed 2026-09-23)

`05_主线纲领与设计\论文_细胞通讯审计_2026-09-23\代码编号对账表_2026-09-23.md`

Reconciles the NC manuscript logical script numbers (4–17) with the physical files across the four cell-line folders; registers the code-52 naming collision (TrackC2 vs chemotaxis protocol-artefact script) and the S5.2/S27 numerical-pin scripts (`05_主线纲领与设计\数值核查_S27_Gamma族_V3.py`, `数值核查_S27_噪声集中界.py`). The package's `04_code\CODE_MAP.csv` is the machine-readable companion to this table.

## 8. Double-recorded failure inventory (SI S10.3 summary)

The programme's failures are archived with the same discipline as its hits:

- **Falsified adjudication clauses**: P2-1, P2-2, P3-3, P4-1, P4-3, P5-1 (SI S7; verdict JSONs in §4).
- **Power-limited intermediates reported, not resolved post hoc**: P3-1, P4-2, P5-2, P6, P7.
- **Two rejected NF-κB repair patches**: codes 12 and 15 (SI S6.2), scripts in `01_细胞线\代码与图\`.
- **ERK static-extrapolation falsification**: v0.8 prediction (amplitude α ≳ 0.5) overturned by source-figure discrimination on Albeck 2013 (α ≈ 0); preserved in version history (SI S6.3).
- **Three double-recorded chemotaxis pipeline errors**: codes 55, 57, 59 (SI S8.3), scripts in `03_细胞线3\代码与图\`.
- **P5 pipeline error** (code 47, voided by erratum 47b) and **P4 v0.1.2 guardrail empty-run history** (abolished before running), both double-recorded (SI S7.8).
- **PTH1R admission-ticket rejection**, retained as boundary-artefact teaching case (SI S9.1).
- **Framework bug lineage**: six internal bugs and sixteen external-review items across cascade-framework versions v2.0–v4.2, plus Bug 7 and the chemotaxis analyst-normalisation Bug 8 (SI S1.8, S1.9); framework documents in `02_细胞线2\`.
- **Code-70 extraction-convention error** (V1 s-field records total MeAsp), caught by file-identity check and repaired, double-recorded (SI S8.6).

## 9. Evidence-strength register S13 (2026-09-23)

S13 evidence-strength register (2026-09-23), including the code 82 mechanism-level stress test M1–M6 on the reconstructed 23-parameter Mönke p53 model. Location: SI section S13 (S13.1 headline-claim register E1–E32 with tier assignments, S13.3 mechanism-level stress test, S13.4 calibrated-wording ledger). Verdict card: `03_细胞线3\结果\判词卡_代码82_p53机制模型_2026-09-23.md` (with companion report `报告_代码82_p53机制模型_2026-09-23.md` and `代码82_机制审计_结果.json` in the same folder). Working-copy Chinese archive of the section: `05_主线纲领与设计\论文_细胞通讯审计_2026-09-23\S13_证据强度登记_2026-09-23.md`.

- 2026-09-25: codes 89-93 verdict cards and audit card (`03_细胞线3/结果/判词卡_代码89..93_*.md`, `审计卡_代码89-91_*.md`) underpin SI S14; result JSONs copied to `03_data_csv/p53_molecular_layer/`.

## 肿瘤预言线（v06 新增，P11–P16）

卡片、脚本与结果已归位投稿包：`07_audit_cards_tumour_chain\`（英文文件名，内容中文为权威工作记录）、`04_code\tumour_chain_P11_P16\`、`03_data_csv\tumour_chain_P11_P16\`；工作副本（中文文件名）保留于 `05_主线纲领与设计\肿瘤预言线\`：

| 编号 | 内容 | 判决 | 预注册卡 / 判词卡 |
|---|---|---|---|
| P11 | p53 计数律跨细胞系 | 一毙一不决一弱反对 | 预注册卡_P11 / 判词卡_P11_p53计数律跨细胞系_2026-09-28 |
| P12 | 增益分布与耐药零方向 | P12a 不成立；P12b 成立 + EGFR 弱成立 | 预注册卡_P12（含 A1–A5）/ 判词卡_P12_增益分布与耐药零方向_2026-09-28 |
| P13 | 信道容量互信息 | 部分成立、分级、复核稳健 | 预注册卡_P13（含 A1）/ 判词卡_P13_信道容量_2026-09-28 |
| P14 | p53 下游误读 | 成立（AUC 校正登记） | 预注册卡_P14（含 A1）/ 判词卡_P14_p53下游误读_2026-09-28 |
| P15 | TCGA 组织解耦 | 不成立（外推边界） | 预注册卡_P15 / 判词卡_P15_TCGA信道解耦_2026-09-28 |
| P16 | G13D 信道预言 | 冻结待验 | 预注册卡_P16_G13D信道预言_2026-09-28 |

结果 JSON：p13_mi_results.json、p13_robustness.json、p14_results.json、
p15_tcga_results.json、p12b_egfr_dms_results.json、p12b_egfr_dms_A5_results.json、
p13_fig6_decay.json。逐样本/逐基因表：data/p15_*_persample.csv、p14_pergene_log2R.csv。
