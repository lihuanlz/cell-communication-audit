# Submission Package — Cell Communication Audit

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23008824.svg)](https://doi.org/10.5281/zenodo.23008824)

**Paper:** *Scale-free statistics are the currency of cellular communication*
**Version:** manuscript v06 / SI v06, 2026-09-28 (v05 archived in ../_archive_旧版本/manuscript_v05/; see 更新清单_v05_to_v06_2026-09-28.md)
**Package assembled:** 2026-09-23, from the local project archive. All files are **copies**; archive originals were never moved, modified, or deleted.

## Package contents

```
submission_package\
  01_manuscript\                  main text + SI, Markdown and docx, v06 only
  02_figures\                     Fig1–6 and EDFig1–16, each as SVG + PNG twins (flat layout)
  03_data_csv\                    key quantitative result tables as clean CSV (UTF-8, header row only)
    tumour_chain_P11_P16\         cancer-chain audit result JSONs + input/digitised CSVs (P11–P15)
  04_code\                        all numbered analysis scripts (ASCII filenames) + CODE_MAP.csv
    tumour_chain_P11_P16\         cancer-chain audit scripts P11–P15 (English filenames, self-contained)
  05_validation_records_index.md  index of pre-registrations, verdict cards, errata, cloud re-runs, prior-art search
  07_audit_cards_tumour_chain\    verdict + pre-registration cards for the cancer chain (P11–P16)
  README.md                       this file
```

### 01_manuscript

| File | Content |
|---|---|
| `Nature_main_v06_2026-09-28.md` / `.docx` | Main text, v06 |
| `Nature_SI_v06_2026-09-28.md` / `.docx` | Supplementary Information, v06 (S0–S19) |

### 02_figures

- Main text: `Fig1_ledger`, `Fig2_static_sensor`, `Fig3_dynamic_encoder`, `Fig4_adjudications`, `Fig5_receptor_rupture`, `Fig6_gpcr_metaaudit` (`.svg` + `.png` each).
- Extended Data: `EDFig1_framework`, `EDFig2_PdPC_full_panels`, `EDFig3_p53_discrimination_grid`, `EDFig4_p53_channel_census_51r3`, `EDFig5_dispersion_noise_colour_51r4_52`, `EDFig6_digital_limit`, `EDFig7_NFkB_two_regime_audit`, `EDFig8_MorrisLecar_duality`, `EDFig9_blind_adjudications_P2_P7`, `EDFig10_chemotaxis_rupture`, `EDFig11_GPCR_fracture_maps`, `EDFig12_atlas_violations_repair`, `EDFig13_workflow` (`.svg` + `.png` each).

### 03_data_csv

| CSV | Source (archive path, quoted verbatim) | SI section |
|---|---|---|
| `P2_msn2_stratification_units.csv` | `01_细胞线\结果\P2P3P4_盲裁决_本机原件\代码44_单元表.csv` | S7.1 |
| `P2_verdict_summary.csv` | `...\P2P3P4_盲裁决_本机原件\代码44_裁决.json` + `代码44b_CI单元表.csv` (44b CI erratum) | S7.1 |
| `P3_verdict_summary.csv` | `...\P2P3P4_盲裁决_本机原件\代码45_裁决.json` | S7.2 |
| `P4_verdict_summary.csv` | `...\P2P3P4_盲裁决_本机原件\代码46_裁决.json` + `代码46_单元表.csv` | S7.3 |
| `P5_verdict_summary.csv` | `01_细胞线\结果\Keshelava2018_P5_47b勘误重裁\代码47b_裁决.json` | S7.4 |
| `P5_population_arm.csv` | `...\Keshelava2018_P5_47b勘误重裁\代码47b_描述臂.csv` | S7.4 |
| `P6_verdict_summary.csv` | `02_细胞线2\结果\P6_Wang2022_FCD归属\代码48_裁决.json` | S7.5 |
| `P7_verdict_summary.csv` | `02_细胞线2\结果\P7_Moore2024_FCD-Weber\代码49_裁决.json` + `代码49b_裁决.json` (two rows: first pass, erratum) | S7.6 |
| `p53_channel_census.csv` | `03_细胞线3\结果\代码51r3_p53_Fisher审计_结果.json` (+ `代码51r4B_Fisher_M4000.json` robustness block) | S4.1–S4.3 |
| `p53_noise_calibration.csv` | `03_细胞线3\结果\代码51r4A_噪声标定与稳健性.json` (+ delay block of `代码51r3_结果.json`) | S4.3/S4.6 |
| `chemotaxis_weber_line.csv` | `05_主线纲领与设计\论文_细胞通讯审计_2026-09-23\figures_svg\ED\scripts\chemo54_rerun.json` (code-54 rerun cache; per-cell K1/2 medians identical to code 62 pipeline, as cross-checked in the code-62 header) | S8.4 |
| `chemotaxis_prospective_arm.csv` | transcribed from `03_细胞线3\代码与图\代码70_混合配体验证臂_P1P2判决_全文.txt` (code-70 archived full-text log) | S8.6 |
| `mor_code80_residuals.csv` | `04_细胞线4\_归档\结果\code80_arm_coordinate_residuals.csv` | S9.8 |
| `atlas_per_paper_medians.csv` | `04_细胞线4\_归档\结果\atlas_audit\atlas_audit_papers.csv` | S9.5 |
| `atlas_pair_violations.csv` | exact per-pair recomputation from `04_细胞线4\_归档\结果\atlas_audit\atlas_audit_entries.csv` (ligands deduplicated per group, all within-group ligand pairs, deviation = \|Δ(log τ) − Δ(logRA)\|; 17,987 rows) | S9.5 |

Verdict-summary CSVs share the column schema `clause, metric, value, CI_lo, CI_hi, frozen_line, verdict`.

**Notes on files that do not exist as flat tables:**

- `atlas_pair_violations.csv` in this package **is** the full 17,987-row ligand-pair violation table, regenerated exactly from `atlas_audit_entries.csv` (code 71's own entry-level output; the archived script `代码71_Atlas全库恒等式流水线.py` did not write the pair table to disk). It reproduces the paper's headline rates exactly: **21.5% of 17,987 pairs > 0.3 dex, 4.2% > 1 dex** (verified 2026-09-23 during package assembly; the same computation underlies Fig. 6d and EDFig. 12). One honest caveat: the quantile fields inside the archived `atlas_audit_summary.json` were computed with a cruder group-median replication shortcut (which would give 13.5% / 1.2% if used for the headline rates) — that shortcut is **not** the basis of any number in the paper; the exact per-pair table shipped here supersedes it. Source artefacts: `04_细胞线4\_归档\结果\atlas_audit\` (`atlas_audit_entries.csv`, `atlas_audit_groups.csv`, `atlas_audit_summary.json`).
- The D2R/AT1R/μOR per-ligand fracture tables of codes 63/65/66 were **never written as CSVs**; their results are embedded in the archived full-text logs `代码63_GPCR功能亲和力_断裂检验_全文.txt`, `代码65_AT1R跨臂亲和力断裂检验_全文.txt`, `代码66_MOR_Gillis内部一致性审计_全文.txt` in `04_细胞线4\_归档\代码与图\`, and summarised in `04_细胞线4\_归档\结果\报告_细胞线4_GPCR偏向性断裂审计_详细报告.md`.

### 04_code

All numbered analysis scripts renamed to ASCII English under the scheme `code<NN>_<short_english_slug>.py`. `CODE_MAP.csv` (columns `package_filename, original_path, NC_logical_number_if_any, si_section, one_line_description`) maps every package file back to its original Chinese filename and absolute archive location, and gives the NC-lineage logical number where the reconciliation table assigns one (codes 16 and 17 live physically under different numbers). Known collisions and notes:

- `code52_TrackC2_discrimination_map.py` and `code52chemo_protocol_artifact_FRET_coordinates.py` are two different scripts that share the number 52 in two line folders (registered in the reconciliation table).
- `code49_P7_blind_adjudication_FCD_Weber.py` was copied from `02_细胞线2\代码与图\` (the `03_细胞线3` folder holds only the erratum 49b); this is recorded in CODE_MAP.
- `code51r1/r2` are the archived earlier revisions of the p53 Fisher audit; `code51r3` is definitive (NC logical number 16). `code51r4_p53_fisher_wrapup.py` is included because `p53_noise_calibration.csv` derives from its outputs.
- `ed_scripts/` contains the Extended Data panel scripts and `chemo54_rerun.json` (data cache consumed by `edfig10.py`).

## How to reproduce

1. Start from `04_code\CODE_MAP.csv` and the frozen reconciliation table `05_主线纲领与设计\论文_细胞通讯审计_2026-09-23\代码编号对账表_2026-09-23.md` (quoted path; SI S0.3). The reconciliation table is the authoritative mapping between NC logical numbers, physical scripts, and their paired data/figures.
2. Every adjudication (P2–P7) runs from the frozen script with fixed seed (SEED = 20260814 for P2, 20260815 for P3–P7) and bootstrap ×2000 against the public datasets below. Verdict JSONs and unit tables are archived in the folders listed in `05_validation_records_index.md` §4; cloud re-executions of P2/P3/P4 are bit-identical (§5 there).
3. Figures: `python make_fig<N>.py` for main figures; `python draw_ed_figures.py` (which drives `ed_scripts/`) for Extended Data.
4. Archive paths referenced throughout use the original Chinese folder names; they are quoted verbatim in CODE_MAP and in `05_validation_records_index.md`.

## Data availability (public datasets)

| Dataset | Used in | Access |
|---|---|---|
| Hansen & Zechner 2021 (Msn2, yeast) | P2 | Zenodo doi:10.5281/zenodo.2755026 |
| Son et al. 2022 (NF-κB spatial gradient) | P3 | Zenodo doi:10.5281/zenodo.6858118 |
| Chavez-Abiega et al. 2022 (ERK-KTR, GPCR dose series) | P4 | publisher deposit (see SI S7.3) |
| Keshelava et al. 2018 (GPCR→Ca²⁺, Nat. Commun. 9:876) | P5 | publisher deposit (see SI S7.4) |
| Wang et al. 2022 (NF-κB sequential stimulation) | P6 | publisher deposit (see SI S7.5) |
| Moore et al. 2024 (chemotaxis single-cell FRET) | P7, S8 | Dryad doi:10.5061/dryad.nvx0k6dzz (CC0) |
| Biased Signaling Atlas | S9.5 | https://biasedsignalingatlas.org/ |

All other data are generated by the scripts in `04_code\`.


---

## Addendum 2026-09-25 (molecular-layer closure, SI S14)

The p53 line was closed at the molecular layer after the 2026-09-23 freeze. Added to this package:

- `04_code/code82` to `code93`: mechanism-level audit stack and molecular-layer decomposition scripts (see CODE_MAP.csv; deterministic seeds 20260923/20260925).
- `03_data_csv/p53_molecular_layer/`: result JSONs for codes 89-93.
- Figures: EDFig. 14 (molecular decomposition) and EDFig. 15 (closed-form cross-validation), PNG + SVG.
- Manuscript: main text carries a molecular-closure sentence in the counting-channel paragraph and a Sixth limitations entry; SI carries the new S14 section and register entry M8. The docx files in `01_manuscript/` were rebuilt on 2026-09-25 from these sources.

Verdict cards for codes 89-93 and the literature-cross-check audit card are archived in the workspace at `03_细胞线3/结果/` (Chinese filenames; English aliases to be assigned at Zenodo snapshot).

---

## Addendum 2026-09-28 (cancer-chain prediction arm, SI S19, main-text "Hijacked, not interrupted")

The cancer-chain audits P11–P15 (plus the frozen forward prediction P16) were executed after the v05 freeze and are integrated in manuscript v06. Added to this package:

- `04_code/tumour_chain_P11_P16/`: twelve self-contained Python scripts (English filenames). Naming: `p11_*` counting-law cross-cell-line audit; `p12a_*` ERK gain-distribution digitisation (violin + Gerosa A375); `p12b_*` resistance zero-direction audits (ABL, clinical expansion, EGFR DMS, EGFR DMS arm A5); `p13_*` mutual-information channel audit + robustness recheck; `p14_*` downstream decoding audit; `p15_*` TCGA decoupling audit.
- `03_data_csv/tumour_chain_P11_P16/`: result JSONs (`p11`–`p15`, including the archived `p12b_egfr_dms_A5_results.json`) and the input/digitised CSVs they consume.
- `07_audit_cards_tumour_chain/`: five verdict cards (`verdict_P11`–`verdict_P15`) and six pre-registration cards (`prereg_P11`–`prereg_P16`), English filenames; card bodies remain in Chinese as the authoritative working record.
- Figure: `EDFig16_cancer_chain` (PNG + SVG) in `02_figures/`; generator script `make_edfig16.py` in `04_code/`.
- Headline numbers: P11-1 slope correlation rho = +0.857 (p = 0.0068); P12b ABL gate 2.94x (p = 0.0326), EGFR DMS 2.2x (p = 0.0022); P13 WT 0.426 bit with oncogenic-mutant deltas down to -0.413 bit (KSG robustness check concordant); P14 AUC ratio 1.65; P15 bulk-tumour decoupling not supported in any of four cancer types (registered negative).
- Working copies (identical content, Chinese filenames) remain at `05_主线纲领与设计\肿瘤预言线\`; errata: P15 verdict card COADREAD/LUAD figure swap corrected in place with an erratum note.
