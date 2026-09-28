# Release Notes Draft — cell-communication-audit v06 (2026-09-28)

Companion release for the manuscript *"Scale-free statistics are the currency of cellular communication: an audit across four signalling systems"* (manuscript v06 / SI v06, preprint candidate).

## New in v06 (vs v05)

- **Cancer-chain prediction arm (P11–P16)** — new main-text Results section "Hijacked, not interrupted" and SI section S19:
  - P11: p53 counting-law cross-cell-line audit (slope correlation rho = +0.857, p = 0.0068)
  - P12: gain-distribution and resistance zero-direction audits (ABL gate 2.94x, p = 0.0326; EGFR DMS 2.2x, p = 0.0022)
  - P13: EGF channel-capacity audit (WT 0.426 bit; oncogenic mutants -0.05 to -0.41 bit; KSG robustness concordant)
  - P14: p53 downstream decoding audit (AUC ratio 1.65)
  - P15: TCGA bulk-tumour decoupling test (negative in all four cancer types; registered negative result)
  - P16: frozen forward prediction for KRAS-G13D (pre-registration card included, awaiting data)
- `07_audit_cards_tumour_chain/`: 5 verdict cards + 6 pre-registration cards (Chinese bodies = authoritative working record).
- `04_code/tumour_chain_P11_P16/`: 12 self-contained audit scripts; `03_data_csv/tumour_chain_P11_P16/`: result JSONs + digitised/input CSVs.
- `EDFig16_cancer_chain` (PNG + SVG) + generator script `04_code/make_edfig16.py`.
- README Addendum 2026-09-28; validation records index extended; references 77–81 added (Wang 2025 page numbers verified).
- Errata: P15 verdict card COADREAD/LUAD figure swap corrected in place with erratum note; P12b A5 arm formally archived.

## Carried over from v05

- Full blind-adjudication stack P2–P7 with frozen seeds and bootstrap x2000; cloud re-executions bit-identical.
- 2026-09-27 campaign: P9 Bicoid anchor, P10 gain confinement, cascade audits, Wnt arm, S27 reconciliation.
- GPCR meta-audit (Biased Signaling Atlas, 17,987 ligand pairs, 21.5% > 0.3 dex).
- p53 Fisher census, molecular-layer closure (S14, codes 82–93).

## Reproducibility

Start from `04_code/CODE_MAP.csv` and `05_validation_records_index.md`. All adjudications run from frozen scripts with fixed seeds; public datasets listed under Data availability in README.

## Citation

Sister repository: `alpha-ion-channel-models` (tag v05). Preprint DOI pending (Zenodo snapshot to be taken from this tag).
