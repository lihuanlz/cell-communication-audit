# Cell Communication Audit

A unified audit program asking a single question — **how do cells actually
encode and transmit information?** — across signalling pathways (NF-κB, p53,
ERK, Msn2, PdPC), GPCR atlases, chemotaxis receptors, retinal ganglion cells,
and human neuron electrophysiology.

Starting from a scale-invariance principle (TCS), we show that naive
parameter fitting in these systems is degenerate (Fisher barriers up to
chi ~ 1e300), and that the information lives in scale-invariant, directly
extractable quantities (time-channel occupancy, pulse counting, structural
constants) rather than in fitted amplitudes.

Companion code-and-data repository of the manuscript
*Nature_main_v05* (see `01_manuscript/`).

## Layout

| Path | Content |
|---|---|
| `01_manuscript/` | Main text and SI (Markdown + Word) |
| `02_figures/` | Figures (PNG + SVG) |
| `03_data_csv/` | Numeric outputs behind figures/tables (CSV) |
| `04_code/` | Audit and analysis pipelines |
| `05_validation_records_index.md` | Registry of validation cards |
| `06_supplement_interactive/` | Interactive supplements |
| `SI_verification_scripts_EN/` | English SI verification scripts (TCS_scripts_EN) |

## Method note

This work was produced under an explicit human–AI collaborative protocol:
pre-registered predictions, sealed verdict cards, alias/archival audits, and
post-seal replications are logged in the validation registry; the
collaboration protocol itself is described in the manuscript Methods.

## License

For academic use; manuscript under review.
