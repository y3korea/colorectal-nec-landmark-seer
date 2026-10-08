# What this repository contains (deposit inventory)

This page lists what the public repository `colorectal-nec-landmark-seer` contains and what it does not. It is the reference for the manuscript's Data Availability Statement and the Supplement section "Code and analysis-output availability".

## Contained

| Item | Location |
|---|---|
| SEER*Stat case-listing specification: database, session settings, site/histology selection of each export, the full list of 101 variables (in order) | `data_spec/SEER_export_specification.md` |
| Cohort-building script (eligibility filters, step counts of Supplementary Table S4) | `R/00_build_cohorts.R`, `python/flow_counts.py` |
| Analysis and figure code (independent re-implementation used for the verification and revision analyses A01–A12 and the manuscript figures) | `python/`, `R/` |
| Random seeds | in the scripts; listed in `docs/seeds.md`; recorded in every result JSON |
| Aggregate outputs of the verified run (every quantity used in Tables 2–3, Figures 1–3, Supplementary Tables S2–S43 and Supplementary Figures S1–S7), as publication-safe copies in which cells based on 1–4 patients are suppressed (SEER DUA item 14) | `reference_results/*.json`; regenerated in `results/` (incl. the 0.25-month hazard grid `A09_hazard_grid.csv`, which is not stored here; monthly values in `reference_results/A09.json`) |
| Small-cell suppression script (publication-safe copies of the result JSON) | `tools/suppress_small_cells.py` |
| Graphical Abstract script | `python/figures/graphical_abstract.py` |
| Map of scripts to tables and figures | `docs/output_map.md` |
| Analysis-sequence table (Supplementary Table S1) | `docs/analysis_sequence.md` |
| Figure legends | `docs/figure_legends.md` |
| STROBE/RECORD checklist | `docs/STROBE_RECORD_checklist.md` |
| Reproduction tests and small-cell check | `tests/check_anchors.py`, `tests/compare_with_reference.py`, `tests/check_small_cell_text.py` |
| Synthetic-data generator (runs the pipeline without SEER access) | `tools/make_synthetic_seer_export.py` |

Software of the verified run: Python 3.13.6 (NumPy 2.2.0, pandas 2.3.3, SciPy 1.17.1, Matplotlib 3.10.8) and R 4.5.3 (survival 3.8.6, sandwich 3.1.1, data.table 1.18.4, cmprsk 2.2.12, muhaz 1.2.6.5, jsonlite 2.0.0).

## Not contained

| Item | Reason |
|---|---|
| SEER data: raw case listings, cohort files, any patient-level row | forbidden by the SEER Research Data Use Agreement; obtain from <https://seer.cancer.gov/data/access.html> |
| SEER*Stat session files (`.sl`) and export dictionaries (`.dic`) | not included. The `.sl` files store a statement library whose entries do not all match the lung and pancreas exports, and an export dictionary records the authors' local output folder. Everything needed to recreate the exports is in `data_spec/SEER_export_specification.md`, which was derived from the `.dic` files. |
| The authors' primary analysis code (Python 3.12.10 / NumPy 1.26.4 / pandas 2.2.3; R 4.6.1 / survival 3.8-9) | not included; this repository is the independent re-implementation, which reproduces every anchor quantity of the primary analysis (`tests/check_anchors.py`) |
| Dated written analysis plans of 16 and 28 September 2026 | not included; available from the corresponding author on request (wschoi@ikw.ac.kr); their content and dates are summarised in `docs/analysis_sequence.md` |
| Table 1 cells marked †, Methods S6 / Table S44 (external published cohorts) | not regenerated (see `docs/output_map.md`) |

## Wording used in the manuscript (version of 8 October 2026)

Data Availability Statement (excerpt): the cohort-building, analysis and figure code with all random seeds, the SEER*Stat export specification, aggregate reference results for every reported estimate (cells based on 1–4 patients suppressed in accordance with the SEER Research Data Use Agreement), synthetic demonstration data and a Jupyter notebook are available in this repository; "The dated written analysis plans of September 16 and 28, 2026 (Supplementary Table S1) are not part of the repository and are available from the corresponding author on request."

Supplement, "Code and analysis-output availability": (i) the SEER*Stat export specification; (ii) the cohort-building script; (iii) the analysis toolkit and scripts (A01–A12 and the figure scripts) with the seeds of Supplementary Methods S2; (iv) the publication-safe aggregate reference results; (v) the synthetic-data generator and the Colab notebook; (vi) the Graphical Abstract script, the figure legends, the analysis-sequence table and the STROBE/RECORD checklist. The dated written analysis plans, the SEER*Stat session files and the authors' primary analysis code are not part of the repository; the plans are available from the corresponding author on request.
