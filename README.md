# colorectal-nec-landmark-seer

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23249641.svg)](https://doi.org/10.5281/zenodo.23249641) [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/y3korea/colorectal-nec-landmark-seer/blob/main/notebooks/reproduce_colab.ipynb)

Analysis code for:

> Wang J, Zhang S, Liang W, Choi W. **Conditional one-year cancer mortality shows no early decline in distant-stage colorectal neuroendocrine carcinoma: a SEER landmark analysis with cross-tumor comparison.** Submitted to *International Journal of Cancer* (Research Article, Cancer Epidemiology), 2026.

Corresponding author and repository maintainer: **Wansuk Choi**, Department of Physical Therapy, College of Health Sciences, Kyungwoon University, Gumi, Republic of Korea; wschoi@ikw.ac.kr; ORCID [0000-0003-1324-6723](https://orcid.org/0000-0003-1324-6723).

## Purpose

The repository regenerates every SEER-based result and figure of the manuscript from the user's own SEER*Stat exports:

- conditional 12-month probability of cancer death after each monthly landmark (Aalen–Johansen, other-cause death as a competing event) in distant-stage colorectal neuroendocrine carcinoma (NEC), with patient-level bootstrap intervals;
- the broken-stick "plateau" description, window-width sensitivity, model comparison and equivalence (TOST) framing;
- piecewise-exponential (Poisson) interval rate ratios with cluster-robust errors, multiplicity adjustment and a monotonicity test; kernel-smoothed hazards;
- strata, direct standardization (per-landmark multi-state Cox), proportional-hazards and Fine–Gray checks, case-definition and calendar-split analyses;
- the cross-tumour comparison with small-cell lung cancer, colorectal and pancreatic adenocarcinoma and lung large-cell NEC;
- Tables 2–3, Figures 1–3, Supplementary Tables S2–S43, Supplementary Figures S1–S7 and the Graphical Abstract (see [`docs/output_map.md`](docs/output_map.md)).

## Citation

Please cite the article above and this repository:

> Wang J, Zhang S, Liang W, Choi W. colorectal-nec-landmark-seer: reproducible analysis code for a SEER landmark study of conditional cancer mortality in distant-stage colorectal neuroendocrine carcinoma. Version 1.0.0. Zenodo; 2026. doi:[10.5281/zenodo.23249641](https://doi.org/10.5281/zenodo.23249641)

The concept DOI [10.5281/zenodo.23249640](https://doi.org/10.5281/zenodo.23249640) always resolves to the latest version. Machine-readable metadata: [`CITATION.cff`](CITATION.cff), [`.zenodo.json`](.zenodo.json).

## SEER data access and Data Use Agreement (read first)

**This repository contains no SEER data and no patient-level records, and none may ever be added.**
SEER data are released only under the *SEER Research Data Use Agreement*, which forbids sharing or redistributing the data (including derived patient-level files). To reproduce the analysis you need your own SEER account:

1. Request access at <https://seer.cancer.gov/data/access.html> and accept the Data Use Agreement.
2. Install SEER*Stat and create the four case-listing exports exactly as specified in [`data_spec/SEER_export_specification.md`](data_spec/SEER_export_specification.md) (database **Incidence - SEER Research Data, 17 Registries, Nov 2025 Sub (2000-2023)**, 101 variables, quoted labels, tab-delimited).
3. Keep the exports and everything in `NEC_WORK_DIR` (patient-level cohort files) on your own machine. `.gitignore` excludes `data/`, `work/`, `*.txt` exports, `cohort_*.csv`, `*.rds` and `results/`; do not override it.

What the repository does and does not contain (e.g. no SEER*Stat `.sl` session files) is listed in [`docs/data_availability.md`](docs/data_availability.md).

Only aggregate results (JSON/Markdown/figures) are produced in `results/` and `figures/output/`. `reference_results/` holds the aggregate JSON outputs of the authors' verified run (absolute paths removed; two per-patient fields of A10 — survival days and diagnosis years of the 10 unknown-cause deaths — were removed and are no longer written by `python/A10_data_quality_censoring.py`).

### Small cells (SEER DUA item 14) and publication-safe copies

The Data Use Agreement requires that "statistics about any geography (county, state, registry) or demographic information (age, sex, race/ethnicity) based on counts 1 to 4 must be suppressed in all publications", and that no information on an individual is published. The analysis scripts write their outputs **unchanged** (they are working files for your own SEER-licensed environment and contain a few such cells). Before publishing or sharing any output, make publication-safe copies:

```bash
python3 tools/suppress_small_cells.py            # results/*.json + figure_data.json -> results/publication_safe/
```

The tool replaces every affected cell by the marker `"suppressed (SEER small-cell rule)"`, applies complementary suppression where a suppressed value could otherwise be recovered by subtraction, records what it did under the key `_suppression` of each changed file, and then audits all copies for remaining 1–4 cells attached to demographic or geographic labels (exit code 2 if any are found). In the authors' run it suppresses 92 fields: the stage (number and % distant), diagnosis-period shares, median age, % male and % rectum of the three morphology codes with one patient each (8042, 8043, 8044; only their patient counts are kept), % male of code 8241 (8 patients, so the share implies 1–4 men or women) and, complementarily, % male of code 8045 (`A04.json`); the race/ethnicity category "Non-Hispanic Unknown Race" (1–4 patients) and, complementarily, "Non-Hispanic American Indian/Alaska Native" (`A08.json`; only their sum, the category "Other", is published), and the landmark-6 cause-specific hazard ratio, CI and P of the race/ethnicity level "Other" (1–4 patients at risk at that landmark; `A08.json` `coef_all_combined_cancer_csh`); the lung-export eligibility step "age ≥ 18" (which excluded 1–4 records) with the following unchanged step (`figure_data.json`); and, in `A10.json`, information on individual patients withheld exactly as in the Supplement — survival-day minimum and maximum, maximum days to first treatment of every group, the categories that the Supplement pools because they hold single patients (days to treatment 91–180 and > 180; 4th and 5th primaries and 5 or 7 tumours, all stages; 3 and 4 involved metastatic sites; the rare systemic-therapy/surgery sequences; the no-surgery reasons "recommended, unknown if performed" and "patient refused") together with their pooled partners and cross-tabulations, every cell outside the distant column of the later-primary rows (2nd to 5th primary) of the sequence × stage table (suppressing only its 1–4 cells would not protect them, because they could be recovered from the published row and stage totals), and the multiple-primary sequence of the few patients present in both the NEC and a comparator cohort. **`reference_results/` contains only these publication-safe copies**, and `tests/compare_with_reference.py` skips the marked fields and compares everything else. The `results/*.md` summaries and the per-step logs are working outputs and are not suppressed; the scripts do not write demographic counts of 1–4 into them (for example, `R/A06_piecewise_poisson_suite.R` reports only the pooled "other or unknown" race/ethnicity level). `tests/check_small_cell_text.py` scans the code, documentation and notebooks (and, with `--also results`, the Markdown reports of a run) for race/ethnicity, sex, age or geography labels followed by a count of 1–4; `run_all.py` runs it after the other tests. The manuscript and Supplement apply the same rule (Table 1 combines "Non-Hispanic American Indian or Alaska Native" and "Unknown"; Supplementary Table S2 combines codes 8042–8044 and suppresses % male of codes 8045 and 8241; Supplementary Table S40(e) suppresses the landmark-6 hazard ratio of the race/ethnicity level "Other or unknown"). Besides labelled counts, the JSON audit also flags any % male whose implied number of men or women is 1–4.

## Repository layout

```
README.md  LICENSE  CITATION.cff  requirements.txt  run_all.sh  run_all.py
R/            config.R (paths) · pp.R (shared helpers) · 00_build_cohorts.R · A04/A05/A06/A07/A08/A09/A12 analyses · install_packages.R
python/       config.py (paths) · lm.py (landmark / Aalen-Johansen toolkit) · flow_counts.py · A01-A11 analyses · FIG_manuscript_final.py (+ FIG_final_figures.py helpers)
              figures/graphical_abstract.py (Graphical Abstract from stored results)
tests/        check_anchors.py (key manuscript numbers) · compare_with_reference.py (every JSON vs reference_results/; suppressed cells skipped)
              check_small_cell_text.py (free text: no demographic/geographic counts of 1-4 in code, docs, notebooks, reports)
reference_results/   aggregate JSON of the verified run (A01-A12, FIG_flow_counts, figure_data), publication-safe (cells of 1-4 suppressed)
data_spec/    SEER_export_specification.md (database, session, selections, 101 variables)
docs/         output_map.md (script -> table / figure / supplementary item) · seeds.md · figure_legends.md · analysis_sequence.md
              STROBE_RECORD_checklist.md · data_availability.md (what is and is not deposited)
tools/        make_synthetic_seer_export.py (synthetic exports for a run without SEER access)
              suppress_small_cells.py (publication-safe copies of the result JSON; SEER DUA item 14)
notebooks/    reproduce_colab.ipynb (step-by-step reproduction in Google Colab or local Jupyter; see notebooks/README.md)
data/raw/     your SEER exports (git-ignored)
```

## Environment

| Component | Verified version | Install |
|---|---|---|
| Python | 3.13.6 | in a virtual environment (see below): numpy 2.2.0, pandas 2.3.3, scipy 1.17.1, matplotlib 3.10.8 |
| R | 4.5.3 | `Rscript R/install_packages.R` (data.table 1.18.4, survival 3.8.6, sandwich 3.1.1, jsonlite 2.0.0, cmprsk 2.2.12, muhaz 1.2.6.5; `NEC_PIN_VERSIONS=1` installs exactly these) |
| Font | Arial | manuscript figures use Arial; without it they fall back to DejaVu Sans with a warning (`NEC_REQUIRE_ARIAL=1` makes this an error) |
| OS | macOS 26 (Apple silicon), 10 cores, 16+ GB RAM recommended (the colorectal adenocarcinoma export is 1.26 GB) | |

Python packages: install into a virtual environment (system Pythons such as Homebrew's are "externally managed" under PEP 668 and refuse a plain `pip install`):

```bash
python3 -m venv .venv && . .venv/bin/activate     # .venv/ is git-ignored
pip install -r requirements.txt
```

With the virtual environment active, `./run_all.sh` and `python3 run_all.py` use its interpreter (or set `PYTHON=/path/to/python` for `run_all.sh`).

`Rscript` must be on the `PATH`. A08 uses forked workers (`NEC_CORES`, default 4; 1 on Windows); results do not depend on the number of cores because the bootstrap indices are drawn before forking.

## Configuration

All paths are set through environment variables. Relative values are resolved against the directory from which you start `run_all.sh` / `run_all.py` (the runner converts them to absolute paths and passes those to every step); unset variables default to folders inside the repository root. Absolute paths are the safest choice. Shared logic: [`python/config.py`](python/config.py) and [`R/config.R`](R/config.R).

| Variable | Default | Meaning |
|---|---|---|
| `NEC_DATA_DIR` | `./data/raw` | folder with the four SEER*Stat exports |
| `NEC_WORK_DIR` | `./work` | derived patient-level cohort files (`cohort_*.csv`) and the aggregate `flow_all.csv` |
| `NEC_RESULTS_DIR` | `./results` | aggregate results (`*.json`, `*.md`, `A09_hazard_grid.csv`, `logs/`) |
| `NEC_FIG_DIR` | `./figures/output` | analysis figures; manuscript figures in `manuscript_final/` |
| `NEC_BOOT` | `2000` (as in the paper) | bootstrap replicates of the main analyses; the secondary bootstraps (500 / 1,000 in the paper) are scaled proportionally |
| `NEC_QUICK` | `0` | `1` = smoke test: `NEC_BOOT` defaults to 40, anchors informational, pipeline keeps going after a failing step |
| `NEC_STRICT` | `1` (`0` in QUICK mode or on synthetic data) | stop when a reproduction anchor of the paper is not met |
| `NEC_FILE_NEN`, `NEC_FILE_AC`, `NEC_FILE_LUNG`, `NEC_FILE_PANC` | the authors' file names (see data spec) | export file names inside `NEC_DATA_DIR` |
| `NEC_CORES` | `4` | forked workers for the A08 bootstrap |

Seeds are fixed in the scripts (20260917 for Table 2; 20260928, 20260929, 20261007, 20261008 and documented offsets elsewhere), are recorded in every result JSON and are listed by analysis in [`docs/seeds.md`](docs/seeds.md).

## One-command run

```bash
export NEC_DATA_DIR=/path/to/your/seer/exports      # or copy the four .txt files to data/raw/
./run_all.sh                                        # = python3 run_all.py
```

`run_all.py` runs the 17 steps in order (`python3 run_all.py --list`), writes one log per step to `results/logs/`, then runs `tests/check_anchors.py` and `tests/compare_with_reference.py`. Useful options: `--quick`, `--boot N`, `--from STEP`, `--only A01,A03`, `--keep-going`, `--no-tests`.

Two scripts run after the pipeline, on its stored results:

```bash
python3 python/figures/graphical_abstract.py     # Graphical Abstract -> figures/output/manuscript_final/GraphicalAbstract.{pdf,png}
python3 tools/suppress_small_cells.py            # publication-safe copies -> results/publication_safe/
```

### Without SEER access (synthetic data)

```bash
python3 tools/make_synthetic_seer_export.py --out data/raw   # 4 synthetic exports with the exact 101-column headers
./run_all.sh --quick                                          # whole pipeline, small bootstrap
```

The synthetic files are random numbers with the exact layout of the SEER exports and are labelled as synthetic (a `SYNTHETIC_DATA_README.txt` marker that the pipeline detects; manuscript figures get a "SYNTHETIC DATA" watermark). The pipeline runs end to end, but **the results are meaningless and do not match the paper**; the anchor check is reported as informational and the reference comparison is skipped.

### Notebook (Google Colab or local Jupyter)

[`notebooks/reproduce_colab.ipynb`](notebooks/reproduce_colab.ipynb) runs the same 17 steps one by one with dependency installation, a count-only data check, per-step timing, and a verification table (manuscript value vs reproduced value). Modes `full`, `quick` and `synthetic`; parameters come from a Colab form or, locally, from `NEC_REPO_DIR`, `NEC_DATA_DIR` and `NEC_MODE`. **Running it in Colab processes SEER data on Google cloud infrastructure; check whether your Data Use Agreement permits this (local Jupyter is recommended).** See [`notebooks/README.md`](notebooks/README.md).

## Expected runtime

Full paper settings (B = 2,000), complete pipeline from the raw exports including both checks, each from a fresh copy of the repository on Apple silicon (10 cores, macOS): **about 12–17 minutes wall-clock** (three verified runs on 2026-10-08: 12.3, 13.8 and 16.6 min). The longest steps are the multi-state Cox standardizations (A05R ≈ 3.3–4.4 min, A08 ≈ 3.0–3.7 min), A11 (≈ 1.6–2.1 min), A06 (≈ 1.3–2.1 min) and building the cohorts from the 1.7 GB of exports (≈ 0.5–1.3 min, depending on disk cache). The synthetic QUICK run takes ≈ 7 min (dominated by A05R/A08, whose model fits do not depend on the bootstrap size). Expect longer times on slower disks or machines with less memory.

## Outputs and the manuscript

See [`docs/output_map.md`](docs/output_map.md) for the full map. In brief:

| Manuscript item | Script(s) | Result |
|---|---|---|
| Table 2 | `python/A01_landmark_main.py` | `results/A01.json` → `table2` |
| Table 3 | `python/A11_cross_tumour_comparison.py`, `R/A06_piecewise_poisson_suite.R`, `R/A12_lcnec_rate_ratios.R` | `A11.json` → `table3`; `A06.json`; `A12.json` |
| Figures 1–3, Suppl. Figures S1–S7 | `python/FIG_manuscript_final.py` | `figures/output/manuscript_final/*.png|pdf`, `figure_data.json` |
| Graphical Abstract | `python/figures/graphical_abstract.py` (reads `A06.json`, `A11.json`, `A12.json`) | `figures/output/manuscript_final/GraphicalAbstract.png|pdf`, `graphical_abstract_data.json` |
| Supplementary Tables S2–S43 | A01–A12 | `results/A*.json` / `A*.md` |
| Table 1 | A04/A05/A07/A08/A10 (re-derived cells) | cells marked † in the manuscript are carried over from the authors' original analysis and are not regenerated here |

## Reproduction anchors

`tests/check_anchors.py` asserts the key numbers of the paper (all must pass on the authors' exports):

| Quantity | Value |
|---|---|
| Distant colorectal NEC, 12-month cancer-death probability at landmarks 0/3/6/12/24 months | 72.29 / 71.62 / 72.33 / 60.03 / 30.51 % |
| Numbers at risk at those landmarks | 908 / 553 / 396 / 169 / 48 |
| Main colorectal NEC cohort / distant stage | n = 1,417 / 908 |
| Broken-stick plateau end (95% CI) | 6.75 months (4.00–9.75) |
| Adjusted rate ratio, months 6–9 vs 3–6 | 1.49 (1.17–1.90) |
| SCLC, landmarks 0/6/12 | 69.98 / 69.49 / 61.93 % |
| Colorectal adenocarcinoma, landmarks 0/6/12 | 40.59 / 33.92 / 34.41 % |
| Pancreatic adenocarcinoma, landmarks 0/6/12 | 77.12 / 69.72 / 65.38 % |
| Lung LCNEC (A12): n distant; RR 6–9 vs 3–6; shape test vs colorectal NEC | 3,185; 0.97 (0.84–1.11); W = 15.35, P = 0.018 (2016–2023: W = 19.77) |
| Morphology × organ (A12): RR 6–9 vs 3–6 colorectal small/large, lung small/large | 1.75 / 1.18 / 1.38 / 0.97; interaction 1.04 |

`tests/compare_with_reference.py` then compares every produced `results/*.json` (and `figure_data.json`) with `reference_results/` — tolerance 1e-6 for estimates, exact equality for seeds; run times, dates, software versions and file locations are not compared, and the 92 fields marked `"suppressed (SEER small-cell rule)"` in the publication-safe reference copies are skipped (counted in the report: A04 23, A08 6, A10 61, figure_data 2). In the verified run all 17 reference files (36,945 values) were reproduced exactly (maximum absolute difference 0), and the rebuilt cohort files were byte-identical to those of the authors' analysis; after the final small-cell suppression of 8 October 2026 the comparison was repeated from a fresh copy (8 October 2026, full pipeline at B = 2,000 from the raw exports, 16.9 min: 17/17 steps, 60/60 anchors, 17/17 reference files reproduced with 36,851 values compared, maximum absolute difference 0, 92 suppressed fields skipped, `tests/check_small_cell_text.py` clean on the code, documentation and the run's Markdown reports; `python3 tools/suppress_small_cells.py` on that run reproduced the stored publication-safe copies, and `python3 python/figures/graphical_abstract.py` regenerated a Graphical Abstract PNG pixel-identical to the submitted file).

Each analysis script also carries its own reproduction check against these anchors and stops in strict mode when they are not met.

## License

Code and documentation: MIT (see [`LICENSE`](LICENSE)), © 2026 Jiang Wang, Shuchun Zhang, Wei Liang and Wansuk Choi. The license does not cover SEER data.
