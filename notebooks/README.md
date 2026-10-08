# notebooks/

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/y3korea/colorectal-nec-landmark-seer/blob/main/notebooks/reproduce_colab.ipynb)

`reproduce_colab.ipynb` reproduces every SEER-based result and figure of the manuscript from **your own** SEER*Stat exports, in Google Colab or in a local Jupyter session. It is a guided front end to the same pipeline as `run_all.py`; nothing in it changes the analysis.

Corresponding author and maintainer: Wansuk Choi, Department of Physical Therapy, College of Health Sciences, Kyungwoon University (wschoi@ikw.ac.kr).

## SEER Data Use Agreement — read first

* No SEER data are included in the repository or the notebook. The SEER Research Data Use Agreement forbids sharing SEER data or derived patient-level files.
* **In Colab, the SEER exports and the derived cohorts are processed on Google cloud infrastructure** (the Colab virtual machine, and Google Drive if you keep the exports there). Each user must decide whether this is compatible with **their own DUA and institutional rules**. If in doubt, run the notebook **locally** (recommended) or use the `synthetic` mode. In Colab, the `full` and `quick` modes only start after you tick `I_CONFIRM_DUA_PERMITS_CLOUD`.
* The notebook never prints or saves record-level rows (the data check reports file size, row count and column count only), and the last cell deletes `work/` (derived cohorts) plus any export copy or synthetic files from the runtime. In Colab, finish with *Runtime ▸ Disconnect and delete runtime*.
* **Small cells.** Outputs in `results/` are unsuppressed working files and include a few cells based on 1–4 patients. The DUA (item 14) requires statistics on demographic characteristics (age, sex, race/ethnicity) or geography based on 1–4 patients to be suppressed in every publication: run `python3 tools/suppress_small_cells.py` before publishing or sharing any output (it writes `results/publication_safe/`). `reference_results/` already holds publication-safe copies; `tests/compare_with_reference.py` skips the suppressed fields.
* Obtain the data as described in [`../data_spec/SEER_export_specification.md`](../data_spec/SEER_export_specification.md) (four case-listing exports, *Incidence – SEER Research Data, 17 Registries, Nov 2025 Sub (2000–2023)*, 101 variables).

## Sections

| # | Section | What happens |
|---|---|---|
| 0 | Purpose and DUA notice | English, with a short Korean summary |
| 1 | Environment and parameters | detects Colab (`'google.colab' in sys.modules`); Colab form or local environment variables; fetches the repository (Colab) |
| 2 | Dependencies | Python packages from `requirements.txt` (pip in Colab); R (apt-get `r-base` in Colab if `Rscript` is missing) and R packages through `R/install_packages.R` (CRAN cloud mirror, `Ncpus`) |
| 3 | Data check | the four exports: size, data rows, columns, header vs specification (no rows shown); `synthetic` mode generates the files with `tools/make_synthetic_seer_export.py` |
| 4 | Pipeline | the 17 steps of `run_all.py` (cohorts, FLOW, A01–A12, figures) as subprocesses, with per-step status, run time, the step's reproduction-check lines and logs in `results/logs/` |
| 5 | Verification | `tests/check_anchors.py` as a table *manuscript value vs reproduced value vs status*; `tests/compare_with_reference.py` per reference file |
| 6 | Outputs | Figures 1–3 inline, list of outputs, optional zip of aggregate outputs, cleanup |

## Modes

| `MODE` | Data | Bootstrap | Expected verification result |
|---|---|---|---|
| `full` | your SEER exports | B = 2,000 (paper) | all anchors PASS; 17/17 reference files reproduced |
| `quick` | your SEER exports | B = 40 | point-estimate anchors pass; bootstrap-based values and the reference comparison differ |
| `synthetic` | random files in the export layout | B = 40 | anchor failures and reference mismatches are **expected** (results are meaningless) |

`BOOT` (Colab) / `NEC_BOOT` (local) overrides the bootstrap size.

## Running in Google Colab

1. Open the notebook: *File ▸ Upload notebook* (select `notebooks/reproduce_colab.ipynb`), or *File ▸ Open notebook ▸ GitHub* and paste the repository URL, or use the badge above.
2. Fill in the form in Section 1:
   * `REPO_SOURCE = "github"` with `REPO_URL` (default: https://github.com/y3korea/colorectal-nec-landmark-seer.git), optional `REPO_REF` (branch, tag or commit); or `REPO_SOURCE = "drive"` with `REPO_DRIVE_PATH` pointing to the repository folder copied to Drive, e.g. `/content/drive/MyDrive/.../02_GitHub_Repository/colorectal-nec-landmark-seer`. The code is always copied to `/content/colorectal-nec-landmark-seer` (VM disk), so no derived file is ever written to Drive.
   * `DATA_DIR`: Drive folder with the four exports, e.g. `/content/drive/MyDrive/SEER_exports`; `COPY_DATA_TO_RUNTIME` copies them to the VM disk for faster reading.
   * `MODE`, `BOOT`, `I_CONFIRM_DUA_PERMITS_CLOUD`, `MAKE_ZIP`, `CLEANUP_AT_END`.
3. *Runtime ▸ Run all*, authorise the Drive mount, and wait.

## Running locally (recommended for SEER data)

```bash
cd colorectal-nec-landmark-seer
python3 -m venv .venv && . .venv/bin/activate && pip install -r requirements.txt jupyter
export NEC_REPO_DIR=$PWD NEC_DATA_DIR=/path/to/seer/exports NEC_MODE=full    # or quick / synthetic
jupyter nbconvert --to notebook --execute --ExecutePreprocessor.timeout=-1 \
    notebooks/reproduce_colab.ipynb --output-dir /some/private/folder --output reproduce_executed.ipynb
```

or open the notebook in Jupyter (`jupyter lab`) after setting the same variables. Optional variables: `NEC_BOOT`, `NEC_MAKE_ZIP=1`, `NEC_CLEANUP=0` (keep `work/`), `NEC_SYNTH_SCALE`, `NEC_HEARTBEAT` (seconds between progress lines). Keep executed copies of the notebook private until you have checked them; the notebook is designed to show aggregate output only, but any cell you add is your responsibility. Commit the notebook only with outputs cleared (`jupyter nbconvert --clear-output --inplace notebooks/reproduce_colab.ipynb`).

## Expected run time

Measured locally on 2026-10-08 (Apple silicon, 10 cores, macOS; notebook executed headless with `jupyter nbconvert --execute`, fresh copy of the repository each time):

| Mode | Pipeline (17 steps) | Whole notebook | Verification |
|---|---|---|---|
| `full` (B = 2,000) | 13.8 and 15.3 min (two runs) | 14.0 and 15.5 min | 60/60 anchors PASS; 17/17 reference files identical (max \|diff\| 0) |
| `synthetic` (B = 40) | 7.4 and 8.9 min | 7.6 and 9.0 min | 0/60 anchors, 0/17 reference files (expected: random data) |

The longest steps are A05R and A08 (multi-state Cox standardization, 3–5 min each), A11 (≈ 2 min) and A06 (≈ 1.4 min); building the cohorts from the 1.7 GB of exports takes ≈ 0.5 min. `quick` mode takes about as long as `synthetic` because the A05R/A08 model fits do not depend on the bootstrap size.

**Colab (not measured):** a free runtime has 2 vCPUs and ≈ 12 GB RAM, so expect several times the local figures — roughly 5–15 min for the first installation (R packages are compiled from source), plus very roughly 30–60 min for `full` and 15–25 min for `quick`/`synthetic`. Copying the exports from Drive adds a few minutes. Use `quick` first to check the set-up.

## Troubleshooting

| Symptom | Fix |
|---|---|
| `Set REPO_URL ...` | set `REPO_URL` to the repository address, or use `REPO_SOURCE = "drive"` |
| `PermissionError: MODE='full'/'quick' processes SEER data on Google cloud` | tick `I_CONFIRM_DUA_PERMITS_CLOUD` only if your DUA allows it; otherwise run locally or choose `synthetic` |
| `SEER exports not found` | check `DATA_DIR` and the four file names (Section 0); set `NEC_FILE_NEN/AC/LUNG/PANC` for other names |
| header = specification shows `NO` | the export does not contain the 101 variables of the specification; re-export from SEER*Stat |
| R package compilation fails in Colab | re-run the R cell; missing system libraries can be added with `apt-get install -y <lib>-dev` |
| pip dependency-conflict warnings in Colab | harmless for this pipeline (it runs in subprocesses); if later cells fail to import, *Runtime ▸ Restart session* and *Run all* |
| Colab runs out of memory while building cohorts (the colorectal adenocarcinoma export is 1.26 GB) | use a high-RAM runtime or run locally |
| Colab session disconnects during `full` | a disconnect resets the VM: run all cells again from Section 1 (try `quick` first; free runtimes disconnect after inactivity). Within a live session, a failed Section 4 cell can be re-run on its own after fixing the cause (cohorts in `work/` are reused) |
| a step fails in `full` mode | read the printed log tail and `results/logs/NN_STEP.log`; strict mode stops when a manuscript anchor is not met |
| figures use DejaVu Sans instead of Arial | expected without Arial (e.g. Colab); numbers are unaffected |
| `compare_with_reference` reports differences in `quick` / `synthetic` mode | expected (B = 40 or random data) |
