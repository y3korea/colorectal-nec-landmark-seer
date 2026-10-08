# SEER*Stat export specification

This document specifies the four SEER*Stat case-listing exports that the pipeline expects in `NEC_DATA_DIR`. It was derived from the header line of each of the authors' exports and from the SEER*Stat export dictionaries (`.dic`) saved with them; no data row of any export was copied. The patient-level exports themselves cannot be shared (SEER Research Data Use Agreement) and must be produced by each user with their own SEER*Stat account.

## 1. Database and session

| Setting | Value |
|---|---|
| Software | SEER*Stat 9.1.0 (version as reported by the authors in Supplementary Methods S1; the version number is not recorded in the export files themselves) |
| Database | **Incidence - SEER Research Data, 17 Registries, Nov 2025 Sub (2000-2023)** |
| Full database title as recorded in the export dictionary | Incidence - SEER Research Data, 17 Registries, Nov 2025 Sub (2000-2023) - Linked To County Attributes - Time Dependent (1990-2024) Income/Rurality, 1969-2024 Counties |
| Session type | Case Listing (one row per tumour record) |
| Registries (selection statement `{SEER registry}`) | the 17 SEER registries of the database (registry codes 1, 2, 21, 22, 23, 25, 26, 27, 29, 31, 35, 37, 41, 42, 43, 44, 47) |
| Year of diagnosis | 2000-2023 |
| Follow-up cut-off | 31 December 2023 (study cutoff of the November 2025 submission) |
| Patient selection beyond site / histology | none in SEER*Stat; all other eligibility rules (behaviour, age, diagnostic confirmation, primary by international rules, reporting source, survival known, one record per patient) are applied by `R/00_build_cohorts.R` |

### Export options (identical for all four files)

| Option | Value |
|---|---|
| Variable format | quoted labels (`quotedlabels`): formatted labels, every field in double quotes |
| File format | DOS/Windows (CRLF line ends) |
| Field delimiter | tab |
| Missing character | "NA" |
| Fields with delimiter in quotes | false |
| Remove thousands separators | true |
| Variable names included | true (first line = header) |
| Column variables as stats | false |
| GZipped | false |
| Variable display type | Formatted (labels) for all variables |

## 2. Case selection per export

Site and histology selections as recorded in the title lines of each export dictionary. Primary site uses `Primary Site - labeled`; histology uses `Histologic Type ICD-O-3` (behaviour is filtered later).

| Default file name (`.txt`) | Override variable | Content | Histology (ICD-O-3) | Primary site |
|---|---|---|---|---|
| `colorectal_NEN_SEER17_2000_2023_full_20260806.txt` | `NEC_FILE_NEN` | Colorectal neuroendocrine neoplasms (source of the colorectal NEC cohorts) | 8013, 8041-8045, 8240-8246, 8249 | C18.0-C18.9, C19.9, C20.9 |
| `colorectal_adenocarcinoma_SEER17_2000_2023_full_20260928.txt` | `NEC_FILE_AC` | Colorectal adenocarcinoma (comparator) | 8140-8389, 8480-8490 | C18.0-C18.9, C19.9, C20.9 |
| `SCLC_lung_SEER17_2000_2023_full_20260928.txt` | `NEC_FILE_LUNG` | Lung small-cell carcinoma and large-cell NEC (comparators) | 8002, 8013, 8040-8045 | C34.0-C34.9 (lung and bronchus) |
| `pancreatic_adenocarcinoma_SEER17_2000_2023_full_20260928.txt` | `NEC_FILE_PANC` | Pancreatic adenocarcinoma (comparator) | 8140-8389, 8480-8490, 8500 (main analysis: all listed codes except islet-cell 8150-8157) | C25.0-C25.9 (pancreas) |

Notes. (i) The colorectal NEN export includes the appendix (C18.1); it is removed in `R/00_build_cohorts.R` (colorectal sites kept: C18.0, C18.2-C18.9, C19.9, C20.9). (ii) Pancreatic adenocarcinoma: codes 8140-8389, 8480-8490 and 8500 are kept and islet-cell codes 8150-8157 are removed in the analysis (`panc_main`). (iii) Lung: SCLC = 8041-8045 (main), 8002/8040 added in a sensitivity definition, large-cell NEC = 8013. (iv) Colorectal adenocarcinoma: 8140-8389 and 8480-8490 excluding 8240-8249. (v) The SEER*Stat session files (`.sl`) saved by the authors contain a stored statement library whose entries do not all match the lung and pancreas exports; the selections above are taken from the export dictionaries, which describe the files actually exported.

To use your own exports under different file names, either rename them to the default names above or set `NEC_FILE_NEN`, `NEC_FILE_AC`, `NEC_FILE_LUNG` and `NEC_FILE_PANC` to the file names inside `NEC_DATA_DIR`.

## 3. Variables (101, in this order, identical in all four exports)

Select these variables in the SEER*Stat case-listing *Table* tab in exactly this order (SEER*Stat variable names; all with display type *Formatted*). Variables read by the pipeline are marked **used**; the remaining variables are exported for completeness and documentation but are not read. The cohort builder selects columns by name, so the column order is not critical, but every used variable must be present with exactly this name.

| # | Variable | Read by the pipeline |
|---|---|---|
| 1 | Age recode with <1 year olds and 90+ |  |
| 2 | Sex | **used** |
| 3 | Year of diagnosis | **used** |
| 4 | Race recode (W, B, AI, API) |  |
| 5 | Origin recode NHIA (Hispanic, Non-Hisp) |  |
| 6 | Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic) | **used** |
| 7 | Site recode ICD-O-3 2023 Revision |  |
| 8 | Primary Site - labeled | **used** |
| 9 | Histologic Type ICD-O-3 | **used** |
| 10 | Behavior recode for analysis | **used** |
| 11 | Grade Recode (thru 2017) | **used** |
| 12 | Derived Summary Grade 2018 (2018+) | **used** |
| 13 | Grade Clinical (2018+) |  |
| 14 | Grade Pathological (2018+) |  |
| 15 | Laterality |  |
| 16 | Combined Summary Stage with Expanded Regional Codes (2004+) | **used** |
| 17 | Summary stage 2000 (1998-2017) | **used** |
| 18 | SEER historic stage A (1973-2015) |  |
| 19 | Derived EOD 2018 T Recode (2018+) |  |
| 20 | Derived EOD 2018 N Recode (2018+) |  |
| 21 | Derived EOD 2018 M Recode (2018+) |  |
| 22 | Derived EOD 2018 Stage Group Recode (2018+) |  |
| 23 | Derived AJCC Stage Group, 7th ed (2010-2015) |  |
| 24 | Derived AJCC T, 7th ed (2010-2015) |  |
| 25 | Derived AJCC N, 7th ed (2010-2015) |  |
| 26 | Derived AJCC M, 7th ed (2010-2015) |  |
| 27 | 7th Edition Stage Group Recode (2016-2017) |  |
| 28 | Derived SEER Cmb Stg Grp (2016-2017) |  |
| 29 | Derived SEER Combined T (2016-2017) |  |
| 30 | Derived SEER Combined N (2016-2017) |  |
| 31 | Derived SEER Combined M (2016-2017) |  |
| 32 | Derived AJCC Stage Group, 6th ed (2004-2015) |  |
| 33 | Derived AJCC T, 6th ed (2004-2015) |  |
| 34 | Derived AJCC N, 6th ed (2004-2015) |  |
| 35 | Derived AJCC M, 6th ed (2004-2015) |  |
| 36 | AJCC stage 3rd edition (1988-2003) |  |
| 37 | SEER modified AJCC stage 3rd (1988-2003) |  |
| 38 | T value - based on AJCC 3rd (1988-2003) |  |
| 39 | N value - based on AJCC 3rd (1988-2003) |  |
| 40 | M value - based on AJCC 3rd (1988-2003) |  |
| 41 | TNM Edition Number (2016-2017) |  |
| 42 | RX Summ--Surg Prim Site (1998-2022) | **used** |
| 43 | RX Summ--Surg Prim Site 2023 (2023+) | **used** |
| 44 | RX Summ--Scope Reg LN Sur (2003+) |  |
| 45 | RX Summ--Surg Oth Reg/Dis (2003+) |  |
| 46 | RX Summ--Surg/Rad Seq |  |
| 47 | Reason no cancer-directed surgery | **used** |
| 48 | Radiation recode | **used** |
| 49 | Chemotherapy recode (yes, no/unk) | **used** |
| 50 | Scope of reg lymph nd surg (1998-2002) |  |
| 51 | RX Summ--Reg LN Examined (1998-2002) |  |
| 52 | Surgery of oth reg/dis sites (1998-2002) |  |
| 53 | RX Summ--Systemic/Sur Seq (2007+) | **used** |
| 54 | Time from diagnosis to treatment in days recode | **used** |
| 55 | EOD Primary Tumor Recode (2018+) |  |
| 56 | EOD Regional Nodes Recode (2018+) |  |
| 57 | EOD Mets Recode (2018+) |  |
| 58 | Tumor Size Over Time Recode (1988+) |  |
| 59 | Tumor Size Summary (2016+) |  |
| 60 | Regional nodes examined (1988+) |  |
| 61 | Regional nodes positive (1988+) |  |
| 62 | SEER Combined Mets at DX-bone (2010+) | **used** |
| 63 | SEER Combined Mets at DX-brain (2010+) | **used** |
| 64 | SEER Combined Mets at DX-liver (2010+) | **used** |
| 65 | SEER Combined Mets at DX-lung (2010+) | **used** |
| 66 | Mets at DX-Distant LN (2016+) | **used** |
| 67 | Mets at DX-Other (2016+) | **used** |
| 68 | Lymph-vascular Invasion (2004+ varying by schema) |  |
| 69 | CS tumor size (2004-2015) |  |
| 70 | CS extension (2004-2015) |  |
| 71 | CS lymph nodes (2004-2015) |  |
| 72 | CS mets at dx (2004-2015) |  |
| 73 | CS Tumor Size/Ext Eval (2004-2015) |  |
| 74 | CS Reg Node Eval (2004-2015) |  |
| 75 | CS Mets Eval (2004-2015) |  |
| 76 | EOD 10 - extent (1988-2003) |  |
| 77 | EOD 10 - nodes (1988-2003) |  |
| 78 | EOD 10 - size (1988-2003) |  |
| 79 | Coding system-EOD (1973-2003) |  |
| 80 | COD to site recode | **used** |
| 81 | SEER cause-specific death classification | **used** |
| 82 | SEER other cause of death classification | **used** |
| 83 | Survival months | **used** |
| 84 | Survival Days | **used** |
| 85 | Survival months flag | **used** |
| 86 | Vital status recode (study cutoff used) | **used** |
| 87 | Sequence number | **used** |
| 88 | Primary by international rules | **used** |
| 89 | Record number recode | **used** |
| 90 | Total number of in situ/malignant tumors for patient | **used** |
| 91 | Age recode with single ages and 90+ | **used** |
| 92 | Race recode (with detailed Asian and Native Hawaiian other PI) |  |
| 93 | Year of follow-up recode | **used** |
| 94 | Patient ID | **used** |
| 95 | Type of Reporting Source | **used** |
| 96 | Marital status at diagnosis |  |
| 97 | Median household income inflation adj to 2024 | **used** |
| 98 | Rural-Urban Continuum Code | **used** |
| 99 | Site rec ICD-O-3/WHO 2008 (individual sites only) |  |
| 100 | Race and origin (recommended by SEER) |  |
| 101 | Diagnostic Confirmation | **used** |

41 of the 101 variables are read (`COLS` in `R/00_build_cohorts.R`).

## 4. Checks after export

- First line = 101 quoted variable names separated by tabs; every following line = one tumour record with 101 quoted fields.
- Run `python3 run_all.py --only COHORTS,FLOW` and compare `results/FIG_flow_counts.json` with `reference_results/FIG_flow_counts.json` (e.g. main colorectal NEC n = 1,417, distant 908; distant SCLC 77,707; distant colorectal adenocarcinoma 147,014; distant pancreatic adenocarcinoma 75,634). Different counts usually mean a different database release, registry set, or variable definition.

