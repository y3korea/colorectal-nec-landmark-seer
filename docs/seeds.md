# Random seeds

Authoritative list of the random seeds used by the code in this repository. Every bootstrap seed is also written to the result JSON (`seed` fields), and `tests/compare_with_reference.py` requires exact equality with `reference_results/`. Bootstrap resampling is at the patient level, with replacement, within the analysis cohort or stratum.

| Analysis (script) | Seed(s) | Notes |
|---|---|---|
| Table 2 and Supplementary Table S8, pre-specified landmarks (`python/A01_landmark_main.py`) | 20260917 | B = 2,000 |
| Monthly landmarks 0–24 and subsets (`python/A01_landmark_main.py`) | 20260928 | |
| Window-width sensitivity and model comparison (`python/A02_window_sensitivity_model_comparison.py`) | 20261007 | identical resamples across windows and endpoints |
| Equivalence / TOST, Supplementary Table S16 (`python/A03_equivalence_tost.py`) | NEC 20261007; **SCLC 20260929; colorectal AC 20260930; pancreatic AC 20260931** | comparator seeds = 20260928 + 1/2/3 in the order SCLC, AC, PANC |
| Case definitions, subtypes, era (`python/A04_expanded_cohort_codes_subtypes.py`) | 20260928 | same seed for every cohort-specific bootstrap |
| Strata and time-updated chemotherapy (`python/A05_heterogeneity_stratified.py`) | 20261007 | one `numpy.random.SeedSequence(20261007)` child stream per stratum (`seed_spawn_index` in `A05.json`) |
| Isotonic (monotonicity) parametric bootstrap (`R/A06_piecewise_poisson_suite.R`) | 20261007 | B = 500 |
| Split-sample (calendar) validation (`python/A07_split_sample_validation.py`) | derivation 20261007; confirmation 20261008; eras 20261017–20261019; period comparators 20261027 | SEED, SEED + 1, SEED + 10…12, SEED + 20 |
| Standardization bootstrap (`R/A08_standardization_extended.R`) | 20261007 | B = 500; indices drawn once before forking, so results do not depend on `NEC_CORES` |
| Comparator reference samples of 5,000 (`R/A08_standardization_extended.R`) | 20260928 | colorectal AC then SCLC from one stream; re-seeded (20260928) for pancreatic AC |
| Kernel-smoothed hazards (`R/A09_smoothed_hazards_muhaz.R`) | 20260929 | B = 500; same resample at every bandwidth |
| Cross-tumour comparison, Table 3 (`python/A11_cross_tumour_comparison.py`) | 20260928 + group index: NEC 20260928; **colorectal AC 20260929; SCLC 20260930**; pancreatic AC 20260931; lung LCNEC 20260932 | further groups: colorectal small-/large-cell 20260933/20260934; SCLC definitions 20260935/20260936; chemotherapy from diagnosis 20260937–20260940 (NEC, SCLC, AC, PANC) and 20260956 (LCNEC); 3-month-conditioned chemotherapy subsets 20260941–20260944 (NEC, SCLC, AC, PANC) and 20260957 (LCNEC); cytology-inclusive 20260945–20260947; regional 20260948–20260951; chemotherapy no/unknown 20260952–20260955 (full list: `A11.json` → `group_index`) |
| Lung LCNEC rate ratios and shape tests (`R/A12_lcnec_rate_ratios.R`) | none | no random numbers are drawn |
| Cohort building, flow counts, A10 data quality, figures (incl. `python/figures/graphical_abstract.py`), small-cell suppression (`tools/suppress_small_cells.py`) | none | deterministic |

Note on the cross-tumour seeds. The TOST analysis (A03) and the cross-tumour comparison (A11) use different comparator orders: in A03 SCLC = 20260929 and colorectal AC = 20260930, whereas in A11 colorectal AC = 20260929 and SCLC = 20260930 (group index NEC 0, AC 1, SCLC 2, PANC 3, LCNEC 4). Pancreatic AC is 20260931 in both. The Table 3 NEC confidence intervals in the manuscript come from the main-analysis seed 20260917 (A01); A11 uses 20260928 for NEC, so its NEC intervals can differ by a few tenths of a percentage point.
