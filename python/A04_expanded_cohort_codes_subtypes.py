# -*- coding: utf-8 -*-
"""A04 - case-definition analyses: per-code table, expanded cohort (main + 8246 G3/4), 8246-G3/4-only cohort,
small-cell vs large-cell subtypes, ICD-O-3.2-era sensitivity, Sun 2026 comparability.
Run: python3 python/A04_expanded_cohort_codes_subtypes.py   (calls the companion R/A04_expanded_cohort_codes_subtypes.R for the Poisson models)
"""
import sys, os, json, time, subprocess
import numpy as np, pandas as pd
t_start = time.time()
import config
import lm
RES = config.RES; FIGS = config.FIGS; SCR = config.R_DIR + os.sep
os.makedirs(RES, exist_ok=True); os.makedirs(FIGS, exist_ok=True)
RAW = config.RAW + config.RAW_FILES["NEN"]   # raw colorectal NEN export (read in place; counts only)
B = config.B(2000); SEED = 20260928; LM = np.arange(0, 25); WINDOW = 12; LM_H = 19   # hinge fitted on landmarks 0-18
KEY = [0, 3, 6, 12, 24]
out = {"id": "A04", "B": B, "seed": SEED, "window_months": WINDOW, "hinge_landmarks": "0-18 (grid 0-15 step 0.25)"}

# ---------------------------------------------------------------- data
nen = lm.load("NEN_all")
G3 = nen.grade17.isin(["Poorly differentiated; Grade III", "Undifferentiated; anaplastic; Grade IV"]) | nen.grade18.isin(["Site-specific grade system category (3)", "Site-specific grade system category (4)"])
nen = nen.assign(g3=G3.values)
nen["grp"] = np.where(nen["hist"].isin(lm.SCLC_CODES), "small", np.where(nen["hist"] == "8013", "large", np.where((nen["hist"] == "8246") & nen.g3, "nos8246_g34", "other")))
main = lm.nec_main(nen); expd = lm.nec_expanded(nen); nos = nen[nen.grp == "nos8246_g34"]
D = lambda d: d[d.stage == "Distant"]
md, ed, nd = D(main), D(expd), D(nos)
if not (len(main) == 1417 and len(md) == 908 and len(expd) == 3179 and len(ed) == 1975 and len(nos) == 1762): config.anchor_fail(f"cohort sizes {(len(main), len(md), len(expd), len(ed), len(nos))}")

# ---------------------------------------------------------------- reproduction check (main distant)
k, ev, KMAX = lm.prep(md)
cif, allc, nrisk = lm.curves(k, ev, KMAX, LM, WINDOW)
plateau_main = float(lm.hinge_many(cif[None, :LM_H])[0])
cr = lm.crude_rates(md)
exp_cif = [72.29, 71.62, 72.33, 60.03, 30.51]; exp_n = [908, 553, 396, 169, 48]
exp_rate = [0.1435, 0.0909, 0.1344, 0.1204, 0.1231, 0.0840, 0.0646]
got_cif = [round(float(cif[LM == s][0]), 2) for s in KEY]; got_n = [int(nrisk[LM == s][0]) for s in KEY]
got_rate = [round(float(r), 4) for r in cr.rate]
chk_ok = (np.allclose(got_cif, exp_cif, atol=0.011) and got_n == exp_n and abs(plateau_main - 6.75) < 1e-9 and np.allclose(got_rate, exp_rate, atol=0.00011)
          and len(md[md.year >= 2000]) == 908 and int(lm.surgery(md).eq(1).sum()) == 354 and int(lm.surgery(md).eq(0).sum()) == 554
          and int(md.flagc.sum()) == 881 and int(main.unkcod.sum()) == 10 and int(md.unkcod.sum()) == 7)
print("REPRODUCTION CHECK:", got_cif, got_n, plateau_main, got_rate, "->", "OK" if chk_ok else "MISMATCH")
if not chk_ok: config.anchor_fail("Reproduction check failed - stopping.")
out["reproduction_check"] = {"ok": bool(chk_ok), "cif_0_3_6_12_24": got_cif, "expected_cif": exp_cif, "n_at_risk": got_n, "expected_n": exp_n, "plateau_end": plateau_main, "expected_plateau": 6.75,
                            "crude_rates": got_rate, "expected_rates": exp_rate, "surgery_yes_no": [354, 554], "complete_date_subset": 881, "unknown_cod_all_distant": [10, 7]}

# ---------------------------------------------------------------- (1) per-code table
rows = []
def code_row(label, d):
    return {"code": label, "n": int(len(d)), "n_distant": int((d.stage == "Distant").sum()), "pct_distant": round(100 * (d.stage == "Distant").mean(), 1) if len(d) else None,
            "pct_2016_2023": round(100 * (d.year >= 2016).mean(), 1) if len(d) else None, "pct_2021_2023": round(100 * (d.year >= 2021).mean(), 1) if len(d) else None,
            "median_age": float(d.age.median()) if len(d) else None, "pct_male": round(100 * d.male.mean(), 1) if len(d) else None,
            "pct_rectum": round(100 * d.rectum.mean(), 1) if len(d) else None}
for c in ["8041", "8042", "8043", "8044", "8045", "8013", "8240", "8241", "8242", "8243", "8244", "8245", "8246", "8249"]:
    rows.append(code_row(c, nen[nen["hist"] == c]))
rows.append(code_row("main NEC (8013/8041-8045)", main)); rows.append(code_row("small-cell family (8041-8045)", main[main.grp == "small"]))
rows.append(code_row("8246 with grade 3/4", nos)); rows.append(code_row("expanded NEC (main + 8246 G3/4)", expd)); rows.append(code_row("all NEN codes", nen))
out["per_code_table"] = rows
g = nen[nen["hist"] == "8246"]
pre = g[g.year <= 2017]; post = g[g.year >= 2018]
out["grade_8246"] = {"n_8246": int(len(g)), "n_2000_2017": int(len(pre)), "n_2018_2023": int(len(post)),
                     "grade17_2000_2017": {str(a): int(b) for a, b in pre.grade17.value_counts().items()},
                     "grade18_2018_2023": {str(a): int(b) for a, b in post.grade18.value_counts().items()},
                     "grade17_distant_2000_2017": {str(a): int(b) for a, b in D(pre).grade17.value_counts().items()},
                     "grade18_distant_2018_2023": {str(a): int(b) for a, b in D(post).grade18.value_counts().items()},
                     "n_g34_total": int(g.g3.sum()), "n_g34_distant": int(D(g).g3.sum()),
                     "pct_g34_among_graded_2000_2017": round(100 * pre.g3.sum() / (~pre.grade17.isin(["Unknown", "Blank(s)"])).sum(), 1),
                     "pct_g34_among_graded_2018_2023": round(100 * post.g3.sum() / (~post.grade18.isin(["Grade cannot be assessed; Unknown", "Blank(s)"])).sum(), 1)}
# grade distribution of the main NEC codes for comparison (how often coded NEC carries grade 3/4)
mm_pre = main[main.year <= 2017]; mm_post = main[main.year >= 2018]
out["grade_mainNEC"] = {"grade17_2000_2017": {str(a): int(b) for a, b in mm_pre.grade17.value_counts().items()}, "grade18_2018_2023": {str(a): int(b) for a, b in mm_post.grade18.value_counts().items()},
                        "pct_g34_among_graded_2000_2017": round(100 * mm_pre.g3.sum() / (~mm_pre.grade17.isin(["Unknown", "Blank(s)"])).sum(), 1),
                        "pct_g34_among_graded_2018_2023": round(100 * mm_post.g3.sum() / (~mm_post.grade18.isin(["Grade cannot be assessed; Unknown", "Blank(s)"])).sum(), 1),
                        "pct_grade_unknown_2000_2017": round(100 * mm_pre.grade17.isin(["Unknown", "Blank(s)"]).mean(), 1)}

# ---------------------------------------------------------------- bootstrap engine
def boot_summary(df, label, seed=SEED, Bn=B):
    k, ev, KMAX = lm.prep(df)
    cif, allc, nrisk = lm.curves(k, ev, KMAX, LM, WINDOW)
    bc, ba = lm.bootstrap(df, B=Bn, seed=seed, landmarks=LM, window=WINDOW)
    ok = np.all(np.isfinite(bc[:, :LM_H]), axis=1)
    pl = np.full(Bn, np.nan); pl[ok] = lm.hinge_many(bc[ok][:, :LM_H])
    plateau = float(lm.hinge_many(cif[None, :LM_H])[0]) if np.all(np.isfinite(cif[:LM_H])) else None
    def ci(v): return [float(x) for x in np.nanpercentile(v, [2.5, 97.5])] if np.isfinite(v).sum() > 0 else [None, None]
    d60 = bc[:, 6] - bc[:, 0]; d120 = bc[:, 12] - bc[:, 0]; d240 = bc[:, 24] - bc[:, 0]
    r = {"label": label, "n": int(len(df)), "deaths_cancer": int((df.ev == 1).sum()), "deaths_other": int((df.ev == 2).sum()), "seed": seed, "B": Bn,
         "landmarks": LM.tolist(), "n_at_risk": nrisk.astype(int).tolist(),
         "cif12": [None if not np.isfinite(x) else round(float(x), 2) for x in cif], "cif12_lo": [], "cif12_hi": [],
         "allcause12": [None if not np.isfinite(x) else round(float(x), 2) for x in allc],
         "plateau_end": plateau, "plateau_ci": ci(pl), "plateau_valid_replicates": int(ok.sum()),
         "diff_6_0": round(float(cif[6] - cif[0]), 2), "diff_6_0_ci": ci(d60), "diff_12_0": round(float(cif[12] - cif[0]), 2), "diff_12_0_ci": ci(d120),
         "diff_24_0": None if not np.isfinite(cif[24]) else round(float(cif[24] - cif[0]), 2), "diff_24_0_ci": ci(d240),
         "key": {str(s): {"cif12": round(float(cif[s]), 2) if np.isfinite(cif[s]) else None, "ci": ci(bc[:, s]), "n_at_risk": int(nrisk[s])} for s in KEY}}
    lo, hi = np.nanpercentile(bc, [2.5, 97.5], axis=0)
    r["cif12_lo"] = [None if not np.isfinite(x) else round(float(x), 2) for x in lo]; r["cif12_hi"] = [None if not np.isfinite(x) else round(float(x), 2) for x in hi]
    r["_bc"] = bc; r["_pl"] = pl; r["_d60"] = d60; r["_d120"] = d120
    r["crude_rates"] = [{kk: (round(float(vv), 5) if isinstance(vv, float) else vv) for kk, vv in row.items()} for row in lm.crude_rates(df).to_dict("records")]
    print(f"{label:34s} n={len(df):5d}  cif 0/3/6/12/24 = {[r['key'][str(s)]['cif12'] for s in KEY]}  plateau={plateau} {r['plateau_ci']}  d6-0={r['diff_6_0']} {np.round(r['diff_6_0_ci'],2)}  [{time.time()-t_start:.0f}s]")
    return r
def strip(r): return {kk: vv for kk, vv in r.items() if not kk.startswith("_")}

# ---------------------------------------------------------------- (2) expanded and 8246-only cohorts (distant)
res_main = boot_summary(md, "main NEC distant")
res_exp = boot_summary(ed, "expanded NEC distant")
res_nos = boot_summary(nd, "8246 G3/4 only distant")
# monotone n at risk check
for r in (res_main, res_exp, res_nos): assert all(np.diff(r["n_at_risk"]) <= 0), "n at risk not monotone"
def pair_contrast(a, b, name):
    d6 = a["_d60"] - b["_d60"]; d12 = a["_d120"] - b["_d120"]; pl = a["_pl"] - b["_pl"]
    c0 = a["_bc"][:, 0] - b["_bc"][:, 0]; c6 = a["_bc"][:, 6] - b["_bc"][:, 6]
    return {"contrast": name, "diff6_0_difference": round(a["diff_6_0"] - b["diff_6_0"], 2), "diff6_0_difference_ci": [float(x) for x in np.percentile(d6, [2.5, 97.5])],
            "diff12_0_difference": round(a["diff_12_0"] - b["diff_12_0"], 2), "diff12_0_difference_ci": [float(x) for x in np.percentile(d12, [2.5, 97.5])],
            "plateau_difference": None if a["plateau_end"] is None or b["plateau_end"] is None else round(a["plateau_end"] - b["plateau_end"], 2), "plateau_difference_ci": [float(x) for x in np.nanpercentile(pl, [2.5, 97.5])],
            "cif0_difference": round(res_diff(a, b, 0), 2), "cif0_difference_ci": [float(x) for x in np.percentile(c0, [2.5, 97.5])],
            "cif6_difference": round(res_diff(a, b, 6), 2), "cif6_difference_ci": [float(x) for x in np.percentile(c6, [2.5, 97.5])],
            "note": "independent bootstraps (same seed, patient-level resampling within each cohort) paired by replicate index"}
def res_diff(a, b, s): return a["key"][str(s)]["cif12"] - b["key"][str(s)]["cif12"]
out["expanded"] = {"main": strip(res_main), "expanded": strip(res_exp), "nos8246_g34": strip(res_nos),
                   "contrasts": [pair_contrast(res_exp, res_main, "expanded minus main"), pair_contrast(res_nos, res_main, "8246 G3/4-only minus main")]}
# composition of expanded distant cohort
out["expanded"]["composition_distant"] = {str(a): int(b) for a, b in ed.grp.value_counts().items()}
out["expanded"]["baseline_distant"] = {lab: {"n": int(len(d)), "median_age": float(d.age.median()), "pct_male": round(100 * d.male.mean(), 1), "pct_rectum": round(100 * d.rectum.mean(), 1),
                                             "pct_2016_2023": round(100 * (d.year >= 2016).mean(), 1), "pct_chemo": round(100 * d.chemo.mean(), 1), "pct_surgery": round(100 * lm.surgery(d).eq(1).mean(), 1),
                                             "median_fu_months_alive": float(d[d.ev == 0].tm.median()) if (d.ev == 0).any() else None, "pct_cancer_death": round(100 * (d.ev == 1).mean(), 1)}
                                       for lab, d in [("main", md), ("nos8246_g34", nd), ("expanded", ed)]}

# ---------------------------------------------------------------- (3) subtypes within distant main NEC
sc = md[md.grp == "small"]; lc = md[md.grp == "large"]
if not (len(sc) == 599 and len(lc) == 309): config.anchor_fail(f"subtype sizes {len(sc)}/{len(lc)}")
res_sc = boot_summary(sc, "small-cell family distant")
res_lc = boot_summary(lc, "large-cell 8013 distant")
out["subtypes"] = {"small_cell": strip(res_sc), "large_cell": strip(res_lc), "contrast": pair_contrast(res_sc, res_lc, "small-cell minus large-cell"),
                   "manuscript": {"small_plateau": "7.25 (4.50-10.00)", "small_diff6_0": 1.9, "large_diff6_0": "-3.8 (-11.1 to 3.9)", "difference": "+5.7 (-3.0 to 14.5)"}}

# ---------------------------------------------------------------- (4) ICD-O-3.2 era sensitivity
era = {}
for lab, lo in [("2021_2023", 2021), ("2018_2023", 2018)]:
    for cl, d in [("main", md), ("expanded", ed), ("nos8246_g34", nd)]:
        sub = d[d.year >= lo]
        r = boot_summary(sub, f"{cl} distant {lab}")
        era[f"{cl}_{lab}"] = strip(r)
        era[f"{cl}_{lab}"]["n_all_stages"] = int(((main if cl == "main" else expd if cl == "expanded" else nos).year >= lo).sum())
        era[f"{cl}_{lab}"]["max_potential_followup_months"] = {"2021_2023": 36, "2018_2023": 72}[lab]
out["era"] = era
out["era"]["note"] = "Follow-up ends 31 Dec 2023 (fu_year max 2023); diagnoses in 2023 contribute <= 12 months, so landmark-6 estimates for 2021-2023 rest on heavy censoring within the 12-month window; the Aalen-Johansen estimator uses the at-risk set at each time, n at risk is reported."

# ---------------------------------------------------------------- (5) Sun 2026 comparability
sun = {"sun2026": {"small_cell": 790, "large_cell": 498, "years": "2000-2022", "sites": "colorectal incl. appendix"}}
m22 = main[main.year <= 2022]
sun["ours_2000_2022_after_all_filters"] = {"small_cell": int((m22.grp == "small").sum()), "large_cell": int((m22.grp == "large").sum()), "total": int(len(m22))}
sun["ours_2000_2023_after_all_filters"] = {"small_cell": int((main.grp == "small").sum()), "large_cell": int((main.grp == "large").sum()), "total": int(len(main))}
sun["ours_2023_only"] = {"small_cell": int(((main.year == 2023) & (main.grp == "small")).sum()), "large_cell": int(((main.year == 2023) & (main.grp == "large")).sum())}
try:
    cols = ["Histologic Type ICD-O-3", "Behavior recode for analysis", "Age recode with single ages and 90+", "Primary Site - labeled", "Diagnostic Confirmation",
            "Primary by international rules", "Type of Reporting Source", "Year of diagnosis", "Survival Days", "Patient ID", "Record number recode"]
    raw = pd.read_csv(RAW, sep="\t", dtype=str, usecols=cols, keep_default_na=False)
    raw["hist"] = raw["Histologic Type ICD-O-3"].str.strip()
    r = raw[raw["hist"].isin(lm.NEC_CODES)].copy()
    r["sub"] = np.where(r["hist"] == "8013", "large", "small"); r["year"] = r["Year of diagnosis"].astype(int); r["site5"] = r["Primary Site - labeled"].str[:5]
    r["age"] = pd.to_numeric(r["Age recode with single ages and 90+"].str.extract(r"^(\d+)")[0], errors="coerce")
    r["row"] = np.arange(len(r)); r["recno"] = pd.to_numeric(r["Record number recode"], errors="coerce")
    flow = []
    def st(name, x, note=""):
        flow.append({"step": name, "records": int(len(x)), "small_cell": int((x["sub"] == "small").sum()), "large_cell": int((x["sub"] == "large").sum()), "note": note})
    st("NEC codes 8013/8041-8045 in export (all sites C18-C20 incl. appendix, 2000-2023)", r)
    r = r[r.year <= 2022]; st("restrict to 2000-2022 (Sun 2026 window)", r)
    st("  of which appendix C18.1", r[r.site5 == "C18.1"], "included by Sun 2026, excluded here")
    r = r[r["Behavior recode for analysis"] == "Malignant"]; st("malignant behaviour", r)
    r = r[r.age >= 18]; st("age >= 18", r)
    r = r[r.site5 != "C18.1"]; st("exclude appendix (C18.1)", r)
    r = r[r["Diagnostic Confirmation"] == "Positive histology"]; st("positive histology only", r)
    r = r[r["Primary by international rules"] == "Yes"]; st("primary by international rules", r)
    r = r[~r["Type of Reporting Source"].isin(["Autopsy only", "Death certificate only"])]; st("not autopsy/DCO", r)
    r["sdays"] = pd.to_numeric(r["Survival Days"], errors="coerce"); r = r[r.sdays.notna() & (r.sdays >= 0)]; st("survival days known", r)
    r = r.sort_values(["Patient ID", "year", "recno", "row"]); r = r[~r.duplicated("Patient ID")]; st("one record per patient (first colorectal NEC record)", r)
    sun["flow_2000_2022_from_raw_export"] = flow
    # also the per-patient count when appendix kept (closest to Sun)
    r2 = raw[raw["hist"].isin(lm.NEC_CODES)].copy(); r2["sub"] = np.where(r2["hist"] == "8013", "large", "small"); r2["year"] = r2["Year of diagnosis"].astype(int)
    r2 = r2[(r2.year <= 2022)]; r2 = r2[~r2.duplicated("Patient ID")]
    sun["ours_2000_2022_appendix_included_no_other_filters_one_per_patient"] = {"small_cell": int((r2["sub"] == "small").sum()), "large_cell": int((r2["sub"] == "large").sum()), "total": int(len(r2))}
    # appendix kept but all other filters applied (closest reconstruction of the Sun 2026 definition)
    r3 = raw[raw["hist"].isin(lm.NEC_CODES)].copy(); r3["sub"] = np.where(r3["hist"] == "8013", "large", "small"); r3["year"] = r3["Year of diagnosis"].astype(int)
    r3["age"] = pd.to_numeric(r3["Age recode with single ages and 90+"].str.extract(r"^(\d+)")[0], errors="coerce"); r3["row"] = np.arange(len(r3)); r3["recno"] = pd.to_numeric(r3["Record number recode"], errors="coerce")
    r3 = r3[(r3.year <= 2022) & (r3["Behavior recode for analysis"] == "Malignant") & (r3.age >= 18) & (r3["Diagnostic Confirmation"] == "Positive histology") & (r3["Primary by international rules"] == "Yes")
            & ~r3["Type of Reporting Source"].isin(["Autopsy only", "Death certificate only"]) & pd.to_numeric(r3["Survival Days"], errors="coerce").ge(0)]
    r3 = r3.sort_values(["Patient ID", "year", "recno", "row"]); r3 = r3[~r3.duplicated("Patient ID")]
    sun["ours_2000_2022_appendix_included_all_other_filters"] = {"small_cell": int((r3["sub"] == "small").sum()), "large_cell": int((r3["sub"] == "large").sum()), "total": int(len(r3))}
    sun["residual_vs_sun2026_appendix_included"] = {"small_cell": int((r3["sub"] == "small").sum()) - 790, "large_cell": int((r3["sub"] == "large").sum()) - 498}
    sun["source"] = "raw SEER export (confidential, read in place; counts only)"
except Exception as e:
    sun["flow_from_raw_export_error"] = str(e)
    fl = pd.read_csv(config.DATA + "flow_all.csv"); fl = fl[fl.cohort == "NEN_all"]
    sun["flow_all_NEN_codes_from_flow_file"] = fl.to_dict("records"); sun["source"] = "flow_all.csv (all NEN codes, not NEC-specific)"
fl = pd.read_csv(config.DATA + "flow_all.csv"); fl = fl[fl.cohort == "NEN_all"].reset_index(drop=True)
fl["records_removed"] = -fl.records.diff().fillna(0).astype(int)
sun["flow_all_NEN_codes_from_flow_file"] = fl.to_dict("records")
out["sun2026_comparability"] = sun

# ---------------------------------------------------------------- R part (adjusted RRs)
rr_path = RES + "A04_rr.json"
rc = subprocess.run(["Rscript", os.path.join(SCR, "A04_expanded_cohort_codes_subtypes.R")], capture_output=True, text=True)
print(rc.stdout[-3000:]);
if rc.returncode != 0: print(rc.stderr[-3000:]); sys.exit("Rscript failed")
out["adjusted_RR"] = json.load(open(rr_path))

# ---------------------------------------------------------------- figures
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
C = {"main": "#2a78d6", "expanded": "#eb6834", "nos": "#1baf7a", "small": "#2a78d6", "large": "#eb6834", "e18": "#eda100", "e21": "#4a3aa7"}
def band(ax, r, col, lab, lmax=24):
    x = np.array(r["landmarks"][:lmax + 1]); y = np.array([np.nan if v is None else v for v in r["cif12"][:lmax + 1]], float)
    lo = np.array([np.nan if v is None else v for v in r["cif12_lo"][:lmax + 1]], float); hi = np.array([np.nan if v is None else v for v in r["cif12_hi"][:lmax + 1]], float)
    ax.fill_between(x, lo, hi, color=col, alpha=0.13, linewidth=0); ax.plot(x, y, color=col, lw=2, label=f"{lab} (n={r['n']})")
    if r["plateau_end"] is not None: ax.axvline(r["plateau_end"], color=col, lw=1, ls=":", alpha=0.8)
fig, axes = plt.subplots(1, 2, figsize=(11, 4.4), dpi=300, sharey=True)
ax = axes[0]; band(ax, out["expanded"]["main"], C["main"], "main NEC"); band(ax, out["expanded"]["expanded"], C["expanded"], "expanded (+8246 G3/4)"); band(ax, out["expanded"]["nos8246_g34"], C["nos"], "8246 G3/4 only")
ax.set_title("a  Case definition (distant stage)", loc="left", fontsize=10); ax.set_xlabel("Landmark (months since diagnosis)"); ax.set_ylabel("12-month cancer-death probability (%)")
ax = axes[1]; band(ax, out["expanded"]["main"], C["main"], "main, 2000-2023"); band(ax, out["era"]["main_2018_2023"], C["e18"], "main, 2018-2023", 18); band(ax, out["era"]["main_2021_2023"], C["e21"], "main, 2021-2023", 12)
ax.set_title("b  Diagnosis era (main NEC, distant)", loc="left", fontsize=10); ax.set_xlabel("Landmark (months since diagnosis)")
for ax in axes:
    ax.set_xlim(0, 24); ax.set_ylim(0, 100); ax.grid(color="#e6e6e3", lw=0.6); ax.spines[["top", "right"]].set_visible(False); ax.legend(frameon=False, fontsize=8, loc="lower left")
fig.text(0.01, 0.005, "Shaded: 95% percentile bootstrap (B=2000). Dotted vertical lines: broken-stick plateau end.", fontsize=7, color="#52514e")
fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(FIGS + "A04_definition_era_curves.png"); plt.close(fig)
fig, ax = plt.subplots(figsize=(6, 4.4), dpi=300)
band(ax, out["subtypes"]["small_cell"], C["small"], "small-cell family 8041-8045", 18); band(ax, out["subtypes"]["large_cell"], C["large"], "large-cell NEC 8013", 18)
ax.set_xlim(0, 18); ax.set_ylim(0, 100); ax.grid(color="#e6e6e3", lw=0.6); ax.spines[["top", "right"]].set_visible(False); ax.legend(frameon=False, fontsize=8, loc="lower left")
ax.set_xlabel("Landmark (months since diagnosis)"); ax.set_ylabel("12-month cancer-death probability (%)"); ax.set_title("Subtypes within distant main NEC", loc="left", fontsize=10)
fig.tight_layout(); fig.savefig(FIGS + "A04_subtype_curves.png"); plt.close(fig)
out["figures"] = [config.relpath(FIGS + "A04_definition_era_curves.png"), config.relpath(FIGS + "A04_subtype_curves.png")]

# ---------------------------------------------------------------- JSON
out["runtime_seconds"] = round(time.time() - t_start, 1)
json.dump(out, open(RES + "A04.json", "w"), indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))

# ---------------------------------------------------------------- Markdown
INTERP = """
1. Per-code composition. The main cohort is essentially two codes: 8041 (small cell carcinoma NOS; 844 of 871 small-cell-family cases) and 8013 (large cell NEC; 546). Codes 8042-8044 contribute 3 cases and 8045 (combined small cell) 24. The NOS code 8246/3 is the largest potential NEC reservoir (4,619 cases, 1,791 distant) but is heterogeneous: among graded cases only 50% (2000-2017) and 39% (2018-2023) are grade 3/4, against 96% for the coded NEC. The 8246 G3/4 cases are older-era (29% diagnosed 2016-2023 vs 51% for the main cohort; 8% vs 23% in 2021-2023), consistent with registrars moving to specific codes (8041/8013) and to 8249 (NET G3, 82% diagnosed 2016+) after WHO 2017/2019 and ICD-O-3.2.
2. Expanded cohort (n = 1,975 distant). The landmark-0 level is the same as in the main cohort (71.8% vs 72.3%; paired difference -0.4 pp, 95% CI -4.1 to 3.2), and the broken-stick plateau end is similar (7.25 months, 95% CI 3.25-9.75 vs 6.75, 4.00-9.75). The 6-0 difference is -1.8 pp (-4.6 to 1.1), i.e. still compatible with no decline, and the adjusted RR for months 6-9 vs 3-6 remains above 1 (1.30, 1.10-1.53, P = 0.002) but is attenuated relative to the main cohort (1.49, 1.17-1.90). The expanded-minus-main contrast of the 6-0 change is -1.8 pp (-6.7 to 3.0): the two definitions do not differ detectably.
3. 8246 G3/4 alone (n = 1,067 distant). These cases share the early lethality of coded NEC (landmark 0: 71.5%, 68.7-74.1; months 0-3 crude rate 0.152 vs 0.144 per person-month) but show a gentle monotone decline instead of a plateau-and-rebound: 6-0 change -3.3 pp (-7.1 to 0.4), RR 6-9 vs 3-6 1.15 (0.92-1.45, P = 0.21), RR 6-15 vs 0-6 0.87 (0.75-1.01), and an uninformative plateau interval (7.50, 0.00-13.75). They also carry more primary-site surgery (48% vs 39%) and a higher 12-month risk at the 24-month landmark (40% vs 31%), a profile consistent with a mixture of true NEC and less aggressive tumours (NET G3, poorly differentiated adenocarcinoma with neuroendocrine features). The direction matches the authors' own trial analysis (-3.2 pp all-cause). Conclusion: 8246 G3/4 cases behave 'NEC-like' in level but only partly in shape; the pre-specified expanded cohort preserves the main conclusion (no fall in 12-month cancer-death risk through 6 months, higher rate in months 6-9 than 3-6) with smaller effect sizes, and the manuscript should report both rather than treat the expanded cohort as a silent sensitivity analysis.
4. Subtypes. Within the main cohort the plateau-and-rebound is a small-cell-family feature: 6-0 change +1.9 pp (-3.2 to 6.5), plateau end 7.25 months (4.25-10.00), RR 6-9 vs 3-6 1.75 (1.29-2.38). Large-cell NEC (n = 309) shows -3.8 pp (-10.9 to 3.7), RR 1.18 (0.79-1.78) and no identifiable plateau (6.25, 0.00-15.00). The paired between-subtype difference is +5.7 pp (-3.7 to 14.3), so the two cannot be separated statistically; the manuscript values (7.25 [4.50-10.00], -3.8 [-11.1 to 3.9], +5.7 [-3.0 to 14.5]) are reproduced to within bootstrap noise.
5. ICD-O-3.2 era. Restricting to diagnoses 2018-2023 (main n = 406 distant) or 2021-2023 (n = 225) leaves the landmark-0 level at 68-71% and the 6-0 change at -0.9 pp (-7.4 to 5.3) and +2.2 pp (-7.2 to 11.2); the adjusted RR 6-9 vs 3-6 is 1.80 (1.25-2.59) for 2018-2023 (expanded: 1.58, 1.18-2.12). The pattern therefore is not an artefact of pre-ICD-O-3.2 coding; the 2021-2023 estimates beyond landmark 6 rest on fewer than 50 patients at risk and are reported for completeness only.
6. Sun 2026 comparability. Restricted to 2000-2022 our definition yields 805 small-cell + 493 large-cell (1,298) vs Sun's 790 + 498 (1,288). Of the 1,355 NEC-coded records diagnosed 2000-2022, our filters remove 23 appendiceal (9 small, 14 large), 17 without positive histology, 14 not primary by international rules and 3 duplicate patient records; none are removed by behaviour, age, reporting source or survival-days rules. Keeping the appendix (as Sun did) gives 814 + 507. The residual difference (about +24 small-cell and +9 large-cell relative to Sun, under 3%) cannot be produced by any of our filters and is most plausibly due to data vintage (our export is the 2000-2023 release with updated histology recodes; Sun used the 2000-2022 submission) and unreported inclusion rules in Sun 2026 (e.g. microscopic confirmation handling, exclusion of unknown stage or of 8042-8045).
Deviations from the plan: the same seed (20260928) was used for every cohort-specific bootstrap and contrasts were taken replicate-by-replicate, so CIs for the subtype contrasts differ from the manuscript by bootstrap noise only; the Sun 2026 flow was computed NEC-specifically from the raw SEER export read in place (counts only; no patient-level data written), because flow_all.csv covers all NEN codes; 2018+ 'Site-specific grade system category (4)' was treated as grade 4 per the pre-specified rule (24 of 4,619 8246 cases); adjusted RR models were additionally run for the 2018-2023 subsets (not requested); B = 2000 throughout (runtime under one minute).
"""
def f1(x): return "NA" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.1f}"
def f2(x): return "NA" if x is None or (isinstance(x, float) and not np.isfinite(x)) else f"{x:.2f}"
def ci_s(c, f=f1): return "NA" if c[0] is None else f"{f(c[0])} to {f(c[1])}"
def rr_s(d): return "NA" if d is None or "RR" not in d else f"{d['RR']:.2f} ({d['lo']:.2f}-{d['hi']:.2f}), P={d['p']:.3g}"
L = []
L.append("# A04 - Case-definition analyses: per-code table, expanded cohort (main + 8246 G3/4), subtypes, ICD-O-3.2 era, Sun 2026 comparability\n")
L.append(f"Generated {time.strftime('%Y-%m-%d %H:%M')}; B = {B} patient-level bootstrap replicates (seed {SEED}, same resample for all landmarks within a cohort; independent bootstraps paired by index for contrasts); 12-month window; broken-stick plateau on landmarks 0-18. Runtime {out['runtime_seconds']} s.\n")
L.append("## Reproduction check (main NEC, distant stage, n = 908)\n")
L.append("| Quantity | Manuscript | This script |\n|---|---|---|")
L.append(f"| 12-month cancer-death % at landmarks 0/3/6/12/24 | 72.29/71.62/72.33/60.03/30.51 | {'/'.join(f'{v:.2f}' for v in got_cif)} |")
L.append(f"| n at risk | 908/553/396/169/48 | {'/'.join(str(v) for v in got_n)} |")
L.append(f"| Plateau end (broken stick) | 6.75 | {plateau_main:.2f} |")
L.append(f"| Crude cancer-death rate /person-month 0-3 ... 18-24 | 0.1435/0.0909/0.1344/0.1204/0.1231/0.0840/0.0646 | {'/'.join(f'{v:.4f}' for v in got_rate)} |")
L.append(f"| Adjusted RR 6-9 vs 3-6 (Poisson, cluster-robust) | 1.49 (1.17-1.90), P=0.0013 | {out['adjusted_RR']['reproduction_check']['observed']} |")
L.append(f"| Surgery yes/no; complete-date subset; unknown COD (all/distant) | 354/554; 881; 10/7 | {int(lm.surgery(md).eq(1).sum())}/{int(lm.surgery(md).eq(0).sum())}; {int(md.flagc.sum())}; {int(main.unkcod.sum())}/{int(md.unkcod.sum())} |")
L.append("\nAll anchors reproduced; the script stops if they do not.\n")
L.append("## 1. Per-code table (all colorectal NEN codes after the common filters; n = 32,620)\n")
L.append("| ICD-O-3 code | n | n distant (%) | % diagnosed 2016-2023 | % diagnosed 2021-2023 | Median age | % male | % rectum/rectosigmoid |\n|---|---|---|---|---|---|---|---|")
for r in rows:
    L.append(f"| {r['code']} | {r['n']:,} | {r['n_distant']:,} ({f1(r['pct_distant'])}) | {f1(r['pct_2016_2023'])} | {f1(r['pct_2021_2023'])} | {f1(r['median_age'])} | {f1(r['pct_male'])} | {f1(r['pct_rectum'])} |")
gg = out["grade_8246"]
L.append(f"\n**Grade distribution of 8246/3 (NEC NOS), n = {gg['n_8246']:,}.** Diagnoses 2000-2017 (Grade Recode thru 2017, n = {gg['n_2000_2017']:,}): " + "; ".join(f"{k}: {v:,}" for k, v in gg["grade17_2000_2017"].items()) + ".")
L.append(f"Diagnoses 2018-2023 (Derived Summary Grade 2018, n = {gg['n_2018_2023']:,}): " + "; ".join(f"{k}: {v:,}" for k, v in gg["grade18_2018_2023"].items()) + ".")
L.append(f"Grade 3/4 = {gg['n_g34_total']:,} ({gg['n_g34_distant']:,} distant); grade 3/4 share among graded cases {gg['pct_g34_among_graded_2000_2017']}% (2000-2017) and {gg['pct_g34_among_graded_2018_2023']}% (2018-2023). Among distant-stage 8246 cases, 2000-2017: " + "; ".join(f"{k}: {v}" for k, v in gg["grade17_distant_2000_2017"].items()) + "; 2018-2023: " + "; ".join(f"{k}: {v}" for k, v in gg["grade18_distant_2018_2023"].items()) + ".")
gm = out["grade_mainNEC"]
L.append(f"For comparison, coded NEC (8013/8041-8045) carries grade 3/4 in {gm['pct_g34_among_graded_2000_2017']}% (2000-2017) and {gm['pct_g34_among_graded_2018_2023']}% (2018-2023) of graded cases; grade is unknown/blank in {gm['pct_grade_unknown_2000_2017']}% of 2000-2017 main-cohort cases. Note: 2018+ 'Site-specific grade system category (4)' is used here as grade 4 per the pre-specified rule (it is rare: {gg['grade18_2018_2023'].get('Site-specific grade system category (4)', 0)} of 8246 cases).\n")
L.append("## 2. Expanded cohort vs main cohort (distant stage)\n")
L.append(f"Composition of expanded distant cohort (n = {res_exp['n']:,}): " + ", ".join(f"{k} = {v:,}" for k, v in out["expanded"]["composition_distant"].items()) + ".\n")
bl = out["expanded"]["baseline_distant"]
L.append("| Cohort (distant) | n | Median age | % male | % rectum | % 2016-2023 | % chemo | % primary surgery | % cancer death |\n|---|---|---|---|---|---|---|---|---|")
for lab, nm in [("main", "Main NEC"), ("nos8246_g34", "8246 G3/4 only"), ("expanded", "Expanded")]:
    b = bl[lab]; L.append(f"| {nm} | {b['n']:,} | {b['median_age']:.0f} | {b['pct_male']} | {b['pct_rectum']} | {b['pct_2016_2023']} | {b['pct_chemo']} | {b['pct_surgery']} | {b['pct_cancer_death']} |")
L.append("\n**Table 2a. Landmark 12-month cancer-death probability (%), 95% percentile bootstrap CI, n at risk.**\n")
L.append("| Landmark | Main NEC (n=908) | Expanded (n=1,975) | 8246 G3/4 only (n=1,067) |\n|---|---|---|---|")
for s in KEY:
    cells = [f"{r['key'][str(s)]['cif12']:.1f} ({ci_s(r['key'][str(s)]['ci'])}); n={r['key'][str(s)]['n_at_risk']}" for r in (res_main, res_exp, res_nos)]
    L.append(f"| {s} | " + " | ".join(cells) + " |")
L.append("\n**Table 2b. Plateau and landmark differences.**\n")
L.append("| Quantity | Main NEC | Expanded | 8246 G3/4 only |\n|---|---|---|---|")
L.append("| Plateau end, months (95% CI) | " + " | ".join(f"{f2(r['plateau_end'])} ({ci_s(r['plateau_ci'], f2)})" for r in (res_main, res_exp, res_nos)) + " |")
L.append("| 6-0 difference, pp (95% CI) | " + " | ".join(f"{r['diff_6_0']:+.1f} ({ci_s(r['diff_6_0_ci'])})" for r in (res_main, res_exp, res_nos)) + " |")
L.append("| 12-0 difference, pp (95% CI) | " + " | ".join(f"{r['diff_12_0']:+.1f} ({ci_s(r['diff_12_0_ci'])})" for r in (res_main, res_exp, res_nos)) + " |")
L.append("| 24-0 difference, pp (95% CI) | " + " | ".join(f"{r['diff_24_0']:+.1f} ({ci_s(r['diff_24_0_ci'])})" for r in (res_main, res_exp, res_nos)) + " |")
for c in out["expanded"]["contrasts"]:
    L.append(f"\nPaired contrast *{c['contrast']}*: 6-0 difference {c['diff6_0_difference']:+.1f} pp ({ci_s(c['diff6_0_difference_ci'])}); 12-0 difference {c['diff12_0_difference']:+.1f} pp ({ci_s(c['diff12_0_difference_ci'])}); landmark-0 level {c['cif0_difference']:+.1f} pp ({ci_s(c['cif0_difference_ci'])}); landmark-6 level {c['cif6_difference']:+.1f} pp ({ci_s(c['cif6_difference_ci'])}); plateau end {f2(c['plateau_difference'])} months ({ci_s(c['plateau_difference_ci'], f2)}).")
L.append("\n**Table 2c. Crude cancer-death rates per person-month (cause-specific; other deaths censor).**\n")
L.append("| Interval (months) | Main NEC deaths / PM / rate | Expanded deaths / PM / rate | 8246 G3/4 only deaths / PM / rate |\n|---|---|---|---|")
for i in range(7):
    cells = [f"{r['crude_rates'][i]['deaths']} / {r['crude_rates'][i]['person_months']:.0f} / {r['crude_rates'][i]['rate']:.4f}" for r in (res_main, res_exp, res_nos)]
    L.append(f"| {res_main['crude_rates'][i]['start']}-{res_main['crude_rates'][i]['end']} | " + " | ".join(cells) + " |")
A = out["adjusted_RR"]
L.append("\n**Table 2d. Adjusted rate ratios (piecewise Poisson, cluster-robust SE). Main: age10+age10sq+male+year10+large+rectum; expanded: age10+age10sq+male+year10+rectum+code group (small/large/8246 G3-4); 8246-only and subtypes: age10+age10sq+male+year10+rectum.**\n")
L.append("| Cohort (distant) | n | RR 6-9 vs 3-6 | RR 9-12 vs 3-6 | RR 3-6 vs 0-3 | RR 6-15 vs 0-6 | RR 15-24 vs 6-15 |\n|---|---|---|---|---|---|---|")
for key, nm in [("main", "Main NEC"), ("expanded", "Expanded"), ("nos8246_g34", "8246 G3/4 only"), ("small_cell", "Small-cell family"), ("large_cell", "Large-cell 8013"), ("main_2018_2023", "Main, 2018-2023"), ("expanded_2018_2023", "Expanded, 2018-2023")]:
    a = A[key]; m3 = a["model_3mo"]; mb = a["model_0_6_15_24"]
    L.append(f"| {nm} | {a['n']:,} | {rr_s(m3.get('rr_6_9_vs_3_6'))} | {rr_s(m3.get('rr_9_12_vs_3_6'))} | {rr_s(m3.get('rr_3_6_vs_0_3'))} | {rr_s(mb.get('rr_6_15_vs_0_6'))} | {rr_s(mb.get('rr_15_24_vs_6_15'))} |")
L.append("\n## 3. Subtypes within distant main NEC\n")
L.append("| Landmark | Small-cell family (n=599) | Large-cell 8013 (n=309) |\n|---|---|---|")
for s in range(0, 19, 3):
    L.append(f"| {s} | {res_sc['cif12'][s]:.1f} ({ci_s(res_sc['cif12_lo'][s:s+1] + res_sc['cif12_hi'][s:s+1])}); n={res_sc['n_at_risk'][s]} | {res_lc['cif12'][s]:.1f} ({ci_s(res_lc['cif12_lo'][s:s+1] + res_lc['cif12_hi'][s:s+1])}); n={res_lc['n_at_risk'][s]} |")
cs = out["subtypes"]["contrast"]
L.append("\n| Quantity | Small-cell | Large-cell | Difference (small minus large), paired bootstrap | Manuscript |\n|---|---|---|---|---|")
L.append(f"| Plateau end, months | {f2(res_sc['plateau_end'])} ({ci_s(res_sc['plateau_ci'], f2)}) | {f2(res_lc['plateau_end'])} ({ci_s(res_lc['plateau_ci'], f2)}) | {f2(cs['plateau_difference'])} ({ci_s(cs['plateau_difference_ci'], f2)}) | small 7.25 (4.50-10.00) |")
L.append(f"| 6-0 difference, pp | {res_sc['diff_6_0']:+.1f} ({ci_s(res_sc['diff_6_0_ci'])}) | {res_lc['diff_6_0']:+.1f} ({ci_s(res_lc['diff_6_0_ci'])}) | {cs['diff6_0_difference']:+.1f} ({ci_s(cs['diff6_0_difference_ci'])}) | +1.9; -3.8 (-11.1 to 3.9); +5.7 (-3.0 to 14.5) |")
L.append(f"| 12-0 difference, pp | {res_sc['diff_12_0']:+.1f} ({ci_s(res_sc['diff_12_0_ci'])}) | {res_lc['diff_12_0']:+.1f} ({ci_s(res_lc['diff_12_0_ci'])}) | {cs['diff12_0_difference']:+.1f} ({ci_s(cs['diff12_0_difference_ci'])}) | - |")
L.append(f"| RR 6-9 vs 3-6 (adjusted) | {rr_s(A['small_cell']['model_3mo'].get('rr_6_9_vs_3_6'))} | {rr_s(A['large_cell']['model_3mo'].get('rr_6_9_vs_3_6'))} | - | - |")
L.append(f"| RR 6-15 vs 0-6 (adjusted) | {rr_s(A['small_cell']['model_0_6_15_24'].get('rr_6_15_vs_0_6'))} | {rr_s(A['large_cell']['model_0_6_15_24'].get('rr_6_15_vs_0_6'))} | - | - |")
L.append("\n## 4. ICD-O-3.2 era sensitivity (distant stage)\n")
L.append("| Cohort, diagnosis years | n distant (all stages) | Landmark 0: 12-mo cancer death % (95% CI); n | Landmark 6 (95% CI); n | 6-0 difference, pp (95% CI) | Landmark 12 (95% CI); n | Plateau end (95% CI) |\n|---|---|---|---|---|---|---|")
for key, nm in [("main_2018_2023", "Main NEC, 2018-2023"), ("expanded_2018_2023", "Expanded, 2018-2023"), ("nos8246_g34_2018_2023", "8246 G3/4 only, 2018-2023"), ("main_2021_2023", "Main NEC, 2021-2023"), ("expanded_2021_2023", "Expanded, 2021-2023"), ("nos8246_g34_2021_2023", "8246 G3/4 only, 2021-2023")]:
    r = era[key]; k0, k6, k12 = r["key"]["0"], r["key"]["6"], r["key"]["12"]
    L.append(f"| {nm} | {r['n']:,} ({r['n_all_stages']:,}) | {f1(k0['cif12'])} ({ci_s(k0['ci'])}); {k0['n_at_risk']} | {f1(k6['cif12'])} ({ci_s(k6['ci'])}); {k6['n_at_risk']} | {r['diff_6_0']:+.1f} ({ci_s(r['diff_6_0_ci'])}) | {f1(k12['cif12'])} ({ci_s(k12['ci'])}); {k12['n_at_risk']} | {f2(r['plateau_end'])} ({ci_s(r['plateau_ci'], f2)}; valid replicates {r['plateau_valid_replicates']}) |")
L.append(f"\n{out['era']['note']}\n")
L.append("## 5. Sun 2026 (Cancers) comparability\n")
s5 = out["sun2026_comparability"]
L.append(f"Sun 2026: 790 small-cell + 498 large-cell = 1,288 (SEER, 2000-2022, colorectal including appendix). Ours after all filters: 2000-2023 {s5['ours_2000_2023_after_all_filters']['small_cell']} + {s5['ours_2000_2023_after_all_filters']['large_cell']} = {s5['ours_2000_2023_after_all_filters']['total']}; restricted to 2000-2022 {s5['ours_2000_2022_after_all_filters']['small_cell']} + {s5['ours_2000_2022_after_all_filters']['large_cell']} = {s5['ours_2000_2022_after_all_filters']['total']} (2023 alone adds {s5['ours_2023_only']['small_cell']} + {s5['ours_2023_only']['large_cell']}).\n")
if "flow_2000_2022_from_raw_export" in s5:
    L.append("NEC-specific flow (records with codes 8013/8041-8045 in the SEER 17 export, diagnoses 2000-2022):\n")
    L.append("| Step | Records | Small-cell | Large-cell | Note |\n|---|---|---|---|---|")
    for f in s5["flow_2000_2022_from_raw_export"]: L.append(f"| {f['step']} | {f['records']:,} | {f['small_cell']:,} | {f['large_cell']:,} | {f['note']} |")
    a = s5["ours_2000_2022_appendix_included_no_other_filters_one_per_patient"]
    b2 = s5["ours_2000_2022_appendix_included_all_other_filters"]; rs = s5["residual_vs_sun2026_appendix_included"]
    L.append(f"\nWith appendix kept and no other filter (one record per patient), 2000-2022: {a['small_cell']} small-cell + {a['large_cell']} large-cell = {a['total']}. With appendix kept and all other filters applied: {b2['small_cell']} + {b2['large_cell']} = {b2['total']}, i.e. {rs['small_cell']:+d} small-cell and {rs['large_cell']:+d} large-cell relative to Sun 2026 (790/498); no combination of our filters reproduces their split exactly.\n")
L.append("All-NEN-code flow from flow_all.csv (not code-specific): " + "; ".join(f"{f['step']} {f['records']:,} (-{f['records_removed']:,})" for f in s5["flow_all_NEN_codes_from_flow_file"]) + ".\n")
L.append("## Interpretation\n")
L.append(INTERP.strip())
open(RES + "A04.md", "w").write("\n".join(L))
print("wrote", RES + "A04.json", RES + "A04.md", "runtime", out["runtime_seconds"], "s")
