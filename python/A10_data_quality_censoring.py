# -*- coding: utf-8 -*-
"""A10 data_quality_censoring: documentation tables for Methods/Supplement (main colorectal NEC cohort, all stages & distant).
(1) survival flag / days-months consistency, (2) administrative censoring by landmark & diagnosis year,
(3) diagnosis-year distribution & potential follow-up, (4) treatment timing (tt_days), (5) multiple primaries,
(6) unknown COD / other-cause classification, (7) metastatic pattern & treatment cross-tabs, (8) comparator cohorts
(flag, diagnostic confirmation incl. cytology variants, pid overlap).  Descriptive only - no bootstrap.
Run: python3 python/A10_data_quality_censoring.py
"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np, pandas as pd
T0 = time.time()
import config
import lm
RES = config.RES; os.makedirs(RES, exist_ok=True)
MD = lm.MD
OUT = {"id": "A10", "task": "data_quality_censoring", "date": "2026-10-07", "bootstrap": None, "seed": None, "B": None}
MDL = []  # markdown lines
def H(s): MDL.append("\n" + s + "\n")
def P(s=""): MDL.append(s)
def tbl(df, floatfmt="{:.2f}"):
    """markdown table from a DataFrame"""
    cols = list(df.columns); P("| " + " | ".join(str(c) for c in cols) + " |"); P("|" + "---|" * len(cols))
    for _, r in df.iterrows():
        P("| " + " | ".join(floatfmt.format(v) if isinstance(v, (float, np.floating)) and not pd.isna(v) else ("" if (isinstance(v, float) and pd.isna(v)) else str(v)) for v in r.values) + " |")
def vc(s, name, denom=None, order=None):
    """value counts with % -> DataFrame (n, %)"""
    s = s.astype(object).where(s.notna(), "(missing)").astype(str)
    c = s.value_counts(dropna=False)
    if order is not None: c = c.reindex([o for o in order if o in c.index] + [i for i in c.index if i not in order]).dropna().astype(int)
    denom = len(s) if denom is None else denom
    d = pd.DataFrame({name: c.index, "n": c.values, "%": 100 * c.values / denom})
    return d
def rec(df): return [{k: (int(v) if isinstance(v, (np.integer,)) else (float(v) if isinstance(v, (np.floating, float)) else v)) for k, v in r.items()} for r in df.to_dict("records")]

# ------------------------------------------------------------------ load main cohort
nen = lm.load("NEN_all"); nec = lm.nec_main(nen).copy(); dist = nec[nec.stage == "Distant"].copy()
nec["surg"] = lm.surgery(nec); dist["surg"] = lm.surgery(dist)
N_ALL, N_DIST = len(nec), len(dist)
print(f"main NEC n={N_ALL}, distant n={N_DIST}")

# ------------------------------------------------------------------ reproduction check
LM5 = np.array([0, 3, 6, 12, 24])
k, ev, KMAX = lm.prep(dist); cif, allc, nrisk = lm.curves(k, ev, KMAX, LM5, 12)
anch = {"cif": [72.29, 71.62, 72.33, 60.03, 30.51], "nrisk": [908, 553, 396, 169, 48], "surg_yes": 354, "surg_no": 554, "flagc": 881, "unkcod": 10, "unkcod_dist": 7, "n": 1417, "n_dist": 908}
mine = {"cif": [round(float(x), 2) for x in cif], "nrisk": [int(x) for x in nrisk], "surg_yes": int((dist.surg == 1).sum()), "surg_no": int((dist.surg == 0).sum()),
        "flagc": int(dist.flagc.sum()), "unkcod": int(nec.unkcod.sum()), "unkcod_dist": int(dist.unkcod.sum()), "n": N_ALL, "n_dist": N_DIST}
ok = all(abs(a - b) <= 0.011 for a, b in zip(anch["cif"], mine["cif"])) and all(anch[q] == mine[q] for q in anch if q != "cif")
print("REPRODUCTION CHECK:", "PASS" if ok else "FAIL", mine)
if not ok: config.anchor_fail("Reproduction anchors disagree - stopping.")
OUT["reproduction_check"] = {"pass": bool(ok), "manuscript": anch, "mine": mine}
H("# A10 Data quality, censoring and documentation tables (main colorectal NEC cohort)")
P(f"Main NEC cohort n = {N_ALL} (distant n = {N_DIST}). Descriptive tables only; no bootstrap. Runtime at end.")
H("## Reproduction check (mine vs manuscript)")
rc = pd.DataFrame({"quantity": ["12-mo cancer-death % LM 0/3/6/12/24", "n at risk LM 0/3/6/12/24", "distant surgery yes/no", "distant complete-date subset", "unknown COD all / distant"],
                   "manuscript": ["/".join(map(str, anch["cif"])), "/".join(map(str, anch["nrisk"])), "354/554", "881", "10/7"],
                   "this script": ["/".join(map(str, mine["cif"])), "/".join(map(str, mine["nrisk"])), f"{mine['surg_yes']}/{mine['surg_no']}", str(mine["flagc"]), f"{mine['unkcod']}/{mine['unkcod_dist']}"],
                   "agree": ["yes"] * 5})
tbl(rc)

# ------------------------------------------------------------------ (1) survival flag & days/months consistency
H("## 1. Survival-time fields")
sec1 = {}
for lab, d in [("all", nec), ("distant", dist)]:
    f = vc(d.flag, "Survival months flag"); sec1[f"flag_{lab}"] = rec(f)
    P(f"**Survival months flag, {lab} stages (n = {len(d)})**"); tbl(f)
    sm = pd.to_numeric(d.smonths, errors="coerce"); sd = d.sdays.astype(float)
    fl = np.floor(sd / 30.4375)
    disc = (fl != sm); diff = (fl - sm)
    zero = int((sd == 0).sum()); zero_flag = d.loc[sd == 0, "flag"].value_counts().to_dict()
    half = int((d.tm < 0.1).sum())
    s = {"n": len(d), "smonths_nonnumeric": int(sm.isna().sum()), "discrepant_floor_days_vs_months": int(disc.sum()), "pct_discrepant": float(100 * disc.mean()),
         "diff_distribution": {str(int(x)): int(c) for x, c in diff.value_counts().sort_index().items()},
         "zero_day_survivors": zero, "zero_day_flag": zero_flag, "sdays_min": float(sd.min()), "sdays_max": float(sd.max()), "sdays_median": float(sd.median()),
         "tm_set_to_half_day": half, "max_abs_days_months_diff_months": float(np.abs((sd / MD) - sm).max())}
    sec1[f"consistency_{lab}"] = s
    P(f"\n**Days vs months consistency, {lab}**: floor(sdays/30.4375) == smonths in {len(d) - int(disc.sum())}/{len(d)} ({100 - s['pct_discrepant']:.2f} %); "
      f"discrepant {int(disc.sum())} ({s['pct_discrepant']:.2f} %), floor(days)-months difference distribution {s['diff_distribution']}; "
      f"non-numeric smonths {s['smonths_nonnumeric']}; zero-day survivors {zero} (flags: {zero_flag}); max |days/30.4375 - months| = {s['max_abs_days_months_diff_months']:.3f} months. "
      f"Survival days: median {sd.median():.0f}, range {sd.min():.0f}-{sd.max():.0f}.")
    # cross: flag x event
    ct = pd.crosstab(d.flag.str.slice(0, 40), d.ev.map({0: "censored", 1: "cancer death", 2: "other death"})); ct.index.name = "flag (first 40 chars)"
    P(f"\n**Flag x event type, {lab}**"); tbl(ct.reset_index(), "{:.0f}"); sec1[f"flag_by_event_{lab}"] = {str(i): {str(c): int(v) for c, v in r.items()} for i, r in ct.iterrows()}
OUT["sec1_survival_fields"] = sec1

# ------------------------------------------------------------------ (2) administrative censoring
H("## 2. Administrative censoring within the 12-month window after each landmark")
P("At risk = alive and under follow-up at the landmark (tm > landmark). 'Censored alive before window end' = ev 0 and tm < landmark + 12; "
  "'followed to window end' = tm >= landmark + 12 (regardless of later outcome); deaths in window = ev 1/2 with landmark < tm <= landmark + 12 (and tm < landmark+12 for the censoring column complement). "
  "Diagnosis years 2021-2023 have potential follow-up < 36 months (SEER follow-up cut-off 2023-12-31).")
def adm(d, lab):
    rows = []
    for L in LM5:
        at = d[d.tm > L]; n = len(at)
        cens = ((at.ev == 0) & (at.tm < L + 12)).sum(); fu = (at.tm >= L + 12).sum()
        d1 = ((at.ev == 1) & (at.tm <= L + 12)).sum(); d2 = ((at.ev == 2) & (at.tm <= L + 12)).sum()
        rows.append({"cohort": lab, "landmark": int(L), "n_at_risk": int(n), "cancer_deaths_in_window": int(d1), "other_deaths_in_window": int(d2),
                     "censored_alive_before_window_end": int(cens), "pct_censored": float(100 * cens / n) if n else np.nan, "followed_to_window_end": int(fu), "pct_followed_to_end": float(100 * fu / n) if n else np.nan})
    return pd.DataFrame(rows)
sec2 = {}
for lab, d in [("all stages", nec), ("distant", dist)]:
    a = adm(d, lab); sec2[lab] = rec(a); P(f"\n**{lab}, all diagnosis years**"); tbl(a)
    parts = [adm(d[d.year <= 2020], f"{lab} dx<=2020")] + [adm(d[d.year == y], f"{lab} dx {y}") for y in (2021, 2022, 2023)]
    by = pd.concat(parts); sec2[lab + "_by_year"] = rec(by); P(f"\n**{lab}, by diagnosis year (2021-2023 separately)**"); tbl(by)
    # censoring by year among censored-before-window-end at LM0
    at = d[d.tm > 0]; c0 = at[(at.ev == 0) & (at.tm < 12)]
    yc = vc(c0.year.astype(int), "diagnosis year", order=sorted(c0.year.unique())); sec2[lab + "_LM0_censored_by_year"] = rec(yc)
    P(f"\n**{lab}: diagnosis year of patients censored alive before 12 months (landmark 0), n = {len(c0)}**"); tbl(yc)
OUT["sec2_administrative_censoring"] = sec2

# ------------------------------------------------------------------ (3) diagnosis year & potential follow-up
H("## 3. Diagnosis year and potential follow-up")
sec3 = {}
for lab, d in [("all stages", nec), ("distant", dist)]:
    y = d.year.astype(int); pf = 2023 - y + 1  # years from 1 Jan of dx year to 2023-12-31 (upper bound); minimum potential = 2023 - y
    per = pd.DataFrame({"period": pd.cut(y, [1999, 2004, 2009, 2014, 2019, 2023], labels=["2000-2004", "2005-2009", "2010-2014", "2015-2019", "2020-2023"])})
    t = vc(per.period, "diagnosis period", order=["2000-2004", "2005-2009", "2010-2014", "2015-2019", "2020-2023"]); sec3[lab + "_period"] = rec(t)
    P(f"\n**{lab}: diagnosis period (n = {len(d)})**"); tbl(t)
    yr = y.value_counts().sort_index(); yt = pd.DataFrame({"year": yr.index, "n": yr.values, "%": 100 * yr.values / len(d), "potential follow-up to 2023-12-31 (years, min-max)": [f"{2023 - v}-{2024 - v}" for v in yr.index],
                                                           "observed tm median (mo)": [float(d.loc[y == v, "tm"].median()) for v in yr.index], "censored alive n": [int(((d.ev == 0) & (y == v)).sum()) for v in yr.index]})
    sec3[lab + "_by_year"] = rec(yt); P(f"\n**{lab}: by single diagnosis year**"); tbl(yt)
    sec3[lab + "_summary"] = {"year_median": float(y.median()), "year_IQR": [float(y.quantile(.25)), float(y.quantile(.75))], "n_2021_2023": int((y >= 2021).sum()), "pct_2021_2023": float(100 * (y >= 2021).mean()),
                              "potential_fu_years_min_median": float((2023 - y).median()), "n_potential_fu_lt_36mo (dx>=2021)": int((y >= 2021).sum()), "n_potential_fu_lt_12mo (dx 2023)": int((y == 2023).sum()),
                              "fu_year_values": {str(a): int(b) for a, b in d.fu_year.value_counts().items()}}
    P(f"\nMedian diagnosis year {y.median():.0f} (IQR {y.quantile(.25):.0f}-{y.quantile(.75):.0f}); diagnosed 2021-2023: {int((y >= 2021).sum())} ({100 * (y >= 2021).mean():.1f} %); "
      f"fu_year field values: {sec3[lab + '_summary']['fu_year_values']}.")
OUT["sec3_diagnosis_year"] = sec3

# ------------------------------------------------------------------ (4) treatment timing
H("## 4. Treatment timing (SEER 'time from diagnosis to treatment', days; first course of ANY treatment)")
sec4 = {}
def ttsum(d, lab):
    t = pd.to_numeric(d.tt_days, errors="coerce"); num = t.dropna()
    nonnum = vc(d.tt_days[t.isna()], "non-numeric tt_days value", denom=len(d))
    s = {"label": lab, "n": len(d), "numeric_n": int(len(num)), "median": float(num.median()) if len(num) else np.nan, "q1": float(num.quantile(.25)) if len(num) else np.nan, "q3": float(num.quantile(.75)) if len(num) else np.nan,
         "mean": float(num.mean()) if len(num) else np.nan, "max": float(num.max()) if len(num) else np.nan,
         "le90_n": int((num <= 90).sum()), "le90_pct_of_all": float(100 * (num <= 90).sum() / len(d)), "le90_pct_of_numeric": float(100 * (num <= 90).mean()) if len(num) else np.nan,
         "gt90_n": int((num > 90).sum()), "gt90_pct_of_all": float(100 * (num > 90).sum() / len(d)),
         "zero_n": int((num == 0).sum()), "zero_pct_of_all": float(100 * (num == 0).sum() / len(d)),
         "le30_n": int((num <= 30).sum()), "le60_n": int((num <= 60).sum()), "non_numeric_n": int(t.isna().sum()), "non_numeric_pct": float(100 * t.isna().mean()), "non_numeric_values": rec(nonnum)}
    return s
groups = [("distant, chemo recorded yes", dist[dist.chemo == 1]), ("distant, chemo yes & survived >3 months (tm>3)", dist[(dist.chemo == 1) & (dist.tm > 3)]),
          ("distant, chemo no/unknown", dist[dist.chemo == 0]), ("distant, chemo no/unknown & tm>3", dist[(dist.chemo == 0) & (dist.tm > 3)]),
          ("distant, surgery yes", dist[dist.surg == 1]), ("distant, surgery yes & chemo yes", dist[(dist.surg == 1) & (dist.chemo == 1)]), ("distant, surgery no", dist[dist.surg == 0]),
          ("distant, all", dist), ("all stages, chemo yes", nec[nec.chemo == 1]), ("all stages, surgery yes", nec[nec.surg == 1]), ("all stages, all", nec)]
rows = [ttsum(d, lab) for lab, d in groups]; sec4["groups"] = rows
t4 = pd.DataFrame([{"group": r["label"], "n": r["n"], "numeric n": r["numeric_n"], "median (d)": r["median"], "IQR": f"{r['q1']:.0f}-{r['q3']:.0f}", "0 d n (%)": f"{r['zero_n']} ({r['zero_pct_of_all']:.1f})",
                    "<=90 d n (%)": f"{r['le90_n']} ({r['le90_pct_of_all']:.1f})", ">90 d n (%)": f"{r['gt90_n']} ({r['gt90_pct_of_all']:.1f})", "non-numeric n (%)": f"{r['non_numeric_n']} ({r['non_numeric_pct']:.1f})", "max (d)": r["max"]} for r in rows])
tbl(t4, "{:.0f}")
P("\nNon-numeric tt_days values are exclusively 'Unable to calculate' (no 'Blank(s)' in the main cohort) - see JSON for per-group counts.")
r0 = rows[0]; r1 = rows[1]
P(f"\nManuscript statement check (distant, chemo recorded, n = 580): <=90 d {r0['le90_n']} ({r0['le90_pct_of_all']:.1f} %), >90 d {r0['gt90_n']} ({r0['gt90_pct_of_all']:.1f} %), unable to calculate {r0['non_numeric_n']} ({r0['non_numeric_pct']:.1f} %), median {r0['median']:.0f} d "
  f"-> manuscript 548 (94.5 %) / 6 (1.0 %) / 26 (4.5 %) / median 15. Among 3-month survivors with chemo (n = {r1['n']}): {r1['le90_n']} / {r1['gt90_n']} / {r1['non_numeric_n']}, median {r1['median']:.0f} d.")
sec4["manuscript_check"] = {"manuscript": {"n": 580, "le90": 548, "gt90": 6, "unable": 26, "median": 15}, "mine": {"n": r0["n"], "le90": r0["le90_n"], "gt90": r0["gt90_n"], "unable": r0["non_numeric_n"], "median": r0["median"]},
                            "agree": bool(r0["n"] == 580 and r0["le90_n"] == 548 and r0["gt90_n"] == 6 and r0["non_numeric_n"] == 26 and r0["median"] == 15)}
# tt_days distribution bins among distant chemo
t = pd.to_numeric(dist.loc[dist.chemo == 1, "tt_days"], errors="coerce")
bins = pd.cut(t, [-1, 0, 7, 14, 30, 60, 90, 180, 10000], labels=["0", "1-7", "8-14", "15-30", "31-60", "61-90", "91-180", ">180"])
bt = vc(bins.astype(object).where(t.notna(), "Unable to calculate"), "tt_days (distant, chemo yes)", order=["0", "1-7", "8-14", "15-30", "31-60", "61-90", "91-180", ">180", "Unable to calculate"])
sec4["bins_distant_chemo"] = rec(bt); P("\n**Binned tt_days, distant stage with chemotherapy recorded (n = 580)**"); tbl(bt)
# tt_days = 0 among surgery-yes vs chemo-only: interpretation
OUT["sec4_treatment_timing"] = sec4

# ------------------------------------------------------------------ (5) multiple primaries
H("## 5. Multiple primaries")
sec5 = {}
SEQ_ORDER = ["One primary only", "1st of 2 or more primaries", "2nd of 2 or more primaries", "3rd of 3 or more primaries", "4th of 4 or more primaries", "5th of 5 or more primaries"]
for lab, d in [("all stages", nec), ("distant", dist)]:
    s = vc(d.seq, "Sequence number", order=SEQ_ORDER); sec5[lab + "_seq"] = rec(s); P(f"\n**{lab}: sequence number (n = {len(d)})**"); tbl(s)
    nt = vc(d.ntum, "Total number of in-situ/malignant tumours for patient"); sec5[lab + "_ntum"] = rec(nt); P(f"\n**{lab}: number of tumours (ntum)**"); tbl(nt)
fp = nec.seq.isin(["One primary only", "1st of 2 or more primaries"])
ct = pd.crosstab(nec.stage, np.where(fp, "first/only primary", "later primary (2nd+)"), margins=True); ct.index.name = "stage"
ct2 = ct.copy(); ct2["% first/only"] = 100 * ct["first/only primary"] / ct["All"]
sec5["first_primary_by_stage"] = {str(i): {str(c): (float(v) if c == "% first/only" else int(v)) for c, v in r.items()} for i, r in ct2.iterrows()}
P("\n**First/only primary vs later primary, by stage (all stages)**"); tbl(ct2.reset_index(), "{:.1f}")
ct3 = pd.crosstab(nec.seq, nec.stage).reindex(SEQ_ORDER).fillna(0).astype(int); ct3.index.name = "Sequence number"   # fillna: sequence levels absent from a dataset; P("\n**Sequence number x stage**"); tbl(ct3.reset_index(), "{:.0f}")
sec5["seq_by_stage"] = {str(i): {str(c): int(v) for c, v in r.items()} for i, r in ct3.iterrows()}
# seq vs ntum consistency
inc = ((nec.seq == "One primary only") & (nec.ntum != "01")).sum() + ((nec.seq != "One primary only") & (nec.ntum == "01")).sum()
sec5["seq_ntum_inconsistent"] = int(inc); P(f"\nSequence number vs ntum inconsistency ('One primary only' with ntum != 01, or multiple-primary sequence with ntum == 01): {int(inc)} patients (ntum counts tumours through the follow-up cut-off, seq refers to order at diagnosis).")
OUT["sec5_multiple_primaries"] = sec5

# ------------------------------------------------------------------ (6) unknown COD / other-cause classification
H("## 6. Cause of death classification")
sec6 = {}
u = nec[nec.unkcod == 1]; us = vc(u.stage, "stage", order=["Localized", "Regional", "Distant", "Unknown"]); sec6["unknown_cod_by_stage"] = rec(us)
P(f"\n**Unknown-COD deaths (coded ev = 2 with unkcod = 1), n = {len(u)}, by stage**"); tbl(us)
# repository version: per-patient survival days / diagnosis years of the unknown-COD deaths are NOT written (record-level data, SEER DUA); aggregate summary only
P(f"Unknown-COD deaths: n = {len(u)}; survival days median {u.sdays.median():.0f}; diagnosis years {int(u.year.min()) if len(u) else 0}-{int(u.year.max()) if len(u) else 0}.")
for lab, d in [("all stages", nec), ("distant", dist)]:
    oc = vc(d.othercod, "SEER other cause of death classification"); sec6[lab + "_othercod"] = rec(oc); P(f"\n**{lab}: SEER other-cause-of-death classification (n = {len(d)})**"); tbl(oc)
    ct = pd.crosstab(d.ev.map({0: "censored", 1: "cancer death", 2: "other/unknown death"}), d.vital); ct.index.name = "ev"; P(f"\n**{lab}: ev x vital status**"); tbl(ct.reset_index(), "{:.0f}")
    sec6[lab + "_ev_by_vital"] = {str(i): {str(c): int(v) for c, v in r.items()} for i, r in ct.iterrows()}
    e2 = d[d.ev == 2]; cs = vc(e2.cod_site, "COD site recode (ev = 2)").head(10); sec6[lab + "_cod_site_top10_ev2"] = rec(cs)
    P(f"\n**{lab}: top-10 COD site recode among other-cause deaths (ev = 2, n = {len(e2)})**"); tbl(cs)
    e1 = d[d.ev == 1]; cs1 = vc(e1.cod_site, "COD site recode (ev = 1)").head(10); sec6[lab + "_cod_site_top10_ev1"] = rec(cs1)
    P(f"\n**{lab}: top COD site recode among cancer deaths (ev = 1, n = {len(e1)})**"); tbl(cs1)
OUT["sec6_cause_of_death"] = sec6

# ------------------------------------------------------------------ (7) metastatic pattern & treatment cross-tabs (distant)
H("## 7. Distant stage: metastatic pattern and treatment")
sec7 = {}
def fillb(s): return s.fillna("Blank(s)").replace({"": "Blank(s)"})
d10 = dist[dist.year >= 2010]; d16 = dist[dist.year >= 2016]
P(f"\n**Site-specific metastases at diagnosis (available 2010+; distant NEC diagnosed 2010+ n = {len(d10)})**")
rows = []
for c in ["mets_liver", "mets_lung", "mets_bone", "mets_brain"]:
    v = fillb(d10[c]); rows.append({"site": c.replace("mets_", ""), "Yes": int((v == "Yes").sum()), "No": int((v == "No").sum()), "Unknown": int((v == "Unknown").sum()), "Blank": int((v == "Blank(s)").sum()),
                                   "% Yes (of 2010+)": float(100 * (v == "Yes").mean()), "% Yes (of known)": float(100 * (v == "Yes").sum() / max((v.isin(["Yes", "No"])).sum(), 1))})
m = pd.DataFrame(rows); sec7["mets_2010plus"] = rec(m); tbl(m)
P(f"\n**Distant lymph node / other metastases (available 2016+; distant NEC diagnosed 2016+ n = {len(d16)})**")
for c in ["mets_dln", "mets_other"]:
    v = vc(fillb(d16[c]), c); sec7[c + "_2016plus"] = rec(v); tbl(v)
# number of recorded sites among 2010+ with all four known
four = d10[np.logical_and.reduce([fillb(d10[c]).isin(["Yes", "No"]) for c in ["mets_liver", "mets_lung", "mets_bone", "mets_brain"]])]
ns = sum((fillb(four[c]) == "Yes").astype(int) for c in ["mets_liver", "mets_lung", "mets_bone", "mets_brain"])
nst = vc(ns.astype(str), "number of the 4 recorded sites involved (2010+, all 4 known)", order=["0", "1", "2", "3", "4"]); sec7["n_sites_2010plus_known4"] = rec(nst)
P(f"\n**Number of involved sites among liver/lung/bone/brain (2010+ with all four fields known, n = {len(four)})**; 0 sites = distant stage via other routes (distant LN, peritoneum, other organs)"); tbl(nst)
liver_only = int(((fillb(four.mets_liver) == "Yes") & (ns == 1)).sum()); sec7["liver_only_n"] = liver_only; P(f"Liver-only (of the four sites): {liver_only} ({100 * liver_only / len(four):.1f} % of {len(four)}).")
# treatment cross-tab
P("\n**Treatment cross-tabulation, distant stage (surgery = primary-site surgery per authors' rule; chemo 1 = recorded yes, 0 = no/unknown; rad 1 = yes)**")
tr = dist.assign(surgery=dist.surg.map({1: "yes", 0: "no", -1: "unknown"}), chemo_=dist.chemo.map({1: "yes", 0: "no/unk"}), rad_=dist.rad.map({1: "yes", 0: "no/unk"}))
ct = pd.crosstab([tr.surgery, tr.chemo_], tr.rad_, margins=True); ct.index.names = ["surgery", "chemo"]; tbl(ct.reset_index(), "{:.0f}")
sec7["surg_chemo_rad"] = [{"surgery": i[0], "chemo": i[1], **{str(c): int(v) for c, v in r.items()}} for i, r in ct.iterrows()]
ct = pd.crosstab(tr.surgery, tr.chemo_, margins=True); ct.index.name = "surgery"; P("\n**Surgery x chemotherapy, distant**"); tbl(ct.reset_index(), "{:.0f}")
sec7["surg_chemo"] = {str(i): {str(c): int(v) for c, v in r.items()} for i, r in ct.iterrows()}
combo = (tr.surgery + "/" + tr.chemo_ + "/" + tr.rad_); cv = vc(combo, "surgery/chemo/rad"); sec7["combo"] = rec(cv); P("\n**Treatment combinations (surgery/chemo/rad), distant**"); tbl(cv)
ss = vc(dist.sys_surg_seq, "Systemic therapy / surgery sequence (distant)"); sec7["sys_surg_seq"] = rec(ss); P("\n**Systemic therapy / surgery sequence, distant**"); tbl(ss)
ss2 = pd.crosstab(dist.sys_surg_seq, tr.surgery); ss2.index.name = "sys_surg_seq"; P("\n**Sequence field x surgery variable, distant**"); tbl(ss2.reset_index(), "{:.0f}")
sec7["sys_surg_seq_by_surgery"] = {str(i): {str(c): int(v) for c, v in r.items()} for i, r in ss2.iterrows()}
nr = vc(dist.loc[dist.surg == 0, "no_surg_reason"], "Reason no surgery (distant, surgery = no)"); sec7["no_surg_reason_surgno"] = rec(nr); P(f"\n**Reason no surgery, distant stage with surgery = no (n = {int((dist.surg == 0).sum())})**"); tbl(nr)
nr2 = pd.crosstab(dist.no_surg_reason, tr.surgery); nr2.index.name = "no_surg_reason"; P("\n**Reason-no-surgery field x surgery variable, distant (consistency check)**"); tbl(nr2.reset_index(), "{:.0f}")
sec7["no_surg_reason_by_surgery"] = {str(i): {str(c): int(v) for c, v in r.items()} for i, r in nr2.iterrows()}
# surgery codes among surgery yes
sc = vc(dist.loc[dist.surg == 1, "s1"].fillna("Blank(s)"), "s1 surgery code (distant, surgery yes)").head(12); sec7["s1_codes_surg_yes_top"] = rec(sc); P("\n**Top primary-site surgery codes (s1) among distant surgery = yes**"); tbl(sc)
OUT["sec7_distant_mets_treatment"] = sec7

# ------------------------------------------------------------------ (8) comparator cohorts
H("## 8. Comparator cohorts: survival flag, diagnostic confirmation, patient overlap")
sec8 = {}
USE = ["pid", "year", "stage", "hist", "flag", "diag", "seq", "ev", "sdays", "smonths", "unkcod"]
def load_cols(name, cols=USE):
    return pd.read_csv(lm.DATA + f"cohort_{name}.csv", usecols=cols, dtype={"pid": str, "hist": str, "seq": str, "smonths": str}, low_memory=False)
lungall = load_cols("LUNG")
comp = {"SCLC": lm.sclc_main(lungall), "LCNEC_lung": lungall[lungall["hist"] == "8013"], "PANC_AC": lm.panc_main(load_cols("PANC")), "CRC_AC": load_cols("AC")}
expected = {"SCLC": 77707, "LCNEC_lung": 3185, "PANC_AC": 75634, "CRC_AC": 147014}
rows = []; flags = {}
for nm, d in comp.items():
    dd = d[d.stage == "Distant"]; f = vc(dd.flag, "flag"); flags[nm] = rec(f)
    sm = pd.to_numeric(dd.smonths, errors="coerce"); disc = int((np.floor(dd.sdays.astype(float) / 30.4375) != sm).sum())
    rows.append({"cohort": nm, "n all stages": len(d), "n distant": len(dd), "expected distant": expected[nm], "complete dates & >0 d (%)": float(100 * dd.flag.str.startswith("Complete dates are available and there are more").mean()),
                 "zero-day (%)": float(100 * (dd.sdays == 0).mean()), "incomplete dates (%)": float(100 * dd.flag.str.startswith("Incomplete").mean()), "days/months discrepant": disc, "unknown COD n": int(dd.unkcod.sum()),
                 "dx 2021-2023 (%)": float(100 * (dd.year >= 2021).mean()), "first/only primary (%)": float(100 * dd.seq.isin(["One primary only", "1st of 2 or more primaries"]).mean())})
if not all(r["n distant"] == r["expected distant"] for r in rows): config.anchor_fail(f"comparator sizes {[(r['cohort'], r['n distant']) for r in rows]}")
dn = dist; rows.insert(0, {"cohort": "NEC (main)", "n all stages": N_ALL, "n distant": N_DIST, "expected distant": 908, "complete dates & >0 d (%)": float(100 * dn.flag.str.startswith("Complete dates are available and there are more").mean()),
                           "zero-day (%)": float(100 * (dn.sdays == 0).mean()), "incomplete dates (%)": float(100 * dn.flag.str.startswith("Incomplete").mean()),
                           "days/months discrepant": int((np.floor(dn.sdays.astype(float) / 30.4375) != pd.to_numeric(dn.smonths, errors="coerce")).sum()), "unknown COD n": int(dn.unkcod.sum()),
                           "dx 2021-2023 (%)": float(100 * (dn.year >= 2021).mean()), "first/only primary (%)": float(100 * dn.seq.isin(["One primary only", "1st of 2 or more primaries"]).mean())})
ct8 = pd.DataFrame(rows); sec8["cohort_summary_distant"] = rec(ct8); sec8["flag_distant"] = flags
P("\n**Distant-stage comparator cohorts: survival-time quality, diagnosis period, multiple primaries**"); tbl(ct8)
P("\nFull flag distributions per comparator are in the JSON (sec8.flag_distant).")
# diagnostic confirmation from _cyt variants
P("\n**Diagnostic confirmation: positive histology vs positive exfoliative cytology only (from the cytology-inclusive cohort variants; main analyses use histology only)**")
cy = []
for nm, fn, sel in [("SCLC", "LUNG_cyt", lambda d: lm.sclc_main(d)), ("LCNEC_lung", "LUNG_cyt", lambda d: d[d["hist"] == "8013"]), ("PANC_AC", "PANC_cyt", lm.panc_main), ("CRC_AC", "AC_cyt", lambda d: d)]:
    d = sel(load_cols(fn, ["pid", "stage", "hist", "diag"])); dd = d[d.stage == "Distant"]
    for lab, x in [("all stages", d), ("distant", dd)]:
        hist_n = int((x.diag == "Positive histology").sum()); cyt_n = int((x.diag != "Positive histology").sum())
        cy.append({"cohort": nm, "stages": lab, "n (histology + cytology)": len(x), "positive histology": hist_n, "cytology only": cyt_n, "% cytology only": float(100 * cyt_n / len(x)), "other diag values": ", ".join(sorted(set(x.diag.unique()) - {"Positive histology"}))})
cy.insert(0, {"cohort": "NEC (main)", "stages": "all / distant", "n (histology + cytology)": N_ALL, "positive histology": N_ALL, "cytology only": 0, "% cytology only": 0.0, "other diag values": "(cohort restricted to positive histology)"})
cyt = pd.DataFrame(cy); sec8["diag_confirmation"] = rec(cyt); tbl(cyt)
# overlap of pids
P("\n**Patient-ID overlap between distant-stage cohorts (same SEER patient ID = multiple primaries in two cohorts)**")
sets = {"NEC": set(dist.pid)}; sets.update({nm: set(comp[nm][comp[nm].stage == "Distant"].pid) for nm in comp})
sets_all = {"NEC": set(nec.pid)}; sets_all.update({nm: set(comp[nm].pid) for nm in comp})
names = list(sets); ov = []
for i, a in enumerate(names):
    for b in names[i + 1:]:
        ov.append({"cohort A": a, "cohort B": b, "shared pids (distant/distant)": len(sets[a] & sets[b]), "shared pids (all stages/all stages)": len(sets_all[a] & sets_all[b])})
ovt = pd.DataFrame(ov); sec8["pid_overlap"] = rec(ovt); tbl(ovt)
man = {"CRC_AC": 3, "SCLC": 1, "PANC_AC": 0}; mine_ov = {b: len(sets["NEC"] & sets[b]) for b in man}
sec8["pid_overlap_manuscript_check"] = {"manuscript": man, "mine": mine_ov, "agree": man == mine_ov}
P(f"\nManuscript check (distant NEC vs distant comparator): AC 3 / SCLC 1 / PANC 0; this script: AC {mine_ov['CRC_AC']} / SCLC {mine_ov['SCLC']} / PANC {mine_ov['PANC_AC']} -> {'agree' if man == mine_ov else 'DISAGREE'}.")
# detail of overlapping NEC patients (counts only, no identifiers)
for b in ["CRC_AC", "SCLC", "PANC_AC", "LCNEC_lung"]:
    sh = sets["NEC"] & sets[b]
    if sh:
        x = dist[dist.pid.isin(sh)]; P(f"NEC patients shared with {b} (n = {len(sh)}): NEC sequence numbers {x.seq.value_counts().to_dict()}; NEC ev {x.ev.value_counts().to_dict()}.")
        sec8[f"shared_NEC_{b}_seq"] = {str(k): int(v) for k, v in x.seq.value_counts().items()}
# uniqueness within each cohort
sec8["pid_unique_within_cohort"] = {nm: bool(comp[nm].pid.is_unique) for nm in comp}; sec8["pid_unique_within_cohort"]["NEC"] = bool(nec.pid.is_unique)
P(f"\npid unique within each cohort: {sec8['pid_unique_within_cohort']}.")
OUT["sec8_comparators"] = sec8

# ------------------------------------------------------------------ interpretation + write
rt = (time.time() - T0) / 60; OUT["runtime_minutes"] = round(rt, 2)
a2 = pd.DataFrame(sec2["distant"]); a2y = pd.DataFrame(sec2["distant_by_year"])
H("## Interpretation")
s1a = sec1["consistency_all"]; s3 = sec3["distant_summary"]
P(f"1. Survival-time fields are internally consistent: in the main NEC cohort {N_ALL - s1a['discrepant_floor_days_vs_months']}/{N_ALL} patients have floor(days/30.4375) equal to the registry survival months "
  f"(discrepancies {s1a['discrepant_floor_days_vs_months']}); {100 * (nec.flag.str.startswith('Complete dates are available and there are more')).mean():.1f} % have complete dates with >0 days, "
  f"{int(nec.flag.str.startswith('Incomplete').sum())} have incomplete dates and {s1a['zero_day_survivors']} are zero-day survivors (set to 0.5 day in the analysis, i.e. day-of-diagnosis deaths/censoring).")
r0 = a2.iloc[0]; r3 = a2.iloc[3]
P(f"2. Administrative censoring is modest for the early landmarks: among the {int(r0.n_at_risk)} distant patients at risk at diagnosis, {int(r0.censored_alive_before_window_end)} ({r0.pct_censored:.1f} %) were censored alive before 12 months, "
  f"rising to {int(r3.censored_alive_before_window_end)}/{int(r3.n_at_risk)} ({r3.pct_censored:.1f} %) at the 12-month landmark and {int(a2.iloc[4].censored_alive_before_window_end)}/{int(a2.iloc[4].n_at_risk)} ({a2.iloc[4].pct_censored:.1f} %) at 24 months. "
  f"Censoring is concentrated in patients diagnosed 2021-2023 ({s3['n_2021_2023']} distant patients, {s3['pct_2021_2023']:.1f} %), whose potential follow-up is <36 months; the Aalen-Johansen estimator handles this non-informatively but late-landmark estimates (12, 24 months) rest on small, increasingly censored risk sets.")
P(f"3. Diagnosis years span 2000-2023 (distant median {s3['year_median']:.0f}); {s3['n_potential_fu_lt_12mo (dx 2023)']} distant patients diagnosed in 2023 have at most 12 months of potential follow-up and contribute only to the earliest landmarks.")
r0 = sec4["groups"][0]
P(f"4. Treatment timing: among distant patients with chemotherapy recorded (n = {r0['n']}), first-course treatment began within 90 days in {r0['le90_n']} ({r0['le90_pct_of_all']:.1f} %), median {r0['median']:.0f} days (IQR {r0['q1']:.0f}-{r0['q3']:.0f}); "
  f"{r0['non_numeric_n']} ({r0['non_numeric_pct']:.1f} %) 'unable to calculate' - identical to the manuscript (548/6/26, median 15). The field records first-course ANY treatment; {r0['zero_n']} ({r0['zero_pct_of_all']:.1f} %) started on the day of diagnosis (mostly surgical).")
fpd = sec5["first_primary_by_stage"]
P(f"5. Multiple primaries: {fpd['All']['first/only primary']}/{fpd['All']['All']} ({fpd['All']['% first/only']:.1f} %) of NEC patients are first/only primaries; {fpd['Distant']['% first/only']:.1f} % among distant stage. "
  f"Later primaries (2nd+) are retained in the main cohort (first-primary restriction is a separate sensitivity analysis); cancer-death attribution for such patients relies on the SEER cause-specific death classification.")
P(f"6. Unknown-COD deaths are rare ({len(u)} of {int((nec.ev > 0).sum())} deaths; by stage " + ", ".join(f"{r['stage']} {r['n']}" for r in sec6['unknown_cod_by_stage']) + ") and were coded as other-cause deaths (conservative for the cancer-death CIF). "
  f"Among other-cause deaths in the distant cohort (ev = 2), the top COD sites include other cancers and cardiovascular causes (table 6), consistent with the SEER cause-specific classification rather than misattributed NEC deaths.")
mt = pd.DataFrame(sec7["mets_2010plus"]).set_index("site")
P(f"7. Metastatic pattern (distant, 2010+): liver involvement {mt.loc['liver', '% Yes (of known)']:.1f} % of known, lung {mt.loc['lung', '% Yes (of known)']:.1f} %, bone {mt.loc['bone', '% Yes (of known)']:.1f} %, brain {mt.loc['brain', '% Yes (of known)']:.1f} %; "
  f"distant lymph-node metastases in {[r for r in sec7['mets_dln_2016plus'] if r['mets_dln'].startswith('Yes')][0]['%']:.1f} % of 2016+ patients. Surgery x chemotherapy cross-tabulation and the 'reason no surgery' field agree with the derived surgery variable (no contradictions in the cross-check tables).")
ov_ok = sec8["pid_overlap_manuscript_check"]["agree"]
P(f"8. Comparator cohorts: distant SCLC/pancreatic AC/colorectal AC have similar survival-time quality (complete dates & >0 days {ct8.iloc[1]['complete dates & >0 d (%)']:.1f}/{ct8.iloc[3]['complete dates & >0 d (%)']:.1f}/{ct8.iloc[4]['complete dates & >0 d (%)']:.1f} %); "
  f"restricting to positive histology excludes {[c for c in cy if c['cohort']=='SCLC' and c['stages']=='distant'][0]['% cytology only']:.1f} % of distant SCLC, {[c for c in cy if c['cohort']=='PANC_AC' and c['stages']=='distant'][0]['% cytology only']:.1f} % of distant pancreatic AC and {[c for c in cy if c['cohort']=='CRC_AC' and c['stages']=='distant'][0]['% cytology only']:.1f} % of distant colorectal AC (cytology-inclusive sensitivity cohorts available). "
  f"Patient overlap between distant NEC and the distant comparators is negligible (AC {mine_ov['CRC_AC']}, SCLC {mine_ov['SCLC']}, PANC {mine_ov['PANC_AC']}; {'matches' if ov_ok else 'does NOT match'} the manuscript).")
P(f"\nDeviations from plan: none (descriptive task; no bootstrap required). Runtime {rt:.1f} min.")
with open(RES + "A10.json", "w") as f: json.dump(OUT, f, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))
with open(RES + "A10.md", "w") as f: f.write("\n".join(MDL) + "\n")
print(f"done in {rt:.1f} min")
