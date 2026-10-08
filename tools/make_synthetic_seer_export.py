#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""SYNTHETIC SEER*Stat case-listing exports for running the pipeline end-to-end WITHOUT SEER data access.

    python3 tools/make_synthetic_seer_export.py --out data/raw            # then: ./run_all.sh --quick

Writes four tab-delimited files with exactly the 101 column headers, file names, quoting ("quotedlabels") and
line endings (CRLF) of the authors' SEER*Stat exports, plus SYNTHETIC_DATA_README.txt (the pipeline detects this
marker, switches to non-strict anchor checking and labels the run as synthetic).

EVERY VALUE IS RANDOMLY GENERATED. No SEER record, and nothing derived from an individual SEER record, is contained
in or used by this script: category labels are the public SEER*Stat variable formats, and survival times are drawn
from simple piecewise-exponential hazards chosen to look plausible. Results computed from these files are
meaningless and will NOT reproduce the paper (the reproduction anchors are expected to fail).
"""
import argparse
import os
import numpy as np

COLUMNS = [
    'Age recode with <1 year olds and 90+', 'Sex', 'Year of diagnosis', 'Race recode (W, B, AI, API)',
    'Origin recode NHIA (Hispanic, Non-Hisp)', 'Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)',
    'Site recode ICD-O-3 2023 Revision', 'Primary Site - labeled', 'Histologic Type ICD-O-3', 'Behavior recode for analysis',
    'Grade Recode (thru 2017)', 'Derived Summary Grade 2018 (2018+)', 'Grade Clinical (2018+)', 'Grade Pathological (2018+)',
    'Laterality', 'Combined Summary Stage with Expanded Regional Codes (2004+)', 'Summary stage 2000 (1998-2017)',
    'SEER historic stage A (1973-2015)', 'Derived EOD 2018 T Recode (2018+)', 'Derived EOD 2018 N Recode (2018+)',
    'Derived EOD 2018 M Recode (2018+)', 'Derived EOD 2018 Stage Group Recode (2018+)', 'Derived AJCC Stage Group, 7th ed (2010-2015)',
    'Derived AJCC T, 7th ed (2010-2015)', 'Derived AJCC N, 7th ed (2010-2015)', 'Derived AJCC M, 7th ed (2010-2015)',
    '7th Edition Stage Group Recode (2016-2017)', 'Derived SEER Cmb Stg Grp (2016-2017)', 'Derived SEER Combined T (2016-2017)',
    'Derived SEER Combined N (2016-2017)', 'Derived SEER Combined M (2016-2017)', 'Derived AJCC Stage Group, 6th ed (2004-2015)',
    'Derived AJCC T, 6th ed (2004-2015)', 'Derived AJCC N, 6th ed (2004-2015)', 'Derived AJCC M, 6th ed (2004-2015)',
    'AJCC stage 3rd edition (1988-2003)', 'SEER modified AJCC stage 3rd (1988-2003)', 'T value - based on AJCC 3rd (1988-2003)',
    'N value - based on AJCC 3rd (1988-2003)', 'M value - based on AJCC 3rd (1988-2003)', 'TNM Edition Number (2016-2017)',
    'RX Summ--Surg Prim Site (1998-2022)', 'RX Summ--Surg Prim Site 2023 (2023+)', 'RX Summ--Scope Reg LN Sur (2003+)',
    'RX Summ--Surg Oth Reg/Dis (2003+)', 'RX Summ--Surg/Rad Seq', 'Reason no cancer-directed surgery', 'Radiation recode',
    'Chemotherapy recode (yes, no/unk)', 'Scope of reg lymph nd surg (1998-2002)', 'RX Summ--Reg LN Examined (1998-2002)',
    'Surgery of oth reg/dis sites (1998-2002)', 'RX Summ--Systemic/Sur Seq (2007+)', 'Time from diagnosis to treatment in days recode',
    'EOD Primary Tumor Recode (2018+)', 'EOD Regional Nodes Recode (2018+)', 'EOD Mets Recode (2018+)',
    'Tumor Size Over Time Recode (1988+)', 'Tumor Size Summary (2016+)', 'Regional nodes examined (1988+)',
    'Regional nodes positive (1988+)', 'SEER Combined Mets at DX-bone (2010+)', 'SEER Combined Mets at DX-brain (2010+)',
    'SEER Combined Mets at DX-liver (2010+)', 'SEER Combined Mets at DX-lung (2010+)', 'Mets at DX-Distant LN (2016+)',
    'Mets at DX-Other (2016+)', 'Lymph-vascular Invasion (2004+ varying by schema)', 'CS tumor size (2004-2015)',
    'CS extension (2004-2015)', 'CS lymph nodes (2004-2015)', 'CS mets at dx (2004-2015)', 'CS Tumor Size/Ext Eval (2004-2015)',
    'CS Reg Node Eval (2004-2015)', 'CS Mets Eval (2004-2015)', 'EOD 10 - extent (1988-2003)', 'EOD 10 - nodes (1988-2003)',
    'EOD 10 - size (1988-2003)', 'Coding system-EOD (1973-2003)', 'COD to site recode', 'SEER cause-specific death classification',
    'SEER other cause of death classification', 'Survival months', 'Survival Days', 'Survival months flag',
    'Vital status recode (study cutoff used)', 'Sequence number', 'Primary by international rules', 'Record number recode',
    'Total number of in situ/malignant tumors for patient', 'Age recode with single ages and 90+',
    'Race recode (with detailed Asian and Native Hawaiian other PI)', 'Year of follow-up recode', 'Patient ID',
    'Type of Reporting Source', 'Marital status at diagnosis', 'Median household income inflation adj to 2024',
    'Rural-Urban Continuum Code', 'Site rec ICD-O-3/WHO 2008 (individual sites only)', 'Race and origin (recommended by SEER)',
    'Diagnostic Confirmation',
]
assert len(COLUMNS) == 101

FILES = {"NEN": "colorectal_NEN_SEER17_2000_2023_full_20260806.txt", "AC": "colorectal_adenocarcinoma_SEER17_2000_2023_full_20260928.txt",
         "LUNG": "SCLC_lung_SEER17_2000_2023_full_20260928.txt", "PANC": "pancreatic_adenocarcinoma_SEER17_2000_2023_full_20260928.txt"}
COLO_SITES = ['C18.0-Cecum', 'C18.1-Appendix', 'C18.2-Ascending colon', 'C18.3-Hepatic flexure of colon', 'C18.4-Transverse colon',
              'C18.5-Splenic flexure of colon', 'C18.6-Descending colon', 'C18.7-Sigmoid colon', 'C18.8-Overlapping lesion of colon',
              'C18.9-Colon, NOS', 'C19.9-Rectosigmoid junction', 'C20.9-Rectum, NOS']
COLO_P = [.10, .05, .09, .03, .06, .02, .03, .14, .02, .04, .07, .35]
LUNG_SITES = ['C34.0-Main bronchus', 'C34.1-Upper lobe, lung', 'C34.2-Middle lobe, lung', 'C34.3-Lower lobe, lung',
              'C34.8-Overlapping lesion of lung', 'C34.9-Lung, NOS']
LUNG_P = [.10, .38, .04, .22, .03, .23]
PANC_SITES = ['C25.0-Head of pancreas', 'C25.1-Body of pancreas', 'C25.2-Tail of pancreas', 'C25.3-Pancreatic duct',
              'C25.4-Islets of Langerhans', 'C25.7-Other specified parts of pancreas', 'C25.8-Overlapping lesion of pancreas', 'C25.9-Pancreas, NOS']
PANC_P = [.45, .1, .1, .01, .01, .03, .1, .2]
# tumour definitions: histology codes / probabilities, stage probabilities, monthly cancer-death hazards (distant stage)
TUMOURS = {
    "NEN": dict(n=6000, sites=COLO_SITES, sp=COLO_P, site_recode="Colon and Rectum", cod="Colon excluding Rectum",
                hist={"8041": .20, "8013": .13, "8045": .012, "8042": .004, "8043": .004, "8044": .004, "8240": .22, "8241": .08,
                      "8243": .02, "8244": .02, "8245": .01, "8246": .22, "8249": .066},
                cyt=.01, haz=[.15, .09, .13, .12, .10, .07]),
    "AC": dict(n=30000, sites=COLO_SITES, sp=COLO_P, site_recode="Colon and Rectum", cod="Colon excluding Rectum",
               hist={"8140": .70, "8480": .10, "8210": .03, "8263": .04, "8260": .03, "8481": .05, "8490": .02, "8255": .02, "8240": .01},
               cyt=.01, haz=[.07, .05, .045, .04, .035, .03]),
    "LUNG": dict(n=20000, sites=LUNG_SITES, sp=LUNG_P, site_recode="Lung and Bronchus", cod="Lung and Bronchus",
                 hist={"8041": .84, "8013": .07, "8042": .03, "8045": .03, "8044": .02, "8043": .006, "8002": .003, "8040": .001},
                 cyt=.15, haz=[.14, .10, .12, .10, .08, .05]),
    "PANC": dict(n=15000, sites=PANC_SITES, sp=PANC_P, site_recode="Pancreas", cod="Pancreas",
                 hist={"8140": .52, "8500": .36, "8480": .03, "8150": .03, "8010": .02, "8246": .01, "8560": .03},
                 cyt=.18, haz=[.16, .11, .09, .08, .07, .05]),
}
NEC_LIKE = {"8013", "8041", "8042", "8043", "8044", "8045", "8002", "8040", "8246"}
BREAKS = np.array([0, 3, 6, 9, 12, 18])                        # months; last hazard applies beyond 18 months
STAGE_MULT = {"Localized": .07, "Regional": .3, "Distant": 1.0, "Unknown": .6}
INCOME = ['< $40,000', '$40,000 - $44,999', '$45,000 - $49,999', '$50,000 - $54,999', '$55,000 - $59,999', '$60,000 - $64,999',
          '$65,000 - $69,999', '$70,000 - $74,999', '$75,000 - $79,999', '$80,000 - $84,999', '$85,000 - $89,999',
          '$90,000 - $94,999', '$95,000 - $99,999', '$100,000 - $109,999', '$110,000 - $119,999', '$120,000+']
RURAL = ['Counties in metropolitan areas ge 1 million pop', 'Counties in metropolitan areas of 250,000 to 1 million pop',
         'Counties in metropolitan areas of lt 250 thousand pop', 'Nonmetropolitan counties adjacent to a metropolitan area',
         'Nonmetropolitan counties not adjacent to a metropolitan area']
RACE = [('Non-Hispanic White', 'White', 'Non-Spanish-Hispanic-Latino', .64), ('Non-Hispanic Black', 'Black', 'Non-Spanish-Hispanic-Latino', .12),
        ('Hispanic (All Races)', 'White', 'Spanish-Hispanic-Latino', .13), ('Non-Hispanic Asian or Pacific Islander', 'Asian or Pacific Islander', 'Non-Spanish-Hispanic-Latino', .09),
        ('Non-Hispanic American Indian/Alaska Native', 'American Indian/Alaska Native', 'Non-Spanish-Hispanic-Latino', .01),
        ('Non-Hispanic Unknown Race', 'Unknown', 'Non-Spanish-Hispanic-Latino', .01)]
AGE5 = ['00 years', '01-04 years', '05-09 years', '10-14 years', '15-19 years', '20-24 years', '25-29 years', '30-34 years', '35-39 years',
        '40-44 years', '45-49 years', '50-54 years', '55-59 years', '60-64 years', '65-69 years', '70-74 years', '75-79 years',
        '80-84 years', '85-89 years', '90+ years']
G17 = ['Well differentiated; Grade I', 'Moderately differentiated; Grade II', 'Poorly differentiated; Grade III', 'Undifferentiated; anaplastic; Grade IV', 'Unknown']
G18 = ['Well differentiated', 'Moderately differentiated', 'Poorly differentiated', 'Undifferentiated and anaplastic',
       'Site-specific grade system category (1)', 'Site-specific grade system category (2)', 'Site-specific grade system category (3)',
       'Site-specific grade system category (4)', 'Grade cannot be assessed; Unknown']
OTHER_COD = ['Diseases of Heart', 'Cerebrovascular Diseases', 'Chronic Obstructive Pulmonary Disease and Allied Cond', 'Other Cause of Death', 'Alzheimers (ICD-9 and 10 only)']
MARITAL = ['Married (including common law)', 'Single (never married)', 'Widowed', 'Divorced', 'Separated', 'Unknown']


def pick(rng, levels, p, n):
    p = np.asarray(p, float); return np.asarray(levels, dtype=object)[rng.choice(len(levels), size=n, p=p / p.sum())]


def pw_times(rng, rates, n):
    """event times (months) from piecewise-constant hazards on BREAKS; rates: (n, len(BREAKS))"""
    e = rng.exponential(1.0, n); widths = np.diff(BREAKS)
    cum = np.concatenate([np.zeros((n, 1)), np.cumsum(rates[:, :-1] * widths, axis=1)], axis=1)   # cumulative hazard at each break
    t = np.full(n, np.inf)
    for j in range(len(BREAKS)):
        last = j == len(BREAKS) - 1
        hit = np.isinf(t) & (e >= cum[:, j]) & (last | (e < (cum[:, j + 1] if not last else np.inf)))
        t[hit] = BREAKS[j] + (e[hit] - cum[hit, j]) / rates[hit, j]
    return t


def make(name, cfg, rng, pid_start, scale, share_pids=None):
    n = max(200, int(cfg["n"] * scale))
    hist = pick(rng, list(cfg["hist"]), list(cfg["hist"].values()), n)
    year = rng.integers(2000, 2024, n)
    nec_like = np.isin(hist, list(NEC_LIKE))
    stage_p_nec = [.10, .25, .61, .04]; stage_p_other = {"NEN": [.55, .25, .15, .05], "AC": [.38, .37, .21, .04], "LUNG": [.05, .2, .7, .05], "PANC": [.12, .28, .55, .05]}[name]
    stage = np.where(nec_like & (name == "NEN"), pick(rng, ["Localized", "Regional", "Distant", "Unknown"], stage_p_nec, n),
                     pick(rng, ["Localized", "Regional", "Distant", "Unknown"], stage_p_other, n))
    age = np.clip(np.rint(rng.normal(66, 12, n)), 16, 99).astype(int)
    male = rng.random(n) < .52
    race = rng.choice(len(RACE), size=n, p=[r[3] for r in RACE])
    site = pick(rng, cfg["sites"], cfg["sp"], n)
    # survival: cancer-specific piecewise hazard x stage multiplier x mild age effect; other-cause hazard rising with age
    mult = np.array([STAGE_MULT[s] for s in stage]) * np.exp(0.01 * (age - 66)) * np.where(nec_like | (name != "NEN"), 1.0, 0.08)
    rates = np.outer(mult, cfg["haz"])
    tc = pw_times(rng, rates, n)
    to = rng.exponential(1 / (0.002 * np.exp(0.07 * (age - 66))), n)
    dx_frac = rng.random(n)
    fu = (2023 + 1 - (year + dx_frac)) * 12.0                     # potential follow-up to 31 Dec 2023, months
    lost = np.where(rng.random(n) < .03, rng.uniform(0, fu), np.inf)
    t = np.minimum.reduce([tc, to, fu, lost])
    ev = np.where(t == tc, 1, np.where(t == to, 2, 0))
    days = np.floor(t * 30.4375).astype(int)
    unk_cod = (ev == 2) & (rng.random(n) < .02)
    chemo = rng.random(n) < np.select([stage == "Distant", stage == "Regional"], [.62, .4], .1)
    surg_p = {"Localized": .85, "Regional": .75, "Distant": .35, "Unknown": .3}
    surg = rng.random(n) < np.array([surg_p[s] for s in stage])
    surg_unk = rng.random(n) < .01
    rad = rng.random(n) < .1
    rows = {c: np.full(n, "Blank(s)", dtype=object) for c in COLUMNS}
    R = rows
    R['Age recode with single ages and 90+'] = np.array([("90+ years" if a >= 90 else f"{a:02d} years") for a in age], dtype=object)
    R['Age recode with <1 year olds and 90+'] = np.array([AGE5[min(19, 1 + (a - 1) // 5) if a > 0 else 0] if a < 90 else '90+ years' for a in age], dtype=object)
    R['Sex'] = np.where(male, "Male", "Female")
    R['Year of diagnosis'] = year.astype(str)
    R['Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)'] = np.array([RACE[i][0] for i in race], dtype=object)
    R['Race recode (W, B, AI, API)'] = np.array([RACE[i][1] for i in race], dtype=object)
    R['Origin recode NHIA (Hispanic, Non-Hisp)'] = np.array([RACE[i][2] for i in race], dtype=object)
    R['Site recode ICD-O-3 2023 Revision'] = cfg["site_recode"]; R['Site rec ICD-O-3/WHO 2008 (individual sites only)'] = cfg["site_recode"]
    R['Primary Site - labeled'] = site
    R['Histologic Type ICD-O-3'] = hist
    R['Behavior recode for analysis'] = np.where(rng.random(n) < .997, "Malignant", "In situ")
    hi_grade = nec_like & (rng.random(n) < (np.where(hist == "8246", .45, .9)))
    R['Grade Recode (thru 2017)'] = np.where(year <= 2017, np.where(hi_grade, pick(rng, G17[2:4], [.7, .3], n), pick(rng, G17, [.25, .3, .15, .02, .28], n)), "Blank(s)")
    R['Derived Summary Grade 2018 (2018+)'] = np.where(year >= 2018, np.where(hi_grade, pick(rng, G18[6:8] + G18[2:3], [.5, .05, .45], n), pick(rng, G18, [.1, .2, .1, .01, .1, .1, .03, .01, .35], n)), "Blank(s)")
    reg_lab = pick(rng, ['Regional by direct extension only', 'Regional lymph nodes involved only', 'Regional by both direct extension and lymph node involvement'], [.3, .3, .4], n)
    cs = np.select([stage == "Localized", stage == "Regional", stage == "Distant"], ["Localized only", reg_lab, "Distant site(s)/node(s) involved"], "Unknown/unstaged/unspecified/DCO")
    ss = np.select([stage == "Localized", stage == "Regional", stage == "Distant"], ["Localized", "Regional", "Distant"], "Unknown/unstaged")
    R['Combined Summary Stage with Expanded Regional Codes (2004+)'] = np.where(year >= 2004, cs, "Blank(s)")
    R['Summary stage 2000 (1998-2017)'] = np.where(year <= 2017, ss, "Blank(s)")
    s1 = np.where(surg_unk, "99", np.where(surg, pick(rng, ["30", "40", "50", "60"], [.4, .3, .2, .1], n), "00"))
    s2 = np.where(surg_unk, "A990", np.where(surg, pick(rng, ["A300", "A400", "A500"], [.5, .3, .2], n), "A000"))
    R['RX Summ--Surg Prim Site (1998-2022)'] = np.where(year <= 2022, s1, "Blank(s)")
    R['RX Summ--Surg Prim Site 2023 (2023+)'] = np.where(year == 2023, s2, "Blank(s)")
    R['Reason no cancer-directed surgery'] = np.where(surg & ~surg_unk, "Surgery performed", np.where(surg_unk, "Recommended, unknown if performed", "Not recommended"))
    R['Radiation recode'] = np.where(rad, "Beam radiation", "None/Unknown")
    R['Chemotherapy recode (yes, no/unk)'] = np.where(chemo, "Yes", "No/Unknown")
    R['RX Summ--Systemic/Sur Seq (2007+)'] = np.where(year >= 2007, np.where(surg & chemo, "Systemic therapy after surgery", np.where(chemo, "No systemic therapy and/or surgical procedures", "No systemic therapy and/or surgical procedures")), "Blank(s)")
    treated = chemo | surg | rad
    tt = np.minimum(np.rint(rng.exponential(20, n)), 730).astype(int)
    tt = np.minimum(tt, np.maximum(days, 0))
    R['Time from diagnosis to treatment in days recode'] = np.where(treated & (rng.random(n) < .95), tt.astype(str), "Unable to calculate")
    for c, p_yes in (('SEER Combined Mets at DX-liver (2010+)', .65), ('SEER Combined Mets at DX-bone (2010+)', .15),
                     ('SEER Combined Mets at DX-brain (2010+)', .06), ('SEER Combined Mets at DX-lung (2010+)', .2)):
        v = np.where(stage == "Distant", pick(rng, ["Yes", "No", "Unknown"], [p_yes, .97 - p_yes, .03], n), "No")
        R[c] = np.where(year >= 2010, v, "Blank(s)")
    R['Mets at DX-Distant LN (2016+)'] = np.where(year >= 2016, np.where(stage == "Distant", pick(rng, ['Yes; distant lymph node metastases', 'None; no lymph node metastases', 'Unknown'], [.3, .65, .05], n), 'None; no lymph node metastases'), "Blank(s)")
    R['Mets at DX-Other (2016+)'] = np.where(year >= 2016, np.where(stage == "Distant", pick(rng, ['Yes; distant mets in known site(s) other than bone, brain, liver, lung, dist LN', 'None; no other metastases', 'Unknown'], [.3, .65, .05], n), 'None; no other metastases'), "Blank(s)")
    R['COD to site recode'] = np.where(ev == 0, "Alive", np.where(ev == 1, cfg["cod"], np.where(unk_cod, "State DC not available or state DC available but no COD", pick(rng, OTHER_COD, [1] * len(OTHER_COD), n))))
    R['SEER cause-specific death classification'] = np.where(ev == 1, "Dead (attributable to this cancer dx)", np.where(unk_cod, "Dead (missing/unknown COD)", "Alive or dead of other cause"))
    R['SEER other cause of death classification'] = np.where(ev == 2, np.where(unk_cod, "Dead (missing/unknown COD)", "Dead (attributable to causes other than this cancer dx)"), "Alive or dead due to cancer")
    sm = np.floor(days / 30.4375).astype(int)
    unknown_surv = rng.random(n) < .002
    R['Survival Days'] = np.where(unknown_surv, "Unknown", days.astype(str))
    R['Survival months'] = np.where(unknown_surv, "Unknown", np.array([f"{m:04d}" for m in sm], dtype=object))
    R['Survival months flag'] = np.where(days == 0, "Complete dates are available and there are 0 days of survival",
                                         np.where(rng.random(n) < .97, "Complete dates are available and there are more than 0 days of survival",
                                                  "Incomplete dates are available and there cannot be zero days of follow-up"))
    R['Vital status recode (study cutoff used)'] = np.where(ev > 0, "Dead", "Alive")
    seq = pick(rng, ["One primary only", "1st of 2 or more primaries", "2nd of 2 or more primaries", "3rd of 3 or more primaries"], [.78, .1, .1, .02], n)
    R['Sequence number'] = seq
    R['Total number of in situ/malignant tumors for patient'] = np.where(seq == "One primary only", "01", np.where(seq == "3rd of 3 or more primaries", "03", "02"))
    R['Record number recode'] = "01"
    R['Primary by international rules'] = np.where(rng.random(n) < .99, "Yes", "No")
    R['Year of follow-up recode'] = np.minimum(year + np.floor((dx_frac * 12 + t) / 12).astype(int), 2023).astype(str)
    R['Type of Reporting Source'] = pick(rng, ["Hospital inpatient/outpatient or clinic", "Physicians office/private medical practitioner (LMD)",
                                               "Laboratory only (hospital or private)", "Death certificate only", "Autopsy only"], [.93, .03, .03, .007, .003], n)
    R['Marital status at diagnosis'] = pick(rng, MARITAL, [.5, .15, .15, .12, .03, .05], n)
    R['Median household income inflation adj to 2024'] = pick(rng, INCOME, [1] * len(INCOME), n)
    R['Rural-Urban Continuum Code'] = pick(rng, RURAL, [.6, .2, .08, .08, .04], n)
    R['Race recode (with detailed Asian and Native Hawaiian other PI)'] = R['Race recode (W, B, AI, API)']
    R['Diagnostic Confirmation'] = pick(rng, ["Positive histology", "Positive exfoliative cytology, no positive histology", "Clinical diagnosis only",
                                              "Radiography without microscopic confirm"], [1 - cfg["cyt"] - .02, cfg["cyt"], .01, .01], n)
    pid = np.arange(pid_start, pid_start + n)
    if share_pids is not None:                                     # a few patients with a second primary in another cohort
        k = min(len(share_pids), n // 50); pid[:k] = share_pids[:k]
    R['Patient ID'] = np.array([f"{p:08d}" for p in pid], dtype=object)
    R = {c: np.array(np.broadcast_to(np.asarray(R[c], dtype=object), (n,))) for c in COLUMNS}   # scalars -> full columns
    # duplicate records (second primary of the same patient in the same export, record number 02)
    dup = rng.choice(n, size=max(1, n // 80), replace=False)
    out = {c: np.concatenate([R[c], R[c][dup]]) for c in COLUMNS}
    out['Record number recode'][n:] = "02"; out['Sequence number'][n:] = "2nd of 2 or more primaries"
    out['Total number of in situ/malignant tumors for patient'][n:] = "02"
    out['Year of diagnosis'][n:] = np.minimum(year[dup] + rng.integers(0, 3, len(dup)), 2023).astype(str)
    return out, pid


def write(path, cols):
    n = len(cols[COLUMNS[0]])
    with open(path, "w", encoding="utf-8", newline="") as f:
        f.write("\t".join(f'"{c}"' for c in COLUMNS) + "\r\n")
        for i in range(n):
            f.write("\t".join(f'"{cols[c][i]}"' for c in COLUMNS) + "\r\n")
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="data/raw", help="output folder (default data/raw); point NEC_DATA_DIR at it")
    ap.add_argument("--seed", type=int, default=20261008)
    ap.add_argument("--scale", type=float, default=1.0, help="multiply the default cohort sizes (NEN 6,000 / AC 30,000 / lung 20,000 / pancreas 15,000)")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    start = 10_000_000; shared = None
    for name in ("NEN", "AC", "LUNG", "PANC"):
        cols, pid = make(name, TUMOURS[name], rng, start, a.scale, shared)
        if shared is None: shared = pid[:400]
        n = write(os.path.join(a.out, FILES[name]), cols)
        start += 5_000_000
        print(f"wrote {FILES[name]}: {n:,} synthetic records")
    with open(os.path.join(a.out, "SYNTHETIC_DATA_README.txt"), "w") as f:
        f.write("SYNTHETIC DATA - generated by tools/make_synthetic_seer_export.py (seed %d, scale %g).\n"
                "Every value is random; these files contain no SEER records. Results computed from them are meaningless\n"
                "and do not reproduce the paper. The pipeline detects this file and runs in non-strict (anchor-informational) mode.\n" % (a.seed, a.scale))
    print(f"synthetic exports written to {os.path.abspath(a.out)}  (marker: SYNTHETIC_DATA_README.txt)")


if __name__ == "__main__":
    main()
