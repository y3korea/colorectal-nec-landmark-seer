# -*- coding: utf-8 -*-
"""Publication-safe copies of the aggregate result JSON files (SEER small-cell rule).

SEER Research Data Use Agreement (November 2025), item 14: "Statistics about any geography (county, state,
registry) or demographic information (age, sex, race/ethnicity) based on counts 1 to 4 must be suppressed in all
publications", and no information on an individual may be published.

The analysis scripts write their outputs unchanged (results/*.json and figures/output/manuscript_final/
figure_data.json); those files are working outputs for the user's own SEER-licensed environment. This tool reads
them and writes copies in which every such cell is replaced by the marker "suppressed (SEER small-cell rule)". Only these copies may be
published; reference_results/ in this repository holds exactly these copies, and tests/compare_with_reference.py
skips the marked fields (and compares everything else).

    python3 tools/suppress_small_cells.py                      # results/ + figure_data.json -> results/publication_safe/
    python3 tools/suppress_small_cells.py --out DIR            # another output folder
    python3 tools/suppress_small_cells.py --results reference_results \\
        --figure-data reference_results/figure_data.json --out reference_results   # (maintainers) re-apply in place

Rules (all cells are aggregate; the rules are specific to the outputs of this pipeline):
  A04.json   per_code_table: a morphology code with 1-4 patients keeps only its patient count n (histology, no
             demographic content); its stage (n_distant, pct_distant), diagnosis-period shares, median age, % male
             and % rectum are suppressed, because with 1-4 patients these describe individual patients. In every
             other row, % male is suppressed when the number of men or of women it implies (n x % / 100, or n minus
             that) is 1-4 (e.g. a histology code whose % male implies 1-4 men or women). Complementary suppression: within an additive
             group whose total row is published (small-cell family 8041-8045; all NEN codes), if the sex split of the
             suppressed members could be recovered by subtraction from the rounded total to within +/-2.5 patients
             and the recovered number of men or women would be 1-4, % male of the smallest unsuppressed member is
             suppressed as well (repeated until nothing can be recovered).
  A06.json   race_income_adjusted.race_levels / A05.json strata (age, sex, race/ethnicity): any level or stratum
             with 1-4 patients (or 1-4 at risk at a landmark) is suppressed (none in the authors' run).
  A08.json   n.race_raw: race/ethnicity categories with 1-4 patients are suppressed; because the collapsed
             category "Other" (n.race5.Other) is published, every category collapsed into it is suppressed
             (complementary suppression). n.income_bands (county attribute): bands with 1-4 patients and, if any,
             the next smallest band. coef_*: a landmark-specific hazard ratio (HR, CI, P) of a race/ethnicity level
             whose number of patients at risk at that landmark is 1-4 is suppressed (at-risk counts are taken from
             the derived cohort NEC_WORK_DIR/cohort_NEN_all.csv when it is available; without it, every landmark > 0
             coefficient of a level with fewer than 10 patients at landmark 0 is suppressed, conservatively).
  A10.json   information on individual patients, withheld exactly as in the Supplement (Tables S5, S34, S35, S37):
             survival-day minimum and maximum (sec1 consistency_*: sdays_min, sdays_max) and the maximum days to
             treatment of every group (sec4 groups[*].max); categories the Supplement pools because they hold
             single patients, together with their pooled partners so that no component can be recovered:
             days-to-treatment bins "91-180" and ">180" (published: > 90 d), sequence "4th"/"5th of ... primaries"
             and tumour counts 05/07 (all stages; published: 4th or later, 5 or more); in seq_by_stage (sequence x
             stage; the Supplement shows only the distant column, Table S35) every cell outside the distant column
             of the later-primary rows (2nd, 3rd, 4th, 5th), including zeros, because single-patient cells there
             could otherwise be recovered from the published row totals and the later-primary totals by stage
             (first_primary_by_stage); 3 and 4 involved metastatic sites (published: 3 or 4), the systemic-therapy/surgery
             sequences "both before and after", "surgery both before and after" and "sequence unknown"
             (published: other or unknown sequence) and the no-surgery reasons "recommended, unknown if performed"
             and "patient refused" (published together); sec8_comparators.shared_*_seq: the multiple-primary sequence
             of the few patients present in both the NEC and a comparator cohort (the overlap counts themselves,
             which carry no demographic information, are kept).
  figure_data.json  figures.SuppFigS1.{NEN_flow, comparator_flow.*}: an eligibility step on age (age >= 18) that
             removes 1-4 records or patients is suppressed together with the following steps whose counts are
             identical (otherwise the removed number could be recovered); the same for A_steps.
Afterwards a heuristic audit scans every copy for remaining 1-4 counts attached to demographic or geographic labels
and exits with code 2 if it finds any (review them, add a rule, re-run). Free text (code, Markdown reports written by
the scripts, documentation) is checked separately by tests/check_small_cell_text.py.
"""
import argparse
import copy
import glob
import json
import math
import os
import re
import shutil
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "python"))

SUPPRESSED = "suppressed (SEER small-cell rule)"
META_KEY = "_suppression"
DUA_NOTE = ("Cells based on 1-4 patients are suppressed in accordance with the SEER Research Data Use Agreement "
            "(item 14: statistics about geography or demographic information based on counts 1 to 4).")
SMALL = range(1, 5)                                     # 1..4 patients
CHAR_FIELDS_A04 = ["n_distant", "pct_distant", "pct_2016_2023", "pct_2021_2023", "median_age", "pct_male", "pct_rectum"]
ADDITIVE_DEMOGRAPHIC_A04 = ["pct_male"]                 # shares that could be recovered by subtraction
A04_GROUPS = {"small-cell family (8041-8045)": ["8041", "8042", "8043", "8044", "8045"]}
DEMOG_STRATA = re.compile(r"^(age|sex|race|race/ethnicity|ethnicity)$", re.I)
DEMOG_LABEL = re.compile(r"\b(age|male|female|sex|race|ethnic|hispanic|white|black|asian|pacific|indian|alaska|aian|"
                         r"api|nhw|registry|county|state|income|rural|urban)\b", re.I)


def is_small(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool) and float(v).is_integer() and int(v) in SMALL


class Log:
    def __init__(self):
        self.items = {}

    def add(self, fname, path, reason):
        self.items.setdefault(fname, []).append({"path": path, "reason": reason})


# ------------------------------------------------------------------------------------------------ rules
def _num(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool) and not (isinstance(v, float) and math.isnan(v))


def _implied_small(n, pct):
    """Number of men (and women) implied by a rounded % male among n patients; True if either can be 1-4."""
    if not (_num(n) and _num(pct)) or n < 1:
        return False
    half = n * 0.0005                                  # pct is rounded to 0.1 percentage points
    lo, hi = n * pct / 100 - half, n * pct / 100 + half
    cand = {k for k in range(int(math.floor(lo)), int(math.ceil(hi)) + 1) if lo - 1e-9 <= k <= hi + 1e-9 and 0 <= k <= n}
    return any(k in SMALL or (n - k) in SMALL for k in cand)


def rule_a04(d, log):
    rows = d.get("per_code_table")
    if not isinstance(rows, list):
        return
    by_code = {str(r.get("code")): r for r in rows if isinstance(r, dict)}
    for i, r in enumerate(rows):
        if not isinstance(r, dict):
            continue
        if is_small(r.get("n")):
            for f in CHAR_FIELDS_A04:
                if f in r and r[f] != SUPPRESSED:
                    r[f] = SUPPRESSED
                    log.add("A04", f"per_code_table[{i}] (code {r.get('code')}, n = {r['n']})/{f}", "1-4 patients")
        elif _implied_small(r.get("n"), r.get("pct_male")):
            r["pct_male"] = SUPPRESSED
            log.add("A04", f"per_code_table[{i}] (code {r.get('code')}, n = {r['n']})/pct_male",
                    "number of men or women implied by % male is 1-4")
    groups = dict(A04_GROUPS)
    groups["all NEN codes"] = [c for c in by_code if re.fullmatch(r"\d{4}", c)]
    for gname, members in groups.items():
        tot = by_code.get(gname)
        if not isinstance(tot, dict) or not _num(tot.get("pct_male")) or not _num(tot.get("n")):
            continue
        while True:
            mem = [by_code[c] for c in members if c in by_code and _num(by_code[c].get("n")) and by_code[c]["n"] > 0]
            hidden = [r for r in mem if r.get("pct_male") == SUPPRESSED]
            shown = [r for r in mem if _num(r.get("pct_male"))]
            if not hidden or not shown:
                break
            n_h = sum(r["n"] for r in hidden)
            est = tot["n"] * tot["pct_male"] / 100 - sum(r["n"] * r["pct_male"] / 100 for r in shown)
            half = 0.0005 * (tot["n"] + sum(r["n"] for r in shown))
            cand = [k for k in range(0, n_h + 1) if est - half - 1e-9 <= k <= est + half + 1e-9]
            if not (half < 2.5 and any(k in SMALL or (n_h - k) in SMALL for k in cand)):
                break
            comp = min((r for r in shown if r["n"] >= 5), key=lambda r: r["n"], default=None)
            if comp is None:
                break
            comp["pct_male"] = SUPPRESSED
            log.add("A04", f"per_code_table[{rows.index(comp)}] (code {comp.get('code')}, n = {comp['n']})/pct_male",
                    f"complementary suppression ({gname}; codes {', '.join(str(r.get('code')) for r in hidden)} suppressed)")


def _suppress_landmarks(obj, where, fname, log):
    """Landmark dict {lm: {est, ci, n_at_risk}} and curve arrays: suppress entries with 1-4 at risk."""
    lm = obj.get("landmark")
    if isinstance(lm, dict):
        for k, v in lm.items():
            if isinstance(v, dict) and is_small(v.get("n_at_risk")):
                for f in list(v):
                    v[f] = SUPPRESSED
                log.add(fname, f"{where}/landmark/{k}", "1-4 at risk in a demographic stratum")
    cur = obj.get("curve_0_18")
    if isinstance(cur, dict) and isinstance(cur.get("n_at_risk"), list):
        for j, nr in enumerate(cur["n_at_risk"]):
            if is_small(nr):
                for f in cur:
                    if isinstance(cur[f], list) and j < len(cur[f]):
                        cur[f][j] = SUPPRESSED
                log.add(fname, f"{where}/curve_0_18[{j}]", "1-4 at risk in a demographic stratum")


def rule_a05(d, log):
    for i, s in enumerate(d.get("strata") or []):
        if not isinstance(s, dict) or not DEMOG_STRATA.match(str(s.get("variable", ""))):
            continue
        where = f"strata[{i}] ({s.get('variable')} {s.get('level')})"
        if is_small(s.get("n")):
            for f in list(s):
                if f not in ("variable", "level", "B", "seed", "seed_spawn_index"):
                    s[f] = SUPPRESSED
            log.add("A05", where, "stratum with 1-4 patients")
            continue
        _suppress_landmarks(s, where, "A05", log)
        c3 = s.get("conditional_on_3mo_survival")
        if isinstance(c3, dict):
            _suppress_landmarks(c3, where + "/conditional_on_3mo_survival", "A05", log)


def rule_a06(d, log):
    rl = (d.get("race_income_adjusted") or {}).get("race_levels")
    for i, r in enumerate(rl or []):
        if isinstance(r, dict) and is_small(r.get("N")):
            r["N"] = SUPPRESSED
            log.add("A06", f"race_income_adjusted/race_levels[{i}] ({r.get('race5')})/N", "1-4 patients")


def rule_a08(d, log):
    n = d.get("n")
    if not isinstance(n, dict):
        return
    raw = n.get("race_raw")
    if isinstance(raw, dict):
        small = [k for k, v in raw.items() if is_small(v)]
        for k in small:
            raw[k] = SUPPRESSED
            log.add("A08", f"n/race_raw/{k}", "1-4 patients")
        collapsed = n.get("race_collapsed_to_other") or []
        if any(k in collapsed for k in small):
            for k in collapsed:
                if k in raw and raw[k] != SUPPRESSED:
                    raw[k] = SUPPRESSED
                    log.add("A08", f"n/race_raw/{k}",
                            "complementary suppression (collapsed into the published category n/race5/Other)")
    bands = n.get("income_bands")
    if isinstance(bands, dict):
        small = [k for k, v in bands.items() if is_small(v)]
        if small:
            for k in small:
                bands[k] = SUPPRESSED
                log.add("A08", f"n/income_bands/{k}", "1-4 patients (county attribute)")
            rest = sorted((v, k) for k, v in bands.items() if v != SUPPRESSED)
            if rest:
                bands[rest[0][1]] = SUPPRESSED
                log.add("A08", f"n/income_bands/{rest[0][1]}", "complementary suppression (tertile totals are published)")

RACE5_MAP = {"Non-Hispanic White": "White", "Non-Hispanic Black": "Black", "Hispanic (All Races)": "Hispanic",
             "Non-Hispanic Asian or Pacific Islander": "API"}          # every other category -> "Other" (as in A08)
NEC_CODES = {"8013", "8041", "8042", "8043", "8044", "8045"}


def _race5_at_risk(cohort_csv):
    """{landmark: {race5 level: n at risk}} for distant main NEC (alive and under follow-up at the landmark), or None."""
    if not cohort_csv or not os.path.exists(cohort_csv):
        return None
    import csv
    rows = []
    with open(cohort_csv, newline="", encoding="utf-8") as f:
        for r in csv.DictReader(f):
            if r.get("hist") in NEC_CODES and r.get("stage") == "Distant":
                rows.append((float(r["tm"]), RACE5_MAP.get(r.get("race"), "Other")))
    def at(s):
        out = {}
        for tm, lev in rows:
            if tm > s:
                out[lev] = out.get(lev, 0) + 1
        return out
    return at


def rule_a08_coef(d, log, cohort_csv=None):
    n5 = (d.get("n") or {}).get("race5") if isinstance(d.get("n"), dict) else None
    at = _race5_at_risk(cohort_csv)
    for key, rows in d.items():
        if not (key.startswith("coef") and isinstance(rows, list)):
            continue
        for i, r in enumerate(rows):
            if not (isinstance(r, dict) and str(r.get("term", "")).startswith("race5")):
                continue
            lev, lm = str(r["term"])[5:], r.get("landmark", 0)
            if at is not None:
                hit = at(lm).get(lev, 0) in SMALL
                why = "race/ethnicity level with 1-4 patients at risk at this landmark"
            else:
                n0 = n5.get(lev) if isinstance(n5, dict) else None
                hit = _num(lm) and lm > 0 and _num(n0) and n0 < 10
                why = "race/ethnicity level with < 10 patients at landmark 0 (at-risk count unavailable; conservative)"
            if hit:
                for f in ("HR", "lo", "hi", "p"):
                    if f in r and r[f] != SUPPRESSED:
                        r[f] = SUPPRESSED
                        log.add("A08", f"{key}[{i}] (landmark {lm}, {r['term']})/{f}", why)


def _flow_steps(steps, where, log, fname="figure_data"):
    """steps: ordered dict step -> [records, patients]. Suppress an age step that removes 1-4 records/patients and the
    following steps with identical counts."""
    keys = list(steps)
    for j in range(1, len(keys)):
        k = keys[j]
        if not k.lower().startswith("age"):
            continue
        prev, cur = steps[keys[j - 1]], steps[k]
        if not (isinstance(prev, list) and isinstance(cur, list)) or SUPPRESSED in prev or SUPPRESSED in cur:
            continue
        removed = [a - b for a, b in zip(prev, cur)]
        if not any(r in SMALL for r in removed):
            continue
        ref = list(cur)
        steps[k] = SUPPRESSED
        log.add(fname, f"{where}/{k}", "age-based exclusion of 1-4 records/patients")
        for k2 in keys[j + 1:]:
            if steps[k2] == ref:
                steps[k2] = SUPPRESSED
                log.add(fname, f"{where}/{k2}", "complementary suppression (same counts as the suppressed age step)")
            else:
                break


def _suppfig_s1(d):
    figs = d.get("figures") if isinstance(d, dict) else None
    s1 = figs.get("SuppFigS1") if isinstance(figs, dict) else None
    return s1 if isinstance(s1, dict) else {}


def rule_figure_data(d, log):
    s1 = _suppfig_s1(d)
    if isinstance(s1.get("NEN_flow"), dict):
        _flow_steps(s1["NEN_flow"], "figures/SuppFigS1/NEN_flow", log)
    for coh, steps in (s1.get("comparator_flow") or {}).items():
        if isinstance(steps, dict):
            _flow_steps(steps, f"figures/SuppFigS1/comparator_flow/{coh}", log)
    for i, st in enumerate(s1.get("A_steps") or []):
        if isinstance(st, dict) and str(st.get("step", "")).lower().startswith("age") and is_small(st.get("excluded")):
            st["excluded"] = SUPPRESSED; st["remaining"] = SUPPRESSED
            log.add("figure_data", f"figures/SuppFigS1/A_steps[{i}]", "age-based exclusion of 1-4")


A10_POOLED = {   # (section, table, label column, categories pooled in the Supplement; all are suppressed)
    ("sec4_treatment_timing", "bins_distant_chemo", "tt_days (distant, chemo yes)", ("91-180", ">180")),
    ("sec5_multiple_primaries", "all stages_seq", "Sequence number", ("4th of 4 or more primaries", "5th of 5 or more primaries")),
    ("sec5_multiple_primaries", "all stages_ntum", "Total number of in-situ/malignant tumours for patient", ("05", "06", "07", "08", "09")),
    ("sec7_distant_mets_treatment", "n_sites_2010plus_known4", "number of the 4 recorded sites involved (2010+, all 4 known)", ("3", "4")),
    ("sec7_distant_mets_treatment", "sys_surg_seq", "Systemic therapy / surgery sequence (distant)",
     ("Systemic therapy both before and after surgery", "Sequence unknown", "Surgery both before and after systemic therapy")),
    ("sec7_distant_mets_treatment", "no_surg_reason_surgno", "Reason no surgery (distant, surgery = no)",
     ("Recommended, unknown if performed", "Recommended but not performed, patient refused")),
}
A10_POOLED_XTAB = {   # cross-tabulations {category: {column: n}} of the same pooled categories
    ("sec7_distant_mets_treatment", "sys_surg_seq_by_surgery"):
        ("Systemic therapy both before and after surgery", "Sequence unknown", "Surgery both before and after systemic therapy"),
    ("sec7_distant_mets_treatment", "no_surg_reason_by_surgery"):
        ("Recommended, unknown if performed", "Recommended but not performed, patient refused"),
}
WHY_IND = "information on an individual patient (withheld in the Supplement)"
WHY_POOL = "category pooled in the Supplement (single patients; components suppressed so none can be recovered)"


def rule_a10(d, log):
    """Information on individual patients, withheld exactly as in the Supplement."""
    def sup(obj, k, path, why):
        if isinstance(obj, dict) and k in obj and obj[k] != SUPPRESSED and obj[k] is not None:
            obj[k] = SUPPRESSED
            log.add("A10", path, why)
    s1 = d.get("sec1_survival_fields")
    if isinstance(s1, dict):
        for c in ("consistency_all", "consistency_distant"):
            for k in ("sdays_min", "sdays_max"):
                sup(s1.get(c), k, f"sec1_survival_fields/{c}/{k}", WHY_IND)
    s4 = d.get("sec4_treatment_timing")
    if isinstance(s4, dict):
        for i, g in enumerate(s4.get("groups") or []):
            sup(g, "max", f"sec4_treatment_timing/groups[{i}]/max", WHY_IND)
    for sec, tab, col, cats in A10_POOLED:
        rows = (d.get(sec) or {}).get(tab) if isinstance(d.get(sec), dict) else None
        for i, r in enumerate(rows or []):
            if isinstance(r, dict) and str(r.get(col)) in cats:
                for f in list(r):
                    if f != col:
                        sup(r, f, f"{sec}/{tab}[{i}] ({r.get(col)})/{f}", WHY_POOL)
    for (sec, tab), cats in A10_POOLED_XTAB.items():
        x = (d.get(sec) or {}).get(tab) if isinstance(d.get(sec), dict) else None
        if isinstance(x, dict):
            for cat in cats:
                if isinstance(x.get(cat), dict):
                    for col in list(x[cat]):
                        if x[cat][col] != 0:
                            sup(x[cat], col, f"{sec}/{tab}/{cat}/{col}", WHY_POOL)
    # sequence number x stage: the Supplement shows the distant column only (Table S35). Outside the distant column,
    # every cell of the later-primary rows (2nd-5th, zeros included) is suppressed: suppressing only the 1-4 cells
    # would not protect them, because they could be recovered from the published row totals (all stages_seq, Table
    # S35) and the later-primary totals by stage (first_primary_by_stage). Only the stage totals of later primaries
    # remain derivable (all >= 5 in the authors' run).
    x = (d.get("sec5_multiple_primaries") or {}).get("seq_by_stage") if isinstance(d.get("sec5_multiple_primaries"), dict) else None
    if isinstance(x, dict):
        for cat, cols in x.items():
            if isinstance(cols, dict) and re.match(r"(2nd|3rd|[4-9]th)\b", str(cat)):
                for col in list(cols):
                    if col != "Distant":
                        sup(cols, col, f"sec5_multiple_primaries/seq_by_stage/{cat}/{col}",
                            "sequence x stage cell outside the distant column (later-primary rows; not shown in the "
                            "Supplement; suppressed together so that single-patient cells cannot be recovered)")
    sec = d.get("sec8_comparators")
    if isinstance(sec, dict):
        for k, v in sec.items():
            if k.startswith("shared_") and k.endswith("_seq") and isinstance(v, dict):
                for cat, cnt in v.items():
                    if is_small(cnt):
                        v[cat] = SUPPRESSED
                        log.add("A10", f"sec8_comparators/{k}/{cat}", "characteristic of 1-4 individual patients (shared between two cohorts)")


RULES = {"A04": rule_a04, "A05": rule_a05, "A06": rule_a06, "A08": rule_a08, "A10": rule_a10, "figure_data": rule_figure_data}


# ------------------------------------------------------------------------------------------------ audit
def audit(name, d):
    """Heuristic scan for remaining 1-4 counts attached to demographic / geographic labels."""
    out = []

    def label_of(x):
        return " ".join(str(v) for k, v in x.items() if isinstance(v, str) and k != META_KEY) if isinstance(x, dict) else ""

    def walk(x, path, ctx):
        if isinstance(x, dict):
            lab = label_of(x)
            nn = next((x[k] for k in ("n", "N") if k in x), None)
            if "pct_male" in x and not is_small(nn) and _implied_small(nn, x.get("pct_male")):
                out.append((f"{path}/pct_male  [n = {nn}; implied number of men or women 1-4]", x["pct_male"]))
            if is_small(nn):   # a row of 1-4 patients that still carries a demographic statistic
                for k, v in x.items():
                    if DEMOG_LABEL.search(k.replace("_", " ")) and v not in (None, SUPPRESSED) and not isinstance(v, str):
                        out.append((f"{path}/{k}  [row n = {nn}]", v))
            for k, v in x.items():
                if k == META_KEY:
                    continue
                kp = f"{path}/{k}"
                if is_small(v) and DEMOG_LABEL.search(k) and not re.search(r"seed|index|df|digits|version|minutes|B$", k, re.I):
                    out.append((kp, v))
                elif is_small(v) and k in ("n", "N", "count") and DEMOG_LABEL.search(lab) and not re.search(r"covariates|formula|label", lab):
                    out.append((kp + f"  [{lab[:60]}]", v))
                walk(v, kp, ctx)
        elif isinstance(x, list):
            for i, v in enumerate(x):
                walk(v, f"{path}[{i}]", ctx)

    walk(d, name, "")
    # flows: age steps
    s1 = _suppfig_s1(d)
    for coh, steps in ([("NEN_flow", s1.get("NEN_flow"))] + list((s1.get("comparator_flow") or {}).items())):
        if isinstance(steps, dict):
            keys = list(steps)
            for j in range(1, len(keys)):
                a, b = steps[keys[j - 1]], steps[keys[j]]
                if keys[j].lower().startswith("age") and isinstance(a, list) and isinstance(b, list) and \
                        SUPPRESSED not in a and SUPPRESSED not in b and any((p - q) in SMALL for p, q in zip(a, b)):
                    out.append((f"{name}/figures/SuppFigS1/{coh}/{keys[j]}", [p - q for p, q in zip(a, b)]))
    return out


# ------------------------------------------------------------------------------------------------ main
def main():
    import config  # noqa: E402  (paths from NEC_* environment variables)
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--results", default=config.RESULTS_DIR, help="folder with the result JSON files (default: NEC_RESULTS_DIR)")
    ap.add_argument("--figure-data", default=os.path.join(config.MANUSCRIPT_FIG_DIR, "figure_data.json"),
                    help="manuscript figure_data.json (default: NEC_FIG_DIR/manuscript_final/figure_data.json)")
    ap.add_argument("--out", default=None, help="output folder (default: <results>/publication_safe)")
    ap.add_argument("--cohort", default=os.path.join(config.WORK_DIR, "cohort_NEN_all.csv"),
                    help="derived NEN cohort (patient-level, read locally only) for landmark at-risk counts by race/ethnicity "
                         "(default: NEC_WORK_DIR/cohort_NEN_all.csv; if absent a conservative rule is used)")
    ap.add_argument("--no-audit", action="store_true", help="do not fail on heuristic audit findings")
    a = ap.parse_args()
    out_dir = a.out or os.path.join(a.results, "publication_safe")
    os.makedirs(out_dir, exist_ok=True)
    files = {os.path.basename(p)[:-5]: p for p in sorted(glob.glob(os.path.join(a.results, "*.json")))}
    if os.path.exists(a.figure_data):
        files["figure_data"] = a.figure_data
    if not files:
        sys.exit(f"no result JSON files in {a.results}; run the pipeline first")
    log, findings = Log(), []
    for name, src in files.items():
        with open(src, encoding="utf-8") as f:
            d = json.load(f)
        before = copy.deepcopy(d)
        if name in RULES:
            RULES[name](d, log)
        if name == "A08":
            rule_a08_coef(d, log, a.cohort)
        changed = d != before
        if log.items.get(name):
            meta = d.get(META_KEY) if isinstance(d.get(META_KEY), dict) else {}
            fields = {x["path"]: x["reason"] for x in meta.get("fields", [])}
            for x in log.items[name]:
                fields[x["path"]] = x["reason"]
            d[META_KEY] = {"marker": SUPPRESSED, "note": DUA_NOTE, "tool": "tools/suppress_small_cells.py",
                           "fields": [{"path": k, "reason": v} for k, v in fields.items()]}
            changed = True
        dst = os.path.join(out_dir, name + ".json")
        if changed:
            with open(dst, "w", encoding="utf-8") as f:
                json.dump(d, f, indent=1, ensure_ascii=False, default=float)
        elif os.path.abspath(src) != os.path.abspath(dst):
            shutil.copyfile(src, dst)
        findings += audit(name, d)
    n = sum(len(v) for v in log.items.values())
    print(f"{len(files)} JSON files -> {out_dir}; {n} fields newly suppressed")
    for name, items in log.items.items():
        for x in items:
            print(f"  {name}: {x['path']}  ({x['reason']})")
    if findings:
        print(f"\nAUDIT: {len(findings)} possible small demographic/geographic cell(s) not covered by a rule:")
        for p, v in findings:
            print(f"  {p} = {v}")
        if not a.no_audit:
            sys.exit(2)
    else:
        print("audit: no remaining 1-4 cells attached to demographic or geographic labels")


if __name__ == "__main__":
    main()
