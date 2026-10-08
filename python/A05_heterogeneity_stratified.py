# -*- coding: utf-8 -*-
"""A05 (part 2/2, revised after audit): population-heterogeneity (mixture) check of the distant main-NEC landmark plateau.
Stratum-specific landmark curves (12-month cancer-death probability, Aalen-Johansen), bootstrap CIs,
6-0 (and 6-3 / 12-3 for treatment strata) differences, plateau end, crude interval rates, forest + multi-panel figures,
plus (revision) a supplementary immortal-time test of the chemotherapy split using time to first treatment (tt_days):
(a) landmark-specific 'chemotherapy initiated by the landmark' classification, (b) time-dependent-exposure crude interval rates.
Merges results/A05_std.json (from A05_heterogeneity_std.R, run first) into results/A05.json and writes results/A05.md.
Run: Rscript R/A05_heterogeneity_std.R && python3 python/A05_heterogeneity_stratified.py
"""
import sys, os, json, time
t_start = time.time()
import config
import numpy as np, pandas as pd, lm
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
RES = config.RES
FIG = config.FIGS
os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)
SEED = 20261007; LMS = np.arange(0, 19); WIN = 12; REP_LM = [0, 3, 6, 12]
# independent bootstrap streams: one SeedSequence child per analysis unit (strata, conditional analyses, pooled, sensitivity, time-dependent)
SS = np.random.SeedSequence(SEED); CHILDREN = SS.spawn(64); _spawn = [0]
def next_stream():
    j = _spawn[0]; _spawn[0] += 1; return CHILDREN[j], j

# ---------------------------------------------------------------- data + reproduction check
nen = lm.load("NEN_all"); nec = lm.nec_main(nen); d = nec[nec.stage == "Distant"].copy().reset_index(drop=True)
d["surg"] = lm.surgery(d).to_numpy()
k0, ev0, K0 = lm.prep(d)
c_chk, _, n_chk = lm.curves(k0, ev0, K0, np.array([0, 3, 6, 12, 24]), WIN)
pl_chk = float(lm.hinge_many(lm.curves(k0, ev0, K0, LMS, WIN)[0][None, :])[0])
cr_chk = lm.crude_rates(d).rate.to_numpy()
exp_c = np.array([72.29, 71.62, 72.33, 60.03, 30.51]); exp_n = np.array([908, 553, 396, 169, 48])
exp_cr = np.array([0.1435, 0.0909, 0.1344, 0.1204, 0.1231, 0.0840, 0.0646])
print("REPRODUCTION: 12-mo cancer-death % at LM 0/3/6/12/24 =", np.round(c_chk, 2), "n at risk", n_chk.astype(int), "plateau", pl_chk, "crude", np.round(cr_chk, 4))
ok = np.allclose(c_chk, exp_c, atol=0.01) and np.array_equal(n_chk.astype(int), exp_n) and pl_chk == 6.75 and np.allclose(cr_chk, exp_cr, atol=5e-5) \
     and (d.surg == 1).sum() == 354 and (d.surg == 0).sum() == 554 and len(d) == 908
if not ok: config.anchor_fail("REPRODUCTION CHECK FAILED - stopping")
repro = {"landmark_12mo_cancer_death_pct": {"ours": np.round(c_chk, 2).tolist(), "manuscript": exp_c.tolist()},
         "n_at_risk": {"ours": n_chk.astype(int).tolist(), "manuscript": exp_n.tolist()},
         "plateau_end": {"ours": pl_chk, "manuscript": 6.75}, "crude_rates": {"ours": np.round(cr_chk, 4).tolist(), "manuscript": exp_cr.tolist()},
         "surgery_yes_no": {"ours": [int((d.surg == 1).sum()), int((d.surg == 0).sum())], "manuscript": [354, 554]}, "ok": bool(ok)}

# ---------------------------------------------------------------- strata definitions
d10 = d[(d.year >= 2010) & d.mets_liver.isin(["Yes", "No"])]
liver_only = (d10.mets_liver == "Yes") & (d10.mets_bone == "No") & (d10.mets_brain == "No") & (d10.mets_lung == "No")
liver_yes_other_yes = (d10.mets_liver == "Yes") & ((d10.mets_bone == "Yes") | (d10.mets_brain == "Yes") | (d10.mets_lung == "Yes"))
liver_yes_unknown_only = (d10.mets_liver == "Yes") & ~liver_only & ~liver_yes_other_yes   # liver Yes, no other site Yes, >=1 site Unknown/blank
n_unknown_only = int(liver_yes_unknown_only.sum())
STRATA = [  # (variable, level, mask over d, treatment flag)
    ("Surgery (primary site)", "Yes", d.surg == 1, True), ("Surgery (primary site)", "No", d.surg == 0, True),
    ("Liver metastasis (2010+)", "Yes", d.index.isin(d10.index[d10.mets_liver == "Yes"]), False),
    ("Liver metastasis (2010+)", "No", d.index.isin(d10.index[d10.mets_liver == "No"]), False),
    ("Liver-only pattern (2010+)", "Liver-only", d.index.isin(d10.index[liver_only]), False),
    ("Liver-only pattern (2010+)", "Other patterns", d.index.isin(d10.index[~liver_only]), False),
    ("Age", "<65", d.age < 65, False), ("Age", ">=65", d.age >= 65, False),
    ("Primary site", "Colon", d.rectum == 0, False), ("Primary site", "Rectum", d.rectum == 1, False),
    ("Sex", "Male", d.male == 1, False), ("Sex", "Female", d.male == 0, False),
    ("Histology", "Small cell (8041-8045)", d["hist"] != "8013", False), ("Histology", "Large cell (8013)", d["hist"] == "8013", False),
    ("Race", "Non-Hispanic White", d.race == "Non-Hispanic White", False), ("Race", "Other", d.race != "Non-Hispanic White", False),
    ("Chemotherapy recorded", "Yes", d.chemo == 1, True), ("Chemotherapy recorded", "No/unknown", d.chemo != 1, True),
]

def boot_curves(df, B, stream, landmarks=LMS):
    k, ev, KMAX = lm.prep(df); n = len(k); rng = np.random.default_rng(stream)
    pt, _, nrisk = lm.curves(k, ev, KMAX, landmarks, WIN)
    bc = np.empty((B, len(landmarks)))
    for b in range(B):
        w = np.bincount(rng.integers(0, n, n), minlength=n).astype(float)
        bc[b] = lm.curves(k, ev, KMAX, landmarks, WIN, w)[0]
    return pt, nrisk, bc

def ci(v): return [round(float(x), 2) for x in np.nanpercentile(v, [2.5, 97.5])]
def plateau(Y):
    """hinge on each row; NaN entries (landmarks with 0 at risk in a bootstrap resample) are filled with that ROW's mean (revision: was global mean)"""
    Y = np.asarray(Y, float); rm = np.nanmean(Y, axis=1, keepdims=True)
    return lm.hinge_many(np.where(np.isnan(Y), rm, Y))
LO3 = int(np.floor(2 * lm.MD * 3)) + 1  # first half-day index counted at landmark 3

def stratum_result(var, lev, sub, trt):
    n = len(sub); B = config.B(2000) if n >= 100 else config.B(1000)
    stream, j = next_stream(); pt, nrisk, bc = boot_curves(sub, B, stream)
    i = {s: int(np.where(LMS == s)[0][0]) for s in LMS}
    pl_pt = float(plateau(pt[None, :])[0]); pl_b = plateau(bc)
    cr = lm.crude_rates(sub, cuts=(0, 3, 6, 9, 12))
    r = {"variable": var, "level": lev, "n": n, "B": B, "seed": SEED, "seed_spawn_index": j,
         "deaths_cancer": int((sub.ev == 1).sum()), "deaths_other": int((sub.ev == 2).sum()),
         "landmark": {str(s): {"est": round(float(pt[i[s]]), 2), "ci": ci(bc[:, i[s]]), "n_at_risk": int(nrisk[i[s]])} for s in REP_LM},
         "curve_0_18": {"est": np.round(pt, 2).tolist(), "lo": np.round(np.nanpercentile(bc, 2.5, axis=0), 2).tolist(),
                        "hi": np.round(np.nanpercentile(bc, 97.5, axis=0), 2).tolist(), "n_at_risk": nrisk.astype(int).tolist()},
         "diff_6_0": {"est": round(float(pt[i[6]] - pt[i[0]]), 2), "ci": ci(bc[:, i[6]] - bc[:, i[0]])},
         "diff_12_0": {"est": round(float(pt[i[12]] - pt[i[0]]), 2), "ci": ci(bc[:, i[12]] - bc[:, i[0]])},
         "plateau_end": {"est": pl_pt, "ci": ci(pl_b), "boot_median": round(float(np.median(pl_b)), 2),
                         "frac_boot_plateau_ge_3": round(float(np.mean(pl_b >= 3)), 3)},
         "n_at_risk_lm6_lt_100": bool(nrisk[i[6]] < 100),
         "crude_rates": [{"interval": f"{int(a)}-{int(b)}", "deaths": int(x), "person_months": round(float(p), 1), "rate": round(float(rt), 4)}
                         for a, b, x, p, rt in zip(cr.start, cr.end, cr.deaths, cr.person_months, cr.rate)]}
    r["crude_rr_6_9_vs_3_6"] = round(float(cr.rate[2] / cr.rate[1]), 2) if cr.rate[1] > 0 else None
    r["crude_dip_3_6_then_rise_6_9"] = bool(cr.rate[1] < cr.rate[0] and cr.rate[2] > cr.rate[1])
    store = {"pt": pt, "bc": bc}
    if trt:  # condition on surviving 3 months: resample among patients still at risk at landmark 3
        k, ev, KMAX = lm.prep(sub); sub3 = sub[k >= LO3]
        stream3, j3 = next_stream(); pt3, nrisk3, bc3 = boot_curves(sub3, B, stream3, landmarks=np.array([3, 6, 9, 12]))
        r["conditional_on_3mo_survival"] = {"n": len(sub3), "B": B, "seed_spawn_index": j3,
            "landmark": {str(s): {"est": round(float(pt3[jj]), 2), "ci": ci(bc3[:, jj]), "n_at_risk": int(nrisk3[jj])} for jj, s in enumerate([3, 6, 9, 12])},
            "diff_6_3": {"est": round(float(pt3[1] - pt3[0]), 2), "ci": ci(bc3[:, 1] - bc3[:, 0])},
            "diff_12_3": {"est": round(float(pt3[3] - pt3[0]), 2), "ci": ci(bc3[:, 3] - bc3[:, 0])},
            "check_same_point_estimates_as_unconditional": bool(np.allclose(pt3, pt[[3, 6, 9, 12]], equal_nan=True))}
        store["bc3"] = bc3; store["pt3"] = pt3
    return r, store

results = []; boot_store = {}
for var, lev, mask, trt in STRATA:
    r, st = stratum_result(var, lev, d[np.asarray(mask)], trt); results.append(r); boot_store[(var, lev)] = st
    pt = st["pt"]; print(f"{var:28s} {lev:24s} n={r['n']:4d} B={r['B']} spawn={r['seed_spawn_index']:2d} LM0/3/6/12 = {pt[0]:.1f}/{pt[3]:.1f}/{pt[6]:.1f}/{pt[12]:.1f}  d6-0={r['diff_6_0']['est']:+.1f} {r['diff_6_0']['ci']}  plateau={r['plateau_end']['est']} {r['plateau_end']['ci']}  nrisk6={r['landmark']['6']['n_at_risk']}")

# sensitivity (audit minor issue): 'other patterns' excluding the liver-Yes patients whose only non-liver information is Unknown
sens_other, _ = stratum_result("Liver-only pattern (2010+), sensitivity", f"Other patterns excl. {n_unknown_only} liver-Yes/Unknown-only", d[np.asarray(d.index.isin(d10.index[~liver_only & ~liver_yes_unknown_only]))], False)
print(f"sensitivity other-patterns excl. unknown-only: n={sens_other['n']} LM0/6 = {sens_other['landmark']['0']['est']}/{sens_other['landmark']['6']['est']} d6-0={sens_other['diff_6_0']['est']} {sens_other['diff_6_0']['ci']} plateau={sens_other['plateau_end']['est']}")

# pooled reference curve (B=2000, own stream)
stream_p, j_p = next_stream(); pt_all, nrisk_all, bc_all = boot_curves(d, config.B(2000), stream_p)
pooled = {"n": int(len(d)), "B": config.B(2000), "seed_spawn_index": j_p, "curve_0_18": {"est": np.round(pt_all, 2).tolist(), "lo": np.round(np.percentile(bc_all, 2.5, axis=0), 2).tolist(),
          "hi": np.round(np.percentile(bc_all, 97.5, axis=0), 2).tolist(), "n_at_risk": nrisk_all.astype(int).tolist()},
          "diff_6_0": {"est": round(float(pt_all[6] - pt_all[0]), 2), "ci": ci(bc_all[:, 6] - bc_all[:, 0])},
          "plateau_end": {"est": float(plateau(pt_all[None, :])[0]), "ci": ci(plateau(bc_all))}}

# between-stratum contrast of the 6-0 difference (independent bootstrap streams, paired by draw index, first min(B) draws)
contrasts = []
vars_ = list(dict.fromkeys(v for v, _, _, _ in STRATA))
for v in vars_:
    (l1, l2) = [lev for vv, lev, _, _ in STRATA if vv == v]
    a, b = boot_store[(v, l1)], boot_store[(v, l2)]; Bm = min(len(a["bc"]), len(b["bc"]))
    da = a["bc"][:Bm, 6] - a["bc"][:Bm, 0]; db = b["bc"][:Bm, 6] - b["bc"][:Bm, 0]
    e = {"variable": v, "contrast": f"{l1} minus {l2}", "B_paired": Bm,
         "diff_of_6_0_differences": {"est": round(float((a["pt"][6] - a["pt"][0]) - (b["pt"][6] - b["pt"][0])), 2), "ci": ci(da - db)},
         "diff_of_LM0": {"est": round(float(a["pt"][0] - b["pt"][0]), 2), "ci": ci(a["bc"][:Bm, 0] - b["bc"][:Bm, 0])},
         "diff_of_LM6": {"est": round(float(a["pt"][6] - b["pt"][6]), 2), "ci": ci(a["bc"][:Bm, 6] - b["bc"][:Bm, 6])}}
    if "bc3" in a:
        d3a = a["bc3"][:Bm, 1] - a["bc3"][:Bm, 0]; d3b = b["bc3"][:Bm, 1] - b["bc3"][:Bm, 0]
        e["diff_of_6_3_differences_conditional"] = {"est": round(float((a["pt3"][1] - a["pt3"][0]) - (b["pt3"][1] - b["pt3"][0])), 2), "ci": ci(d3a - d3b)}
    contrasts.append(e)

# ---------------------------------------------------------------- supplementary (revision): immortal-time test of the chemotherapy split using tt_days
# tt_days = days from diagnosis to first-course treatment of any modality. In the subset with no primary-site surgery and no radiotherapy,
# the first treatment of a chemotherapy-recorded patient IS the chemotherapy, so tt_days is the chemotherapy start date exactly (primary analysis);
# the full cohort is a secondary analysis in which tt_days bounds the chemotherapy start from below for surgery-first patients.
TD_LM = np.array([0, 1, 2, 3, 6, 9, 12])
def td_setup(sub):
    tt = pd.to_numeric(sub.tt_days, errors="coerce").to_numpy(float)
    tm = sub.tm.to_numpy(float); chemo = (sub.chemo == 1).to_numpy()
    known = chemo & ~np.isnan(tt); med = float(np.median(tt[known]))
    n_imp = int((chemo & np.isnan(tt)).sum()); n_clip = int((chemo & (tt > sub.sdays.to_numpy(float))).sum())
    tstart = np.full(len(sub), np.inf)
    tstart[chemo] = np.where(np.isnan(tt[chemo]), med, tt[chemo]) / lm.MD        # months; unknown tt among chemo-yes imputed at the median
    tstart[chemo] = np.minimum(tstart[chemo], tm[chemo])                        # treatment recorded after the survival end date -> clip to it (data inconsistency); chemo-no stay at +inf
    return tstart, {"n": int(len(sub)), "n_chemo_yes": int(chemo.sum()), "n_chemo_yes_tt_unknown_imputed": n_imp, "imputed_tt_days": med,
                    "n_chemo_yes_tt_after_survival_end_clipped": n_clip,
                    "share_chemo_yes_initiated_by_LM": {str(s): round(float(np.mean(tstart[chemo] <= s)), 3) for s in TD_LM}}

def td_rates(sub, tstart, cuts=(0, 3, 6, 9, 12)):
    tm = sub.tm.to_numpy(float); ev = sub.ev.to_numpy(int); rows = []
    for a, b in zip(cuts[:-1], cuts[1:]):
        at = tm > a; end = np.minimum(tm, b)
        un = np.clip(np.minimum(end, tstart) - a, 0, None) * at; ex = np.clip(end - np.maximum(a, tstart), 0, None) * at
        dead = (ev == 1) & (tm > a) & (tm <= b); d_ex = int((dead & (tm >= tstart)).sum()); d_un = int(dead.sum()) - d_ex
        rows.append({"interval": f"{a}-{b}", "unexposed": {"deaths": d_un, "person_months": round(float(un.sum()), 1), "rate": round(d_un / un.sum(), 4) if un.sum() > 0 else None},
                     "exposed": {"deaths": d_ex, "person_months": round(float(ex.sum()), 1), "rate": round(d_ex / ex.sum(), 4) if ex.sum() > 0 else None}})
    def rr(g, i1, i0):
        d1, d0 = rows[i1][g]["deaths"], rows[i0][g]["deaths"]; r1, r0 = rows[i1][g]["rate"], rows[i0][g]["rate"]
        if not r1 or not r0 or d1 == 0 or d0 == 0: return None
        se = np.sqrt(1 / d1 + 1 / d0); e = np.log(r1 / r0)
        return {"RR": round(float(np.exp(e)), 2), "ci": [round(float(np.exp(e - 1.96 * se)), 2), round(float(np.exp(e + 1.96 * se)), 2)]}
    return rows, {"unexposed_6_9_vs_3_6": rr("unexposed", 2, 1), "exposed_6_9_vs_3_6": rr("exposed", 2, 1),
                  "unexposed_3_6_vs_0_3": rr("unexposed", 1, 0), "exposed_3_6_vs_0_3": rr("exposed", 1, 0)}

def td_landmark(sub, tstart, B, stream):
    """Landmark-specific exposure: at landmark s, 'initiated' = chemotherapy start <= s among patients at risk at s; one patient-level resample per draw shared by all landmarks/groups."""
    k, ev, KMAX = lm.prep(sub); n = len(k); rng = np.random.default_rng(stream)
    masks = {s: (tstart <= s) for s in TD_LM}
    def one(w):
        out = np.full((2, len(TD_LM)), np.nan); nr = np.zeros((2, len(TD_LM)))
        for jj, s in enumerate(TD_LM):
            for g, m in enumerate([masks[s], ~masks[s]]):
                c, _, nn = lm.curves(k, ev, KMAX, np.array([s]), WIN, w * m); out[g, jj] = c[0]; nr[g, jj] = nn[0]
        return out, nr
    pt, nr = one(np.ones(n)); bc = np.empty((B, 2, len(TD_LM)))
    for b in range(B):
        w = np.bincount(rng.integers(0, n, n), minlength=n).astype(float); bc[b] = one(w)[0]
    L = {str(s): jj for jj, s in enumerate(TD_LM)}
    def grp(g, name):
        return {"group": name, "landmark": {str(s): {"est": round(float(pt[g, L[str(s)]]), 2), "ci": ci(bc[:, g, L[str(s)]]), "n_at_risk": int(nr[g, L[str(s)]])} for s in TD_LM},
                "diff_6_1": {"est": round(float(pt[g, L["6"]] - pt[g, L["1"]]), 2), "ci": ci(bc[:, g, L["6"]] - bc[:, g, L["1"]])},
                "diff_6_3": {"est": round(float(pt[g, L["6"]] - pt[g, L["3"]]), 2), "ci": ci(bc[:, g, L["6"]] - bc[:, g, L["3"]])},
                "diff_12_3": {"est": round(float(pt[g, L["12"]] - pt[g, L["3"]]), 2), "ci": ci(bc[:, g, L["12"]] - bc[:, g, L["3"]])}}
    gap = {str(s): {"est": round(float(pt[0, L[str(s)]] - pt[1, L[str(s)]]), 2), "ci": ci(bc[:, 0, L[str(s)]] - bc[:, 1, L[str(s)]])} for s in TD_LM}
    return {"initiated": grp(0, "chemotherapy initiated by the landmark"), "not_initiated": grp(1, "not initiated by the landmark (never, or not yet)"),
            "gap_initiated_minus_not": gap, "B": B}, pt, bc, nr

clean_mask = (d.surg == 0) & (d.rad == 0)
td = {}
for key, sub, label in [("clean_no_surgery_no_rt", d[np.asarray(clean_mask)], "no primary-site surgery and no radiotherapy (tt_days = chemotherapy start)"),
                        ("full_cohort", d, "full distant NEC cohort (tt_days = first-course treatment start, any modality)")]:
    sub = sub.reset_index(drop=True); tstart, meta = td_setup(sub); rates, rr = td_rates(sub, tstart)
    naive = {}
    for lev, m in [("chemo_yes", sub.chemo == 1), ("chemo_no", sub.chemo != 1)]:
        cr = lm.crude_rates(sub[m], cuts=(0, 3, 6, 9, 12)); naive[lev] = {"n": int(m.sum()), "rates": np.round(cr.rate.to_numpy(), 4).tolist(), "rr_6_9_vs_3_6": round(float(cr.rate[2] / cr.rate[1]), 2)}
    stream_t, j_t = next_stream(); lmres, pt_td, bc_td, nr_td = td_landmark(sub, tstart, config.B(2000), stream_t)
    # naive (baseline-classified) chemo strata at the same landmarks, for comparison in this subset (point estimates only)
    k, ev, KMAX = lm.prep(sub); naive_lm = {}
    for lev, m in [("chemo_yes", (sub.chemo == 1).to_numpy()), ("chemo_no", (sub.chemo != 1).to_numpy())]:
        c, _, nn = lm.curves(k, ev, KMAX, TD_LM, WIN, m.astype(float)); naive_lm[lev] = {"est": np.round(c, 2).tolist(), "n_at_risk": nn.astype(int).tolist()}
    td[key] = {"label": label, "setup": meta, "time_dependent_rates": rates, "rate_ratios_poisson_ci": rr, "naive_stratum_rates": naive,
               "landmark_classification": lmres, "naive_stratum_landmarks": naive_lm, "seed_spawn_index": j_t, "landmarks": TD_LM.tolist()}
    td[key]["_pt"] = pt_td; td[key]["_bc"] = bc_td; td[key]["_nr"] = nr_td
    print(f"[time-dependent] {key}: n={meta['n']} chemo-yes {meta['n_chemo_yes']} (imputed {meta['n_chemo_yes_tt_unknown_imputed']}); initiated-by-LM share {meta['share_chemo_yes_initiated_by_LM']}")
    print("   rates unexposed:", [x['unexposed']['rate'] for x in rates], " exposed:", [x['exposed']['rate'] for x in rates], " RR6-9/3-6", rr["unexposed_6_9_vs_3_6"], rr["exposed_6_9_vs_3_6"])
    print("   LM initiated   :", [lmres['initiated']['landmark'][str(s)]['est'] for s in TD_LM], [lmres['initiated']['landmark'][str(s)]['n_at_risk'] for s in TD_LM])
    print("   LM not-initiat.:", [lmres['not_initiated']['landmark'][str(s)]['est'] for s in TD_LM], [lmres['not_initiated']['landmark'][str(s)]['n_at_risk'] for s in TD_LM])

# ---------------------------------------------------------------- figures
labels, ests, los, his, cols = [], [], [], [], []
def add(lbl, est, c, col): labels.append(lbl); ests.append(est); los.append(c[0]); his.append(c[1]); cols.append(col)
add("Pooled distant NEC (n=908)", pooled["diff_6_0"]["est"], pooled["diff_6_0"]["ci"], "black")
for r in results:
    add(f"{r['variable']}: {r['level']} (n={r['n']})", r["diff_6_0"]["est"], r["diff_6_0"]["ci"], "tab:blue")
    if "conditional_on_3mo_survival" in r:
        c3 = r["conditional_on_3mo_survival"]; add(f"   same, 6-3 among 3-mo survivors (n={c3['n']})", c3["diff_6_3"]["est"], c3["diff_6_3"]["ci"], "tab:orange")
tdc = td["clean_no_surgery_no_rt"]["landmark_classification"]
add(f"Time-updated: chemo initiated by LM, no surgery/RT subset: 6-3 (n at LM3={tdc['initiated']['landmark']['3']['n_at_risk']})", tdc["initiated"]["diff_6_3"]["est"], tdc["initiated"]["diff_6_3"]["ci"], "tab:green")
add(f"Time-updated: not initiated by LM, same subset: 6-3 (n at LM3={tdc['not_initiated']['landmark']['3']['n_at_risk']})", tdc["not_initiated"]["diff_6_3"]["est"], tdc["not_initiated"]["diff_6_3"]["ci"], "tab:green")
fig, ax = plt.subplots(figsize=(11, 0.32 * len(labels) + 1.4))
y = np.arange(len(labels))[::-1]
for yi, e, lo, hi, c in zip(y, ests, los, his, cols):
    ax.plot([lo, hi], [yi, yi], color=c, lw=1.4); ax.plot(e, yi, "s", color=c, ms=5)
ax.axvline(0, color="grey", lw=0.8, ls="--"); ax.set_yticks(y); ax.set_yticklabels(labels, fontsize=7.5)
ax.set_xlabel("Change in 12-month cancer-death probability (pp, 95% bootstrap CI)\nblue = LM 6 minus LM 0; orange/green = LM 6 minus LM 3", fontsize=9)
ax.set_title("Distant colorectal NEC: landmark change by stratum", fontsize=10)
plt.tight_layout(); plt.savefig(FIG + "A05_forest_6_0.png", dpi=300); plt.close()

fig, axes = plt.subplots(3, 3, figsize=(13, 10), sharex=True, sharey=True); axes = axes.ravel()
for ax, v in zip(axes, vars_):
    ax.plot(LMS, pt_all, color="black", lw=1, ls="--", label="pooled (n=908)")
    for (vv, lev, _, _), col in zip([s for s in STRATA if s[0] == v], ["tab:blue", "tab:red"]):
        r = next(x for x in results if x["variable"] == vv and x["level"] == lev); cv = r["curve_0_18"]
        ax.plot(LMS, cv["est"], color=col, lw=1.6, label=f"{lev} (n={r['n']}, plateau {r['plateau_end']['est']} mo)")
        ax.fill_between(LMS, cv["lo"], cv["hi"], color=col, alpha=0.15)
        if r["n_at_risk_lm6_lt_100"]: ax.plot(6, cv["est"][6], "x", color=col, ms=9, mew=2)
    ax.set_title(v, fontsize=10); ax.legend(fontsize=7, loc="lower left"); ax.set_ylim(0, 100); ax.grid(alpha=0.3)
for ax in axes[6:]: ax.set_xlabel("Landmark (months since diagnosis)")
for ax in axes[::3]: ax.set_ylabel("12-month cancer-death probability (%)")
fig.suptitle("Distant colorectal NEC: stratum-specific landmark curves (Aalen-Johansen, 12-month window; shaded = 95% bootstrap CI; x = <100 at risk at landmark 6)", fontsize=10)
plt.tight_layout(); plt.savefig(FIG + "A05_curves_by_stratum.png", dpi=300); plt.close()

fig, axes = plt.subplots(2, 2, figsize=(12, 9))
for col_i, key in enumerate(["clean_no_surgery_no_rt", "full_cohort"]):
    t = td[key]; pt_td = t["_pt"]; bc_td = t["_bc"]; nr_td = t["_nr"]; ax = axes[0, col_i]
    for g, (name, c) in enumerate([("chemotherapy initiated by landmark", "tab:green"), ("not initiated by landmark", "tab:purple")]):
        okm = nr_td[g] >= 20; pm = np.where(okm, pt_td[g], np.nan)   # cells with <20 at risk are not drawn
        ax.plot(TD_LM, pm, "-o", color=c, lw=1.6, ms=4, label=name + " (cells with <20 at risk omitted)")
        ax.fill_between(TD_LM, np.where(okm, np.nanpercentile(bc_td[:, g, :], 2.5, axis=0), np.nan), np.where(okm, np.nanpercentile(bc_td[:, g, :], 97.5, axis=0), np.nan), color=c, alpha=0.15)
    ax.plot(TD_LM, t["naive_stratum_landmarks"]["chemo_yes"]["est"], "--", color="tab:blue", lw=1, label="baseline-classified: chemo yes")
    ax.plot(TD_LM, t["naive_stratum_landmarks"]["chemo_no"]["est"], "--", color="tab:red", lw=1, label="baseline-classified: chemo no/unknown")
    ax.set_title(f"{'A' if col_i == 0 else 'B'}. Landmark-specific exposure, {t['label'].split(' (')[0]} (n={t['setup']['n']})", fontsize=9)
    ax.set_ylim(0, 100); ax.set_xlabel("Landmark (months)"); ax.set_ylabel("12-month cancer-death probability (%)"); ax.grid(alpha=0.3); ax.legend(fontsize=7, loc="lower left")
    ax = axes[1, col_i]; x = np.arange(4); wdt = 0.2
    ax.bar(x - 1.5 * wdt, [r_["unexposed"]["rate"] or 0 for r_ in t["time_dependent_rates"]], wdt, color="tab:purple", label="time-dependent: not (yet) on chemotherapy")
    ax.bar(x - 0.5 * wdt, [r_["exposed"]["rate"] or 0 for r_ in t["time_dependent_rates"]], wdt, color="tab:green", label="time-dependent: after chemotherapy start")
    ax.bar(x + 0.5 * wdt, t["naive_stratum_rates"]["chemo_no"]["rates"], wdt, color="tab:red", alpha=0.6, label="baseline-classified: chemo no/unknown")
    ax.bar(x + 1.5 * wdt, t["naive_stratum_rates"]["chemo_yes"]["rates"], wdt, color="tab:blue", alpha=0.6, label="baseline-classified: chemo yes")
    ax.set_xticks(x); ax.set_xticklabels(["0-3", "3-6", "6-9", "9-12"]); ax.set_xlabel("Months since diagnosis"); ax.set_ylabel("Cancer deaths per person-month")
    ax.set_title(f"{'C' if col_i == 0 else 'D'}. Crude cancer-death rates by exposure classification", fontsize=9); ax.legend(fontsize=7); ax.grid(alpha=0.3, axis="y")
fig.suptitle("Distant colorectal NEC: chemotherapy split with time-updated exposure (tt_days), vs baseline classification", fontsize=10)
plt.tight_layout(); plt.savefig(FIG + "A05_chemo_timedep.png", dpi=300); plt.close()
for key in td: td[key].pop("_pt"); td[key].pop("_bc"); td[key].pop("_nr")

# ---------------------------------------------------------------- merge standardized (R) results and write JSON
std = None
if os.path.exists(RES + "A05_std.json"):
    with open(RES + "A05_std.json") as f: std = json.load(f)
else: print("WARNING: results/A05_std.json not found - run A05_heterogeneity_std.R first; standardized section omitted")
runtime_min = (time.time() - t_start) / 60
baseline = [r for r in results if not r["variable"].startswith("Chemo")]
chemo = [r for r in results if r["variable"].startswith("Chemo")]
summary_flags = {"n_strata": len(results), "strata_with_ci_excluded_0_6_decline": [f"{r['variable']}: {r['level']}" for r in results if r["diff_6_0"]["ci"][1] < 0],
                 "strata_with_ci_excluded_0_6_rise": [f"{r['variable']}: {r['level']}" for r in results if r["diff_6_0"]["ci"][0] > 0],
                 "baseline_strata_with_crude_dip_3_6_then_rise_6_9": int(sum(r["crude_dip_3_6_then_rise_6_9"] for r in baseline)), "n_baseline_strata": len(baseline),
                 "chemo_strata_with_crude_dip_3_6_then_rise_6_9": int(sum(r["crude_dip_3_6_then_rise_6_9"] for r in chemo))}
out = {"id": "A05", "title": "Population-heterogeneity (mixture) check of the distant main-NEC landmark plateau (revised after audit)",
       "date": "2026-10-07", "cohort": "distant main NEC (hist 8013/8041-8045), n=908", "seed": SEED, "window_months": WIN,
       "B_rule": "2000 if stratum n>=100 else 1000; patient-level resampling, same resample for all landmarks; 95% percentile CI; each analysis unit uses its own SeedSequence(20261007).spawn() child stream (seed_spawn_index)",
       "reproduction_check": repro, "pooled": pooled, "strata": results, "sensitivity_other_patterns_excl_unknown_only": sens_other,
       "between_stratum_contrasts": contrasts, "summary_flags": summary_flags,
       "chemotherapy_time_dependent_supplement": td, "standardized_curves_R": std,
       "definitions": {"surgery": "lm.surgery(): no = s1 '00' or s2 'A000'; unknown = 99/98/A990/A980/blank; else yes (no unknowns in distant NEC)",
                       "liver": "mets_liver Yes vs No among 2010+ diagnoses (Blank pre-2010 and Unknown/NA excluded, n=4 excluded)",
                       "liver_only": f"mets_liver Yes and mets_bone/brain/lung all 'No' (n=389); 'other patterns' (n=302) = liver No (n=132), or liver Yes with any of bone/brain/lung Yes (n=157) or Unknown-only (n={n_unknown_only}); sensitivity excludes the {n_unknown_only} Unknown-only cases",
                       "conditional_on_3mo": "patients still at risk at landmark 3 (half-day index >= floor(2*MD*3)+1); bootstrap resamples drawn within this subset",
                       "plateau": "broken-stick hinge on landmark 0-18 point estimates, grid 0-15 by 0.25; CI = percentile of hinge over bootstrap curves; bootstrap curves with NaN at late landmarks (0 at risk in that resample) have the NaN filled with that curve's own mean before the hinge fit (revision: previously the global mean)",
                       "time_dependent_chemo": "tt_days (days from diagnosis to first-course treatment) defines the exposure start for chemotherapy-recorded patients (unknown tt imputed at the median of known chemo-yes tt in the same cohort; tt later than the survival end date clipped to it); chemo-no patients are never exposed. Primary cohort = no primary-site surgery and no radiotherapy (first treatment = chemotherapy); secondary = full cohort. Landmark classification: at landmark s, 'initiated' = exposure start <= s among patients at risk at s. Rates: person-time in each interval split at the exposure start; deaths assigned to the exposure state at death. Rate-ratio CIs are Poisson (log-scale, 1/d1+1/d2)."},
       "runtime_minutes_python": round(runtime_min, 2)}
with open(RES + "A05.json", "w") as f: json.dump(out, f, indent=1)

# ---------------------------------------------------------------- markdown
L = []
L.append("# A05 - Population heterogeneity (mixture) check of the distant-NEC plateau (revised after audit)\n")
L.append("Cohort: distant main colorectal NEC (ICD-O-3 8013/8041-8045), n = 908. Landmark estimand: Aalen-Johansen 12-month cancer-death probability among patients alive and under follow-up at the landmark. Bootstrap: patient-level, seed 20261007 (one SeedSequence child stream per analysis unit), B = 2000 (n >= 100) or 1000 (n < 100), 95% percentile CIs. Plateau end: broken-stick hinge on landmarks 0-18.\n")
L.append("## Reproduction check (ours vs manuscript)\n")
L.append("| Quantity | Ours | Manuscript |\n|---|---|---|")
L.append(f"| 12-mo cancer-death % at LM 0/3/6/12/24 | {'/'.join(f'{x:.2f}' for x in c_chk)} | 72.29/71.62/72.33/60.03/30.51 |")
L.append(f"| n at risk | {'/'.join(str(int(x)) for x in n_chk)} | 908/553/396/169/48 |")
L.append(f"| Plateau end (months) | {pl_chk} | 6.75 |")
L.append(f"| Crude rates 0-3 ... 18-24 | {'/'.join(f'{x:.4f}' for x in cr_chk)} | 0.1435/0.0909/0.1344/0.1204/0.1231/0.0840/0.0646 |")
L.append(f"| Surgery yes/no | {(d.surg==1).sum()}/{(d.surg==0).sum()} | 354/554 |")
if std: L.append(f"| Standardized curve 0/6/12 (R, original covariates) | {'/'.join(f'{x:.2f}' for x in std['original']['at_0_6_12'])}, plateau {std['original']['plateau_end']} | 72.11/71.49/58.31, plateau 7.5 |")
L.append("\n## Table 1. Stratum-specific landmark 12-month cancer-death probability (%), 95% bootstrap CI\n")
L.append("| Stratum | n | B | LM 0 | LM 3 | LM 6 | LM 12 | 6 minus 0 (pp) | 12 minus 0 (pp) | Plateau end (mo) [CI] | n at risk LM 6 |\n|---|---|---|---|---|---|---|---|---|---|---|")
def fm(x): return f"{x['est']:.1f} [{x['ci'][0]:.1f}, {x['ci'][1]:.1f}]"
L.append(f"| Pooled | 908 | 2000 | {pt_all[0]:.1f} | {pt_all[3]:.1f} | {pt_all[6]:.1f} | {pt_all[12]:.1f} | {fm(pooled['diff_6_0'])} | - | {pooled['plateau_end']['est']} [{pooled['plateau_end']['ci'][0]}, {pooled['plateau_end']['ci'][1]}] | {int(nrisk_all[6])} |")
for r in results + [sens_other]:
    lmk = r["landmark"]; flag = " *" if r["n_at_risk_lm6_lt_100"] else ""
    L.append(f"| {r['variable']}: {r['level']} | {r['n']} | {r['B']} | {fm(lmk['0'])} | {fm(lmk['3'])} | {fm(lmk['6'])} | {fm(lmk['12'])} | {fm(r['diff_6_0'])} | {fm(r['diff_12_0'])} | {r['plateau_end']['est']} [{r['plateau_end']['ci'][0]}, {r['plateau_end']['ci'][1]}] | {lmk['6']['n_at_risk']}{flag} |")
L.append("\n\\* fewer than 100 patients at risk at landmark 6. The last row is a sensitivity analysis (audit): 'other patterns' without the liver-Yes patients whose only non-liver information is Unknown.\n")
L.append("## Table 2. Treatment strata conditioned on surviving 3 months (immortal-time guard)\n")
L.append("| Stratum | n (alive at 3 mo) | LM 3 | LM 6 | LM 9 | LM 12 | 6 minus 3 (pp) | 12 minus 3 (pp) |\n|---|---|---|---|---|---|---|---|")
for r in results:
    if "conditional_on_3mo_survival" in r:
        c3 = r["conditional_on_3mo_survival"]; lmk = c3["landmark"]
        L.append(f"| {r['variable']}: {r['level']} | {c3['n']} | {fm(lmk['3'])} | {fm(lmk['6'])} | {fm(lmk['9'])} | {fm(lmk['12'])} | {fm(c3['diff_6_3'])} | {fm(c3['diff_12_3'])} |")
L.append("\n## Table 3. Crude cancer-death rates per person-month by interval\n")
L.append("| Stratum | 0-3 | 3-6 | 6-9 | 9-12 | crude RR 6-9 vs 3-6 | dip 3-6 then rise 6-9 |\n|---|---|---|---|---|---|---|")
for r in results:
    cr = r["crude_rates"]; L.append(f"| {r['variable']}: {r['level']} | " + " | ".join(f"{x['rate']:.4f} ({x['deaths']}/{x['person_months']:.0f})" for x in cr) + f" | {r['crude_rr_6_9_vs_3_6']} | {'yes' if r['crude_dip_3_6_then_rise_6_9'] else 'no'} |")
L.append("\n## Table 4. Between-stratum contrasts (independent bootstrap streams paired by draw index)\n")
L.append("| Variable | Contrast | LM 0 difference (pp) | LM 6 difference (pp) | Difference of (6-0) changes (pp) | Difference of (6-3) changes, 3-mo survivors |\n|---|---|---|---|---|---|")
for e in contrasts:
    c3 = e.get("diff_of_6_3_differences_conditional"); L.append(f"| {e['variable']} | {e['contrast']} | {fm(e['diff_of_LM0'])} | {fm(e['diff_of_LM6'])} | {fm(e['diff_of_6_0_differences'])} | {fm(c3) if c3 else '-'} |")
if std:
    L.append("\n## Table 5. Standardized landmark curves (R, pp.R std_curve: per-landmark multi-state Cox, direct standardization to the cohort)\n")
    L.append("| Covariate set | n | LM 0 | LM 6 | LM 12 | Plateau end (mo) |\n|---|---|---|---|---|---|")
    for key in ["original", "plus_surg", "plus_nhw", "plus_surg_nhw", "sub2010_original", "sub2010_plus_liver", "sub2010_plus_liver_surg_nhw"]:
        s = std[key]; L.append(f"| {s['label']} | {s['n']} | {s['at_0_6_12'][0]:.2f} | {s['at_0_6_12'][1]:.2f} | {s['at_0_6_12'][2]:.2f} | {s['plateau_end']} |")
L.append("\n## Table 6. Supplementary (revision): chemotherapy with time-updated exposure from time to first treatment (tt_days)\n")
L.append("Landmark-specific classification: at each landmark, patients at risk are split into 'chemotherapy initiated by the landmark' vs 'not initiated' (never, or not yet). This removes immortal time from the treated group because a patient contributes to the treated curve only from landmarks after the treatment start. Primary cohort: no primary-site surgery and no radiotherapy, where the first treatment of a chemotherapy-recorded patient is the chemotherapy itself; secondary: full cohort (tt_days = first treatment of any modality, a lower bound for the chemotherapy start in surgery-first patients).\n")
L.append("| Cohort | Group | LM 0 | LM 1 | LM 2 | LM 3 | LM 6 | LM 9 | LM 12 | 6 minus 1 (pp) | 6 minus 3 (pp) | 12 minus 3 (pp) |\n|---|---|---|---|---|---|---|---|---|---|---|---|")
def fmn(x): return (f"{x['est']:.1f} [{x['ci'][0]:.1f}, {x['ci'][1]:.1f}] (n={x['n_at_risk']})" if not np.isnan(x['est']) else f"- (n={x['n_at_risk']})") + (" \u2020" if x['n_at_risk'] < 20 else "")
for key, t in td.items():
    lc = t["landmark_classification"]
    for g in ["initiated", "not_initiated"]:
        G = lc[g]; L.append(f"| {t['label'].split(' (')[0]} (n={t['setup']['n']}) | {G['group']} | " + " | ".join(fmn(G['landmark'][str(s)]) for s in TD_LM) + f" | {fm(G['diff_6_1'])} | {fm(G['diff_6_3'])} | {fm(G['diff_12_3'])} |")
    L.append(f"| same | gap: initiated minus not initiated | " + " | ".join(fm(lc['gap_initiated_minus_not'][str(s)]) for s in TD_LM) + " | - | - | - |")
    L.append(f"| same | baseline-classified chemo yes (n={t['naive_stratum_rates']['chemo_yes']['n']}), point estimates | " + " | ".join(f"{x:.1f} (n={n_})" for x, n_ in zip(t['naive_stratum_landmarks']['chemo_yes']['est'], t['naive_stratum_landmarks']['chemo_yes']['n_at_risk'])) + " | - | - | - |")
    L.append(f"| same | baseline-classified chemo no/unknown (n={t['naive_stratum_rates']['chemo_no']['n']}), point estimates | " + " | ".join(f"{x:.1f} (n={n_})" for x, n_ in zip(t['naive_stratum_landmarks']['chemo_no']['est'], t['naive_stratum_landmarks']['chemo_no']['n_at_risk'])) + " | - | - | - |")
L.append("\n\u2020 fewer than 20 patients at risk: not interpretable (shown for completeness). Differences 6 minus 1 / 6 minus 3 / 12 minus 3 are between the landmark-specific groups, which have different membership at each landmark.\n")
L.append("\nCrude cancer-death rates per person-month with person-time split at the chemotherapy start (time-dependent exposure) vs baseline classification:\n")
L.append("| Cohort | Classification | 0-3 | 3-6 | 6-9 | 9-12 | RR 6-9 vs 3-6 [Poisson CI] | RR 3-6 vs 0-3 |\n|---|---|---|---|---|---|---|---|")
def fr(x): return f"{x['RR']} [{x['ci'][0]}, {x['ci'][1]}]" if x else "-"
for key, t in td.items():
    lab = t["label"].split(" (")[0]; rr = t["rate_ratios_poisson_ci"]
    for g, nm in [("unexposed", "time-dependent: not (yet) on chemotherapy"), ("exposed", "time-dependent: after chemotherapy start")]:
        L.append(f"| {lab} | {nm} | " + " | ".join(f"{x[g]['rate']:.4f} ({x[g]['deaths']}/{x[g]['person_months']:.0f})" if x[g]['rate'] is not None else "-" for x in t["time_dependent_rates"]) + f" | {fr(rr[g + '_6_9_vs_3_6'])} | {fr(rr[g + '_3_6_vs_0_3'])} |")
    for g, nm in [("chemo_no", "baseline: chemo no/unknown"), ("chemo_yes", "baseline: chemo yes")]:
        L.append(f"| {lab} | {nm} (n={t['naive_stratum_rates'][g]['n']}) | " + " | ".join(f"{x:.4f}" for x in t["naive_stratum_rates"][g]["rates"]) + f" | {t['naive_stratum_rates'][g]['rr_6_9_vs_3_6']} | - |")
L.append(f"\nSetup: primary cohort n = {td['clean_no_surgery_no_rt']['setup']['n']} (chemotherapy recorded {td['clean_no_surgery_no_rt']['setup']['n_chemo_yes']}, of whom {td['clean_no_surgery_no_rt']['setup']['n_chemo_yes_tt_unknown_imputed']} had unknown tt_days imputed at {td['clean_no_surgery_no_rt']['setup']['imputed_tt_days']:.0f} days; share initiated by landmark 1/2/3 = {td['clean_no_surgery_no_rt']['setup']['share_chemo_yes_initiated_by_LM']['1']:.2f}/{td['clean_no_surgery_no_rt']['setup']['share_chemo_yes_initiated_by_LM']['2']:.2f}/{td['clean_no_surgery_no_rt']['setup']['share_chemo_yes_initiated_by_LM']['3']:.2f}); full cohort n = 908 (chemotherapy recorded {td['full_cohort']['setup']['n_chemo_yes']}, {td['full_cohort']['setup']['n_chemo_yes_tt_unknown_imputed']} imputed at {td['full_cohort']['setup']['imputed_tt_days']:.0f} days, {td['full_cohort']['setup']['n_chemo_yes_tt_after_survival_end_clipped']} with treatment date after the survival end date clipped).\n")
L.append("## Figures\n\n- figures/output/A05_forest_6_0.png - forest plot of the landmark 6 minus landmark 0 change by stratum (orange: 6 minus 3 among 3-month survivors for treatment strata; green: 6 minus 3 with time-updated chemotherapy classification).\n- figures/output/A05_curves_by_stratum.png - stratum-specific landmark curves (landmarks 0-18) with bootstrap bands and the pooled curve.\n- figures/output/A05_chemo_timedep.png - chemotherapy split with time-updated exposure (landmark classification, top; time-dependent crude rates, bottom) vs baseline classification.\n")

# ---------------------------------------------------------------- interpretation (data-driven text, reworded after audit)
early_decline = [r for r in results if r["diff_6_0"]["ci"][1] < 0]; early_rise = [r for r in results if r["diff_6_0"]["ci"][0] > 0]
plateau_lt3 = [r for r in results if r["plateau_end"]["est"] < 3]
sy = next(r for r in results if r["variable"].startswith("Surgery") and r["level"] == "Yes"); sn = next(r for r in results if r["variable"].startswith("Surgery") and r["level"] == "No")
cy = next(r for r in results if r["variable"].startswith("Chemo") and r["level"] == "Yes"); cn = next(r for r in results if r["variable"].startswith("Chemo") and r["level"] != "Yes")
lv = [r for r in results if r["variable"].startswith("Liver metastasis")]; lo_ = [r for r in results if r["variable"].startswith("Liver-only")]
gap = {e["variable"]: e for e in contrasts}
gaps_sorted = sorted(contrasts, key=lambda e: -abs(e["diff_of_LM0"]["est"]))
tdc = td["clean_no_surgery_no_rt"]; tdf = td["full_cohort"]; lc, lf = tdc["landmark_classification"], tdf["landmark_classification"]
L.append("## Interpretation\n")
L.append(f"1. The pooled plateau (landmark 6 minus 0 = {pooled['diff_6_0']['est']:+.1f} pp, 95% CI {pooled['diff_6_0']['ci']}; plateau end {pooled['plateau_end']['est']} months) is reproduced within all {len(baseline)} baseline strata (age, sex, primary site, histology, race, surgery, liver metastasis and metastatic pattern): the 6-minus-0 change ranges from {min(r['diff_6_0']['est'] for r in baseline):+.1f} to {max(r['diff_6_0']['est'] for r in baseline):+.1f} pp, every CI includes zero" + (" or is positive" if any(r in early_rise for r in baseline) else "") + f", and point-estimate plateau ends are {min(r['plateau_end']['est'] for r in baseline)}-{max(r['plateau_end']['est'] for r in baseline)} months. Overall, {len(results) - len(early_decline)} of {len(results)} strata show no CI-excluded early decline; the exception is " + "; ".join(f"{r['variable']} {r['level']} ({r['diff_6_0']['est']:+.1f} pp, CI {r['diff_6_0']['ci']})" for r in early_decline) + ".")
L.append(f"2. Point-estimate plateau ends are >= 3 months in {len(results) - len(plateau_lt3)} of {len(results)} strata" + ("; strata with a plateau end < 3 months: " + "; ".join(f"{r['variable']} {r['level']} ({r['plateau_end']['est']} mo, n at risk at LM 6 = {r['landmark']['6']['n_at_risk']})" for r in plateau_lt3) if plateau_lt3 else "") + ". The hinge location is imprecise within strata (wide bootstrap CIs), so the plateau end should be interpreted as a pooled quantity.")
L.append(f"3. Surgery is the treatment stratifier of primary interest (LM 0: {sy['landmark']['0']['est']:.1f}% resected vs {sn['landmark']['0']['est']:.1f}% unresected, gap {gap['Surgery (primary site)']['diff_of_LM0']['est']:+.1f} pp, CI {gap['Surgery (primary site)']['diff_of_LM0']['ci']}), but it is not the largest baseline gap: liver metastasis ({gap['Liver metastasis (2010+)']['diff_of_LM0']['est']:+.1f} pp, CI {gap['Liver metastasis (2010+)']['diff_of_LM0']['ci']}) and the chemotherapy split ({gap['Chemotherapy recorded']['diff_of_LM0']['est']:+.1f} pp, CI {gap['Chemotherapy recorded']['diff_of_LM0']['ci']}) separate the LM 0 estimate far more. Both surgery components are flat from 0 to 6 months (resected {sy['diff_6_0']['est']:+.1f} pp, unresected {sn['diff_6_0']['est']:+.1f} pp) and, among 3-month survivors, 6 minus 3 = {sy['conditional_on_3mo_survival']['diff_6_3']['est']:+.1f} vs {sn['conditional_on_3mo_survival']['diff_6_3']['est']:+.1f} pp; a resected/unresected mixture therefore does not generate the plateau.")
L.append(f"4. The chemotherapy split is different: it DOES decompose the flat pooled 0-6 month segment into one falling and one rising component. The no/unknown-chemotherapy stratum starts at {cn['landmark']['0']['est']:.1f}% and falls to {cn['landmark']['6']['est']:.1f}% ({cn['diff_6_0']['est']:+.1f} pp, CI {cn['diff_6_0']['ci']}; plateau end {cn['plateau_end']['est']}; only {cn['landmark']['6']['n_at_risk']} at risk at LM 6) with crude rates that decline monotonically ({'/'.join(f'{x['rate']:.3f}' for x in cn['crude_rates'])}; RR 6-9 vs 3-6 = {cn['crude_rr_6_9_vs_3_6']}), whereas the chemotherapy stratum starts at {cy['landmark']['0']['est']:.1f}% and rises to {cy['landmark']['6']['est']:.1f}% ({cy['diff_6_0']['est']:+.1f} pp, CI {cy['diff_6_0']['ci']}) with crude rates that rise monotonically to 9 months ({'/'.join(f'{x['rate']:.3f}' for x in cy['crude_rates'])}; no 3-6 month dip). This is the mixture signature the analysis was designed to detect. It is consistent with immortal-time/treatment-selection bias of a treatment variable recorded at diagnosis (patients who die early cannot receive chemotherapy, so the untreated stratum is enriched for imminent deaths at LM 0) and/or with genuine treatment-related heterogeneity (treatment response followed by progression in the treated, declining frailty in the untreated); the baseline-classified design cannot separate these. Within the treatment strata the plateau claim is therefore restricted to landmarks >= 3: among 3-month survivors, 6 minus 3 = {cy['conditional_on_3mo_survival']['diff_6_3']['est']:+.1f} pp (CI {cy['conditional_on_3mo_survival']['diff_6_3']['ci']}) with chemotherapy and {cn['conditional_on_3mo_survival']['diff_6_3']['est']:+.1f} pp (CI {cn['conditional_on_3mo_survival']['diff_6_3']['ci']}) without, both CIs including zero, and LM 6 is {cy['landmark']['6']['est']:.1f}% vs {cn['landmark']['6']['est']:.1f}%.")
rrc, rrf = tdc["rate_ratios_poisson_ci"], tdf["rate_ratios_poisson_ci"]
L.append(f"5. Supplementary immortal-time test (Table 6, figure A05_chemo_timedep.png). (a) Landmark probabilities. With the landmark-specific classification in the no-surgery/no-radiotherapy cohort (n = {tdc['setup']['n']}; {tdc['setup']['share_chemo_yes_initiated_by_LM']['1']:.0%}/{tdc['setup']['share_chemo_yes_initiated_by_LM']['2']:.0%}/{tdc['setup']['share_chemo_yes_initiated_by_LM']['3']:.0%} of chemotherapy patients had started by LM 1/2/3), the 'not initiated by the landmark' curve is flat: {lc['not_initiated']['landmark']['0']['est']:.1f}% (LM 0; this is essentially the whole cohort), {lc['not_initiated']['landmark']['1']['est']:.1f}% (LM 1, n = {lc['not_initiated']['landmark']['1']['n_at_risk']}), {lc['not_initiated']['landmark']['2']['est']:.1f}% (LM 2, n = {lc['not_initiated']['landmark']['2']['n_at_risk']}), {lc['not_initiated']['landmark']['3']['est']:.1f}% (LM 3, n = {lc['not_initiated']['landmark']['3']['n_at_risk']}), whereas the baseline-classified no/unknown stratum in the same cohort falls from {tdc['naive_stratum_landmarks']['chemo_no']['est'][0]:.1f}% to {tdc['naive_stratum_landmarks']['chemo_no']['est'][1]:.1f}/{tdc['naive_stratum_landmarks']['chemo_no']['est'][2]:.1f}/{tdc['naive_stratum_landmarks']['chemo_no']['est'][3]:.1f}% at LM 1/2/3; the 'initiated' curve is {lc['initiated']['landmark']['1']['est']:.1f}/{lc['initiated']['landmark']['2']['est']:.1f}/{lc['initiated']['landmark']['3']['est']:.1f}/{lc['initiated']['landmark']['6']['est']:.1f}% at LM 1/2/3/6 (6 minus 1 = {lc['initiated']['diff_6_1']['est']:+.1f} pp, CI {lc['initiated']['diff_6_1']['ci']}), and the treated-minus-untreated gap is {lc['gap_initiated_minus_not']['1']['est']:+.1f} (CI {lc['gap_initiated_minus_not']['1']['ci']}), {lc['gap_initiated_minus_not']['2']['est']:+.1f} and {lc['gap_initiated_minus_not']['3']['est']:+.1f} pp at LM 1/2/3. The full cohort behaves the same way (not initiated {lf['not_initiated']['landmark']['0']['est']:.1f}/{lf['not_initiated']['landmark']['1']['est']:.1f}/{lf['not_initiated']['landmark']['2']['est']:.1f}/{lf['not_initiated']['landmark']['3']['est']:.1f}% vs baseline no/unknown {tdf['naive_stratum_landmarks']['chemo_no']['est'][0]:.1f}/{tdf['naive_stratum_landmarks']['chemo_no']['est'][1]:.1f}/{tdf['naive_stratum_landmarks']['chemo_no']['est'][2]:.1f}/{tdf['naive_stratum_landmarks']['chemo_no']['est'][3]:.1f}% at LM 0/1/2/3; initiated {lf['initiated']['landmark']['1']['est']:.1f}/{lf['initiated']['landmark']['3']['est']:.1f}/{lf['initiated']['landmark']['6']['est']:.1f}% at LM 1/3/6, 6 minus 1 = {lf['initiated']['diff_6_1']['est']:+.1f} pp, CI {lf['initiated']['diff_6_1']['ci']}, vs the baseline-classified chemotherapy stratum's CI-excluded rise of {cy['diff_6_0']['est']:+.1f} pp). At the level of 12-month landmark probabilities the falling-untreated/rising-treated decomposition therefore largely disappears once exposure is classified at the landmark, which supports immortal-time/classification bias as its main source. (b) Interval hazards. Splitting person-time at the chemotherapy start, the full-cohort rate after the start is {'/'.join(f'{x['exposed']['rate']:.3f}' for x in tdf['time_dependent_rates'])} per person-month over 0-3/3-6/6-9/9-12 (RR 6-9 vs 3-6 = {fr(rrf['exposed_6_9_vs_3_6'])}) and the rate while not (yet) treated is {'/'.join(f'{x['unexposed']['rate']:.3f}' for x in tdf['time_dependent_rates'])} (RR 3-6 vs 0-3 = {fr(rrf['unexposed_3_6_vs_0_3'])}; RR 6-9 vs 3-6 = {fr(rrf['unexposed_6_9_vs_3_6'])}, {sum(x['unexposed']['deaths'] for x in tdf['time_dependent_rates'][1:])} deaths after 3 months); the clean cohort gives {'/'.join(f'{x['exposed']['rate']:.3f}' for x in tdc['time_dependent_rates'])} after the start (RR 6-9 vs 3-6 = {fr(rrc['exposed_6_9_vs_3_6'])}) and {'/'.join(f'{x['unexposed']['rate']:.3f}' for x in tdc['time_dependent_rates'])} while untreated (RR 6-9 vs 3-6 = {fr(rrc['unexposed_6_9_vs_3_6'])}, {sum(x['unexposed']['deaths'] for x in tdc['time_dependent_rates'][1:])} deaths after 3 months). A post-initiation rise of the cancer-death hazard from 3-6 to 6-9 months is thus present with time-updated exposure, while the untreated hazard is highest in 0-3 months and sparse afterwards. The pooled dip-then-rise hazard is therefore at least partly compositional (early person-time dominated by untreated high-hazard patients, later person-time by treated patients whose hazard rises after 6 months); this hazard-level heterogeneity is not removed by the time-updated classification and is compatible with treatment-related heterogeneity (response followed by progression). In short: the landmark-probability decomposition is mostly attributable to immortal-time/classification bias, the hazard-level decomposition is not, and the present design cannot separate treatment effect from treatment selection in the latter.")
L.append(f"6. Liver metastasis (2010+): LM 0 {lv[0]['landmark']['0']['est']:.1f}% (yes) vs {lv[1]['landmark']['0']['est']:.1f}% (no), the largest prognostic baseline gap; 6 minus 0 = {lv[0]['diff_6_0']['est']:+.1f} vs {lv[1]['diff_6_0']['est']:+.1f} pp. The liver-negative stratum is small (n = {lv[1]['n']}, {lv[1]['landmark']['6']['n_at_risk']} at risk at LM 6) and its curve is noisy. Liver-only (n = {lo_[0]['n']}) vs other patterns (n = {lo_[1]['n']}; liver No, or liver Yes with another site Yes or Unknown-only): 6 minus 0 = {lo_[0]['diff_6_0']['est']:+.1f} vs {lo_[1]['diff_6_0']['est']:+.1f} pp; excluding the {n_unknown_only} Unknown-only cases from 'other patterns' gives {sens_other['diff_6_0']['est']:+.1f} pp (CI {sens_other['diff_6_0']['ci']}, n = {sens_other['n']}). No measured metastatic pattern isolates a declining-hazard subgroup.")
L.append(f"7. Crude cancer-death rates show the 3-6 month dip followed by a 6-9 month rise in {summary_flags['baseline_strata_with_crude_dip_3_6_then_rise_6_9']} of {summary_flags['n_baseline_strata']} baseline strata (crude RR 6-9 vs 3-6: surgery yes {sy['crude_rr_6_9_vs_3_6']}, no {sn['crude_rr_6_9_vs_3_6']}; liver yes {lv[0]['crude_rr_6_9_vs_3_6']}, no {lv[1]['crude_rr_6_9_vs_3_6']}; age <65 {next(r for r in results if r['level']=='<65')['crude_rr_6_9_vs_3_6']}, >=65 {next(r for r in results if r['level']=='>=65')['crude_rr_6_9_vs_3_6']}; colon {next(r for r in results if r['level']=='Colon')['crude_rr_6_9_vs_3_6']}, rectum {next(r for r in results if r['level']=='Rectum')['crude_rr_6_9_vs_3_6']}; small cell {next(r for r in results if r['level'].startswith('Small'))['crude_rr_6_9_vs_3_6']}, large cell {next(r for r in results if r['level'].startswith('Large'))['crude_rr_6_9_vs_3_6']}; NHW {next(r for r in results if r['level']=='Non-Hispanic White')['crude_rr_6_9_vs_3_6']}, other {next(r for r in results if r['variable']=='Race' and r['level']=='Other')['crude_rr_6_9_vs_3_6']}), but within neither chemotherapy stratum ({summary_flags['chemo_strata_with_crude_dip_3_6_then_rise_6_9']} of 2: no/unknown RR {cn['crude_rr_6_9_vs_3_6']} with declining rates, yes RR {cy['crude_rr_6_9_vs_3_6']} with rising rates and no dip).")
if std:
    L.append(f"8. Direct standardization (R) leaves the curve essentially unchanged when surgery and race are added to the original covariates (LM 0/6/12: original {std['original']['at_0_6_12']}; + surg {std['plus_surg']['at_0_6_12']}; + NHW {std['plus_nhw']['at_0_6_12']}; + surg + NHW {std['plus_surg_nhw']['at_0_6_12']}; plateau ends {std['original']['plateau_end']}, {std['plus_surg']['plateau_end']}, {std['plus_nhw']['plateau_end']}, {std['plus_surg_nhw']['plateau_end']} months). In the 2010+ subset with known liver status (n = {std['sub2010_original']['n']}), adding the liver indicator gives {std['sub2010_plus_liver']['at_0_6_12']} (plateau {std['sub2010_plus_liver']['plateau_end']}) vs {std['sub2010_original']['at_0_6_12']} (plateau {std['sub2010_original']['plateau_end']}) without it.")
L.append("9. Conclusion: the plateau persists within every baseline (non-treatment) stratum, so a mixture of measured prognostic subgroups (age, sex, site, histology, race, surgery, metastatic pattern) does not explain it. The chemotherapy split does decompose the pooled 0-6 month plateau into a rising treated and a falling untreated component; with landmark-specific (time-updated) exposure classification the decomposition of the landmark probabilities largely disappears, so at that level it is mostly an immortal-time/classification artefact; at the hazard level a post-treatment rise and an early untreated peak remain, which is compatible with genuine treatment-related heterogeneity that the present observational design cannot separate from treatment selection. Within the baseline-classified treatment strata the plateau is supported only from landmark 3 onward. Unmeasured frailty (performance status, treatment response) cannot be excluded by this analysis and is addressed by the frailty/selection sensitivity tasks.\n")
L.append("## Revision after audit\n")
L.append(f"- Interpretation points 1, 3, 4, 7 and 9 (previously 1, 3, 4, 5, 7) reworded: the plateau claim is now stated for the {len(baseline)} baseline strata, the chemotherapy split is reported as a genuine decomposition of the 0-6 month segment (rising treated / falling untreated component) that is consistent with immortal-time bias and/or treatment-related heterogeneity, and the within-treatment-stratum plateau claim is restricted to landmarks >= 3. 'Surgery is the strongest prognostic stratifier' replaced by 'treatment stratifier of primary interest', with liver metastasis and chemotherapy named as the largest baseline gaps. 'Dip-then-rise in essentially every stratum' replaced by the counted statement (baseline strata vs chemotherapy strata; Table 3 now has a dip-then-rise column).")
L.append("- New supplementary analysis (Table 6, figure A05_chemo_timedep.png, interpretation point 5): landmark-specific 'chemotherapy initiated by the landmark' classification and time-dependent-exposure crude rates using tt_days, primary in the no-surgery/no-radiotherapy cohort and secondary in the full cohort.")
L.append(f"- Definition of 'other patterns' corrected (liver No, or liver Yes with another site Yes or Unknown-only; {n_unknown_only} Unknown-only cases) and a sensitivity analysis excluding those cases added (Table 1, last row).")
L.append("- plateau(): NaN entries of bootstrap curves are now filled with the row's own mean (previously the global mean of the bootstrap matrix); this changes only the plateau-end CIs of sparse strata.")
L.append("- Bootstrap streams: every analysis unit now uses its own SeedSequence(20261007) child (seed_spawn_index in the JSON), so between-stratum contrasts are genuinely independent bootstraps paired by draw index; all bootstrap CIs therefore differ slightly from the first version (point estimates unchanged).")
L.append("- key_findings: 'colon 1.56-ish' corrected to the JSON value; the R script and the standardized curves are unchanged (re-run from scratch).\n")
L.append(f"Deviations from plan: none in the stratum set. The 'conditional on surviving 3 months' analysis uses the patients still at risk at landmark 3 (identical point estimates to the unconditional landmarks 3-12; the bootstrap is drawn within survivors). Liver strata exclude 4 patients with Unknown/missing liver status. Late-landmark multi-state Cox fits in R emitted 'coefficient may be infinite' warnings for sparse other-death transitions (few events at landmarks >= 12); these affect only the other-death hazard and the late tail of the standardized curve. Time-dependent chemotherapy: tt_days unknown for {td['full_cohort']['setup']['n_chemo_yes_tt_unknown_imputed']} chemotherapy-recorded patients (imputed at the median); Poisson CIs for the time-dependent rate ratios instead of bootstrap. Runtime (Python): {runtime_min:.1f} min" + (f"; R standardization: {std['meta']['runtime_sec']/60:.1f} min." if std else "."))
with open(RES + "A05.md", "w") as f: f.write("\n".join(L))
print(f"wrote {RES}A05.json and A05.md; runtime {runtime_min:.1f} min")
