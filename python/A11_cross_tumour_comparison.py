# -*- coding: utf-8 -*-
"""A11 cross_tumour_comparison: distant-stage cross-tumour comparison (Table 3 and extensions).
(1) Table 3 reproduction (NEC, SCLC, AC, PANC): landmark 0/6/12 estimates, 6-0 / 12-6 / 12-0 changes, plateau end,
    paired NEC-minus-comparator differences; cancer death (Aalen-Johansen CIF) and all-cause death.
(2) Lung LCNEC (8013) as fourth comparator + 2x2 (colorectal SC / colorectal LC / lung SC / lung LC).
(3) SCLC definitions: Breitling (8002/8040-8045) and 8041 only.
(4) Immortal-time demonstration: chemo-recorded-yes from diagnosis vs conditioned on 3-month survival (+ chemo no/unknown).
(5) Cytology-inclusive comparators (AC_cyt, LUNG_cyt small-cell codes, PANC_cyt excl. islet).
(6) Regional-stage auxiliary comparison (NEC vs AC, PANC, SCLC).
(7) Figures: four-tumour monthly curves + paired-difference forest; 2x2 figure.
Run: python3 python/A11_cross_tumour_comparison.py   (single-threaded, non-interactive; paths from python/config.py)
"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1"); os.environ.setdefault("MKL_NUM_THREADS", "1")
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
T0 = time.time()
import config
import lm
RES = config.RES; FIGS = config.FIGS
os.makedirs(RES, exist_ok=True); os.makedirs(FIGS, exist_ok=True)
B = config.B(2000); SEED0 = 20260928
LM25 = np.arange(0, 25); M19 = np.arange(0, 19).astype(float); WINDOW = 12
# group index -> seed = SEED0 + index (authors' convention: NEC 0, AC 1, SCLC 2, PANC 3, LCNEC-lung 4, others 5+)
GI = {"NEC": 0, "AC": 1, "SCLC": 2, "PANC": 3, "LCNEC": 4,
      "CRC_SC": 5, "CRC_LC": 6, "SCLC_breitling": 7, "SCLC_8041": 8,
      "NEC_chemo": 9, "SCLC_chemo": 10, "AC_chemo": 11, "PANC_chemo": 12,
      "NEC_chemo_c3": 13, "SCLC_chemo_c3": 14, "AC_chemo_c3": 15, "PANC_chemo_c3": 16,
      "AC_cyt": 17, "SCLC_cyt": 18, "PANC_cyt": 19,
      "NEC_regional": 20, "AC_regional": 21, "PANC_regional": 22, "SCLC_regional": 23,
      "NEC_nochemo": 24, "SCLC_nochemo": 25, "AC_nochemo": 26, "PANC_nochemo": 27,
      "LCNEC_chemo": 28, "LCNEC_chemo_c3": 29}
LABEL = {"NEC": "Colorectal NEC", "SCLC": "Small-cell lung cancer", "AC": "Colorectal adenocarcinoma", "PANC": "Pancreatic adenocarcinoma",
         "LCNEC": "Lung LCNEC (8013)", "CRC_SC": "Colorectal small-cell NEC (8041-8045)", "CRC_LC": "Colorectal large-cell NEC (8013)"}
# colours: dataviz default categorical slots 1-5 (fixed order by entity)
COL = {"NEC": "#2a78d6", "SCLC": "#eb6834", "AC": "#1baf7a", "PANC": "#eda100", "LCNEC": "#e87ba4"}
USECOLS = ["pid", "year", "stage", "hist", "ev", "unkcod", "tm", "sdays", "age", "male", "chemo", "diag", "grade17", "grade18"]

# ------------------------------------------------------------------ helpers
def load_cols(name):
    """same as lm.load but only the columns this task needs (memory: ~11 agents share the machine)"""
    return pd.read_csv(lm.DATA + f"cohort_{name}.csv", usecols=USECOLS,
                       dtype={"pid": str, "hist": str, "grade17": str, "grade18": str, "diag": str, "stage": str}, low_memory=False)

def deaths_in_window(k, ev, landmarks, window, KMAX):
    out = []
    for s in landmarks:
        lo = int(np.floor(2 * lm.MD * s)) + 1; hi = min(int(np.floor(2 * lm.MD * (s + window))), KMAX + 1)
        out.append(int(((ev == 1) & (k >= lo) & (k <= hi)).sum()))
    return out

def point_and_boot(df, seed, landmarks=LM25, window=WINDOW, B=B, kmax_months=48):
    """point estimates + bootstrap matrices (patient-level resampling, same resample for every landmark)"""
    k, ev, KMAX = lm.prep(df, kmax_months); n = len(k)
    cif, allc, nrisk = lm.curves(k, ev, KMAX, landmarks, window)
    d1 = deaths_in_window(k, ev, landmarks, window, KMAX)
    rng = np.random.default_rng(seed)
    bc = np.empty((B, len(landmarks))); ba = np.empty((B, len(landmarks)))
    for b in range(B):
        w = np.bincount(rng.integers(0, n, n), minlength=n).astype(float)
        bc[b], ba[b], _ = lm.curves(k, ev, KMAX, landmarks, window, w)
    return dict(cif=cif, allc=allc, nrisk=nrisk, d1=d1, bc=bc, ba=ba, n=n, seed=int(seed))

def ci(v): return lm.pct(v)
def r(x, d=2):
    if x is None: return None
    x = float(x)
    return None if np.isnan(x) else round(x, d)

def hinge_point(y, m=M19): return float(lm.hinge_many(np.asarray(y)[None, :len(m)], m)[0])
def hinge_boot(bmat, m=M19):
    Y = bmat[:, :len(m)]; ok = ~np.isnan(Y).any(axis=1)
    return lm.hinge_many(Y[ok], m), int(ok.sum())

def summarize(res, key="cif", bkey="bc", landmarks=LM25, lms=(0, 6, 12), plateau=True):
    """landmark estimates with CI, changes with CI (6-0, 12-6, 12-0), plateau end with CI"""
    pt = res[key]; bm = res[bkey]; idx = {s: int(np.where(landmarks == s)[0][0]) for s in lms}
    out = {"n": int(res["n"]), "seed": res["seed"], "B": int(bm.shape[0]), "landmark": {}, "change": {}}
    for s in lms:
        j = idx[s]
        out["landmark"][str(s)] = {"n_at_risk": int(res["nrisk"][j]), "cancer_deaths_window": int(res["d1"][j]) if key == "cif" else None,
                                  "pct": r(pt[j]), "ci": [r(x) for x in ci(bm[:, j])]}
    for a, b_ in ((0, 6), (6, 12), (0, 12)):
        d = pt[idx[b_]] - pt[idx[a]]; bd = bm[:, idx[b_]] - bm[:, idx[a]]
        out["change"][f"{b_}-{a}"] = {"pct_points": r(d), "ci": [r(x) for x in ci(bd)]}
    if plateau:
        c = hinge_point(pt); cs, nb = hinge_boot(bm)
        out["plateau_end"] = {"months": c, "ci": [r(x) for x in ci(cs)], "boot_rows": nb,
                              "share_zero": r(float(np.mean(cs == 0)), 3), "median": r(float(np.median(cs)))}
    out["monthly"] = {"landmarks": landmarks.tolist(), "pct": [r(x) for x in pt], "lo": [r(x) for x in np.nanpercentile(bm, 2.5, axis=0)],
                      "hi": [r(x) for x in np.nanpercentile(bm, 97.5, axis=0)], "n_at_risk": res["nrisk"].astype(int).tolist()}
    return out

def paired(resA, resB, bkey="bc", ptkey="cif", landmarks=LM25):
    """A minus B differences of the landmark changes, independent bootstraps paired by index"""
    idx = {s: int(np.where(landmarks == s)[0][0]) for s in (0, 6, 12)}; out = {}
    for a, b_ in ((0, 6), (6, 12), (0, 12)):
        dA = resA[ptkey][idx[b_]] - resA[ptkey][idx[a]]; dB = resB[ptkey][idx[b_]] - resB[ptkey][idx[a]]
        bA = resA[bkey][:, idx[b_]] - resA[bkey][:, idx[a]]; bB = resB[bkey][:, idx[b_]] - resB[bkey][:, idx[a]]
        bd = bA - bB
        out[f"{b_}-{a}"] = {"diff_pct_points": r(dA - dB), "ci": [r(x) for x in ci(bd)], "boot_p_two_sided": r(2 * min(np.mean(bd <= 0), np.mean(bd >= 0)), 4)}
    # difference in the landmark-0 level as well (lethality descriptor)
    bd0 = resA[bkey][:, idx[0]] - resB[bkey][:, idx[0]]
    out["level_0"] = {"diff_pct_points": r(resA[ptkey][idx[0]] - resB[ptkey][idx[0]]), "ci": [r(x) for x in ci(bd0)]}
    cA, _ = hinge_boot(resA[bkey]); cB, _ = hinge_boot(resB[bkey]); m = min(len(cA), len(cB))
    out["plateau_end_diff"] = {"diff_months": hinge_point(resA[ptkey]) - hinge_point(resB[ptkey]), "ci": [r(x) for x in ci(cA[:m] - cB[:m])]}
    return out

def at_risk_at(df, s):
    """patients at risk at landmark s on the half-day grid (exactly the landmark-s risk set used by curves())"""
    k, ev, KMAX = lm.prep(df); lo = int(np.floor(2 * lm.MD * s)) + 1
    return df[k >= lo].reset_index(drop=True)

def fmt(x, d=1): return "NA" if x is None else f"{x:.{d}f}"
def fci(e, d=1): return f"{fmt(e['pct'] if 'pct' in e else e.get('pct_points', e.get('diff_pct_points')), d)} ({fmt(e['ci'][0], d)} to {fmt(e['ci'][1], d)})"

def check_monotone(res, tag):
    nr = res["nrisk"]; ok = bool(np.all(np.diff(nr) <= 0))
    rng_ok = bool(np.nanmin(res["cif"]) >= 0 and np.nanmax(res["cif"]) <= 100 and np.all(res["allc"][~np.isnan(res["allc"])] >= res["cif"][~np.isnan(res["cif"])] - 1e-9))
    if not (ok and rng_ok): print(f"  SANITY WARNING {tag}: monotone={ok} range={rng_ok}")
    return ok and rng_ok

# ------------------------------------------------------------------ data
print("loading cohorts ...")
nen = load_cols("NEN_all"); nec = lm.nec_main(nen).reset_index(drop=True)
lung = load_cols("LUNG"); ac = load_cols("AC"); panc_all = load_cols("PANC"); panc = lm.panc_main(panc_all).reset_index(drop=True)
D = lambda d: d[d.stage == "Distant"].reset_index(drop=True)
R = lambda d: d[d.stage == "Regional"].reset_index(drop=True)
coh = {"NEC": D(nec), "SCLC": D(lm.sclc_main(lung)), "AC": D(ac), "PANC": D(panc), "LCNEC": D(lung[lung["hist"] == "8013"])}
coh["CRC_SC"] = coh["NEC"][coh["NEC"]["hist"].isin(lm.SCLC_CODES)].reset_index(drop=True)
coh["CRC_LC"] = coh["NEC"][coh["NEC"]["hist"] == "8013"].reset_index(drop=True)
coh["SCLC_breitling"] = D(lung[lung["hist"].isin(["8002", "8040"] + lm.SCLC_CODES)])
coh["SCLC_8041"] = D(lung[lung["hist"] == "8041"])
for t in ("NEC", "SCLC", "AC", "PANC", "LCNEC"):
    coh[t + "_chemo"] = coh[t][coh[t].chemo == 1].reset_index(drop=True)
    coh[t + "_chemo_c3"] = at_risk_at(coh[t + "_chemo"], 3)
    if t != "LCNEC": coh[t + "_nochemo"] = coh[t][coh[t].chemo == 0].reset_index(drop=True)
coh["NEC_regional"] = R(nec); coh["AC_regional"] = R(ac); coh["PANC_regional"] = R(panc); coh["SCLC_regional"] = R(lm.sclc_main(lung))
del ac, panc_all, panc, lung
ac_cyt = load_cols("AC_cyt"); coh["AC_cyt"] = D(ac_cyt); del ac_cyt
lung_cyt = load_cols("LUNG_cyt"); coh["SCLC_cyt"] = D(lm.sclc_main(lung_cyt)); del lung_cyt
panc_cyt = load_cols("PANC_cyt"); coh["PANC_cyt"] = D(lm.panc_main(panc_cyt)); del panc_cyt
lung_hist_distant = {}
print({k: len(v) for k, v in coh.items()}); print(f"load time {time.time() - T0:.0f}s")
# overlap between NEC and comparators (shared pids), descriptive
overlap = {t: int(coh["NEC"].pid.isin(coh[t].pid).sum()) for t in ("SCLC", "AC", "PANC", "LCNEC")}
print("distant NEC patients shared with comparators:", overlap)

# ------------------------------------------------------------------ bootstrap everything
print("\nbootstrapping (B=%d) ..." % B)
res = {}
for t in GI:
    t1 = time.time(); res[t] = point_and_boot(coh[t], SEED0 + GI[t]); check_monotone(res[t], t)
    print(f"  {t:16s} n={res[t]['n']:7d} seed={SEED0 + GI[t]}  CIF0={res[t]['cif'][0]:.2f} CIF6={res[t]['cif'][6]:.2f} CIF12={res[t]['cif'][12]:.2f}  plateau={hinge_point(res[t]['cif']):.2f}  ({time.time() - t1:.0f}s)")

# ------------------------------------------------------------------ reproduction check
LM5 = [0, 3, 6, 12, 24]
nec5 = [round(float(res["NEC"]["cif"][s]), 2) for s in LM5]; nr5 = [int(res["NEC"]["nrisk"][s]) for s in LM5]
def three(t): return [round(float(res[t]["cif"][s]), 2) for s in (0, 6, 12)]
checks = {
    "NEC_cif_0_3_6_12_24": (nec5, [72.29, 71.62, 72.33, 60.03, 30.51]), "NEC_nrisk": (nr5, [908, 553, 396, 169, 48]),
    "NEC_plateau": (hinge_point(res["NEC"]["cif"]), 6.75),
    "SCLC_cif_0_6_12": (three("SCLC"), [69.98, 69.49, 61.93]), "SCLC_plateau": (hinge_point(res["SCLC"]["cif"]), 8.5),
    "AC_cif_0_6_12": (three("AC"), [40.59, 33.92, 34.41]), "AC_plateau": (hinge_point(res["AC"]["cif"]), 0.0),
    "PANC_cif_0_6_12": (three("PANC"), [77.12, 69.72, 65.38]), "PANC_plateau": (hinge_point(res["PANC"]["cif"]), 0.0),
    "n_distant": ([len(coh["NEC"]), len(coh["SCLC"]), len(coh["AC"]), len(coh["PANC"]), len(coh["LCNEC"])], [908, 77707, 147014, 75634, 3185]),
    "n_main_nec": (len(nec), 1417), "n_2x2_crc": ([len(coh["CRC_SC"]), len(coh["CRC_LC"])], [599, 309]),
    "nrisk_6_12_comparators": ([int(res[t]["nrisk"][s]) for t in ("SCLC", "AC", "PANC") for s in (6, 12)], [37152, 17826, 98550, 75885, 26125, 12373]),
}
print("\nREPRODUCTION CHECK"); allok = True
for kk, (mine, ref) in checks.items():
    ok = bool(np.allclose(np.asarray(mine, float), np.asarray(ref, float), atol=0.011))
    allok &= ok; print(f"  {kk}: mine={mine} anchor={ref} -> {'OK' if ok else 'MISMATCH'}")
if not allok: config.anchor_fail("Reproduction anchors not reproduced; stopping.")

out = {"id": "A11", "B": B, "seed_rule": "20260928 + group index", "group_index": GI, "window_months": WINDOW, "landmarks": LM25.tolist(),
       "plateau_grid": "0-18 month landmarks, broken-stick, c in [0,15] step 0.25", "n": {t: int(len(coh[t])) for t in coh},
       "overlap_NEC_with_comparators_distant": overlap,
       "reproduction_check": {kk: {"mine": v[0], "anchor": v[1]} for kk, v in checks.items()}, "reproduction_ok": allok}

# ------------------------------------------------------------------ (1) Table 3
print("\n(1) Table 3 (distant, 12-month window): cancer death and all-cause death")
MAN3 = {"NEC": {"lm": {0: (72.3, 69.2, 75.2), 6: (72.3, 67.6, 76.8), 12: (60.0, 52.4, 67.5)}, "ch60": (0.03, -4.0, 4.1), "plateau": (6.75, 4.00, 9.75), "diff": None},
        "SCLC": {"lm": {0: (70.0, 69.7, 70.3), 6: (69.5, 69.0, 70.0), 12: (61.9, 61.3, 62.7)}, "ch60": (-0.5, -0.9, -0.04), "plateau": (8.50, 8.00, 8.75), "diff": (0.5, -3.5, 4.7)},
        "AC": {"lm": {0: (40.6, 40.3, 40.8), 6: (33.9, 33.6, 34.2), 12: (34.4, 34.1, 34.7)}, "ch60": (-6.7, -6.9, -6.4), "plateau": (0.0, 0.0, 0.0), "diff": (6.7, 2.6, 10.8)},
        "PANC": {"lm": {0: (77.1, 76.8, 77.4), 6: (69.7, 69.2, 70.3), 12: (65.4, 64.6, 66.2)}, "ch60": (-7.4, -7.9, -6.9), "plateau": (0.0, 0.0, 0.0), "diff": (7.4, 3.3, 11.5)},
        "LCNEC": {"lm": {0: (66.8, None, None)}, "ch60": (-5.6, -7.8, -3.5), "plateau": None, "diff": None}}
table3 = {}
for t in ("NEC", "SCLC", "AC", "PANC"):
    table3[t] = {"cancer": summarize(res[t]), "allcause": summarize(res[t], "allc", "ba"), "manuscript": MAN3[t]}
    if t != "NEC":
        table3[t]["nec_minus_comparator_cancer"] = paired(res["NEC"], res[t])
        table3[t]["nec_minus_comparator_allcause"] = paired(res["NEC"], res[t], "ba", "allc")
    c = table3[t]["cancer"]; a = table3[t]["allcause"]
    print(f"  {t:5s} cancer: LM0 {fci(c['landmark']['0'])} LM6 {fci(c['landmark']['6'])} LM12 {fci(c['landmark']['12'])} | 6-0 {fci(c['change']['6-0'])} 12-6 {fci(c['change']['12-6'])} 12-0 {fci(c['change']['12-0'])} | plateau {c['plateau_end']['months']} [{c['plateau_end']['ci'][0]}-{c['plateau_end']['ci'][1]}]")
    print(f"        allcause: LM0 {fci(a['landmark']['0'])} LM6 {fci(a['landmark']['6'])} LM12 {fci(a['landmark']['12'])} | 6-0 {fci(a['change']['6-0'])} 12-6 {fci(a['change']['12-6'])} 12-0 {fci(a['change']['12-0'])} | plateau {a['plateau_end']['months']} [{a['plateau_end']['ci'][0]}-{a['plateau_end']['ci'][1]}]")
    if t != "NEC":
        p = table3[t]["nec_minus_comparator_cancer"]; q = table3[t]["nec_minus_comparator_allcause"]
        print(f"        NEC minus {t}: cancer 6-0 {fci(p['6-0'])} 12-6 {fci(p['12-6'])} 12-0 {fci(p['12-0'])} | allcause 6-0 {fci(q['6-0'])} 12-6 {fci(q['12-6'])} 12-0 {fci(q['12-0'])}")
out["table3"] = table3

# ------------------------------------------------------------------ (2) LCNEC + 2x2
print("\n(2) Lung LCNEC as fourth comparator; 2x2 morphology x organ")
lc = {"cancer": summarize(res["LCNEC"]), "allcause": summarize(res["LCNEC"], "allc", "ba"), "manuscript": MAN3["LCNEC"],
      "nec_minus_comparator_cancer": paired(res["NEC"], res["LCNEC"]), "nec_minus_comparator_allcause": paired(res["NEC"], res["LCNEC"], "ba", "allc")}
print(f"  LCNEC cancer: LM0 {fci(lc['cancer']['landmark']['0'])} LM6 {fci(lc['cancer']['landmark']['6'])} | 6-0 {fci(lc['cancer']['change']['6-0'])} plateau {lc['cancer']['plateau_end']['months']} {lc['cancer']['plateau_end']['ci']} | NEC minus LCNEC 6-0 {fci(lc['nec_minus_comparator_cancer']['6-0'])}")
cells = {"colorectal_SC": "CRC_SC", "colorectal_LC": "CRC_LC", "lung_SC": "SCLC", "lung_LC": "LCNEC"}
two = {"cells": {}, "paired": {}}
for nm, t in cells.items():
    two["cells"][nm] = {"group": t, "cancer": summarize(res[t]), "allcause": summarize(res[t], "allc", "ba")}
    c = two["cells"][nm]["cancer"]
    print(f"  {nm:14s} n={c['n']:6d} LM0 {fci(c['landmark']['0'])} LM6 {fci(c['landmark']['6'])} 6-0 {fci(c['change']['6-0'])} 12-0 {fci(c['change']['12-0'])} plateau {c['plateau_end']['months']} [{c['plateau_end']['ci'][0]}-{c['plateau_end']['ci'][1]}] share0={c['plateau_end']['share_zero']}")
for nm, (a, b_) in {"colorectalSC_minus_lungSC": ("CRC_SC", "SCLC"), "colorectalLC_minus_lungLC": ("CRC_LC", "LCNEC"),
                    "colorectalSC_minus_colorectalLC": ("CRC_SC", "CRC_LC"), "lungSC_minus_lungLC": ("SCLC", "LCNEC")}.items():
    two["paired"][nm] = {"cancer": paired(res[a], res[b_]), "allcause": paired(res[a], res[b_], "ba", "allc")}
    print(f"  {nm:32s} 6-0 diff {fci(two['paired'][nm]['cancer']['6-0'])}  12-0 diff {fci(two['paired'][nm]['cancer']['12-0'])}")
# morphology x organ interaction on the 6-0 change: (CRC_SC - CRC_LC) - (SCLC - LCNEC)
i6 = lambda t: res[t]["bc"][:, 6] - res[t]["bc"][:, 0]; p6 = lambda t: res[t]["cif"][6] - res[t]["cif"][0]
inter_b = (i6("CRC_SC") - i6("CRC_LC")) - (i6("SCLC") - i6("LCNEC")); inter = (p6("CRC_SC") - p6("CRC_LC")) - (p6("SCLC") - p6("LCNEC"))
two["interaction_6_0"] = {"definition": "(colorectal SC - colorectal LC) - (lung SC - lung LC), 6-0 change, pct points", "estimate": r(inter), "ci": [r(x) for x in ci(inter_b)]}
# morphology main effect pooled over organs (average of SC-minus-LC within organ)
me_b = 0.5 * ((i6("CRC_SC") - i6("CRC_LC")) + (i6("SCLC") - i6("LCNEC"))); me = 0.5 * ((p6("CRC_SC") - p6("CRC_LC")) + (p6("SCLC") - p6("LCNEC")))
two["morphology_effect_6_0"] = {"definition": "mean over organs of (SC - LC) 6-0 change, pct points", "estimate": r(me), "ci": [r(x) for x in ci(me_b)]}
oe_b = 0.5 * ((i6("CRC_SC") - i6("SCLC")) + (i6("CRC_LC") - i6("LCNEC"))); oe = 0.5 * ((p6("CRC_SC") - p6("SCLC")) + (p6("CRC_LC") - p6("LCNEC")))
two["organ_effect_6_0"] = {"definition": "mean over morphologies of (colorectal - lung) 6-0 change, pct points", "estimate": r(oe), "ci": [r(x) for x in ci(oe_b)]}
print(f"  interaction (SC-LC | colorectal) - (SC-LC | lung): {fmt(two['interaction_6_0']['estimate'])} {two['interaction_6_0']['ci']}; morphology effect {fmt(me)} {two['morphology_effect_6_0']['ci']}; organ effect {fmt(oe)} {two['organ_effect_6_0']['ci']}")
out["lcnec"] = lc; out["two_by_two"] = two

# ------------------------------------------------------------------ (3) SCLC definitions
print("\n(3) SCLC definitions")
defs = {}
for t, desc in (("SCLC", "main: 8041-8045"), ("SCLC_breitling", "Breitling: 8002/8040-8045"), ("SCLC_8041", "8041 only")):
    defs[t] = {"definition": desc, "cancer": summarize(res[t]), "nec_minus_comparator_cancer": paired(res["NEC"], res[t])}
    c = defs[t]["cancer"]; p = defs[t]["nec_minus_comparator_cancer"]
    print(f"  {desc:28s} n={c['n']:6d} LM0 {fci(c['landmark']['0'])} LM6 {fci(c['landmark']['6'])} 6-0 {fci(c['change']['6-0'])} plateau {c['plateau_end']['months']} [{c['plateau_end']['ci'][0]}-{c['plateau_end']['ci'][1]}] | NEC minus: 6-0 {fci(p['6-0'])}")
out["sclc_definitions"] = defs

# ------------------------------------------------------------------ (4) immortal time
print("\n(4) Immortal-time demonstration (chemo recorded yes)")
MAN_CH = {"NEC": 7.1, "SCLC": 7.0, "AC": 4.7, "PANC": 1.5}
imm = {}
for t in ("NEC", "SCLC", "AC", "PANC", "LCNEC"):
    fy = summarize(res[t + "_chemo"]); c3 = summarize(res[t + "_chemo_c3"], plateau=False)
    # conditioned on 3-month survival: landmarks 3/6/12 from the subset at risk at 3 (0-month entry is the same risk set)
    j3, j6, j12 = 3, 6, 12; bm = res[t + "_chemo_c3"]["bc"]; pt = res[t + "_chemo_c3"]["cif"]
    cond = {"n_alive_3m": int(res[t + "_chemo_c3"]["n"]), "seed": res[t + "_chemo_c3"]["seed"],
            "lm3": {"n_at_risk": int(res[t + "_chemo_c3"]["nrisk"][j3]), "pct": r(pt[j3]), "ci": [r(x) for x in ci(bm[:, j3])]},
            "lm6": {"n_at_risk": int(res[t + "_chemo_c3"]["nrisk"][j6]), "pct": r(pt[j6]), "ci": [r(x) for x in ci(bm[:, j6])]},
            "lm12": {"n_at_risk": int(res[t + "_chemo_c3"]["nrisk"][j12]), "pct": r(pt[j12]), "ci": [r(x) for x in ci(bm[:, j12])]},
            "change_6-3": {"pct_points": r(pt[j6] - pt[j3]), "ci": [r(x) for x in ci(bm[:, j6] - bm[:, j3])]},
            "change_12-3": {"pct_points": r(pt[j12] - pt[j3]), "ci": [r(x) for x in ci(bm[:, j12] - bm[:, j3])]}}
    # same contrasts in the full-cohort (all patients, unselected by chemo) for reference
    full = res[t]; cond_all = {"change_6-3": {"pct_points": r(full["cif"][6] - full["cif"][3]), "ci": [r(x) for x in ci(full["bc"][:, 6] - full["bc"][:, 3])]},
                               "change_12-3": {"pct_points": r(full["cif"][12] - full["cif"][3]), "ci": [r(x) for x in ci(full["bc"][:, 12] - full["bc"][:, 3])]},
                               "lm3_pct": r(full["cif"][3]), "n_at_risk_3": int(full["nrisk"][3])}
    imm[t] = {"chemo_yes_from_diagnosis": fy, "manuscript_6_0_change": MAN_CH.get(t), "chemo_yes_conditioned_3m": cond, "all_patients_reference": cond_all,
              "chemo_share_at_diagnosis": r(len(coh[t + "_chemo"]) / len(coh[t]), 4)}
    if t != "LCNEC":
        imm[t]["chemo_no_or_unknown_from_diagnosis"] = summarize(res[t + "_nochemo"])
        nc = imm[t]["chemo_no_or_unknown_from_diagnosis"]
    print(f"  {t:5s} chemo-yes n={fy['n']:6d}: LM0 {fci(fy['landmark']['0'])} LM6 {fci(fy['landmark']['6'])} 6-0 {fci(fy['change']['6-0'])} (ms {MAN_CH.get(t)}) plateau {fy['plateau_end']['months']} | conditioned on 3m (n={cond['n_alive_3m']}): LM3 {fci(cond['lm3'])} 6-3 {fci(cond['change_6-3'])} 12-3 {fci(cond['change_12-3'])} | all pts 6-3 {fci(cond_all['change_6-3'])}"
          + (f" | chemo-no n={nc['n']}: LM0 {fmt(nc['landmark']['0']['pct'])} 6-0 {fci(nc['change']['6-0'])}" if t != "LCNEC" else ""))
out["immortal_time"] = imm

# ------------------------------------------------------------------ (5) cytology-inclusive comparators
print("\n(5) Cytology-inclusive comparators")
cyt = {}
for t, base, desc in (("AC_cyt", "AC", "colorectal AC, histology or cytology"), ("SCLC_cyt", "SCLC", "SCLC 8041-8045, histology or cytology"), ("PANC_cyt", "PANC", "pancreatic AC excl. islet, histology or cytology")):
    cyt[t] = {"definition": desc, "n_histology_only": int(len(coh[base])), "cancer": summarize(res[t]), "nec_minus_comparator_cancer": paired(res["NEC"], res[t])}
    c = cyt[t]["cancer"]; p = cyt[t]["nec_minus_comparator_cancer"]
    print(f"  {desc:45s} n={c['n']:6d} (hist-only {len(coh[base])}) LM0 {fci(c['landmark']['0'])} 6-0 {fci(c['change']['6-0'])} plateau {c['plateau_end']['months']} [{c['plateau_end']['ci'][0]}-{c['plateau_end']['ci'][1]}] | NEC minus: 6-0 {fci(p['6-0'])}")
out["cytology_inclusive"] = cyt

# ------------------------------------------------------------------ (6) regional
print("\n(6) Regional stage auxiliary comparison")
reg = {}
for t in ("NEC_regional", "AC_regional", "PANC_regional", "SCLC_regional"):
    reg[t] = {"cancer": summarize(res[t]), "allcause": summarize(res[t], "allc", "ba")}
    if t != "NEC_regional": reg[t]["nec_minus_comparator_cancer"] = paired(res["NEC_regional"], res[t])
    c = reg[t]["cancer"]
    print(f"  {t:14s} n={c['n']:6d} LM0 {fci(c['landmark']['0'])} LM6 {fci(c['landmark']['6'])} 6-0 {fci(c['change']['6-0'])} plateau {c['plateau_end']['months']} [{c['plateau_end']['ci'][0]}-{c['plateau_end']['ci'][1]}]"
          + (f" | NEC minus: 6-0 {fci(reg[t]['nec_minus_comparator_cancer']['6-0'])}" if t != "NEC_regional" else ""))
reg["manuscript"] = {"NEC_minus_AC_6_0": (2.7, -2.7, 8.2), "PANC_regional_6_0": 4.0, "PANC_regional_plateau": 12.5}
out["regional"] = reg

# ------------------------------------------------------------------ (7) figures
print("\n(7) figures")
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": "#888", "axes.labelcolor": "#222", "xtick.color": "#444", "ytick.color": "#444"})
def band(ax, t, label, color, lms=LM25, key="cif", bkey="bc"):
    rr = res[t]; y = rr[key]; lo = np.nanpercentile(rr[bkey], 2.5, axis=0); hi = np.nanpercentile(rr[bkey], 97.5, axis=0)
    ax.fill_between(lms, lo, hi, color=color, alpha=0.18, linewidth=0)
    ax.plot(lms, y, color=color, lw=2, label=label)
    small = rr["nrisk"] < 100
    ax.plot(lms[~small], y[~small], "o", color=color, ms=4, mec="white", mew=0.8)
    if small.any(): ax.plot(lms[small], y[small], "o", color="white", ms=4, mec=color, mew=1.2)

fig, axes = plt.subplots(1, 2, figsize=(12, 4.8), gridspec_kw={"width_ratios": [1.45, 1]})
ax = axes[0]
for t in ("NEC", "SCLC", "AC", "PANC"):
    band(ax, t, f"{LABEL[t]} (n={res[t]['n']:,})", COL[t])
ax.axvspan(0, 6, color="#999", alpha=0.08, lw=0)
ax.set_xlabel("Landmark: months survived since diagnosis"); ax.set_ylabel("Cancer death within the next 12 months, % (Aalen-Johansen)")
ax.set_xticks(range(0, 25, 3)); ax.set_ylim(0, 100); ax.grid(axis="y", color="#eee"); ax.legend(frameon=False, loc="lower left", fontsize=8)
ax.set_title("a  Distant stage, monthly landmarks (bands: 95% bootstrap CI; open markers: n at risk < 100)", loc="left", fontsize=8.5)
# n at risk text
txt = "n at risk at 0/6/12/24 months: " + "; ".join(f"{t} {int(res[t]['nrisk'][0]):,}/{int(res[t]['nrisk'][6]):,}/{int(res[t]['nrisk'][12]):,}/{int(res[t]['nrisk'][24]):,}" for t in ("NEC", "SCLC", "AC", "PANC"))
ax.text(0, -0.17, txt, transform=ax.transAxes, fontsize=6.5, color="#555")
ax = axes[1]
rows = []
for t in ("SCLC", "AC", "PANC", "LCNEC"):
    p = (table3[t] if t in table3 else lc)["nec_minus_comparator_cancer"]
    for ch in ("6-0", "12-6", "12-0"): rows.append((t, ch, p[ch]["diff_pct_points"], p[ch]["ci"]))
ypos = np.arange(len(rows))[::-1]
for yy, (t, ch, d, (lo, hi)) in zip(ypos, rows):
    ax.plot([lo, hi], [yy, yy], color=COL[t], lw=2); ax.plot(d, yy, "o", color=COL[t], ms=6, mec="white", mew=0.8)
ax.axvline(0, color="#999", lw=1, ls="--")
SHORT = {"SCLC": "SCLC", "AC": "Colorectal AC", "PANC": "Pancreatic AC", "LCNEC": "Lung LCNEC"}
ax.set_yticks(ypos); ax.set_yticklabels([f"{SHORT[t]}, {ch} mo" for t, ch, _, _ in rows], fontsize=7.5)
ax.set_xlabel("Percentage points (positive = NEC declines less)", fontsize=8)
ax.set_title("b  NEC minus comparator: change in risk (95% CI)", loc="left", fontsize=8.5)
ax.grid(axis="x", color="#eee")
fig.tight_layout(rect=(0, 0.02, 0.99, 1)); f1 = FIGS + "A11_four_tumour_curves_forest.png"; fig.savefig(f1, dpi=300); plt.close(fig)

fig, axes = plt.subplots(2, 2, figsize=(9, 7), sharex=True, sharey=True)
order = [("colorectal_SC", "CRC_SC", COL["NEC"]), ("lung_SC", "SCLC", COL["SCLC"]), ("colorectal_LC", "CRC_LC", COL["NEC"]), ("lung_LC", "LCNEC", COL["LCNEC"])]
for ax, (nm, t, col) in zip(axes.ravel(), order):
    band(ax, t, LABEL[t], col); ax.axvspan(0, 6, color="#999", alpha=0.08, lw=0)
    c = two["cells"][nm]["cancer"]
    ax.set_title(f"{LABEL[t]}  n={c['n']:,}", loc="left", fontsize=9)
    ax.text(0.02, 0.06, f"6-0 change {fci(c['change']['6-0'])} pp\nplateau end {c['plateau_end']['months']:.2f} mo ({c['plateau_end']['ci'][0]:.2f}-{c['plateau_end']['ci'][1]:.2f})",
            transform=ax.transAxes, fontsize=7.5, color="#333", va="bottom")
    ax.set_ylim(0, 100); ax.grid(axis="y", color="#eee"); ax.set_xticks(range(0, 25, 3))
for ax in axes[1]: ax.set_xlabel("Landmark: months survived since diagnosis")
for ax in axes[:, 0]: ax.set_ylabel("Cancer death within next 12 months, %")
fig.suptitle("Distant stage: morphology (rows: small-cell / large-cell) x organ (columns: colorectal / lung)", fontsize=10, x=0.02, ha="left")
fig.tight_layout(); f2 = FIGS + "A11_two_by_two.png"; fig.savefig(f2, dpi=300); plt.close(fig)
out["figures"] = [config.relpath(f1), config.relpath(f2)]

# ------------------------------------------------------------------ interpretation flags
flags = {t: {"6-0": table3[t]["cancer"]["change"]["6-0"], "plateau_end": table3[t]["cancer"]["plateau_end"]["months"]} for t in table3}
flags["LCNEC"] = {"6-0": lc["cancer"]["change"]["6-0"], "plateau_end": lc["cancer"]["plateau_end"]["months"]}
flags["CRC_SC"] = {"6-0": two["cells"]["colorectal_SC"]["cancer"]["change"]["6-0"], "plateau_end": two["cells"]["colorectal_SC"]["cancer"]["plateau_end"]["months"]}
flags["CRC_LC"] = {"6-0": two["cells"]["colorectal_LC"]["cancer"]["change"]["6-0"], "plateau_end": two["cells"]["colorectal_LC"]["cancer"]["plateau_end"]["months"]}
out["explanation_flags"] = flags
out["runtime_minutes"] = round((time.time() - T0) / 60, 2)
with open(RES + "A11.json", "w") as f: json.dump(out, f, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))

# ------------------------------------------------------------------ markdown
L = []
L.append("# A11 Cross-tumour comparison (distant stage, 12-month window, B = 2,000, seeds 20260928 + group index)\n")
L.append("## Reproduction check (mine vs manuscript / orchestrator anchors)\n")
L.append("| Quantity | This script | Anchor | OK |\n|---|---|---|---|")
for kk, (mine, ref) in checks.items(): L.append(f"| {kk} | {mine} | {ref} | {'yes' if np.allclose(np.asarray(mine, float), np.asarray(ref, float), atol=0.011) else 'NO'} |")
L.append("")
L.append("## 1. Table 3: distant-stage NEC vs three comparators (cancer death = Aalen-Johansen cumulative incidence within 12 months after each landmark)\n")
L.append("| Quantity | Colorectal NEC | SCLC | Colorectal AC | Pancreatic AC |\n|---|---|---|---|---|")
T = table3
def rowlm(s, key):
    cells_ = []
    for t in ("NEC", "SCLC", "AC", "PANC"):
        e = T[t][key]["landmark"][str(s)]; cells_.append(f"{fmt(e['pct'])} ({fmt(e['ci'][0])}-{fmt(e['ci'][1])}); n={e['n_at_risk']:,}")
    return " | ".join(cells_)
def rowch(ch, key):
    return " | ".join(fci(T[t][key]["change"][ch]) for t in ("NEC", "SCLC", "AC", "PANC"))
def rowpl(key):
    return " | ".join(f"{T[t][key]['plateau_end']['months']:.2f} ({T[t][key]['plateau_end']['ci'][0]:.2f}-{T[t][key]['plateau_end']['ci'][1]:.2f})" for t in ("NEC", "SCLC", "AC", "PANC"))
for s in (0, 6, 12): L.append(f"| Cancer death from landmark {s} mo, % (95% CI); n at risk | {rowlm(s, 'cancer')} |")
for ch in ("6-0", "12-6", "12-0"): L.append(f"| Change {ch} mo, pp (95% CI) | {rowch(ch, 'cancer')} |")
L.append(f"| Plateau end, months (95% CI) | {rowpl('cancer')} |")
for ch in ("6-0", "12-6", "12-0"): L.append(f"| NEC minus comparator, change {ch}, pp (95% paired bootstrap CI) | reference | " + " | ".join(fci(T[t]['nec_minus_comparator_cancer'][ch]) for t in ("SCLC", "AC", "PANC")) + " |")
L.append(f"| Manuscript Table 3: 6-0 change | " + " | ".join(f"{MAN3[t]['ch60'][0]} ({MAN3[t]['ch60'][1]} to {MAN3[t]['ch60'][2]})" for t in ("NEC", "SCLC", "AC", "PANC")) + " |")
L.append(f"| Manuscript Table 3: NEC minus comparator 6-0 | reference | " + " | ".join(f"{MAN3[t]['diff'][0]} ({MAN3[t]['diff'][1]} to {MAN3[t]['diff'][2]})" for t in ("SCLC", "AC", "PANC")) + " |")
L.append(f"| Manuscript Table 3: plateau end | " + " | ".join(f"{MAN3[t]['plateau'][0]} ({MAN3[t]['plateau'][1]}-{MAN3[t]['plateau'][2]})" for t in ("NEC", "SCLC", "AC", "PANC")) + " |")
L.append("\n### 1b. Same quantities for all-cause death (Kaplan-Meier within the 12-month window)\n")
L.append("| Quantity | Colorectal NEC | SCLC | Colorectal AC | Pancreatic AC |\n|---|---|---|---|---|")
for s in (0, 6, 12): L.append(f"| All-cause death from landmark {s} mo, % (95% CI); n at risk | {rowlm(s, 'allcause')} |")
for ch in ("6-0", "12-6", "12-0"): L.append(f"| Change {ch} mo, pp (95% CI) | {rowch(ch, 'allcause')} |")
L.append(f"| Plateau end, months (95% CI) | {rowpl('allcause')} |")
for ch in ("6-0", "12-6", "12-0"): L.append(f"| NEC minus comparator, change {ch}, pp (95% CI) | reference | " + " | ".join(fci(T[t]['nec_minus_comparator_allcause'][ch]) for t in ("SCLC", "AC", "PANC")) + " |")
L.append(f"\nDistant NEC patients also present in a comparator cohort (kept, as in the manuscript): {overlap}.\n")

L.append("## 2. Lung LCNEC (8013) as fourth comparator and the 2x2 morphology x organ layout\n")
c = lc["cancer"]; p = lc["nec_minus_comparator_cancer"]
L.append(f"Lung LCNEC distant n = {c['n']:,} (seed {c['seed']}): cancer death from 0 mo {fci(c['landmark']['0'])}, 6 mo {fci(c['landmark']['6'])}, 12 mo {fci(c['landmark']['12'])}; "
         f"change 6-0 {fci(c['change']['6-0'])} pp (manuscript -5.6, -7.8 to -3.5), 12-6 {fci(c['change']['12-6'])}, 12-0 {fci(c['change']['12-0'])}; plateau end {c['plateau_end']['months']:.2f} mo ({c['plateau_end']['ci'][0]:.2f}-{c['plateau_end']['ci'][1]:.2f}). "
         f"NEC minus LCNEC: 6-0 {fci(p['6-0'])}, 12-6 {fci(p['12-6'])}, 12-0 {fci(p['12-0'])} pp; all-cause 6-0 {fci(lc['nec_minus_comparator_allcause']['6-0'])}.\n")
L.append("| Cell | n | Cancer death from 0 mo, % (95% CI) | From 6 mo | Change 6-0, pp (95% CI) | Change 12-0, pp (95% CI) | Plateau end, mo (95% CI) | Bootstrap share with plateau end = 0 |\n|---|---|---|---|---|---|---|---|")
for nm in ("colorectal_SC", "colorectal_LC", "lung_SC", "lung_LC"):
    c = two["cells"][nm]["cancer"]
    L.append(f"| {nm.replace('_', ' ')} ({LABEL[two['cells'][nm]['group']]}) | {c['n']:,} | {fci(c['landmark']['0'])} | {fci(c['landmark']['6'])} | {fci(c['change']['6-0'])} | {fci(c['change']['12-0'])} | {c['plateau_end']['months']:.2f} ({c['plateau_end']['ci'][0]:.2f}-{c['plateau_end']['ci'][1]:.2f}) | {c['plateau_end']['share_zero']} |")
L.append("\n| Paired contrast (independent bootstraps paired by index) | 6-0 change difference, pp (95% CI) | 12-0 change difference, pp (95% CI) | Level at 0 mo difference, pp (95% CI) |\n|---|---|---|---|")
for nm, pr in two["paired"].items():
    q = pr["cancer"]; L.append(f"| {nm.replace('_', ' ')} | {fci(q['6-0'])} | {fci(q['12-0'])} | {fci(q['level_0'])} |")
L.append(f"| morphology x organ interaction (SC-LC in colorectal) - (SC-LC in lung) | {fmt(two['interaction_6_0']['estimate'])} ({fmt(two['interaction_6_0']['ci'][0])} to {fmt(two['interaction_6_0']['ci'][1])}) | | |")
L.append(f"| morphology main effect (mean SC-LC over organs) | {fmt(two['morphology_effect_6_0']['estimate'])} ({fmt(two['morphology_effect_6_0']['ci'][0])} to {fmt(two['morphology_effect_6_0']['ci'][1])}) | | |")
L.append(f"| organ main effect (mean colorectal-lung over morphologies) | {fmt(two['organ_effect_6_0']['estimate'])} ({fmt(two['organ_effect_6_0']['ci'][0])} to {fmt(two['organ_effect_6_0']['ci'][1])}) | | |\n")

L.append("## 3. SCLC definition sensitivity\n")
L.append("| Definition | n | From 0 mo, % (95% CI) | From 6 mo | Change 6-0, pp (95% CI) | Plateau end, mo (95% CI) | NEC minus SCLC, 6-0, pp (95% CI) |\n|---|---|---|---|---|---|---|")
for t in ("SCLC", "SCLC_breitling", "SCLC_8041"):
    c = defs[t]["cancer"]; p = defs[t]["nec_minus_comparator_cancer"]
    L.append(f"| {defs[t]['definition']} | {c['n']:,} | {fci(c['landmark']['0'])} | {fci(c['landmark']['6'])} | {fci(c['change']['6-0'])} | {c['plateau_end']['months']:.2f} ({c['plateau_end']['ci'][0]:.2f}-{c['plateau_end']['ci'][1]:.2f}) | {fci(p['6-0'])} |")
L.append("\nNo distant-stage 8040 records exist in the lung cohort; the Breitling definition adds only the 8002 records.\n")

L.append("## 4. Immortal-time demonstration: chemotherapy recorded 'yes'\n")
L.append("| Tumour | Chemo-yes n (share) | From diagnosis: 0 mo % | 6 mo % | Change 6-0, pp (95% CI) | Manuscript 6-0 | Plateau end, mo | Conditioned on 3-mo survival: n | 3 mo % (95% CI) | Change 6-3, pp (95% CI) | Change 12-3, pp (95% CI) | All patients (unselected): change 6-3 | Chemo no/unknown: 0 mo %; change 6-0 |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for t in ("NEC", "SCLC", "AC", "PANC", "LCNEC"):
    e = imm[t]; fy = e["chemo_yes_from_diagnosis"]; cd = e["chemo_yes_conditioned_3m"]; al = e["all_patients_reference"]
    nc = e.get("chemo_no_or_unknown_from_diagnosis")
    L.append(f"| {LABEL[t]} | {fy['n']:,} ({100 * e['chemo_share_at_diagnosis']:.1f}%) | {fmt(fy['landmark']['0']['pct'])} | {fmt(fy['landmark']['6']['pct'])} | {fci(fy['change']['6-0'])} | {e['manuscript_6_0_change']} | {fy['plateau_end']['months']:.2f} | {cd['n_alive_3m']:,} | {fci(cd['lm3'])} | {fci(cd['change_6-3'])} | {fci(cd['change_12-3'])} | {fci(al['change_6-3'])} | "
             + (f"{fmt(nc['landmark']['0']['pct'])}; {fci(nc['change']['6-0'])}" if nc else "-") + " |")
L.append("")
L.append("## 5. Cytology-inclusive comparators (diagnostic confirmation = positive histology OR positive exfoliative cytology)\n")
L.append("| Comparator | n (histology-only n) | From 0 mo, % (95% CI) | Change 6-0, pp (95% CI) | Plateau end, mo (95% CI) | NEC minus comparator, 6-0, pp (95% CI) |\n|---|---|---|---|---|---|")
for t in ("SCLC_cyt", "AC_cyt", "PANC_cyt"):
    c = cyt[t]["cancer"]; p = cyt[t]["nec_minus_comparator_cancer"]
    L.append(f"| {cyt[t]['definition']} | {c['n']:,} ({cyt[t]['n_histology_only']:,}) | {fci(c['landmark']['0'])} | {fci(c['change']['6-0'])} | {c['plateau_end']['months']:.2f} ({c['plateau_end']['ci'][0]:.2f}-{c['plateau_end']['ci'][1]:.2f}) | {fci(p['6-0'])} |")
L.append("")
L.append("## 6. Regional stage (auxiliary)\n")
L.append("| Cohort | n | From 0 mo, % (95% CI) | From 6 mo | Change 6-0, pp (95% CI) | Change 12-0 | Plateau end, mo (95% CI) | NEC minus comparator, 6-0, pp (95% CI) |\n|---|---|---|---|---|---|---|---|")
for t in ("NEC_regional", "AC_regional", "PANC_regional", "SCLC_regional"):
    c = reg[t]["cancer"]; p = reg[t].get("nec_minus_comparator_cancer")
    L.append(f"| {t.replace('_regional', '')} regional | {c['n']:,} | {fci(c['landmark']['0'])} | {fci(c['landmark']['6'])} | {fci(c['change']['6-0'])} | {fci(c['change']['12-0'])} | {c['plateau_end']['months']:.2f} ({c['plateau_end']['ci'][0]:.2f}-{c['plateau_end']['ci'][1]:.2f}) | {fci(p['6-0']) if p else 'reference'} |")
L.append("\nManuscript: NEC minus AC +2.7 (-2.7 to 8.2); pancreatic AC regional +4.0 pp, plateau end 12.5 mo.\n")
L.append("## 7. Figures\n")
L.append(f"- `{f1}`: (a) four-tumour monthly landmark curves with 95% bootstrap bands and (b) forest plot of the paired NEC-minus-comparator differences (6-0, 12-6, 12-0), including lung LCNEC.")
L.append(f"- `{f2}`: 2x2 layout (rows: small-cell / large-cell morphology; columns: colorectal / lung).\n")

# ---------------- interpretation (written from the computed numbers)
def g(t, key="cancer"): return (table3[t] if t in table3 else {"cancer": lc["cancer"]} if t == "LCNEC" else two["cells"][t])[key]
cs_, cl_, ls_, ll_ = (two["cells"][n]["cancer"] for n in ("colorectal_SC", "colorectal_LC", "lung_SC", "lung_LC"))
L.append("## Interpretation\n")
L.append(f"1. Table 3 reproduces: distant colorectal NEC shows no decline in the 12-month cancer-death probability over the first 6 months ({fci(T['NEC']['cancer']['change']['6-0'])} pp; plateau end {T['NEC']['cancer']['plateau_end']['months']:.2f} mo), SCLC declines negligibly ({fci(T['SCLC']['cancer']['change']['6-0'])}; plateau end {T['SCLC']['cancer']['plateau_end']['months']:.2f} mo), whereas colorectal AC ({fci(T['AC']['cancer']['change']['6-0'])}) and pancreatic AC ({fci(T['PANC']['cancer']['change']['6-0'])}) fall immediately (plateau end 0 in every bootstrap resample). The pre-specified paired contrasts are NEC minus SCLC {fci(T['SCLC']['nec_minus_comparator_cancer']['6-0'])}, minus AC {fci(T['AC']['nec_minus_comparator_cancer']['6-0'])}, minus pancreatic AC {fci(T['PANC']['nec_minus_comparator_cancer']['6-0'])} pp; all-cause death gives the same pattern (NEC minus AC {fci(T['AC']['nec_minus_comparator_allcause']['6-0'])}, minus pancreatic AC {fci(T['PANC']['nec_minus_comparator_allcause']['6-0'])}, minus SCLC {fci(T['SCLC']['nec_minus_comparator_allcause']['6-0'])}).")
L.append(f"2. The 12-6 and 12-0 contrasts separate the tumours less: between 6 and 12 months NEC ({fci(T['NEC']['cancer']['change']['12-6'])}) and SCLC ({fci(T['SCLC']['cancer']['change']['12-6'])}) both fall, so the 12-0 NEC-minus-SCLC difference is {fci(T['SCLC']['nec_minus_comparator_cancer']['12-0'])} and NEC-minus-AC {fci(T['AC']['nec_minus_comparator_cancer']['12-0'])}; the distinctive feature of NEC and SCLC is the shape in the first 6-9 months, not a permanent difference in conditional risk.")
L.append(f"3. Lethality explanation: pancreatic AC starts higher than NEC (level difference NEC minus PANC {fci(T['PANC']['nec_minus_comparator_cancer']['level_0'])} pp) yet drops {fmt(T['PANC']['cancer']['change']['6-0']['pct_points'])} pp by 6 months with plateau end 0 in all resamples. Under the pre-written rule this contrast does NOT support 'plateau because the disease is lethal'. Likewise lung LCNEC is lethal (0-mo {fmt(lc['cancer']['landmark']['0']['pct'])}%) but declines ({fci(lc['cancer']['change']['6-0'])}).")
L.append(f"4. NEC-cell-type explanation (plateau in every high-grade neuroendocrine carcinoma regardless of small/large-cell morphology): lung LCNEC shows no plateau ({fci(lc['cancer']['change']['6-0'])} pp; plateau end {lc['cancer']['plateau_end']['months']:.2f} mo, share of resamples with plateau end 0 = {lc['cancer']['plateau_end']['share_zero']}), and NEC minus LCNEC is {fci(lc['nec_minus_comparator_cancer']['6-0'])} pp. Colorectal large-cell NEC gives {fci(cl_['change']['6-0'])} pp (n = {cl_['n']}, wide). The 'NEC cell type' explanation is therefore not supported in the lung and cannot be assessed with precision in the colorectum.")
L.append(f"5. Small-cell-morphology explanation: both small-cell cells show a plateau (colorectal SC {fci(cs_['change']['6-0'])} pp, plateau end {cs_['plateau_end']['months']:.2f} mo; lung SC {fci(ls_['change']['6-0'])}, plateau end {ls_['plateau_end']['months']:.2f} mo), both large-cell cells decline (colorectal LC {fci(cl_['change']['6-0'])}; lung LC {fci(ll_['change']['6-0'])}). Paired: colorectal SC minus lung SC {fci(two['paired']['colorectalSC_minus_lungSC']['cancer']['6-0'])}; colorectal LC minus lung LC {fci(two['paired']['colorectalLC_minus_lungLC']['cancer']['6-0'])}; lung SC minus lung LC {fci(two['paired']['lungSC_minus_lungLC']['cancer']['6-0'])}; colorectal SC minus LC {fci(two['paired']['colorectalSC_minus_colorectalLC']['cancer']['6-0'])}. The pooled morphology effect (SC minus LC, averaged over organs) is {fmt(two['morphology_effect_6_0']['estimate'])} ({fmt(two['morphology_effect_6_0']['ci'][0])} to {fmt(two['morphology_effect_6_0']['ci'][1])}) pp, the organ effect {fmt(two['organ_effect_6_0']['estimate'])} ({fmt(two['organ_effect_6_0']['ci'][0])} to {fmt(two['organ_effect_6_0']['ci'][1])}) pp, and the morphology x organ interaction {fmt(two['interaction_6_0']['estimate'])} ({fmt(two['interaction_6_0']['ci'][0])} to {fmt(two['interaction_6_0']['ci'][1])}). Of the three explanations, 'small-cell morphology' is the only one consistent with all four cells, with the colorectal large-cell cell being the imprecise one; the lung contrast (SC vs LC, both n > 3,000) is the decisive and precise one.")
L.append(f"6. SCLC definition does not matter: Breitling (8002/8040-8045) and 8041-only give 6-0 changes of {fci(defs['SCLC_breitling']['cancer']['change']['6-0'])} and {fci(defs['SCLC_8041']['cancer']['change']['6-0'])} pp (plateau end {defs['SCLC_breitling']['cancer']['plateau_end']['months']:.2f} and {defs['SCLC_8041']['cancer']['plateau_end']['months']:.2f} mo). Cytology-inclusive comparators change nothing material (AC {fci(cyt['AC_cyt']['cancer']['change']['6-0'])}, SCLC {fci(cyt['SCLC_cyt']['cancer']['change']['6-0'])}, pancreatic AC {fci(cyt['PANC_cyt']['cancer']['change']['6-0'])}).")
L.append(f"7. Immortal time: restricting to chemotherapy-recorded patients from diagnosis manufactures a rise in all four tumours (6-0: NEC {fmt(imm['NEC']['chemo_yes_from_diagnosis']['change']['6-0']['pct_points'])}, SCLC {fmt(imm['SCLC']['chemo_yes_from_diagnosis']['change']['6-0']['pct_points'])}, AC {fmt(imm['AC']['chemo_yes_from_diagnosis']['change']['6-0']['pct_points'])}, pancreatic AC {fmt(imm['PANC']['chemo_yes_from_diagnosis']['change']['6-0']['pct_points'])} pp; manuscript +7.1/+7.0/+4.7/+1.5), and the chemo no/unknown group shows the mirror image (steep early decline). Conditioning on 3-month survival removes most of the artefact: 6-3 changes in chemo-yes patients are NEC {fci(imm['NEC']['chemo_yes_conditioned_3m']['change_6-3'])}, SCLC {fci(imm['SCLC']['chemo_yes_conditioned_3m']['change_6-3'])}, AC {fci(imm['AC']['chemo_yes_conditioned_3m']['change_6-3'])}, pancreatic AC {fci(imm['PANC']['chemo_yes_conditioned_3m']['change_6-3'])} pp, which now track the unselected cohorts (NEC {fci(imm['NEC']['all_patients_reference']['change_6-3'])}, AC {fci(imm['AC']['all_patients_reference']['change_6-3'])}). The cross-tumour plateau comparison must therefore be made in all distant patients, as the manuscript does.")
L.append(f"8. Regional stage does not discriminate: NEC minus AC {fci(reg['AC_regional']['nec_minus_comparator_cancer']['6-0'])} pp (manuscript +2.7, -2.7 to 8.2); pancreatic AC regional rises in the first 6 months ({fci(reg['PANC_regional']['cancer']['change']['6-0'])}, plateau end {reg['PANC_regional']['cancer']['plateau_end']['months']:.2f} mo), SCLC regional {fci(reg['SCLC_regional']['cancer']['change']['6-0'])}.")
L.append("\nPlateau-end CIs are percentiles of the bootstrap distribution of the grid-searched break point (0.25-month grid); non-grid values such as 4.24 arise from percentile interpolation between grid points.")
L.append("\nDeviations from the plan: none in method. Additions beyond the task list: chemo no/unknown groups, lung LCNEC chemo groups, SCLC regional, a 12-0 / level-0 / plateau-end paired contrast, and pooled morphology / organ / interaction contrasts for the 2x2 (all exploratory, same bootstrap). The conditioned-on-3-month chemo subsets were re-bootstrapped as their own cohorts (seeds 20260941-20260944, 20260957); the regional cohorts use seeds 20260948-20260951, so their CIs are not expected to match the manuscript's third decimal. Manuscript NEC CIs in Table 3 come from the main-analysis seed 20260917; here the monthly-program seed 20260928 is used, so NEC CIs may differ by up to a few tenths of a point.")
L.append(f"\nRuntime: {out['runtime_minutes']} min, single-threaded.")
with open(RES + "A11.md", "w") as f: f.write("\n".join(L) + "\n")
print(f"\nDONE  runtime {out['runtime_minutes']} min; JSON/MD written; figures: {f1}, {f2}")
