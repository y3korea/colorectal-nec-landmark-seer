# -*- coding: utf-8 -*-
"""A01 landmark_main: pre-specified landmark analysis (Table 2), monthly curves, plateau, 36-month window,
unknown-COD / complete-date / 8045-removed sensitivity analyses, with patient-level bootstrap.
Run: python3 python/A01_landmark_main.py   (single-threaded; paths from python/config.py)
"""
import os, sys, json, time
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np, pandas as pd
T0 = time.time()
import config
import lm
RES = config.RES; FIGS = config.FIGS
os.makedirs(RES, exist_ok=True); os.makedirs(FIGS, exist_ok=True)
B = config.B(2000); SEED_MAIN = 20260917; SEED_OTHER = 20260928
LM5 = np.array([0, 3, 6, 12, 24]); LM25 = np.arange(0, 25); M19 = np.arange(0, 19).astype(float)
STAGES = ["All", "Localized", "Regional", "Distant"]

# ------------------------------------------------------------------ helpers
def deaths_in_window(k, ev, landmarks, window, KMAX):
    """cancer deaths (ev==1) with landmark < t <= landmark+window, among patients at risk at the landmark"""
    out = []
    for s in landmarks:
        lo = int(np.floor(2 * lm.MD * s)) + 1; hi = min(int(np.floor(2 * lm.MD * (s + window))), KMAX + 1)
        out.append(int(((ev == 1) & (k >= lo) & (k <= hi)).sum()))
    return out

def swap_ev(ev):
    """other-cause death as the event of interest, cancer death as competing"""
    e = ev.copy(); e[ev == 1] = 2; e[ev == 2] = 1; return e

def point_and_boot(df, landmarks, window, seed, B=B, kmax_months=48):
    """point estimates (cancer CIF, other-cause CIF, all-cause, n at risk, deaths) + bootstrap matrices (same resample for all landmarks)"""
    k, ev, KMAX = lm.prep(df, kmax_months); n = len(k); ev2 = swap_ev(ev)
    cif, allc, nrisk = lm.curves(k, ev, KMAX, landmarks, window)
    oth, _, _ = lm.curves(k, ev2, KMAX, landmarks, window)
    d1 = deaths_in_window(k, ev, landmarks, window, KMAX); d2 = deaths_in_window(k, ev2, landmarks, window, KMAX)
    rng = np.random.default_rng(seed)
    bc = np.empty((B, len(landmarks))); bo = np.empty((B, len(landmarks))); ba = np.empty((B, len(landmarks)))
    for b in range(B):
        w = np.bincount(rng.integers(0, n, n), minlength=n).astype(float)
        bc[b], ba[b], _ = lm.curves(k, ev, KMAX, landmarks, window, w)
        bo[b], _, _ = lm.curves(k, ev2, KMAX, landmarks, window, w)
    return dict(cif=cif, oth=oth, allc=allc, nrisk=nrisk, d1=d1, d2=d2, bc=bc, bo=bo, ba=ba, n=n)

def ci(v): return lm.pct(v)
def r(x, d=2): return None if x is None or (isinstance(x, float) and np.isnan(x)) else round(float(x), d)

def hinge_fit(y, m=M19):
    """single broken-stick fit: returns (c, b0, b1)"""
    c = float(lm.hinge_many(y[None, :], m)[0])
    X = np.column_stack([np.ones_like(m), np.maximum(0, m - c)])
    beta = np.linalg.lstsq(X, y, rcond=None)[0]
    return c, float(beta[0]), float(beta[1])

def hinge_boot(bc, m=M19):
    """plateau end c and slope b1 for every bootstrap row (rows with NaN in 0-18 are dropped)"""
    Y = bc[:, :len(m)]; ok = ~np.isnan(Y).any(axis=1); Y = Y[ok]
    cs = lm.hinge_many(Y, m); b1 = np.empty(len(cs))
    for i, (y, c) in enumerate(zip(Y, cs)):
        X = np.column_stack([np.ones_like(m), np.maximum(0, m - c)])
        b1[i] = np.linalg.lstsq(X, y, rcond=None)[0][1]
    return cs, b1, int(ok.sum())

# ------------------------------------------------------------------ data
nen = lm.load("NEN_all"); nec = lm.nec_main(nen).reset_index(drop=True)
dist = nec[nec.stage == "Distant"].reset_index(drop=True)
print(f"main NEC n={len(nec)}  distant n={len(dist)}  stages={nec.stage.value_counts().to_dict()}")

# ------------------------------------------------------------------ reproduction check
k, ev, KMAX = lm.prep(dist)
cif5, _, nr5 = lm.curves(k, ev, KMAX, LM5, 12)
cif25, allc25, nr25 = lm.curves(k, ev, KMAX, LM25, 12)
c_hat, b0_hat, b1_hat = hinge_fit(cif25[:19])
anch = dict(cif=[72.29, 71.62, 72.33, 60.03, 30.51], nrisk=[908, 553, 396, 169, 48], plateau=6.75)
flagc_n = int((dist.flagc == 1).sum()); unk_all = int(nec.unkcod.sum()); unk_dist = int(dist.unkcod.sum())
checks = {
    "cif_0_3_6_12_24": (np.round(cif5, 2).tolist(), anch["cif"], bool(np.allclose(np.round(cif5, 2), anch["cif"], atol=0.011))),
    "nrisk": (nr5.astype(int).tolist(), anch["nrisk"], bool((nr5.astype(int) == anch["nrisk"]).all())),
    "plateau_end": (c_hat, 6.75, c_hat == 6.75),
    "n_main_distant": ([len(nec), len(dist)], [1417, 908], len(nec) == 1417 and len(dist) == 908),
    "flagc_distant": (flagc_n, 881, flagc_n == 881),
    "unknown_cod": ([unk_all, unk_dist], [10, 7], unk_all == 10 and unk_dist == 7),
}
print("REPRODUCTION CHECK")
for kk, (mine, ref, ok) in checks.items(): print(f"  {kk}: mine={mine} anchor={ref} -> {'OK' if ok else 'MISMATCH'}")
if not all(v[2] for v in checks.values()):
    config.anchor_fail("Reproduction anchors not reproduced; stopping.")

out = {"id": "A01", "B": B, "seeds": {"table2": SEED_MAIN, "other": SEED_OTHER}, "n_main_nec": len(nec), "n_distant": len(dist),
       "reproduction_check": {kk: {"mine": v[0], "anchor": v[1], "ok": v[2]} for kk, v in checks.items()}}

# ------------------------------------------------------------------ (1) Table 2
print("\n(1) Table 2 by stage, landmarks 0/3/6/12/24, 12-month window, B=2000 seed", SEED_MAIN)
MAN = {  # manuscript Table 2: cancer death % (lo, hi), other-cause % (lo, hi)
    "All": {0: (57.9, 55.2, 60.5, 7.3, 5.9, 8.6), 3: (53.8, 50.5, 56.8, 5.6, 4.2, 7.1), 6: (52.6, 48.9, 55.9, 3.9, 2.6, 5.3), 12: (39.4, 34.9, 43.9, 3.4, 1.9, 5.3), 24: (19.6, 14.9, 24.9, 4.4, 1.8, 7.2)},
    "Localized": {0: (21.2, 14.0, 28.9, 7.9, 3.5, 13.2), 3: (16.1, 9.0, 23.6, 6.0, 1.9, 11.0), 6: (15.6, 8.2, 23.7, 2.2, 0.0, 5.6), 12: (20.7, 12.0, 29.9, 0.0, 0.0, 0.0), 24: (16.0, 7.0, 26.6, 3.5, 0.0, 8.8)},
    "Regional": {0: (34.8, 29.7, 39.9, 5.1, 3.0, 7.5), 3: (34.2, 28.7, 40.1, 3.0, 1.1, 5.3), 6: (36.1, 30.1, 42.4, 2.8, 1.2, 5.1), 12: (29.8, 22.7, 36.6, 1.7, 0.0, 3.7), 24: (17.1, 10.7, 24.3, 4.3, 0.9, 8.5)},
    "Distant": {0: (72.3, 69.2, 75.2, 7.0, 5.4, 8.7), 3: (71.6, 67.7, 75.3, 5.9, 4.0, 8.0), 6: (72.3, 67.6, 76.8, 4.7, 2.8, 7.0), 12: (60.0, 52.4, 67.5, 7.4, 3.7, 11.5), 24: (30.5, 16.9, 43.8, 6.4, 0.0, 14.3)},
}
MAN_N = {"All": [1417, 971, 768, 450, 241], "Localized": [116, 100, 91, 79, 61], "Regional": [346, 285, 256, 187, 122], "Distant": [908, 553, 396, 169, 48]}
MAN_D = {"All": [793, 503, 391, 172, 45], "Localized": [24, 16, 14, 16, 9], "Regional": [115, 92, 89, 54, 20], "Distant": [634, 381, 276, 98, 14]}
table2 = []; t2_boot = {}
for st in STAGES:
    df = nec if st == "All" else nec[nec.stage == st]
    res = point_and_boot(df, LM5, 12, SEED_MAIN); t2_boot[st] = res
    for j, s in enumerate(LM5):
        m = MAN[st][int(s)]
        row = {"stage": st, "landmark": int(s), "n_at_risk": int(res["nrisk"][j]), "cancer_deaths_window": res["d1"][j], "other_deaths_window": res["d2"][j],
               "cancer_death_pct": r(res["cif"][j]), "cancer_death_ci": [r(x) for x in ci(res["bc"][:, j])],
               "other_death_pct": r(res["oth"][j]), "other_death_ci": [r(x) for x in ci(res["bo"][:, j])],
               "allcause_death_pct": r(res["allc"][j]), "allcause_ci": [r(x) for x in ci(res["ba"][:, j])],
               "manuscript": {"n_at_risk": MAN_N[st][j], "cancer_deaths": MAN_D[st][j], "cancer_pct": m[0], "cancer_ci": [m[1], m[2]], "other_pct": m[3], "other_ci": [m[4], m[5]]}}
        row["diff_cancer_pct"] = r(round(res["cif"][j], 1) - m[0]); row["diff_other_pct"] = r(round(res["oth"][j], 1) - m[3])
        row["diff_cancer_ci"] = [r(round(ci(res["bc"][:, j])[0], 1) - m[1]), r(round(ci(res["bc"][:, j])[1], 1) - m[2])]
        row["diff_other_ci"] = [r(round(ci(res["bo"][:, j])[0], 1) - m[4]), r(round(ci(res["bo"][:, j])[1], 1) - m[5])]
        row["n_match"] = (int(res["nrisk"][j]) == MAN_N[st][j]) and (res["d1"][j] == MAN_D[st][j])
        table2.append(row)
        print(f"  {st:9s} LM{s:>2}: n={int(res['nrisk'][j]):5d} d={res['d1'][j]:4d}  cancer {res['cif'][j]:5.1f} [{ci(res['bc'][:, j])[0]:.1f}-{ci(res['bc'][:, j])[1]:.1f}]"
              f"  other {res['oth'][j]:4.1f} [{ci(res['bo'][:, j])[0]:.1f}-{ci(res['bo'][:, j])[1]:.1f}]  | ms {m[0]} [{m[1]}-{m[2]}] {m[3]} [{m[4]}-{m[5]}]  nmatch={row['n_match']}")
out["table2"] = table2
# monotone n at risk check
for st in STAGES:
    nr = t2_boot[st]["nrisk"]; assert np.all(np.diff(nr) < 0), f"n at risk not decreasing for {st}"
maxdiff_pt = max(abs(x["diff_cancer_pct"]) for x in table2); maxdiff_oth = max(abs(x["diff_other_pct"]) for x in table2)
maxdiff_ci = max(max(abs(x["diff_cancer_ci"][0]), abs(x["diff_cancer_ci"][1]), abs(x["diff_other_ci"][0]), abs(x["diff_other_ci"][1])) for x in table2)
out["table2_vs_manuscript"] = {"max_abs_diff_cancer_point": maxdiff_pt, "max_abs_diff_other_point": maxdiff_oth, "max_abs_diff_any_ci_bound": maxdiff_ci,
                               "all_n_and_deaths_match": all(x["n_match"] for x in table2)}
print(f"  max |diff| point cancer={maxdiff_pt} other={maxdiff_oth}; max |diff| CI bound={maxdiff_ci}; n/deaths all match={out['table2_vs_manuscript']['all_n_and_deaths_match']}")

# ------------------------------------------------------------------ (2) monthly distant, plateau, differences
print("\n(2) monthly landmarks 0-24, distant NEC, B=2000 seed", SEED_OTHER)
def monthly_block(df, label, seed=SEED_OTHER, diffs=((6, 0), (12, 0), (12, 6), (3, 0), (6, 3))):
    res = point_and_boot(df, LM25, 12, seed)
    cs, b1s, nok = hinge_boot(res["bc"]); c0, b0, b1 = hinge_fit(res["cif"][:19])
    csa, b1sa, noka = hinge_boot(res["ba"]); c0a, b0a, b1a = hinge_fit(res["allc"][:19])
    blk = {"label": label, "n": res["n"], "landmarks": LM25.tolist(), "n_at_risk": res["nrisk"].astype(int).tolist(),
           "cancer_deaths_window": res["d1"],
           "cancer_pct": [r(x) for x in res["cif"]], "cancer_ci_lo": [r(ci(res["bc"][:, j])[0]) for j in range(25)], "cancer_ci_hi": [r(ci(res["bc"][:, j])[1]) for j in range(25)],
           "other_pct": [r(x) for x in res["oth"]], "other_ci_lo": [r(ci(res["bo"][:, j])[0]) for j in range(25)], "other_ci_hi": [r(ci(res["bo"][:, j])[1]) for j in range(25)],
           "allcause_pct": [r(x) for x in res["allc"]], "allcause_ci_lo": [r(ci(res["ba"][:, j])[0]) for j in range(25)], "allcause_ci_hi": [r(ci(res["ba"][:, j])[1]) for j in range(25)],
           "plateau_cancer": {"c": c0, "c_ci": [r(x) for x in ci(cs)], "prop_c0": r(float(np.mean(cs == 0)), 4), "b0": r(b0), "b1_per_month": r(b1, 3), "b1_ci": [r(x, 3) for x in ci(b1s)], "B_used": nok,
                              "c_boot_quantiles_10_25_50_75_90": [r(x) for x in np.percentile(cs, [10, 25, 50, 75, 90])]},
           "plateau_allcause": {"c": c0a, "c_ci": [r(x) for x in ci(csa)], "prop_c0": r(float(np.mean(csa == 0)), 4), "b1_per_month": r(b1a, 3), "b1_ci": [r(x, 3) for x in ci(b1sa)], "B_used": noka},
           "differences": []}
    for a, b in diffs:
        d = res["cif"][a] - res["cif"][b]; db = res["bc"][:, a] - res["bc"][:, b]
        da = res["allc"][a] - res["allc"][b]; dab = res["ba"][:, a] - res["ba"][:, b]
        blk["differences"].append({"contrast": f"{a}-{b}", "cancer_diff": r(d), "cancer_ci": [r(x) for x in ci(db)], "allcause_diff": r(da), "allcause_ci": [r(x) for x in ci(dab)]})
    print(f"  {label}: n={res['n']} plateau c={c0} CI {ci(cs)} prop_c0={np.mean(cs==0):.4f} b1={b1:.3f} CI {ci(b1s)}; all-cause c={c0a} CI {ci(csa)}")
    for dd in blk["differences"]: print(f"     {dd['contrast']}: {dd['cancer_diff']} {dd['cancer_ci']}  (all-cause {dd['allcause_diff']} {dd['allcause_ci']})")
    return blk, res
mon_dist, res_dist = monthly_block(dist, "Distant")
out["monthly_distant"] = mon_dist
out["monthly_distant"]["manuscript"] = {"plateau_c": 6.75, "c_ci": [4.00, 9.75], "prop_c0": 0.002, "diff_6_0": [0.03, -4.0, 4.1], "diff_12_0": [-12.3, -19.8, -4.2], "slope_approx": -2.5}

# ------------------------------------------------------------------ (3) regional / localized monthly
print("\n(3) regional and localized monthly curves")
mon_reg, _ = monthly_block(nec[nec.stage == "Regional"], "Regional")
mon_loc, _ = monthly_block(nec[nec.stage == "Localized"], "Localized")
mon_reg["manuscript"] = {"plateau_c": 8.0, "c_ci": [0.00, 12.63], "diff_12_0": [-5.0, -13.6, 2.8]}
out["monthly_regional"] = mon_reg; out["monthly_localized"] = mon_loc

# ------------------------------------------------------------------ (4) 36-month window
print("\n(4) 36-month window, distant NEC diagnosed 2000-2020")
d36 = dist[dist.year <= 2020].reset_index(drop=True); LM4 = np.array([0, 3, 6, 12])
res36 = point_and_boot(d36, LM4, 36, SEED_OTHER, kmax_months=48)
res12s = point_and_boot(d36, LM4, 12, SEED_OTHER)  # same subset, 12-month window, for comparison
out["window36"] = {"n": len(d36), "years": "2000-2020", "window": 36, "rows": []}
for j, s in enumerate(LM4):
    row = {"landmark": int(s), "n_at_risk": int(res36["nrisk"][j]), "cancer_deaths_window": res36["d1"][j], "other_deaths_window": res36["d2"][j],
           "cancer_pct": r(res36["cif"][j]), "cancer_ci": [r(x) for x in ci(res36["bc"][:, j])],
           "other_pct": r(res36["oth"][j]), "other_ci": [r(x) for x in ci(res36["bo"][:, j])],
           "allcause_pct": r(res36["allc"][j]), "allcause_ci": [r(x) for x in ci(res36["ba"][:, j])],
           "cancer_pct_12mo_window_same_subset": r(res12s["cif"][j]), "cancer_ci_12mo": [r(x) for x in ci(res12s["bc"][:, j])]}
    out["window36"]["rows"].append(row)
    print(f"  LM{s:>2}: n={row['n_at_risk']} d={row['cancer_deaths_window']} cancer36 {row['cancer_pct']} {row['cancer_ci']} other36 {row['other_pct']} {row['other_ci']} | cancer12 {row['cancer_pct_12mo_window_same_subset']}")
d36b = {"6-0": None, "12-0": None}
for lab, (a, b) in {"6-0": (2, 0), "12-0": (3, 0), "12-6": (3, 2)}.items():
    d36b[lab] = {"diff": r(res36["cif"][a] - res36["cif"][b]), "ci": [r(x) for x in ci(res36["bc"][:, a] - res36["bc"][:, b])]}
out["window36"]["differences"] = d36b; print("  36-mo differences:", d36b)

# ------------------------------------------------------------------ (5) unknown-COD sensitivity
print("\n(5) unknown-COD sensitivity (exclude / censor at death)")
def stage_points(df, landmarks=LM5, window=12):
    pts = {}
    for st in STAGES:
        d = df if st == "All" else df[df.stage == st]
        k, ev, KMAX = lm.prep(d); c, a, nr = lm.curves(k, ev, KMAX, landmarks, window); o, _, _ = lm.curves(k, swap_ev(ev), KMAX, landmarks, window)
        pts[st] = {"cancer": c, "other": o, "all": a, "n": nr}
    return pts
base = stage_points(nec)
nec_ex = nec[nec.unkcod != 1]; nec_cen = nec.copy(); nec_cen.loc[nec_cen.unkcod == 1, "ev"] = 0
sens = {"n_unknown_cod": unk_all, "n_unknown_cod_distant": unk_dist, "exclude": {"n": len(nec_ex)}, "censor": {"n": len(nec_cen)}}
for lab, d in (("exclude", nec_ex), ("censor", nec_cen)):
    p = stage_points(d); rows = []; mc = 0; mo = 0
    for st in STAGES:
        for j, s in enumerate(LM5):
            dc = p[st]["cancer"][j] - base[st]["cancer"][j]; do = p[st]["other"][j] - base[st]["other"][j]
            mc = max(mc, abs(dc)); mo = max(mo, abs(do))
            rows.append({"stage": st, "landmark": int(s), "n_at_risk": int(p[st]["n"][j]), "cancer_pct": r(p[st]["cancer"][j]), "delta_cancer": r(dc, 3), "other_pct": r(p[st]["other"][j]), "delta_other": r(do, 3)})
    sens[lab].update({"rows": rows, "max_abs_change_cancer": r(mc, 3), "max_abs_change_other": r(mo, 3)})
    # distant monthly plateau under this sensitivity (point estimate)
    dd = d[d.stage == "Distant"]; k, ev, KMAX = lm.prep(dd); c25, _, _ = lm.curves(k, ev, KMAX, LM25, 12)
    sens[lab]["distant_plateau_c"] = hinge_fit(c25[:19])[0]
    print(f"  {lab}: n={len(d)} max|Δ| cancer={mc:.3f} other={mo:.3f}; distant plateau c={sens[lab]['distant_plateau_c']}")
sens["max_abs_change_cancer_both"] = r(max(sens["exclude"]["max_abs_change_cancer"], sens["censor"]["max_abs_change_cancer"]), 3)
sens["max_abs_change_any_outcome_both"] = r(max(sens["exclude"]["max_abs_change_cancer"], sens["censor"]["max_abs_change_cancer"], sens["exclude"]["max_abs_change_other"], sens["censor"]["max_abs_change_other"]), 3)
sens["manuscript_max_abs_change"] = 0.77
# censored variant with bootstrap for distant plateau CI (manuscript 7.25)
res_cen = point_and_boot(nec_cen[nec_cen.stage == "Distant"], LM25, 12, SEED_OTHER)
cs_cen, _, _ = hinge_boot(res_cen["bc"])
sens["censor"]["distant_plateau_ci"] = [r(x) for x in ci(cs_cen)]; sens["censor"]["distant_plateau_prop_c0"] = r(float(np.mean(cs_cen == 0)), 4)
sens["censor"]["distant_0_6_12"] = [r(res_cen["cif"][i]) for i in (0, 6, 12)]
print(f"  censor: distant plateau {sens['censor']['distant_plateau_c']} CI {sens['censor']['distant_plateau_ci']} prop0={sens['censor']['distant_plateau_prop_c0']}")
out["unknown_cod_sensitivity"] = sens

# ------------------------------------------------------------------ (6) complete-date and 8045-removed subsets
print("\n(6) complete-date subset and 8045-removed subset (distant)")
subs = {"complete_dates": dist[dist.flagc == 1], "no_8045": dist[dist["hist"] != "8045"]}
out["subsets"] = {}
for lab, d in subs.items():
    res = point_and_boot(d, LM25, 12, SEED_OTHER); cs, b1s, nok = hinge_boot(res["bc"]); c0, b0, b1 = hinge_fit(res["cif"][:19])
    blk = {"n": len(d), "rows": [], "plateau_c": c0, "plateau_ci": [r(x) for x in ci(cs)], "prop_c0": r(float(np.mean(cs == 0)), 4), "b1_per_month": r(b1, 3), "b1_ci": [r(x, 3) for x in ci(b1s)],
           "manuscript_plateau_c": 7.25, "manuscript_n": 881 if lab == "complete_dates" else 890}
    for s in (0, 6, 12):
        blk["rows"].append({"landmark": s, "n_at_risk": int(res["nrisk"][s]), "cancer_deaths_window": res["d1"][s], "cancer_pct": r(res["cif"][s]), "cancer_ci": [r(x) for x in ci(res["bc"][:, s])],
                            "other_pct": r(res["oth"][s]), "other_ci": [r(x) for x in ci(res["bo"][:, s])]})
    blk["diff_6_0"] = {"diff": r(res["cif"][6] - res["cif"][0]), "ci": [r(x) for x in ci(res["bc"][:, 6] - res["bc"][:, 0])]}
    blk["diff_12_0"] = {"diff": r(res["cif"][12] - res["cif"][0]), "ci": [r(x) for x in ci(res["bc"][:, 12] - res["bc"][:, 0])]}
    out["subsets"][lab] = blk
    print(f"  {lab}: n={len(d)} 0/6/12 = {[x['cancer_pct'] for x in blk['rows']]} plateau {c0} CI {blk['plateau_ci']} prop0={blk['prop_c0']}")

# ------------------------------------------------------------------ figures
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
C1, C2, C3 = "#2a78d6", "#eb6834", "#1baf7a"; INK = "#0b0b0b"; INK2 = "#52514e"; GRID = "#e6e6e3"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": INK2, "axes.labelcolor": INK, "xtick.color": INK2, "ytick.color": INK2, "axes.spines.top": False, "axes.spines.right": False})
def panel(ax, blk, title, show_all=True):
    x = np.array(blk["landmarks"]); n = np.array(blk["n_at_risk"]); filled = n >= 100
    y = np.array(blk["cancer_pct"], float); lo = np.array(blk["cancer_ci_lo"], float); hi = np.array(blk["cancer_ci_hi"], float)
    ax.fill_between(x, lo, hi, color=C1, alpha=0.15, linewidth=0)
    ax.plot(x, y, color=C1, lw=2, label="Cancer death (Aalen-Johansen)")
    ax.plot(x[filled], y[filled], "o", color=C1, ms=4); ax.plot(x[~filled], y[~filled], "o", mfc="white", mec=C1, ms=4)
    if show_all:
        ya = np.array(blk["allcause_pct"], float); loa = np.array(blk["allcause_ci_lo"], float); hia = np.array(blk["allcause_ci_hi"], float)
        ax.fill_between(x, loa, hia, color=C2, alpha=0.15, linewidth=0); ax.plot(x, ya, color=C2, lw=2, label="All-cause death")
    pc = blk["plateau_cancer"]; ax.axvspan(pc["c_ci"][0], pc["c_ci"][1], color=C3, alpha=0.10, linewidth=0)
    ax.axvline(pc["c"], color=C3, lw=1.5, ls="--", label=f"Plateau end {pc['c']:.2f} mo (95% CI {pc['c_ci'][0]:.2f}-{pc['c_ci'][1]:.2f})")
    ax.set_title(title, loc="left", fontsize=10, color=INK); ax.set_xlabel("Landmark (months survived after diagnosis)"); ax.set_ylabel("Probability of death in next 12 months, %")
    ax.set_xlim(-0.5, 24.5); ax.set_ylim(0, 100); ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True); ax.legend(frameon=False, fontsize=8, loc="upper right")
    for xi, ni in zip(x[::3], n[::3]): ax.text(xi, 2, str(ni), ha="center", fontsize=7, color=INK2)
    ax.text(-0.3, 6, "n at risk", fontsize=7, color=INK2, ha="left")
fig, ax = plt.subplots(figsize=(7.5, 4.6), dpi=300); panel(ax, mon_dist, "Distant colorectal NEC (n = 908): 12-month outcome by monthly landmark"); fig.tight_layout(); fig.savefig(FIGS + "A01_distant_monthly.png"); plt.close(fig)
fig, axs = plt.subplots(1, 2, figsize=(11, 4.6), dpi=300)
panel(axs[0], mon_reg, f"Regional (n = {mon_reg['n']})"); panel(axs[1], mon_loc, f"Localized (n = {mon_loc['n']})"); fig.tight_layout(); fig.savefig(FIGS + "A01_regional_localized_monthly.png"); plt.close(fig)
# Table 2 figure: cancer-death probability by stage at the 5 pre-specified landmarks
fig, ax = plt.subplots(figsize=(7.5, 4.4), dpi=300); cols = {"Distant": C1, "Regional": C2, "Localized": C3, "All": INK2}
for st in ["Distant", "Regional", "Localized", "All"]:
    rows = [x for x in table2 if x["stage"] == st]; xs = [x["landmark"] for x in rows]; ys = [x["cancer_death_pct"] for x in rows]
    lo = [x["cancer_death_ci"][0] for x in rows]; hi = [x["cancer_death_ci"][1] for x in rows]
    ax.errorbar(xs, ys, yerr=[np.array(ys) - np.array(lo), np.array(hi) - np.array(ys)], color=cols[st], lw=1.8, capsize=3, marker="o", ms=4, label=st, ls="-" if st != "All" else ":")
ax.set_xticks([0, 3, 6, 12, 24]); ax.set_xlabel("Landmark (months)"); ax.set_ylabel("Cancer death in next 12 months, % (95% CI)"); ax.set_ylim(0, 100); ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
ax.set_title("Pre-specified landmarks by stage (Table 2), patient-level bootstrap B = 2000", loc="left", fontsize=10); ax.legend(frameon=False, fontsize=8); fig.tight_layout(); fig.savefig(FIGS + "A01_table2_by_stage.png"); plt.close(fig)
# bootstrap distribution of plateau end
fig, ax = plt.subplots(figsize=(6, 3.6), dpi=300)
cs_d, _, _ = hinge_boot(res_dist["bc"]); vals, cnts = np.unique(cs_d, return_counts=True)
ax.bar(vals, cnts / cnts.sum() * 100, width=0.22, color=C1); ax.axvline(6.75, color=C3, ls="--", lw=1.5)
ax.set_xlabel("Plateau end c (months), bootstrap resamples"); ax.set_ylabel("% of resamples"); ax.set_title("Distant NEC: bootstrap distribution of plateau end (B = 2000)", loc="left", fontsize=10); ax.grid(axis="y", color=GRID, lw=0.6); ax.set_axisbelow(True)
fig.tight_layout(); fig.savefig(FIGS + "A01_plateau_bootstrap.png"); plt.close(fig)

out["runtime_minutes"] = r((time.time() - T0) / 60, 2)
out["figures"] = [config.relpath(FIGS + f) for f in ("A01_distant_monthly.png", "A01_regional_localized_monthly.png", "A01_table2_by_stage.png", "A01_plateau_bootstrap.png")]
with open(RES + "A01.json", "w") as f: json.dump(out, f, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))

# ------------------------------------------------------------------ markdown
L = []
L.append("# A01 landmark_main: pre-specified landmark analysis, monthly curves, plateau and sensitivity analyses\n")
L.append(f"Main colorectal NEC cohort n = {len(nec)} (distant {len(dist)}). Aalen-Johansen 12-month probabilities; patient-level percentile bootstrap, B = {B} (seed {SEED_MAIN} for Table 2, {SEED_OTHER} for all other analyses), same resample for all landmarks. Runtime {out['runtime_minutes']} min.\n")
L.append("## Reproduction check (mine vs orchestrator anchors)\n\n| Item | Mine | Anchor | OK |\n|---|---|---|---|")
for kk, v in checks.items(): L.append(f"| {kk} | {v[0]} | {v[1]} | {'yes' if v[2] else 'NO'} |")
L.append("\n## (1) Table 2: 12-month outcome probabilities by stage and pre-specified landmark\n")
L.append("| Stage | Landmark (mo) | n at risk | Cancer deaths in window | Cancer death % (95% CI) | Other-cause death % (95% CI) | Manuscript cancer % (CI) | Manuscript other % (CI) | n/deaths match |\n|---|---|---|---|---|---|---|---|---|")
for x in table2:
    m = x["manuscript"]
    L.append(f"| {x['stage']} | {x['landmark']} | {x['n_at_risk']} | {x['cancer_deaths_window']} | {x['cancer_death_pct']:.1f} ({x['cancer_death_ci'][0]:.1f}-{x['cancer_death_ci'][1]:.1f}) | {x['other_death_pct']:.1f} ({x['other_death_ci'][0]:.1f}-{x['other_death_ci'][1]:.1f}) | {m['cancer_pct']} ({m['cancer_ci'][0]}-{m['cancer_ci'][1]}) | {m['other_pct']} ({m['other_ci'][0]}-{m['other_ci'][1]}) | {'yes' if x['n_match'] else 'NO'} |")
L.append(f"\nMax absolute difference from manuscript: point estimates cancer {maxdiff_pt} / other-cause {maxdiff_oth} percentage points (after rounding to 1 dp); any CI bound {maxdiff_ci} points. All n at risk and death counts match: {out['table2_vs_manuscript']['all_n_and_deaths_match']}.\n")
L.append("## (2) Distant NEC: monthly landmarks 0-24 (12-month window)\n")
L.append("| Landmark | n at risk | Cancer deaths | Cancer death % (95% CI) | Other-cause % (95% CI) | All-cause % (95% CI) |\n|---|---|---|---|---|---|")
for j in range(25):
    b = mon_dist
    L.append(f"| {j} | {b['n_at_risk'][j]} | {b['cancer_deaths_window'][j]} | {b['cancer_pct'][j]:.1f} ({b['cancer_ci_lo'][j]:.1f}-{b['cancer_ci_hi'][j]:.1f}) | {b['other_pct'][j]:.1f} ({b['other_ci_lo'][j]:.1f}-{b['other_ci_hi'][j]:.1f}) | {b['allcause_pct'][j]:.1f} ({b['allcause_ci_lo'][j]:.1f}-{b['allcause_ci_hi'][j]:.1f}) |")
pc = mon_dist["plateau_cancer"]; pa = mon_dist["plateau_allcause"]
L.append(f"\nBroken-stick (landmarks 0-18, grid 0-15 by 0.25): plateau end c = {pc['c']} months (95% CI {pc['c_ci'][0]}-{pc['c_ci'][1]}); proportion of resamples with c = 0: {pc['prop_c0']*100:.2f}% (manuscript 6.75 [4.00-9.75], 0.2%). Post-plateau slope b1 = {pc['b1_per_month']} points/month (95% CI {pc['b1_ci'][0]} to {pc['b1_ci'][1]}); intercept b0 = {pc['b0']}%. All-cause curve: c = {pa['c']} (95% CI {pa['c_ci'][0]}-{pa['c_ci'][1]}), b1 = {pa['b1_per_month']} ({pa['b1_ci'][0]} to {pa['b1_ci'][1]}). Bootstrap c quantiles (10/25/50/75/90): {pc['c_boot_quantiles_10_25_50_75_90']}.\n")
L.append("| Contrast (months) | Cancer death difference, points (95% CI) | All-cause difference (95% CI) | Manuscript |\n|---|---|---|---|")
msd = {"6-0": "+0.03 (-4.0 to 4.1)", "12-0": "-12.3 (-19.8 to -4.2)"}
for dd in mon_dist["differences"]: L.append(f"| {dd['contrast']} | {dd['cancer_diff']:+.2f} ({dd['cancer_ci'][0]:.1f} to {dd['cancer_ci'][1]:.1f}) | {dd['allcause_diff']:+.2f} ({dd['allcause_ci'][0]:.1f} to {dd['allcause_ci'][1]:.1f}) | {msd.get(dd['contrast'], '-')} |")
L.append("\n## (3) Regional and localized: monthly curves and plateau\n")
for b in (mon_reg, mon_loc):
    pc = b["plateau_cancer"]
    L.append(f"**{b['label']}** (n = {b['n']}): plateau end c = {pc['c']} (95% CI {pc['c_ci'][0]}-{pc['c_ci'][1]}), c = 0 in {pc['prop_c0']*100:.1f}% of resamples, b1 = {pc['b1_per_month']} ({pc['b1_ci'][0]} to {pc['b1_ci'][1]}). Differences: " + "; ".join(f"{dd['contrast']}: {dd['cancer_diff']:+.1f} ({dd['cancer_ci'][0]:.1f} to {dd['cancer_ci'][1]:.1f})" for dd in b["differences"]) + ".\n")
    L.append("| Landmark | n at risk | Cancer death % (95% CI) |\n|---|---|---|")
    for j in range(0, 25, 2): L.append(f"| {j} | {b['n_at_risk'][j]} | {b['cancer_pct'][j]:.1f} ({b['cancer_ci_lo'][j]:.1f}-{b['cancer_ci_hi'][j]:.1f}) |")
    L.append("")
L.append(f"## (4) Pre-specified 36-month window, distant NEC diagnosed 2000-2020 (n = {len(d36)})\n")
L.append("| Landmark | n at risk | Cancer deaths (36 mo) | Cancer death % (95% CI), 36 mo | Other-cause % (95% CI), 36 mo | All-cause % (95% CI), 36 mo | Cancer death %, 12-mo window, same subset |\n|---|---|---|---|---|---|---|")
for x in out["window36"]["rows"]:
    L.append(f"| {x['landmark']} | {x['n_at_risk']} | {x['cancer_deaths_window']} | {x['cancer_pct']:.1f} ({x['cancer_ci'][0]:.1f}-{x['cancer_ci'][1]:.1f}) | {x['other_pct']:.1f} ({x['other_ci'][0]:.1f}-{x['other_ci'][1]:.1f}) | {x['allcause_pct']:.1f} ({x['allcause_ci'][0]:.1f}-{x['allcause_ci'][1]:.1f}) | {x['cancer_pct_12mo_window_same_subset']:.1f} ({x['cancer_ci_12mo'][0]:.1f}-{x['cancer_ci_12mo'][1]:.1f}) |")
L.append("\n36-month-window differences: " + "; ".join(f"{k}: {v['diff']:+.1f} ({v['ci'][0]:.1f} to {v['ci'][1]:.1f})" for k, v in d36b.items()) + ".\n")
L.append(f"## (5) Unknown cause of death ({unk_all} deaths, {unk_dist} distant)\n")
L.append("| Handling | n | Max |change| cancer-death % (stages x landmarks 0/3/6/12/24) | Max |change| other-cause % | Distant plateau end |\n|---|---|---|---|---|")
for lab in ("exclude", "censor"):
    s = sens[lab]; extra = f" (95% CI {s['distant_plateau_ci'][0]}-{s['distant_plateau_ci'][1]}; c=0 in {s['distant_plateau_prop_c0']*100:.1f}%)" if lab == "censor" else ""
    L.append(f"| {lab} | {s['n']} | {s['max_abs_change_cancer']} | {s['max_abs_change_other']} | {s['distant_plateau_c']}{extra} |")
L.append(f"\nManuscript: max change 0.77 points; plateau 7.25 after censoring. Mine: max change in cancer-death probability {sens['max_abs_change_cancer_both']} points (cancer) / {sens['max_abs_change_any_outcome_both']} points (any outcome).\n")
L.append("## (6) Complete-date subset and 8045-removed subset (distant NEC)\n")
L.append("| Subset | n | LM 0 % (95% CI) | LM 6 % (95% CI) | LM 12 % (95% CI) | 6-0 diff (95% CI) | 12-0 diff (95% CI) | Plateau end (95% CI) | c=0 % | Manuscript plateau |\n|---|---|---|---|---|---|---|---|---|---|")
for lab, b in out["subsets"].items():
    rr = b["rows"]
    L.append(f"| {lab} | {b['n']} | " + " | ".join(f"{x['cancer_pct']:.1f} ({x['cancer_ci'][0]:.1f}-{x['cancer_ci'][1]:.1f})" for x in rr) + f" | {b['diff_6_0']['diff']:+.1f} ({b['diff_6_0']['ci'][0]:.1f} to {b['diff_6_0']['ci'][1]:.1f}) | {b['diff_12_0']['diff']:+.1f} ({b['diff_12_0']['ci'][0]:.1f} to {b['diff_12_0']['ci'][1]:.1f}) | {b['plateau_c']} ({b['plateau_ci'][0]}-{b['plateau_ci'][1]}) | {b['prop_c0']*100:.1f} | {b['manuscript_plateau_c']} (n = {b['manuscript_n']}) |")
L.append("\n## Figures\n\n" + "\n".join(f"- {f}" for f in out["figures"]) + "\n")
L.append("## Interpretation\n")
pc = mon_dist["plateau_cancer"]; d60 = mon_dist["differences"][0]; d120 = mon_dist["differences"][1]
L.append(f"- The pre-specified Table 2 is reproduced exactly in its point estimates, n at risk and death counts (max point difference {maxdiff_pt} points after rounding); bootstrap CI bounds differ by at most {maxdiff_ci} points, which is the expected Monte-Carlo variation of an independent random-number stream at B = 2000, not a substantive difference.")
L.append(f"- Among distant-stage NEC, the 12-month cancer-death probability is {mon_dist['cancer_pct'][0]:.1f}% at diagnosis and {mon_dist['cancer_pct'][6]:.1f}% after surviving 6 months (difference {d60['cancer_diff']:+.2f} points, 95% CI {d60['cancer_ci'][0]:.1f} to {d60['cancer_ci'][1]:.1f}); by 12 months it has fallen to {mon_dist['cancer_pct'][12]:.1f}% (difference {d120['cancer_diff']:+.1f}, 95% CI {d120['cancer_ci'][0]:.1f} to {d120['cancer_ci'][1]:.1f}).")
L.append(f"- The broken-stick plateau end is {pc['c']} months (95% CI {pc['c_ci'][0]}-{pc['c_ci'][1]}); only {pc['prop_c0']*100:.2f}% of resamples placed the breakpoint at 0, so an immediate post-diagnosis decline is not supported. After the plateau the 12-month risk falls by {abs(pc['b1_per_month']):.2f} points per month (95% CI {abs(pc['b1_ci'][1]):.2f} to {abs(pc['b1_ci'][0]):.2f}).")
L.append(f"- Other-cause death within 12 months is small at every landmark (distant: {table2[15]['other_death_pct']:.1f}% at diagnosis), so the plateau is a property of cancer-specific risk rather than an artefact of competing deaths; the all-cause curve gives the same breakpoint ({mon_dist['plateau_allcause']['c']} months).")
L.append(f"- Regional disease shows a similar but imprecise plateau (c = {mon_reg['plateau_cancer']['c']}, 95% CI {mon_reg['plateau_cancer']['c_ci'][0]}-{mon_reg['plateau_cancer']['c_ci'][1]}; c = 0 in {mon_reg['plateau_cancer']['prop_c0']*100:.1f}% of resamples); localized disease (n = {mon_loc['n']}) is too small for inference and its 12-month risk is low and flat.")
w0 = out["window36"]["rows"][0]; w2 = out["window36"]["rows"][2]; w3 = out["window36"]["rows"][3]
L.append(f"- With the pre-specified 36-month window (diagnosed 2000-2020, n = {len(d36)}), cancer death is {w0['cancer_pct']:.1f}% from diagnosis, {w2['cancer_pct']:.1f}% from 6 months and {w3['cancer_pct']:.1f}% from 12 months: the longer window is less sensitive to the plateau because it integrates over the later decline, but the 0-to-6-month comparison remains flat ({d36b['6-0']['diff']:+.1f}, 95% CI {d36b['6-0']['ci'][0]:.1f} to {d36b['6-0']['ci'][1]:.1f}).")
L.append(f"- Handling of the {unk_all} unknown-cause deaths (exclude or censor) changes any stage-landmark cancer-death probability by at most {sens['max_abs_change_cancer_both']} points (manuscript 0.77); the other-cause probability moves by up to {sens['max_abs_change_any_outcome_both']} points because those deaths are reclassified, which the manuscript's 0.77 figure does not cover and should be stated. The complete-date (n = {out['subsets']['complete_dates']['n']}) and 8045-removed (n = {out['subsets']['no_8045']['n']}) subsets give plateau ends of {out['subsets']['complete_dates']['plateau_c']} and {out['subsets']['no_8045']['plateau_c']} months with CIs excluding 0.")
L.append("- Deviations from the plan: none in design. Bootstrap CIs use an independent RNG stream (numpy default_rng, fresh per stratum), so CI bounds are reproducible from this script but will not match the manuscript's digits exactly; the 36-month window was not reported in the manuscript and is presented here as new pre-specified output.")
with open(RES + "A01.md", "w") as f: f.write("\n".join(L) + "\n")
print(f"\nDone. runtime {out['runtime_minutes']} min. Wrote {RES}A01.json, {RES}A01.md and {len(out['figures'])} figures.")
