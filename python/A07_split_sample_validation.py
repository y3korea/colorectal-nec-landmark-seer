# -*- coding: utf-8 -*-
"""A07 split-sample (calendar) validation, distant main colorectal NEC.
Derivation 2000-2015 vs confirmation 2016-2023: landmark 12-month cancer-death probabilities (AJ), bootstrap CIs,
6-0 change, plateau end (broken stick), crude interval rates, adjusted RR 6-9 vs 3-6 (R, pp.R), three-era
stratification, SCLC/AC/PANC 6-0 change within 2016-2023 and the NEC-minus-comparator contrast.
Run: python3 python/A07_split_sample_validation.py (calls Rscript R/A07_split_sample_validation.R)."""
import os, sys, time, json, subprocess, numpy as np, pandas as pd
T0 = time.time()
import config
import lm
RES, FIG, SCR = config.RES, config.FIGS, config.R_DIR + os.sep
B_MAIN, B_CMP, SEED = config.B(2000), config.B(500), 20261007
LM4 = np.array([0, 3, 6, 12]); LM19 = np.arange(0, 19); LM25 = np.arange(0, 25)

# ---------------- reproduction check (full distant NEC) ----------------
nen = lm.load("NEN_all"); nec = lm.nec_main(nen); D = nec[nec.stage == "Distant"].copy()
k, ev, K = lm.prep(D); c5, _, n5 = lm.curves(k, ev, K, np.array([0, 3, 6, 12, 24]))
c19, _, _ = lm.curves(k, ev, K, LM19); plateau_all = float(lm.hinge_many(c19[None, :])[0])
rates_all = lm.crude_rates(D).rate.to_numpy()
anchor = {"n_distant": int(len(D)), "cif_0_3_6_12_24": [round(float(x), 2) for x in c5], "nrisk": [int(x) for x in n5],
          "plateau_end": plateau_all, "crude_rates": [round(float(x), 4) for x in rates_all]}
exp_c = [72.29, 71.62, 72.33, 60.03, 30.51]; exp_n = [908, 553, 396, 169, 48]; exp_r = [0.1435, 0.0909, 0.1344, 0.1204, 0.1231, 0.0840, 0.0646]
ok = (len(D) == 908 and np.allclose(np.round(c5, 2), exp_c, atol=0.011) and list(n5.astype(int)) == exp_n
      and plateau_all == 6.75 and np.allclose(np.round(rates_all, 4), exp_r, atol=1e-4)
      and len(D[D.year <= 2015]) == 431 and len(D[D.year >= 2016]) == 477)
print("REPRODUCTION CHECK:", anchor, "->", "PASS" if ok else "FAIL")
if not ok: config.anchor_fail("Reproduction anchors do not match; stopping.")

# ---------------- helpers ----------------
def boot_curves(df, B, seed, landmarks=LM25):
    """patient-level resampling; same resample for all landmarks; returns (point curve, n at risk, boot matrix)"""
    k, ev, K = lm.prep(df); n = len(k)
    pt, _, nr = lm.curves(k, ev, K, landmarks)
    rng = np.random.default_rng(seed); bc = np.empty((B, len(landmarks)))
    for b in range(B):
        w = np.bincount(rng.integers(0, n, n), minlength=n).astype(float)
        bc[b] = lm.curves(k, ev, K, landmarks, 12, w)[0]
    return pt, nr, bc
def ci(v): return lm.pct(v)
def boot_p(v):  # two-sided percentile-bootstrap p for H0: quantity = 0
    v = v[~np.isnan(v)]; return float(2 * min((v <= 0).mean(), (v >= 0).mean()))
idx = {int(s): int(np.where(LM25 == s)[0][0]) for s in [0, 3, 6, 12]}

def analyse(df, label, B, seed, full=True):
    pt, nr, bc = boot_curves(df, B, seed)
    out = {"label": label, "n": int(len(df)), "years": [int(df.year.min()), int(df.year.max())], "B": B, "seed": seed,
           "deaths_cancer": int((df.ev == 1).sum()), "deaths_other": int((df.ev == 2).sum()),
           "landmarks": {}, "curve_0_24": [None if np.isnan(x) else round(float(x), 2) for x in pt],
           "nrisk_0_24": [int(x) for x in nr]}
    for s in [0, 3, 6, 12]:
        j = idx[s]; out["landmarks"][str(s)] = {"cif12": round(float(pt[j]), 2), "ci": [round(x, 2) for x in ci(bc[:, j])], "nrisk": int(nr[j])}
    d60 = bc[:, idx[6]] - bc[:, idx[0]]
    out["change_6_minus_0"] = {"est": round(float(pt[idx[6]] - pt[idx[0]]), 2), "ci": [round(x, 2) for x in ci(d60)], "boot_p": boot_p(d60)}
    d120 = bc[:, idx[12]] - bc[:, idx[0]]
    out["change_12_minus_0"] = {"est": round(float(pt[idx[12]] - pt[idx[0]]), 2), "ci": [round(x, 2) for x in ci(d120)]}
    if full:
        pl = lm.hinge_many(bc[:, :19]); plp = float(lm.hinge_many(pt[None, :19])[0])
        out["plateau_end"] = {"est": plp, "ci": [round(x, 2) for x in ci(pl)], "boot_median": float(np.median(pl)),
                              "boot_share_ge_3": float((pl >= 3).mean()), "boot_share_ge_4_5": float((pl >= 4.5).mean())}
        cr = lm.crude_rates(df)
        out["crude_rates"] = [{"interval": f"{int(a)}-{int(b)}", "deaths": int(d), "person_months": round(float(pm), 1), "rate": round(float(r), 4)}
                              for a, b, d, pm, r in cr[["start", "end", "deaths", "person_months", "rate"]].itertuples(index=False)]
        r = cr.rate.to_numpy(); out["crude_rate_ratio_6to9_vs_3to6"] = round(float(r[2] / r[1]), 3)
    return out, bc, pt, nr

results = {"id": "A07", "task": "split_sample_validation", "reproduction_check": anchor, "reproduction_pass": bool(ok),
           "method": "Aalen-Johansen 12-month cancer-death probability conditional on being alive/under follow-up at each landmark; "
                     "patient-level bootstrap, same resample for all landmarks; 95% percentile CI; plateau end = broken-stick hinge on landmarks 0-18 "
                     "(grid 0-15 by 0.25); cross-period/cross-tumour contrasts use independent bootstraps paired by replicate index."}

# ---------------- derivation vs confirmation ----------------
early, late = D[D.year <= 2015], D[D.year >= 2016]
r_early, bc_e, pt_e, nr_e = analyse(early, "derivation 2000-2015", B_MAIN, SEED)
r_late, bc_l, pt_l, nr_l = analyse(late, "confirmation 2016-2023", B_MAIN, SEED + 1)
results["derivation_2000_2015"] = r_early; results["confirmation_2016_2023"] = r_late
dd = (bc_l[:, idx[6]] - bc_l[:, idx[0]]) - (bc_e[:, idx[6]] - bc_e[:, idx[0]])
results["period_difference_in_6_minus_0"] = {"est": round(r_late["change_6_minus_0"]["est"] - r_early["change_6_minus_0"]["est"], 2),
                                              "ci": [round(x, 2) for x in ci(dd)], "boot_p": boot_p(dd),
                                              "note": "late minus early; independent bootstraps paired by replicate index"}
dpl = lm.hinge_many(bc_l[:, :19]) - lm.hinge_many(bc_e[:, :19])
results["period_difference_in_plateau_end"] = {"est": r_late["plateau_end"]["est"] - r_early["plateau_end"]["est"], "ci": [round(x, 2) for x in ci(dpl)], "boot_p": boot_p(dpl)}
for s in [0, 3, 6, 12]:
    dl = bc_l[:, idx[s]] - bc_e[:, idx[s]]
    results[f"period_difference_cif_landmark_{s}"] = {"est": round(float(pt_l[idx[s]] - pt_e[idx[s]]), 2), "ci": [round(x, 2) for x in ci(dl)], "boot_p": boot_p(dl)}

# ---------------- three eras ----------------
results["eras"] = {}
for lab, lo, hi, sd in [("2000-2007", 2000, 2007, SEED + 10), ("2008-2015", 2008, 2015, SEED + 11), ("2016-2023", 2016, 2023, SEED + 12)]:
    sub = D[(D.year >= lo) & (D.year <= hi)]
    r, _, _, _ = analyse(sub, lab, B_MAIN, sd); results["eras"][lab] = r

# ---------------- comparators within 2016-2023 ----------------
comp = {}
ac = lm.load("AC"); sclc = lm.sclc_main(lm.load("LUNG")); panc = lm.panc_main(lm.load("PANC"))
bcl500 = bc_l[:B_CMP]
for name, df in [("SCLC", sclc), ("AC", ac), ("PANC", panc)]:
    df = df[df.stage == "Distant"]
    sub = df[df.year >= 2016]
    r, bc, pt, nr = analyse(sub, f"{name} distant 2016-2023", B_CMP, SEED + 20, full=False)
    pl = lm.hinge_many(pt[None, :19])[0]; r["plateau_end_point"] = float(pl)
    # NEC (2016-2023) minus comparator (2016-2023) contrast in the 6-0 change, paired by replicate index (first 500 NEC replicates)
    dc = (bcl500[:, idx[6]] - bcl500[:, idx[0]]) - (bc[:, idx[6]] - bc[:, idx[0]])
    r["NEC_minus_comparator_change_6_minus_0"] = {"est": round(r_late["change_6_minus_0"]["est"] - r["change_6_minus_0"]["est"], 2),
                                                 "ci": [round(x, 2) for x in ci(dc)], "boot_p": boot_p(dc), "B_pairs": B_CMP}
    # point estimates in 2000-2015 for context (no bootstrap)
    e = df[df.year <= 2015]; ke, eve, Ke = lm.prep(e); pe, _, ne = lm.curves(ke, eve, Ke, LM25)
    r["early_2000_2015_point"] = {"n": int(len(e)), "cif12_0": round(float(pe[0]), 2), "cif12_6": round(float(pe[6]), 2),
                                  "change_6_minus_0": round(float(pe[6] - pe[0]), 2), "plateau_end_point": float(lm.hinge_many(pe[None, :19])[0])}
    comp[name] = r
results["comparators_2016_2023"] = comp
# NEC early-vs-comparator contrast for completeness (point only)
results["cross_tumour_summary"] = {nm: {"NEC_minus_comp_2016_2023": comp[nm]["NEC_minus_comparator_change_6_minus_0"]["est"],
                                        "NEC_minus_comp_2000_2015_point": round(r_early["change_6_minus_0"]["est"] - comp[nm]["early_2000_2015_point"]["change_6_minus_0"], 2)} for nm in comp}

# ---------------- adjusted RR (R) ----------------
rp = subprocess.run(["Rscript", os.path.join(SCR, "A07_split_sample_validation.R")], capture_output=True, text=True)
print(rp.stdout[-2000:]); 
if rp.returncode != 0: sys.exit("R step failed:\n" + rp.stderr[-3000:])
with open(RES + "A07_rr.json") as f: rr = json.load(f)
results["adjusted_RR_R"] = rr
rra = rr["anchor"]["RR_6to9_vs_3to6"]
if not abs(rra["RR"] - 1.49) < 0.01: config.anchor_fail("R anchor mismatch")

# ---------------- figure ----------------
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
fig, axes = plt.subplots(1, 2, figsize=(9, 4.0), sharey=True)
cols = {"early": "#1f5fa8", "late": "#d4731c"}
for ax, (r, bc, pt, nr, key, ttl) in zip(axes, [(r_early, bc_e, pt_e, nr_e, "early", "A. Derivation 2000-2015 (n=431)"),
                                              (r_late, bc_l, pt_l, nr_l, "late", "B. Confirmation 2016-2023 (n=477)")]):
    lo = np.nanpercentile(bc, 2.5, axis=0); hi = np.nanpercentile(bc, 97.5, axis=0)
    m = nr >= 20
    ax.fill_between(LM25[m], lo[m], hi[m], color=cols[key], alpha=0.15, lw=0)
    ax.plot(LM25[m], pt[m], color=cols[key], lw=2, marker="o", ms=3.5)
    pe = r["plateau_end"]["est"]; ax.axvline(pe, color="#777", lw=1, ls="--")
    ax.text(pe + 0.3, 8, f"plateau end {pe:g} m\n(95% CI {r['plateau_end']['ci'][0]:g}-{r['plateau_end']['ci'][1]:g})", color="#444", fontsize=8)
    ax.set_title(ttl, loc="left", fontsize=10); ax.set_xlabel("Landmark (months since diagnosis)")
    ax.set_xticks(range(0, 25, 3)); ax.set_ylim(0, 100); ax.grid(axis="y", color="#e5e5e5", lw=0.6); ax.set_axisbelow(True)
    for s in [0, 6, 12]:
        j = idx[s]; ax.annotate(f"{pt[j]:.1f}%", (s, pt[j]), textcoords="offset points", xytext=(0, 7), ha="center", fontsize=8, color="#222")
    ax.text(0.02, 0.03, "n at risk: " + " ".join(f"{int(nr[j])}" for j in [0, 3, 6, 9, 12, 15, 18, 24]) + "  (m 0,3,...,24)", transform=ax.transAxes, fontsize=7, color="#555")
axes[0].set_ylabel("12-month cancer-death probability (%)")
fig.suptitle("Distant colorectal NEC: landmark 12-month cancer-death probability by calendar period", x=0.01, ha="left", fontsize=10.5)
fig.tight_layout(); fig.savefig(FIG + "A07_derivation_vs_confirmation.png", dpi=300); plt.close(fig)

# ---------------- outputs ----------------
results["runtime_minutes"] = round((time.time() - T0) / 60, 2)
with open(RES + "A07.json", "w") as f: json.dump(results, f, indent=1)

def lmrow(r):
    L = r["landmarks"]; return " | ".join(f"{L[s]['cif12']:.2f} ({L[s]['ci'][0]:.1f}-{L[s]['ci'][1]:.1f}) [n={L[s]['nrisk']}]" for s in ["0", "3", "6", "12"])
def chg(r, key="change_6_minus_0"): c = r[key]; return f"{c['est']:+.2f} ({c['ci'][0]:+.1f} to {c['ci'][1]:+.1f})"
def rrs(x): return f"{x['RR']:.2f} ({x['lo']:.2f}-{x['hi']:.2f}), P={x['p']:.4f}"
md = []
md.append("# A07 Split-sample (calendar) validation - distant main colorectal NEC\n")
md.append("## Reproduction check (this run vs manuscript anchors)\n")
md.append("| Quantity | This run | Manuscript |\n|---|---|---|")
md.append(f"| Distant NEC n | {anchor['n_distant']} | 908 |")
md.append(f"| 12-m cancer-death % at landmarks 0/3/6/12/24 | {'/'.join(f'{x:.2f}' for x in anchor['cif_0_3_6_12_24'])} | 72.29/71.62/72.33/60.03/30.51 |")
md.append(f"| n at risk | {'/'.join(map(str, anchor['nrisk']))} | 908/553/396/169/48 |")
md.append(f"| Plateau end (months) | {anchor['plateau_end']} | 6.75 |")
md.append(f"| Crude cancer-death rate/person-month 0-3..18-24 | {'/'.join(f'{x:.4f}' for x in anchor['crude_rates'])} | 0.1435/0.0909/0.1344/0.1204/0.1231/0.0840/0.0646 |")
md.append(f"| Adjusted RR 6-9 vs 3-6 (R) | {rrs(rra)} | 1.49 (1.17-1.90), P=0.0013 |")
md.append(f"| Derivation / confirmation n | {r_early['n']} / {r_late['n']} | 431 / 477 |\n")
md.append(f"Method: {results['method']} B={B_MAIN} (seed {SEED}; confirmation seed {SEED+1}; eras {SEED+10}-{SEED+12}); comparators B={B_CMP} (seed {SEED+20}).\n")
md.append("## Table 1. Landmark 12-month cancer-death probability (%, 95% bootstrap CI) by period\n")
md.append("| Cohort | n | LM 0 | LM 3 | LM 6 | LM 12 | Change 6-0 | Change 12-0 | Plateau end (95% CI) |\n|---|---|---|---|---|---|---|---|---|")
for r in [r_early, r_late] + [results["eras"][e] for e in ["2000-2007", "2008-2015"]]:
    md.append(f"| {r['label']} | {r['n']} | {lmrow(r)} | {chg(r)} | {chg(r, 'change_12_minus_0')} | {r['plateau_end']['est']:g} ({r['plateau_end']['ci'][0]:g}-{r['plateau_end']['ci'][1]:g}) |")
pd_ = results["period_difference_in_6_minus_0"]; pp_ = results["period_difference_in_plateau_end"]
md.append(f"\nDifference in the 6-0 change, confirmation minus derivation: {pd_['est']:+.2f} points (95% CI {pd_['ci'][0]:+.1f} to {pd_['ci'][1]:+.1f}; bootstrap P = {pd_['boot_p']:.2f}). "
          f"Difference in plateau end: {pp_['est']:+g} months (95% CI {pp_['ci'][0]:g} to {pp_['ci'][1]:g}; P = {pp_['boot_p']:.2f}). "
          + "; ".join(f"LM{s}: {results[f'period_difference_cif_landmark_{s}']['est']:+.1f} ({results[f'period_difference_cif_landmark_{s}']['ci'][0]:+.1f} to {results[f'period_difference_cif_landmark_{s}']['ci'][1]:+.1f})" for s in [0,3,6,12]) + ".\n")
md.append("## Table 2. Crude cancer-death rates per person-month and adjusted RR by period\n")
md.append("| Period | n | 0-3 | 3-6 | 6-9 | 9-12 | 12-15 | 15-18 | 18-24 | crude ratio 6-9/3-6 | adjusted RR 6-9 vs 3-6 | adjusted RR 0-3 vs 3-6 |\n|---|---|---|---|---|---|---|---|---|---|---|---|")
for r, key in [(r_early, "derivation_2000_2015"), (r_late, "confirmation_2016_2023"), (results["eras"]["2000-2007"], "era_2000_2007"), (results["eras"]["2008-2015"], "era_2008_2015")]:
    rates = " | ".join(f"{x['rate']:.4f} ({x['deaths']})" for x in r["crude_rates"])
    md.append(f"| {r['label']} | {r['n']} | {rates} | {r['crude_rate_ratio_6to9_vs_3to6']:.2f} | {rrs(rr[key]['RR_6to9_vs_3to6'])} | {rrs(rr[key]['RR_0to3_vs_3to6'])} |")
rt = rr["ratio_of_RR_late_vs_early"]
md.append(f"\nRates shown as rate (cancer deaths). Adjusted RR: piecewise-exponential Poisson (0-24 months, cuts 3/6/9/12/15/18), covariates age10+age10sq+male+year10+large+rectum, cluster-robust SE by patient. Ratio of adjusted RRs (2016-2023 / 2000-2015) from a pooled period x interval interaction model: {rt['RR']:.2f} ({rt['lo']:.2f}-{rt['hi']:.2f}), P = {rt['p']:.3f}.\n")
md.append("## Table 3. Comparator tumours within 2016-2023 (B=500) and NEC-minus-comparator contrast in the 6-0 change\n")
md.append("| Cohort (distant, 2016-2023) | n | LM 0 | LM 6 | Change 6-0 (95% CI) | Plateau end (point) | NEC minus comparator, 6-0 change (95% CI; boot P) | 2000-2015 change 6-0 (point, n) |\n|---|---|---|---|---|---|---|---|")
md.append(f"| Colorectal NEC | {r_late['n']} | {r_late['landmarks']['0']['cif12']:.2f} | {r_late['landmarks']['6']['cif12']:.2f} | {chg(r_late)} | {r_late['plateau_end']['est']:g} | - | {r_early['change_6_minus_0']['est']:+.2f} (n={r_early['n']}) |")
for nm in ["SCLC", "AC", "PANC"]:
    r = comp[nm]; c = r["NEC_minus_comparator_change_6_minus_0"]; e = r["early_2000_2015_point"]
    md.append(f"| {nm} | {r['n']} | {r['landmarks']['0']['cif12']:.2f} | {r['landmarks']['6']['cif12']:.2f} | {chg(r)} | {r['plateau_end_point']:g} | {c['est']:+.2f} ({c['ci'][0]:+.1f} to {c['ci'][1]:+.1f}; P={c['boot_p']:.3f}) | {e['change_6_minus_0']:+.2f} (n={e['n']}) |")
md.append(f"\nFigure: `{config.relpath(FIG + 'A07_derivation_vs_confirmation.png')}` (landmark curves 0-24 months with 95% bootstrap bands; points with <20 at risk suppressed).\n")
# interpretation
e6, l6 = r_early["change_6_minus_0"], r_late["change_6_minus_0"]
md.append("## Interpretation\n")
md.append(f"1. The conditional-risk plateau replicates in both calendar halves: the 12-month cancer-death probability from landmark 6 is essentially unchanged from landmark 0 in 2000-2015 ({e6['est']:+.2f} points, 95% CI {e6['ci'][0]:+.1f} to {e6['ci'][1]:+.1f}) and in 2016-2023 ({l6['est']:+.2f}, {l6['ci'][0]:+.1f} to {l6['ci'][1]:+.1f}); the estimated plateau end is {r_early['plateau_end']['est']:g} and {r_late['plateau_end']['est']:g} months respectively, and the between-period difference in the 6-0 change is {pd_['est']:+.2f} points ({pd_['ci'][0]:+.1f} to {pd_['ci'][1]:+.1f}; P = {pd_['boot_p']:.2f}), i.e. no evidence that the pattern is confined to one era.")
md.append(f"2. The 6-9 month hazard rise is clearer in the later period: adjusted RR 6-9 vs 3-6 = {rrs(rr['confirmation_2016_2023']['RR_6to9_vs_3to6'])} in 2016-2023 versus {rrs(rr['derivation_2000_2015']['RR_6to9_vs_3to6'])} in 2000-2015 (ratio {rt['RR']:.2f}, {rt['lo']:.2f}-{rt['hi']:.2f}, P = {rt['p']:.3f}). The direction is the same in both halves but the derivation-period estimate is imprecise; the full-cohort RR of 1.49 is therefore driven mainly by 2016-2023 and should be presented as replicated in direction, not magnitude.")
md.append(f"3. Three-era stratification: plateau ends {results['eras']['2000-2007']['plateau_end']['est']:g} / {results['eras']['2008-2015']['plateau_end']['est']:g} / {results['eras']['2016-2023']['plateau_end']['est']:g} months and 6-0 changes {results['eras']['2000-2007']['change_6_minus_0']['est']:+.1f} / {results['eras']['2008-2015']['change_6_minus_0']['est']:+.1f} / {results['eras']['2016-2023']['change_6_minus_0']['est']:+.1f} points; the plateau-end CIs in the two early eras are wide (n = 168 and 263) and the hinge estimate is unstable there, so the three-era table is descriptive only.")
md.append(f"4. The cross-tumour contrast holds within 2016-2023: NEC minus AC = {comp['AC']['NEC_minus_comparator_change_6_minus_0']['est']:+.1f} points ({comp['AC']['NEC_minus_comparator_change_6_minus_0']['ci'][0]:+.1f} to {comp['AC']['NEC_minus_comparator_change_6_minus_0']['ci'][1]:+.1f}), NEC minus PANC = {comp['PANC']['NEC_minus_comparator_change_6_minus_0']['est']:+.1f} ({comp['PANC']['NEC_minus_comparator_change_6_minus_0']['ci'][0]:+.1f} to {comp['PANC']['NEC_minus_comparator_change_6_minus_0']['ci'][1]:+.1f}), NEC minus SCLC = {comp['SCLC']['NEC_minus_comparator_change_6_minus_0']['est']:+.1f} ({comp['SCLC']['NEC_minus_comparator_change_6_minus_0']['ci'][0]:+.1f} to {comp['SCLC']['NEC_minus_comparator_change_6_minus_0']['ci'][1]:+.1f}). AC and PANC show the expected monotone decline in conditional risk from diagnosis onward (plateau end 0) in the recent period, whereas NEC and SCLC do not.")
md.append("5. Caveats: this is a temporal (not geographic or external) validation; the two halves share the same registry and coding rules, and 2016-2023 has shorter follow-up (late landmarks rely on fewer patients). The derivation/confirmation labels follow the task plan; no model was fitted on one half and applied to the other beyond the pre-specified landmark/hinge procedure.")
md.append("\nDeviations from plan: none in the main analyses. Additions: bootstrap CIs (B=2000) for the three eras, a pooled period x interval interaction test for the adjusted RR, 12-0 changes, and 2000-2015 point estimates for the comparators (no bootstrap) for context. Comparator contrasts pair the first 500 of the NEC 2016-2023 replicates with the 500 comparator replicates.")
md.append(f"\nRuntime: {results['runtime_minutes']} min (single thread).")
open(RES + "A07.md", "w").write("\n".join(md))
print("\n".join(md)); print("DONE", results["runtime_minutes"], "min")
