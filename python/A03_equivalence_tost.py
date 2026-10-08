# -*- coding: utf-8 -*-
"""A03 equivalence (TOST) framing of 'no decline' in the 12-month cancer-death probability.
Distant colorectal NEC (main definition) landmark differences 0->3, 0->6, 3->6, 0->12 months, patient-level
bootstrap (B=2000, seed 20261007), 90% percentile CI vs margins +-3/4/5/6.7 points; bootstrap one-sided TOST p-values.
Paired cross-tumour contrasts (NEC 6-0 change minus SCLC / AC / PANC 6-0 change), independent bootstraps paired by index.
Revision after audit (2026-10-07): wording labels margins post hoc and anchors on +-5; P < 1/B shown as '<0.001';
figure titles fixed and +-4 margin drawn; cross-tumour 'largest decline excluded' replaced by lower limits;
referenced anchors (plateau end) stored in JSON; Holm-adjusted TOST P added as supplementary column.
Run: python3 python/A03_equivalence_tost.py
"""
import sys, time, json, numpy as np, pandas as pd
import config
import lm
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

T0 = time.time()
B = config.B(2000)
SEEDS = {"NEC": 20261007, "SCLC": 20260928 + 1, "AC": 20260928 + 2, "PANC": 20260928 + 3}  # AC/PANC seeds: not specified in task, chosen by the same +k convention
LM = np.arange(0, 25); WIN = 12
IDX = {0: 0, 3: 3, 6: 6, 12: 12}
PRIMARY_MARGIN = 5          # clinically motivated post hoc margin (~ the observed 6.7-point colorectal-AC decline)
REFERENCED = {"plateau_end_months": 6.75, "source": "task reproduction anchors / A01 broken-stick fit (not computed in A03)"}

def fp(p, B=B):
    """format a bootstrap P-value; resolution is 1/B so values below 0.5/B are reported as '<0.001'"""
    return "<0.001" if p < 0.0005 else f"{p:.3f}"

# ---------------- cohorts ----------------
nen = lm.load("NEN_all"); nec = lm.nec_main(nen); nec_d = nec[nec.stage == "Distant"]
sclc_d = lm.sclc_main(lm.load("LUNG")).query("stage == 'Distant'")
ac_d = lm.load("AC").query("stage == 'Distant'")
panc_d = lm.panc_main(lm.load("PANC")).query("stage == 'Distant'")
cohorts = {"NEC": nec_d, "SCLC": sclc_d, "AC": ac_d, "PANC": panc_d}
print({k: len(v) for k, v in cohorts.items()})

# ---------------- reproduction anchors ----------------
anch = {"NEC": {"cif": {0: 72.29, 3: 71.62, 6: 72.33, 12: 60.03, 24: 30.51}, "nrisk": {0: 908, 3: 553, 6: 396, 12: 169, 24: 48}},
        "SCLC": {"cif": {0: 69.98, 6: 69.49, 12: 61.93}}, "AC": {"cif": {0: 40.59, 6: 33.92, 12: 34.41}}, "PANC": {"cif": {0: 77.12, 6: 69.72, 12: 65.38}}}
point = {}; repro = []
for nm, df in cohorts.items():
    k, ev, KMAX = lm.prep(df); cif, allc, nr = lm.curves(k, ev, KMAX, LM, WIN)
    point[nm] = {"cif": cif, "nrisk": nr, "n": int(len(df))}
    assert np.all(np.diff(nr) <= 0), "n at risk not monotone"
    for s, v in anch[nm]["cif"].items():
        ok = abs(cif[s] - v) < 0.01; repro.append({"cohort": nm, "quantity": f"CIF12 at landmark {s}", "mine": round(float(cif[s]), 2), "manuscript": v, "ok": bool(ok)})
        if not ok: config.anchor_fail(f"REPRODUCTION FAILURE {nm} landmark {s}: {cif[s]:.2f} vs {v}")
    for s, v in anch[nm].get("nrisk", {}).items():
        ok = int(nr[s]) == v; repro.append({"cohort": nm, "quantity": f"n at risk landmark {s}", "mine": int(nr[s]), "manuscript": v, "ok": bool(ok)})
        if not ok: config.anchor_fail(f"REPRODUCTION FAILURE {nm} n at risk {s}: {nr[s]} vs {v}")
print("Reproduction check passed:", all(r["ok"] for r in repro))

# ---------------- bootstraps ----------------
boot = {}
for nm, df in cohorts.items():
    t = time.time(); bc, _ = lm.bootstrap(df, B=B, seed=SEEDS[nm], landmarks=LM, window=WIN); boot[nm] = bc
    print(f"bootstrap {nm}: n={len(df)} B={B} seed={SEEDS[nm]} {time.time()-t:.1f}s; NaN rows at LM 0/6: {np.isnan(bc[:, [0, 6]]).any(axis=1).sum()}")

def diff(bc, a, b): return bc[:, IDX[b]] - bc[:, IDX[a]]          # later minus earlier landmark: negative = decline
def tost(d, margins, label, est, within=True):
    ci90 = lm.pct(d, (5, 95)); ci95 = lm.pct(d, (2.5, 97.5))
    out = {"contrast": label, "estimate": float(est), "boot_mean": float(np.mean(d)), "boot_sd": float(np.std(d, ddof=1)),
           "ci90": ci90, "ci95": ci95, "B": int(len(d)), "p_resolution": 1.0 / len(d), "margins": []}
    if within:
        out["largest_decline_excluded_95"] = float(-ci95[0]); out["largest_decline_excluded_90"] = float(-ci90[0])
    else:   # cross-tumour: a positive lower limit means the NEC change exceeds the comparator change; 'excluded decline' is not meaningful
        out["lower_limit_95"] = float(ci95[0]); out["lower_limit_90"] = float(ci90[0])
        out["lower_limit_note"] = "positive lower limit = NEC 6-0 change exceeds the comparator 6-0 change by at least this amount"
    for m in margins:
        p_lo = float(np.mean(d <= -m))   # H0a: true diff <= -m  (decline at least m)  -> bootstrap p ~ P*(diff* <= -m)
        p_hi = float(np.mean(d >= m))    # H0b: true diff >= +m
        out["margins"].append({"margin": m, "ci90_within": bool(ci90[0] > -m and ci90[1] < m), "p_lower": p_lo, "p_upper": p_hi,
                               "p_tost": max(p_lo, p_hi), "equivalent_alpha05": bool(max(p_lo, p_hi) < 0.05)})
    return out

def holm(ps):
    """Holm step-down adjusted P-values for a family of P-values (list), returned in the original order"""
    ps = np.asarray(ps, float); n = len(ps); order = np.argsort(ps); adj = np.empty(n)
    running = 0.0
    for rank, i in enumerate(order):
        running = max(running, (n - rank) * ps[i]); adj[i] = min(1.0, running)
    return adj.tolist()

MARG = [3, 4, 5, 6.7]
res = {"id": "A03", "B": B, "seeds": SEEDS, "window_months": WIN, "n": {k: v["n"] for k, v in point.items()},
       "margins_status": "post hoc (not pre-specified in the 2026-09-28 protocol); primary clinically anchored margin +-5 points, approx. the observed 6.7-point colorectal-AC decline",
       "primary_margin": PRIMARY_MARGIN, "referenced_anchors": REFERENCED,
       "point_cif12": {k: {str(s): round(float(v["cif"][s]), 2) for s in (0, 3, 6, 12, 24)} for k, v in point.items()},
       "nrisk": {k: {str(s): int(v["nrisk"][s]) for s in (0, 3, 6, 12, 24)} for k, v in point.items()},
       "reproduction": repro, "within_NEC": [], "cross_tumour": []}
pairs = [(0, 6), (0, 3), (3, 6), (0, 12)]
for a, b in pairs:
    d = diff(boot["NEC"], a, b); est = point["NEC"]["cif"][b] - point["NEC"]["cif"][a]
    res["within_NEC"].append(tost(d, MARG, f"NEC {b} minus {a} months", est))
# Holm adjustment across the four within-NEC contrasts, separately for each margin (family = 4 contrasts); supplementary only
for j, m in enumerate(MARG):
    adj = holm([r["margins"][j]["p_tost"] for r in res["within_NEC"]])
    for r, pa in zip(res["within_NEC"], adj): r["margins"][j]["p_tost_holm_across_contrasts"] = float(pa)
res["holm_note"] = "Holm step-down applied to the bootstrap TOST P-values of the four within-NEC contrasts within each margin (family size 4); supplementary, the primary criterion remains the 90% percentile CI inside the margin"

dn = diff(boot["NEC"], 0, 6); en = point["NEC"]["cif"][6] - point["NEC"]["cif"][0]
for nm in ("SCLC", "AC", "PANC"):
    dc = diff(boot[nm], 0, 6); ec = point[nm]["cif"][6] - point[nm]["cif"][0]
    r = tost(dn - dc, [3, 5, 6.7], f"(NEC 6-0 change) minus ({nm} 6-0 change)", en - ec, within=False)
    r["comparator"] = nm
    r["comparator_change"] = {"estimate": float(ec), "ci95": lm.pct(dc)}; r["nec_change"] = {"estimate": float(en), "ci95": lm.pct(dn)}
    res["cross_tumour"].append(r)

# ---------------- wording ----------------
def pm(r, m): return next(x for x in r["margins"] if x["margin"] == m)
def word_within(r):
    L = r["largest_decline_excluded_95"]
    if r["ci95"][1] < 0:
        return (f"{r['contrast']}: the 12-month cancer-death probability declined by {-r['estimate']:.1f} points (95% CI {r['ci95'][0]:.1f} to {r['ci95'][1]:.1f}); "
                f"a decline was detected (upper 95% limit below 0) and equivalence to 'no change' is rejected at every post hoc margin tested (90% CI {r['ci90'][0]:.1f} to {r['ci90'][1]:.1f}).")
    s = (f"{r['contrast']}: the 12-month cancer-death probability changed by {r['estimate']:+.2f} points (95% CI {r['ci95'][0]:.1f} to {r['ci95'][1]:.1f}); "
         f"no decline was detected, and a decline larger than {L:.1f} points is excluded by the lower 95% limit. ")
    p5 = pm(r, PRIMARY_MARGIN); p4 = pm(r, 4); p3 = pm(r, 3)
    if p5["ci90_within"]:
        s += (f"Under a post hoc equivalence margin of ±{PRIMARY_MARGIN} points (approximately the observed 6.7-point colorectal-adenocarcinoma decline), "
              f"the change was equivalent to zero (90% CI {r['ci90'][0]:.1f} to {r['ci90'][1]:.1f}; bootstrap TOST P {'= ' + fp(p5['p_tost']) if p5['p_tost'] >= 0.0005 else fp(p5['p_tost'])}); "
              f"equivalence {'also held' if p4['ci90_within'] else 'did not hold'} at ±4 points (P {'= ' + fp(p4['p_tost']) if p4['p_tost'] >= 0.0005 else fp(p4['p_tost'])}) "
              f"and {'held' if p3['ci90_within'] else 'did not hold'} at ±3 points (P {'= ' + fp(p3['p_tost']) if p3['p_tost'] >= 0.0005 else fp(p3['p_tost'])}).")
    else:
        s += (f"The 90% CI ({r['ci90'][0]:.1f} to {r['ci90'][1]:.1f}) did not lie within the post hoc margin of ±{PRIMARY_MARGIN} points, so equivalence to 'no change' could not be established.")
    return s
def word_cross(r):
    m_eq = [m["margin"] for m in r["margins"] if m["ci90_within"]]
    s = (f"{r['contrast']}: {r['estimate']:+.2f} points (95% CI {r['ci95'][0]:.1f} to {r['ci95'][1]:.1f}). ")
    if m_eq:
        mm = PRIMARY_MARGIN if PRIMARY_MARGIN in m_eq else min(m_eq); x = pm(r, mm)
        s += (f"No difference between the 6-month changes was detected, and under a post hoc margin of ±{mm:g} points they were equivalent "
              f"(90% CI {r['ci90'][0]:.1f} to {r['ci90'][1]:.1f}; bootstrap TOST P {'= ' + fp(x['p_tost']) if x['p_tost'] >= 0.0005 else fp(x['p_tost'])}); "
              f"equivalence {'also held' if pm(r, 3)['ci90_within'] else 'did not hold'} at ±3 points (P {'= ' + fp(pm(r, 3)['p_tost']) if pm(r, 3)['p_tost'] >= 0.0005 else fp(pm(r, 3)['p_tost'])}). "
              f"'Similar plateau' should therefore be worded as 'no difference detected; equivalent within ±{mm:g} points'.")
    elif r["ci95"][0] > 0:
        s += (f"The NEC change exceeded the comparator change (lower 95% limit {r['ci95'][0]:.1f} > 0); the two were not equivalent at any post hoc margin tested "
              f"(bootstrap TOST P = {fp(pm(r, 6.7)['p_tost'])} even at ±6.7 points).")
    else:
        s += f"No difference was detected but equivalence within ±{min(x['margin'] for x in r['margins']):g} points was not established (90% CI {r['ci90'][0]:.1f} to {r['ci90'][1]:.1f}); 'similar plateau' should be worded as 'no difference detected'."
    return s
for r in res["within_NEC"]: r["wording"] = word_within(r)
for r in res["cross_tumour"]: r["wording"] = word_cross(r)
res["methods_sentence"] = ("Equivalence margins (±3, ±4, ±5 and ±6.7 percentage points) were chosen after the analysis protocol of 28 September 2026 and are reported as post hoc "
                           "sensitivity framings; ±5 points, approximately the observed 6.7-point decline in colorectal adenocarcinoma, is used as the primary clinically anchored margin. "
                           "TOST P-values are bootstrap approximations (proportion of 2000 patient-level resamples at or beyond each bound; resolution 1/2000, values below this are reported as <0.001).")

# ---------------- figure ----------------
fig, axes = plt.subplots(2, 4, figsize=(16, 7))
items = [(r["contrast"], diff(boot["NEC"], a, b)) for (a, b), r in zip(pairs, res["within_NEC"])] + \
        [(f"NEC minus {r['comparator']} (6-0 change)", dn - diff(boot[r["comparator"]], 0, 6)) for r in res["cross_tumour"]]
MCOL = {3: "#DD8452", 4: "#55A868", 5: "#C44E52", 6.7: "#8172B2"}
for ax, (lab, d) in zip(axes.ravel(), items):
    ax.hist(d, bins=60, color="#4C72B0", alpha=.8); ax.axvline(0, color="k", lw=.8)
    for m, c in MCOL.items(): ax.axvline(-m, color=c, ls="--", lw=.9); ax.axvline(m, color=c, ls="--", lw=.9)
    lo, hi = np.percentile(d, (5, 95)); ax.axvspan(lo, hi, color="grey", alpha=.25)
    ax.set_title(lab, fontsize=9); ax.set_xlabel("percentage points")
axes.ravel()[-1].axis("off")
axes.ravel()[-1].text(0, .5, "grey band = 90% percentile CI\ndashed margins: " + ", ".join(f"±{m:g}" for m in MCOL) + "\n(orange / green / red / purple)\nB = 2000 patient-level resamples", fontsize=10)
fig.suptitle("A03 Bootstrap distributions of landmark differences in 12-month cancer-death probability (distant stage)")
fig.tight_layout(); fig.savefig(config.FIGS + "A03_tost_distributions.png", dpi=300); plt.close(fig)

# ---------------- outputs ----------------
res["runtime_minutes"] = round((time.time() - T0) / 60, 2)
json.dump(res, open(config.RES + "A03.json", "w"), indent=1)

md = ["# A03 Equivalence (TOST) framing of 'no decline' — distant-stage 12-month cancer-death probability", "",
      f"Patient-level bootstrap, B = {B}, percentile intervals; seeds: " + ", ".join(f"{k} {v}" for k, v in SEEDS.items()) + f". Runtime {res['runtime_minutes']} min.", "",
      "Margins ±3, ±4, ±5 and ±6.7 points are post hoc (not in the 2026-09-28 protocol); ±5 (approximately the observed 6.7-point colorectal-AC decline) is the primary clinically anchored margin.", "",
      "## Reproduction check (mine vs manuscript anchors)", "", "| Cohort | Quantity | Mine | Manuscript | OK |", "|---|---|---|---|---|"]
md += [f"| {r['cohort']} | {r['quantity']} | {r['mine']} | {r['manuscript']} | {'yes' if r['ok'] else 'NO'} |" for r in repro]
md += ["", "## Table 1. Within-NEC landmark differences (later minus earlier landmark; negative = decline), n = 908 distant NEC", "",
       "| Contrast | Estimate (pp) | 95% CI | 90% CI | Largest decline excluded (lower 95% limit) | " + " | ".join(f"±{m:g}: within / P_TOST" for m in MARG) + " |",
       "|---|---|---|---|---|" + "---|" * len(MARG)]
for r in res["within_NEC"]:
    md.append(f"| {r['contrast']} | {r['estimate']:+.2f} | {r['ci95'][0]:.2f} to {r['ci95'][1]:.2f} | {r['ci90'][0]:.2f} to {r['ci90'][1]:.2f} | {r['largest_decline_excluded_95']:.1f} | " +
              " | ".join(f"{'yes' if m['ci90_within'] else 'no'} / {fp(m['p_tost'])}" for m in r["margins"]) + " |")
md += ["", f"P_TOST = max of the two one-sided bootstrap P-values (proportion of resamples at or beyond the bound); resolution 1/B = {1/B:.4f}, values below 0.0005 shown as <0.001. 'within' = 90% percentile CI strictly inside ±margin (TOST at alpha 0.05).", "",
       "One-sided bootstrap P-values and Holm-adjusted TOST P (family = the four contrasts within each margin; supplementary):", "",
       "| Contrast | Margin | P(decline ≥ margin) | P(increase ≥ margin) | P_TOST | P_TOST Holm (4 contrasts) |", "|---|---|---|---|---|---|"]
for r in res["within_NEC"]:
    for m in r["margins"]: md.append(f"| {r['contrast']} | ±{m['margin']:g} | {fp(m['p_lower'])} | {fp(m['p_upper'])} | {fp(m['p_tost'])} | {fp(m['p_tost_holm_across_contrasts'])} |")
md += ["", "## Table 2. Paired cross-tumour contrasts: (NEC 6-0 change) minus (comparator 6-0 change), independent bootstraps paired by index", "",
       "| Contrast | NEC change (95% CI) | Comparator change (95% CI) | Difference (pp) | 95% CI | 90% CI | ±3: within / P | ±5: within / P | ±6.7: within / P |", "|---|---|---|---|---|---|---|---|---|"]
for r in res["cross_tumour"]:
    md.append(f"| {r['contrast']} | {r['nec_change']['estimate']:+.2f} ({r['nec_change']['ci95'][0]:.2f} to {r['nec_change']['ci95'][1]:.2f}) | {r['comparator_change']['estimate']:+.2f} ({r['comparator_change']['ci95'][0]:.2f} to {r['comparator_change']['ci95'][1]:.2f}) | {r['estimate']:+.2f} | {r['ci95'][0]:.2f} to {r['ci95'][1]:.2f} | {r['ci90'][0]:.2f} to {r['ci90'][1]:.2f} | " +
              " | ".join(f"{'yes' if m['ci90_within'] else 'no'} / {fp(m['p_tost'])}" for m in r["margins"]) + " |")
md += ["", "A positive lower 95% limit (AC, PANC) means the NEC change exceeds the comparator change; 'largest decline excluded' is not reported for these contrasts.", "",
       "## Recommended manuscript wording", ""] + [f"- {r['wording']}" for r in res["within_NEC"] + res["cross_tumour"]]
md += ["", "Methods sentence: " + res["methods_sentence"], "", "## Interpretation", ""]
w06 = res["within_NEC"][0]; w012 = res["within_NEC"][3]; cs = res["cross_tumour"][0]
m_eq06 = [m['margin'] for m in w06['margins'] if m['ci90_within']]
md += [f"1. The 0-to-6-month change in distant NEC is {w06['estimate']:+.2f} points (95% CI {w06['ci95'][0]:.1f} to {w06['ci95'][1]:.1f}; 90% CI {w06['ci90'][0]:.1f} to {w06['ci90'][1]:.1f}). Under the primary post hoc margin of ±{PRIMARY_MARGIN} points equivalence to 'no change' is " + ("established" if pm(w06, PRIMARY_MARGIN)['ci90_within'] else "not established") + f" (bootstrap TOST P = {fp(pm(w06, PRIMARY_MARGIN)['p_tost'])}); it also holds at ±{(min(m_eq06) if m_eq06 else float('nan')):g} points and wider, and at ±3 points it is " + ("also established." if w06['margins'][0]['ci90_within'] else "not established (TOST P = %s)." % fp(w06['margins'][0]['p_tost'])),
       f"2. The lower 95% limit excludes a decline larger than {w06['largest_decline_excluded_95']:.1f} points, in line with the manuscript's 'about 4 points'. The 0-to-3 and 3-to-6 contrasts behave similarly; the 0-to-12 contrast ({w012['estimate']:+.1f} points) is a clear decline and is non-equivalent at every margin, as expected after the plateau ends ({REFERENCED['plateau_end_months']} months, broken-stick fit from A01 / task anchors; not computed here).",
       f"3. Against SCLC the paired difference in 6-month changes is {cs['estimate']:+.2f} points (90% CI {cs['ci90'][0]:.1f} to {cs['ci90'][1]:.1f}): " + ("equivalent within ±5 but not ±3" if cs['margins'][1]['ci90_within'] and not cs['margins'][0]['ci90_within'] else ("equivalent within ±3 and ±5" if cs['margins'][0]['ci90_within'] else "not equivalent within ±5")) + ". The 'similar plateau' claim should therefore be stated with the margin at which equivalence holds.",
       "4. Against colorectal AC and pancreatic AC the paired differences are entirely above zero and non-equivalent at every margin, so the NEC plateau differs from the immediate decline of both adenocarcinomas.",
       f"5. TOST P-values are bootstrap approximations (proportion of resamples beyond each bound; resolution 1/B), not analytic tests; the primary criterion is the 90% percentile CI lying inside the margin. The margins were not pre-specified in the 2026-09-28 protocol and are presented as post hoc sensitivity framings with ±{PRIMARY_MARGIN} as the clinically anchored primary margin. Holm adjustment across the four within-NEC contrasts (supplementary table) leaves the conclusions unchanged: at ±{PRIMARY_MARGIN} the three plateau-phase contrasts remain below 0.05 after adjustment (" + ", ".join(f"{fp(r['margins'][2]['p_tost_holm_across_contrasts'])}" for r in res["within_NEC"][:3]) + ").",
       "", "Deviations: AC and PANC seeds were not given in the task; 20260930 and 20260931 (seed 20260928 + 2/+3) were used. Bootstrap P-values use the resampling distribution directly (percentile/'proportion beyond bound'), which is the approximation requested. Holm adjustment was not requested in the A03 specification and is supplied as a supplementary column only.",
       "", "## Revision after audit (2026-10-07)", "",
       "- Wording: 'pre-specified' replaced by 'post hoc'; the primary statement is anchored on the clinically motivated ±5 margin (approx. the 6.7-point colorectal-AC decline) instead of the narrowest passing margin, with ±4 and ±3 reported alongside; a Methods sentence on the post hoc margins was added; estimates in the wording are given with two decimals (e.g. +0.03, not +0.0).",
       "- P-values below the bootstrap resolution (0.5/B) are printed as '<0.001' in all tables and wording; the 1/B resolution is stated in the Table 1 footnote.",
       "- Cross-tumour entries: the meaningless 'largest_decline_excluded' fields were replaced by lower_limit_95 / lower_limit_90 with an explanatory note.",
       "- Figure: panel titles for the cross-tumour contrasts fixed ('NEC minus SCLC (6-0 change)' etc.); the ±4 margin added to the dashed lines and legend.",
       "- JSON now carries a referenced_anchors block (plateau end 6.75 months, from A01 / task anchors) and the margins_status / primary_margin fields.",
       "- Holm-adjusted TOST P-values across the four within-NEC contrasts (per margin) added as a supplementary column; no change to any estimate, CI or unadjusted P-value."]
open(config.RES + "A03.md", "w").write("\n".join(md))
print("\n".join(md[-22:])); print("done", res["runtime_minutes"], "min")
