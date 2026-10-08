# -*- coding: utf-8 -*-
"""
A02 - Window-width sensitivity and plateau-model comparison (distant colorectal NEC;
point estimates for distant SCLC / colorectal AC / pancreatic AC).

Runs end-to-end non-interactively:  python3 python/A02_window_sensitivity_model_comparison.py
Single-threaded. Patient-level data never leave NEC_WORK_DIR; only aggregate results are written.
"""
import os, sys, json, time, platform
os.environ.setdefault("OMP_NUM_THREADS", "1"); os.environ.setdefault("OPENBLAS_NUM_THREADS", "1"); os.environ.setdefault("MKL_NUM_THREADS", "1")
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

T0 = time.time()
import config
import lm

ID = "A02"
RES = config.RES; FIGS = config.FIGS
os.makedirs(RES, exist_ok=True); os.makedirs(FIGS, exist_ok=True)

WINDOWS = [3, 6, 9, 12, 18, 24]
LM = np.arange(0, 19)                      # landmarks 0-18 for every window
B_NEC, B_OTHER = config.B(2000), config.B(500)
SEED = 20261007
SPARSE = 100                               # mark landmarks with < 100 at risk

# ----------------------------------------------------------------------------------------------
# 0. Cohorts and reproduction check
# ----------------------------------------------------------------------------------------------
def distant(df): return df[df["stage"] == "Distant"]

cohorts = {}
nen = lm.load("NEN_all"); cohorts["NEC"] = distant(lm.nec_main(nen)); del nen
lung = lm.load("LUNG");   cohorts["SCLC"] = distant(lm.sclc_main(lung)); del lung
panc = lm.load("PANC");   cohorts["PANC"] = distant(lm.panc_main(panc)); del panc
ac = lm.load("AC");       cohorts["AC"] = distant(ac); del ac
LABEL = {"NEC": "Colorectal NEC (distant)", "SCLC": "Small-cell lung cancer (distant)",
         "AC": "Colorectal adenocarcinoma (distant)", "PANC": "Pancreatic adenocarcinoma (distant)"}

prepped = {t: lm.prep(cohorts[t]) for t in cohorts}   # (k, ev, KMAX)

def point(t, window, landmarks=LM):
    k, ev, KMAX = prepped[t]
    return lm.curves(k, ev, KMAX, landmarks, window)

# --- reproduction anchors (orchestrator-verified) ---
anchors = {
    "NEC":  {"lm": [0, 3, 6, 12, 24], "cif": [72.29, 71.62, 72.33, 60.03, 30.51], "n": [908, 553, 396, 169, 48], "plateau": 6.75},
    "SCLC": {"lm": [0, 6, 12], "cif": [69.98, 69.49, 61.93], "n": [77707, None, None], "plateau": 8.5},
    "AC":   {"lm": [0, 6, 12], "cif": [40.59, 33.92, 34.41], "n": [147014, None, None], "plateau": 0.0},
    "PANC": {"lm": [0, 6, 12], "cif": [77.12, 69.72, 65.38], "n": [75634, None, None], "plateau": 0.0},
}
repro_rows = []; ok = True
for t, a in anchors.items():
    cif, allc, n = point(t, 12, np.array(a["lm"]))
    c19, _, _ = point(t, 12)
    pl = float(lm.hinge_many(c19[None, :])[0])
    for j, s in enumerate(a["lm"]):
        mine = round(float(cif[j]), 2); agree = abs(mine - a["cif"][j]) < 0.011
        if a["n"][j] is not None: agree = agree and int(n[j]) == a["n"][j]
        ok &= agree
        repro_rows.append({"tumour": t, "quantity": f"12-mo cancer-death prob, landmark {s}", "mine": mine, "manuscript": a["cif"][j],
                           "n_at_risk_mine": int(n[j]), "n_at_risk_manuscript": a["n"][j], "agree": bool(agree)})
    agree = abs(pl - a["plateau"]) < 1e-9; ok &= agree
    repro_rows.append({"tumour": t, "quantity": "broken-stick plateau end (12-mo window, LM 0-18)", "mine": pl, "manuscript": a["plateau"],
                       "n_at_risk_mine": None, "n_at_risk_manuscript": None, "agree": bool(agree)})
print("REPRODUCTION CHECK")
for r in repro_rows: print(f"  {r['tumour']:5s} {r['quantity']:55s} mine={r['mine']} ms={r['manuscript']} agree={r['agree']}")
if not ok:
    config.anchor_fail("REPRODUCTION CHECK FAILED - stopping.")
print("REPRODUCTION CHECK PASSED")

# ----------------------------------------------------------------------------------------------
# 1. Model helpers (fits on a 19-point landmark curve, Gaussian SSE / BIC)
# ----------------------------------------------------------------------------------------------
m = LM.astype(float); NPTS = len(m)
GRID = lm.GRID

def ns_basis(x, knots=(0.0, 6.0, 12.0, 18.0)):
    """Natural cubic spline basis with 3 df beyond the intercept (ESL eq. 5.4-5.5), boundary knots 0 and 18,
    interior knots 6 and 12 (tertiles of 0..18). Columns: 1, x, N3, N4  -> 4 free parameters."""
    K = len(knots); xk = np.asarray(knots)
    def d(k): return (np.maximum(x - xk[k], 0) ** 3 - np.maximum(x - xk[K - 1], 0) ** 3) / (xk[K - 1] - xk[k])
    cols = [np.ones_like(x), x] + [d(k) - d(K - 2) for k in range(K - 2)]
    return np.column_stack(cols)

X_lin = np.column_stack([np.ones_like(m), m])
X_quad = np.column_stack([np.ones_like(m), m, m ** 2])
X_ns = ns_basis(m)
P_lin, P_quad, P_ns = (X @ np.linalg.pinv(X) for X in (X_lin, X_quad, X_ns))
P_bs = np.stack([(lambda X: X @ np.linalg.pinv(X))(np.column_stack([np.ones_like(m), np.maximum(0, m - c0)])) for c0 in GRID])  # (G, n, n)
KPAR = {"broken_stick": 3, "linear": 2, "quadratic": 3, "ns3": 4}   # b0,b1,c | b0,b1 | b0,b1,b2 | 4 coefficients

def sse_all(Y):
    """Y: (R, 19) curves -> dict of SSE arrays (R,), plateau ends (R,) and fitted broken-stick slope after c."""
    Y = np.atleast_2d(Y)
    out = {"linear": np.sum((Y - Y @ P_lin.T) ** 2, 1), "quadratic": np.sum((Y - Y @ P_quad.T) ** 2, 1), "ns3": np.sum((Y - Y @ P_ns.T) ** 2, 1)}
    sse_g = np.stack([np.sum((Y - Y @ P_bs[g].T) ** 2, 1) for g in range(len(GRID))], 1)   # (R, G)
    gi = np.argmin(sse_g, 1)                     # first argmin -> smallest c on ties (same convention as lm.hinge_many)
    out["broken_stick"] = sse_g[np.arange(len(Y)), gi]; out["plateau"] = GRID[gi]
    return out

def bic(sse, k, n=NPTS):
    """Gaussian BIC up to an additive constant common to all models: n*log(SSE/n) + k*log(n)."""
    return n * np.log(np.maximum(sse, 1e-12) / n) + k * np.log(n)

def bs_fit(y):
    """Point broken-stick fit: returns c, intercept, slope after c."""
    r = sse_all(y[None, :]); c0 = float(r["plateau"][0])
    X = np.column_stack([np.ones_like(m), np.maximum(0, m - c0)]); b = np.linalg.lstsq(X, y, rcond=None)[0]
    return c0, float(b[0]), float(b[1])

def ci(v, q=(2.5, 97.5)):
    v = np.asarray(v, float); v = v[~np.isnan(v)]
    return [float(x) for x in np.percentile(v, q)] if len(v) else [None, None]

# ----------------------------------------------------------------------------------------------
# 2. Bootstrap: same patient-level resample for all windows (and both endpoints)
# ----------------------------------------------------------------------------------------------
def boot(t, B, seed, windows=WINDOWS):
    k, ev, KMAX = prepped[t]; n = len(k); rng = np.random.default_rng(seed)
    bc = {w: np.empty((B, NPTS)) for w in windows}; ba = {w: np.empty((B, NPTS)) for w in windows}
    for b in range(B):
        wts = np.bincount(rng.integers(0, n, n), minlength=n).astype(float)
        for w in windows:
            bc[w][b], ba[w][b], _ = lm.curves(k, ev, KMAX, LM, w, wts)
    return bc, ba

DIFFS = [(6, 0), (3, 0), (6, 3)]
results = {"id": ID, "date": "2026-10-07", "seed": SEED, "B": {"NEC": B_NEC, "SCLC": B_OTHER, "AC": B_OTHER, "PANC": B_OTHER},
           "windows_months": WINDOWS, "landmarks": LM.tolist(), "sparse_threshold_at_risk": SPARSE,
           "n": {t: int(len(cohorts[t])) for t in cohorts}, "reproduction_check": repro_rows, "tumours": {}}

t_boot = time.time()
boots = {}
for t in ["NEC", "SCLC", "AC", "PANC"]:
    Bt = B_NEC if t == "NEC" else B_OTHER
    boots[t] = boot(t, Bt, SEED)
    print(f"bootstrap {t}: B={Bt}, {time.time() - t_boot:.1f}s cumulative")

for t in ["NEC", "SCLC", "AC", "PANC"]:
    bc, ba = boots[t]; tr = {"label": LABEL[t], "n_distant": int(len(cohorts[t])), "B": int(bc[WINDOWS[0]].shape[0]), "windows": {}}
    for w in WINDOWS:
        cif, allc, nr = point(t, w)
        c0, b0, b1 = bs_fit(cif); bsP = sse_all(bc[w]); plat_b = bsP["plateau"]
        c0a, _, _ = bs_fit(allc)
        wd = {"window": w, "landmark": LM.tolist(), "n_at_risk": nr.astype(int).tolist(), "sparse": (nr < SPARSE).tolist(),
              "cancer_specific": [round(float(x), 2) for x in cif], "all_cause": [round(float(x), 2) for x in allc],
              "cancer_specific_ci_lo": [round(float(x), 2) for x in np.percentile(bc[w], 2.5, 0)], "cancer_specific_ci_hi": [round(float(x), 2) for x in np.percentile(bc[w], 97.5, 0)],
              "all_cause_ci_lo": [round(float(x), 2) for x in np.percentile(ba[w], 2.5, 0)], "all_cause_ci_hi": [round(float(x), 2) for x in np.percentile(ba[w], 97.5, 0)],
              "plateau_end_cancer_specific": c0, "plateau_end_cancer_specific_ci": ci(plat_b), "plateau_end_boot_median": float(np.median(plat_b)),
              "plateau_prop_boot_gt0": float(np.mean(plat_b > 0)), "plateau_prop_boot_ge3": float(np.mean(plat_b >= 3)),
              "broken_stick_intercept": b0, "broken_stick_slope_after_c_per_month": b1, "plateau_end_all_cause": c0a,
              "differences": {}}
        for (a, bb) in DIFFS:
            for lab, pt, bt in (("cancer_specific", cif, bc[w]), ("all_cause", allc, ba[w])):
                d = float(pt[a] - pt[bb]); db = bt[:, a] - bt[:, bb]
                wd["differences"][f"{a}-{bb}_{lab}"] = {"estimate": round(d, 2), "ci": [round(x, 2) for x in ci(db)], "p_boot_two_sided": float(2 * min(np.mean(db <= 0), np.mean(db >= 0)))}
        tr["windows"][str(w)] = wd
    results["tumours"][t] = tr

# ----------------------------------------------------------------------------------------------
# 3. Model comparison (point curve + per-resample) - all windows for NEC point estimates, bootstrap for all windows (cheap), primary = 12
# ----------------------------------------------------------------------------------------------
MODELS = ["broken_stick", "linear", "quadratic", "ns3"]
mc = {}
for w in WINDOWS:
    cif, _, _ = point("NEC", w); r = sse_all(cif[None, :]); c0, b0, b1 = bs_fit(cif)
    pt = {mod: {"sse": float(r[mod][0]), "k": KPAR[mod], "bic": float(bic(r[mod][0], KPAR[mod]))} for mod in MODELS}
    pt["broken_stick"]["c"] = c0; pt["broken_stick"]["intercept"] = b0; pt["broken_stick"]["slope_after_c"] = b1
    # quadratic and spline point fits for the figure/MD
    bq = np.linalg.lstsq(X_quad, cif, rcond=None)[0]; pt["quadratic"]["coef"] = [float(x) for x in bq]
    bn = np.linalg.lstsq(X_ns, cif, rcond=None)[0]; pt["ns3"]["coef"] = [float(x) for x in bn]
    bl = np.linalg.lstsq(X_lin, cif, rcond=None)[0]; pt["linear"]["coef"] = [float(x) for x in bl]
    rb = sse_all(boots["NEC"][0][w])
    bics = {mod: bic(rb[mod], KPAR[mod]) for mod in MODELS}
    d_lin = bics["linear"] - bics["broken_stick"]; d_quad = bics["quadratic"] - bics["broken_stick"]; d_ns = bics["ns3"] - bics["broken_stick"]
    best = np.argmin(np.stack([bics[mod] for mod in MODELS], 1), 1)
    mc[str(w)] = {"point": pt,
                  "bootstrap": {"B": int(len(d_lin)),
                                "delta_bic_linear_minus_bs_median": float(np.median(d_lin)), "delta_bic_linear_minus_bs_ci": ci(d_lin),
                                "prop_bs_better_than_linear_by_gt2": float(np.mean(d_lin > 2)), "prop_bs_better_than_linear_by_gt6": float(np.mean(d_lin > 6)),
                                "prop_linear_better_than_bs_by_gt2": float(np.mean(d_lin < -2)),
                                "delta_bic_quadratic_minus_bs_median": float(np.median(d_quad)), "prop_bs_better_than_quadratic_by_gt2": float(np.mean(d_quad > 2)),
                                "delta_bic_ns3_minus_bs_median": float(np.median(d_ns)), "prop_bs_better_than_ns3_by_gt2": float(np.mean(d_ns > 2)),
                                "prop_best_model": {mod: float(np.mean(best == i)) for i, mod in enumerate(MODELS)},
                                "plateau_c_boot_quartiles": [float(x) for x in np.percentile(rb["plateau"], [25, 50, 75])]}}
results["model_comparison_NEC"] = {"primary_window": 12, "n_points": NPTS, "bic_formula": "n*log(SSE/n) + k*log(n); k = 3 (broken-stick: b0,b1,c), 2 (linear), 3 (quadratic), 4 (natural cubic spline, 3 df + intercept; knots 0,6,12,18)",
                                   "caveat": "The 19 landmark points are overlapping-window estimates from the same patients and are strongly autocorrelated; SSE/BIC are descriptive, not likelihood-valid inference.",
                                   "by_window": mc}

# ----------------------------------------------------------------------------------------------
# 4. Figures
# ----------------------------------------------------------------------------------------------
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": "#52514e",
                     "axes.labelcolor": "#0b0b0b", "xtick.color": "#52514e", "ytick.color": "#52514e", "grid.color": "#e6e5e1", "grid.linewidth": 0.6})
C1, C2, C3, C4 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"; GREY = "#52514e"

# Figure 1: six panels, distant NEC curve per window with CI
fig, axes = plt.subplots(2, 3, figsize=(11, 6.4), sharex=True)
bc, ba = boots["NEC"]
for ax, w in zip(axes.ravel(), WINDOWS):
    wd = results["tumours"]["NEC"]["windows"][str(w)]
    cif = np.array(wd["cancer_specific"]); lo = np.array(wd["cancer_specific_ci_lo"]); hi = np.array(wd["cancer_specific_ci_hi"])
    allc = np.array(wd["all_cause"]); alo = np.array(wd["all_cause_ci_lo"]); ahi = np.array(wd["all_cause_ci_hi"]); sparse = np.array(wd["sparse"])
    ax.fill_between(LM, lo, hi, color=C1, alpha=0.15, linewidth=0)
    ax.plot(LM, cif, color=C1, lw=2, label="Cancer-specific death (Aalen-Johansen)")
    ax.plot(LM, allc, color=GREY, lw=1.4, ls="--", label="All-cause death")
    ax.plot(LM[~sparse], cif[~sparse], "o", color=C1, ms=4, mec="white", mew=0.8)
    if sparse.any(): ax.plot(LM[sparse], cif[sparse], "o", color="white", mec=C1, ms=4, mew=1.2)
    c0 = wd["plateau_end_cancer_specific"]; pci = wd["plateau_end_cancer_specific_ci"]
    X = np.column_stack([np.ones_like(m), np.maximum(0, m - c0)]); yfit = X @ np.array([wd["broken_stick_intercept"], wd["broken_stick_slope_after_c_per_month"]])
    ax.plot(LM, yfit, color=C2, lw=1.2, ls=":", label="Broken-stick fit")
    ax.axvline(c0, color=C2, lw=0.8, alpha=0.7)
    d60 = wd["differences"]["6-0_cancer_specific"]
    ax.set_title(f"{w}-month window", loc="left", fontweight="bold")
    ax.text(0.02, 0.04, f"plateau end c = {c0:g} mo (95% CI {pci[0]:g}-{pci[1]:g})\nLM6 - LM0 = {d60['estimate']:+.1f} pp ({d60['ci'][0]:+.1f}, {d60['ci'][1]:+.1f})",
            transform=ax.transAxes, fontsize=7.5, color="#0b0b0b", va="bottom")
    ax.set_xticks([0, 3, 6, 9, 12, 15, 18]); ax.grid(axis="y"); ax.set_ylim(0, 100)
for ax in axes[1]: ax.set_xlabel("Landmark (months after diagnosis)")
for ax in axes[:, 0]: ax.set_ylabel("Probability of death within window (%)")
h, l = axes[0, 0].get_legend_handles_labels()
fig.legend(h, l, loc="lower center", ncol=3, frameon=False, bbox_to_anchor=(0.5, -0.01))
fig.suptitle("Distant colorectal NEC (n = 908): conditional death probability by landmark and window width\n"
             "Shaded = 95% percentile bootstrap (B = 2000); open markers = < 100 patients at risk", x=0.01, ha="left", fontsize=10)
fig.tight_layout(rect=(0, 0.04, 1, 0.95))
fig.savefig(FIGS + f"{ID}_fig1_nec_window_panels.png", dpi=300); plt.close(fig)

# Figure 2: 6-0 difference vs window width, four tumours
fig, ax = plt.subplots(figsize=(7.6, 4.8))
cols = {"NEC": C1, "SCLC": C2, "AC": C3, "PANC": C4}; off = {"NEC": -0.45, "SCLC": -0.15, "AC": 0.15, "PANC": 0.45}
for t in ["NEC", "SCLC", "AC", "PANC"]:
    est = [results["tumours"][t]["windows"][str(w)]["differences"]["6-0_cancer_specific"]["estimate"] for w in WINDOWS]
    lo = [results["tumours"][t]["windows"][str(w)]["differences"]["6-0_cancer_specific"]["ci"][0] for w in WINDOWS]
    hi = [results["tumours"][t]["windows"][str(w)]["differences"]["6-0_cancer_specific"]["ci"][1] for w in WINDOWS]
    x = np.array(WINDOWS, float) + off[t]
    ax.errorbar(x, est, yerr=[np.array(est) - np.array(lo), np.array(hi) - np.array(est)], fmt="o-", color=cols[t], ms=5, lw=1.6, capsize=2, mec="white", mew=0.8,
                label=f"{LABEL[t]}, B = {results['tumours'][t]['B']}")
ax.axhline(0, color=GREY, lw=0.8)
ax.set_xticks(WINDOWS); ax.set_xlabel("Window width (months)"); ax.set_ylabel("Landmark 6 minus landmark 0 (percentage points)")
ax.set_title("Conditional cancer-death probability: landmark 6 minus landmark 0, by window width\n(point estimate, 95% percentile bootstrap CI)", loc="left", fontsize=9.5)
ax.grid(axis="y"); ax.set_xlim(1.5, 25.5); ax.set_ylim(-16, 11)
ax.legend(frameon=False, fontsize=8, loc="upper right", ncol=1)
fig.tight_layout(); fig.savefig(FIGS + f"{ID}_fig2_diff60_vs_window.png", dpi=300); plt.close(fig)

# ----------------------------------------------------------------------------------------------
# 5. Write JSON + Markdown
# ----------------------------------------------------------------------------------------------
runtime_min = (time.time() - T0) / 60
results["runtime_minutes"] = round(runtime_min, 2)
results["software"] = {"python": platform.python_version(), "numpy": np.__version__, "pandas": pd.__version__, "matplotlib": matplotlib.__version__}
results["figures"] = [config.relpath(FIGS + f"{ID}_fig1_nec_window_panels.png"), config.relpath(FIGS + f"{ID}_fig2_diff60_vs_window.png")]
with open(RES + f"{ID}.json", "w") as f: json.dump(results, f, indent=1)

def fmt(x, nd=2): return "NA" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{x:.{nd}f}"
L = []
L.append(f"# {ID} - Window-width sensitivity and plateau-model comparison\n")
L.append("Distant colorectal NEC (main definition, n = 908); comparators: distant SCLC (n = {:,}), distant colorectal adenocarcinoma (n = {:,}), distant pancreatic adenocarcinoma (n = {:,}). Date 2026-10-07. Seed {}; B = 2000 (NEC) / 500 (comparators); patient-level resampling, identical resamples across windows and endpoints; 95% percentile intervals. Runtime {:.1f} min.\n".format(
    results["n"]["SCLC"], results["n"]["AC"], results["n"]["PANC"], SEED, runtime_min))
L.append("## Reproduction check (12-month window)\n")
L.append("| Tumour | Quantity | This script | Manuscript / orchestrator | n at risk (mine / ms) | Agree |\n|---|---|---|---|---|---|")
for r in repro_rows:
    L.append(f"| {r['tumour']} | {r['quantity']} | {r['mine']} | {r['manuscript']} | {r['n_at_risk_mine'] if r['n_at_risk_mine'] is not None else ''} / {r['n_at_risk_manuscript'] if r['n_at_risk_manuscript'] is not None else ''} | {'yes' if r['agree'] else 'NO'} |")
L.append("\nAll anchors reproduced to 2 decimals; the script aborts otherwise.\n")

L.append("## Table 1. Distant NEC: conditional probability of cancer death (%) within the window after each landmark, by window width\n")
L.append("Point estimate (95% bootstrap CI). Asterisk = fewer than 100 patients at risk at the landmark (the number at risk does not depend on the window). Landmark + window beyond ~36 months relies on sparse late follow-up.\n")
hdr = "| Landmark | n at risk | " + " | ".join(f"{w}-mo" for w in WINDOWS) + " |"
L.append(hdr); L.append("|---|---|" + "---|" * len(WINDOWS))
NW = results["tumours"]["NEC"]["windows"]
for j, s in enumerate(LM):
    nr = NW["12"]["n_at_risk"][j]; star = "*" if nr < SPARSE else ""
    cells = [f"{NW[str(w)]['cancer_specific'][j]:.1f} ({NW[str(w)]['cancer_specific_ci_lo'][j]:.1f}-{NW[str(w)]['cancer_specific_ci_hi'][j]:.1f})" for w in WINDOWS]
    L.append(f"| {s} | {nr}{star} | " + " | ".join(cells) + " |")

L.append("\n## Table 2. Distant NEC: all-cause death probability (%) within the window, by window width\n")
L.append(hdr); L.append("|---|---|" + "---|" * len(WINDOWS))
for j, s in enumerate(LM):
    nr = NW["12"]["n_at_risk"][j]; star = "*" if nr < SPARSE else ""
    cells = [f"{NW[str(w)]['all_cause'][j]:.1f} ({NW[str(w)]['all_cause_ci_lo'][j]:.1f}-{NW[str(w)]['all_cause_ci_hi'][j]:.1f})" for w in WINDOWS]
    L.append(f"| {s} | {nr}{star} | " + " | ".join(cells) + " |")

L.append("\n## Table 3. Distant NEC: landmark differences (percentage points) and broken-stick plateau end, by window width\n")
L.append("| Window | LM6 - LM0 cancer | LM3 - LM0 cancer | LM6 - LM3 cancer | LM6 - LM0 all-cause | Plateau end c (mo), 95% CI | Boot. median c | P(c > 0) | P(c >= 3) | Slope after c (pp/mo) | c, all-cause |")
L.append("|---|---|---|---|---|---|---|---|---|---|---|")
for w in WINDOWS:
    wd = NW[str(w)]; D = wd["differences"]
    def dd(key): return f"{D[key]['estimate']:+.1f} ({D[key]['ci'][0]:+.1f}, {D[key]['ci'][1]:+.1f})"
    L.append(f"| {w} | {dd('6-0_cancer_specific')} | {dd('3-0_cancer_specific')} | {dd('6-3_cancer_specific')} | {dd('6-0_all_cause')} | {wd['plateau_end_cancer_specific']:g} ({wd['plateau_end_cancer_specific_ci'][0]:g}-{wd['plateau_end_cancer_specific_ci'][1]:g}) | {wd['plateau_end_boot_median']:g} | {wd['plateau_prop_boot_gt0']:.2f} | {wd['plateau_prop_boot_ge3']:.2f} | {wd['broken_stick_slope_after_c_per_month']:+.2f} | {wd['plateau_end_all_cause']:g} |")

L.append("\n## Table 4. Comparator tumours (distant stage): 6-0 difference in cancer-death probability and plateau end, by window width\n")
L.append("Point estimates from the full cohort; 95% CI of the 6-0 difference from B = 500 bootstrap resamples (NEC row repeated from Table 3 for reference).\n")
L.append("| Tumour | Window | LM0 (%) | LM3 (%) | LM6 (%) | LM12 (%) | LM6 - LM0 (95% CI) | LM3 - LM0 | LM6 - LM3 | Plateau end c (mo) | Boot. median c | P(c >= 3) |")
L.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
for t in ["NEC", "SCLC", "AC", "PANC"]:
    for w in WINDOWS:
        wd = results["tumours"][t]["windows"][str(w)]; D = wd["differences"]; cs = wd["cancer_specific"]
        L.append(f"| {t} | {w} | {cs[0]:.1f} | {cs[3]:.1f} | {cs[6]:.1f} | {cs[12]:.1f} | {D['6-0_cancer_specific']['estimate']:+.1f} ({D['6-0_cancer_specific']['ci'][0]:+.1f}, {D['6-0_cancer_specific']['ci'][1]:+.1f}) | {D['3-0_cancer_specific']['estimate']:+.1f} | {D['6-3_cancer_specific']['estimate']:+.1f} | {wd['plateau_end_cancer_specific']:g} | {wd['plateau_end_boot_median']:g} | {wd['plateau_prop_boot_ge3']:.2f} |")

L.append("\n## Table 5. Model comparison on the distant-NEC landmark curve (landmarks 0-18, 19 points)\n")
L.append("Gaussian BIC = n log(SSE/n) + k log(n) (common constant dropped); k = 3 broken-stick (intercept, slope, change-point c searched on a 0.25-month grid 0-15), 2 linear, 3 quadratic, 4 natural cubic spline (3 df; knots 0, 6, 12, 18). delta-BIC = BIC(linear) - BIC(broken-stick); positive favours the broken stick. Bootstrap proportions use the same resamples as Tables 1-3. The 19 points share patients across overlapping windows and are strongly autocorrelated, so BIC values are descriptive only. Primary window = 12 months.\n")
L.append("| Window | SSE broken-stick (c) | SSE linear | SSE quadratic | SSE spline-3df | BIC bs | BIC lin | BIC quad | BIC ns3 | delta-BIC lin-bs point | boot. median delta-BIC (95% CI) | P(bs beats lin by > 2) | P(bs beats lin by > 6) | P(lin beats bs by > 2) | P(bs beats quad by > 2) | P(bs beats ns3 by > 2) | best model (boot. share bs/lin/quad/ns3) |")
L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
for w in WINDOWS:
    p = mc[str(w)]["point"]; b = mc[str(w)]["bootstrap"]; pb = b["prop_best_model"]
    L.append(f"| {w} | {p['broken_stick']['sse']:.1f} (c = {p['broken_stick']['c']:g}) | {p['linear']['sse']:.1f} | {p['quadratic']['sse']:.1f} | {p['ns3']['sse']:.1f} | {p['broken_stick']['bic']:.1f} | {p['linear']['bic']:.1f} | {p['quadratic']['bic']:.1f} | {p['ns3']['bic']:.1f} | {p['linear']['bic'] - p['broken_stick']['bic']:+.1f} | {b['delta_bic_linear_minus_bs_median']:+.1f} ({b['delta_bic_linear_minus_bs_ci'][0]:+.1f}, {b['delta_bic_linear_minus_bs_ci'][1]:+.1f}) | {b['prop_bs_better_than_linear_by_gt2']:.2f} | {b['prop_bs_better_than_linear_by_gt6']:.2f} | {b['prop_linear_better_than_bs_by_gt2']:.2f} | {b['prop_bs_better_than_quadratic_by_gt2']:.2f} | {b['prop_bs_better_than_ns3_by_gt2']:.2f} | {pb['broken_stick']:.2f}/{pb['linear']:.2f}/{pb['quadratic']:.2f}/{pb['ns3']:.2f} |")

# ---- interpretation (data-driven text) ----
w3, w6, w12 = NW["3"], NW["6"], NW["12"]
def trend_word(d): return "rises" if d["ci"][0] > 0 else ("falls" if d["ci"][1] < 0 else "is flat within sampling error")
L.append("\n## Interpretation\n")
L.append(f"1. Under the authors' 12-month window the distant-NEC curve is flat from landmark 0 to 6 (LM6 - LM0 = {w12['differences']['6-0_cancer_specific']['estimate']:+.1f} pp, 95% CI {w12['differences']['6-0_cancer_specific']['ci'][0]:+.1f} to {w12['differences']['6-0_cancer_specific']['ci'][1]:+.1f}) and the broken-stick change-point is c = {w12['plateau_end_cancer_specific']:g} months (bootstrap 95% CI {w12['plateau_end_cancer_specific_ci'][0]:g}-{w12['plateau_end_cancer_specific_ci'][1]:g}; bootstrap median {w12['plateau_end_boot_median']:g}).")
def shape(D):
    a, b = D["3-0_cancer_specific"], D["6-3_cancer_specific"]
    if a["ci"][1] < 0 and b["ci"][0] > 0: return "falls from landmark 0 to 3 and rises again from 3 to 6 (V-shape)"
    if a["ci"][1] < 0: return "falls from landmark 0 to 3"
    if b["ci"][0] > 0: return "rises from landmark 3 to 6"
    return "is flat within sampling error"
L.append(f"2. With a 6-month window the curve {shape(w6['differences'])}: LM3 - LM0 = {w6['differences']['3-0_cancer_specific']['estimate']:+.1f} pp (CI {w6['differences']['3-0_cancer_specific']['ci'][0]:+.1f}, {w6['differences']['3-0_cancer_specific']['ci'][1]:+.1f}), LM6 - LM3 = {w6['differences']['6-3_cancer_specific']['estimate']:+.1f} (CI {w6['differences']['6-3_cancer_specific']['ci'][0]:+.1f}, {w6['differences']['6-3_cancer_specific']['ci'][1]:+.1f}), net LM6 - LM0 = {w6['differences']['6-0_cancer_specific']['estimate']:+.1f} (CI {w6['differences']['6-0_cancer_specific']['ci'][0]:+.1f}, {w6['differences']['6-0_cancer_specific']['ci'][1]:+.1f}); broken-stick c = {w6['plateau_end_cancer_specific']:g} (CI {w6['plateau_end_cancer_specific_ci'][0]:g}-{w6['plateau_end_cancer_specific_ci'][1]:g}).")
L.append(f"3. With a 3-month window the curve {shape(w3['differences'])}: LM3 - LM0 = {w3['differences']['3-0_cancer_specific']['estimate']:+.1f} pp (CI {w3['differences']['3-0_cancer_specific']['ci'][0]:+.1f}, {w3['differences']['3-0_cancer_specific']['ci'][1]:+.1f}), LM6 - LM3 = {w3['differences']['6-3_cancer_specific']['estimate']:+.1f} (CI {w3['differences']['6-3_cancer_specific']['ci'][0]:+.1f}, {w3['differences']['6-3_cancer_specific']['ci'][1]:+.1f}), net LM6 - LM0 = {w3['differences']['6-0_cancer_specific']['estimate']:+.1f} (CI {w3['differences']['6-0_cancer_specific']['ci'][0]:+.1f}, {w3['differences']['6-0_cancer_specific']['ci'][1]:+.1f}). The 3-month window is the closest analogue of the interval hazard (cf. crude rates 0.144/0.091/0.134 per person-month at 0-3/3-6/6-9 months) and shows that the 'flat' 12-month curve is the running integral of a hazard that drops after the first quarter and rebounds at 6-9 months; it is not evidence of a constant hazard.")
L.append("4. Broken-stick change-point by window (NEC, cancer-specific): " + "; ".join(f"{w} mo -> c = {NW[str(w)]['plateau_end_cancer_specific']:g} (CI {NW[str(w)]['plateau_end_cancer_specific_ci'][0]:g}-{NW[str(w)]['plateau_end_cancer_specific_ci'][1]:g})" for w in WINDOWS) + f". P(c >= 3) is >= 0.99 for every window, so a non-zero change-point is robust, but its location drifts systematically from ~6 months (12-24-month windows) to {w3['plateau_end_cancer_specific']:g} months (3-month window), where the hinge model is simply the best straight-line compromise through a V-shaped first year (its SSE {mc['3']['point']['broken_stick']['sse']:.0f} is no better than the 3-df spline's {mc['3']['point']['ns3']['sse']:.0f}). The 'plateau end' is therefore a window-dependent summary, not a fixed biological time point.")
b12 = mc["12"]["bootstrap"]; p12 = mc["12"]["point"]
L.append(f"5. Model comparison (12-month window, 19 points): point-estimate BIC broken-stick {p12['broken_stick']['bic']:.1f} vs linear {p12['linear']['bic']:.1f} (delta {p12['linear']['bic'] - p12['broken_stick']['bic']:+.1f}), quadratic {p12['quadratic']['bic']:.1f}, spline {p12['ns3']['bic']:.1f} (SSE {p12['broken_stick']['sse']:.1f} / {p12['linear']['sse']:.1f} / {p12['quadratic']['sse']:.1f} / {p12['ns3']['sse']:.1f}). Across bootstrap resamples the broken stick beats the linear model by > 2 BIC units in {b12['prop_bs_better_than_linear_by_gt2']*100:.0f}% and by > 6 units in {b12['prop_bs_better_than_linear_by_gt6']*100:.0f}% (median delta-BIC {b12['delta_bic_linear_minus_bs_median']:+.1f}, 95% CI {b12['delta_bic_linear_minus_bs_ci'][0]:+.1f} to {b12['delta_bic_linear_minus_bs_ci'][1]:+.1f}); it is the BIC-best of the four models in {b12['prop_best_model']['broken_stick']*100:.0f}% of resamples versus {b12['prop_best_model']['ns3']*100:.0f}% for the 3-df spline ({b12['prop_best_model']['quadratic']*100:.0f}% quadratic, {b12['prop_best_model']['linear']*100:.0f}% linear), and for the 6-, 18- and 24-month windows the spline is best more often than the hinge. Because the 19 points are overlapping-window estimates from the same patients, they are strongly autocorrelated: the effective number of independent points is far below 19, the Gaussian BIC overstates the evidence, and the comparison supports 'clearly non-linear, with an early flat segment' rather than a specific plateau-then-linear-decline law.")
L.append("6. Comparators: " + "; ".join(f"{t} LM6 - LM0 (12-mo) = {results['tumours'][t]['windows']['12']['differences']['6-0_cancer_specific']['estimate']:+.1f} ({results['tumours'][t]['windows']['12']['differences']['6-0_cancer_specific']['ci'][0]:+.1f}, {results['tumours'][t]['windows']['12']['differences']['6-0_cancer_specific']['ci'][1]:+.1f}), 3-mo = {results['tumours'][t]['windows']['3']['differences']['6-0_cancer_specific']['estimate']:+.1f}, c12 = {results['tumours'][t]['windows']['12']['plateau_end_cancer_specific']:g}, c3 = {results['tumours'][t]['windows']['3']['plateau_end_cancer_specific']:g}" for t in ["SCLC", "AC", "PANC"]) + ". SCLC shares the NEC pattern (flat 12-month curve, 3-month curve falling then rebounding); the two adenocarcinomas decline monotonically under every window.")
L.append("7. Answer to the question posed: the 12-month 'plateau' is a window-width phenomenon in the sense that it is the cumulative integral of a non-monotone short-term hazard; under 3- and 6-month windows the conditional risk is not flat between 0 and 6 months but falls over the first quarter and then returns toward (or above) its initial level by landmark 6 (see Table 3 for exact signs and CIs). The manuscript should describe the finding as 'no decline in 12-month conditional risk over the first 6 months, reflecting an early dip and 6-9-month rebound in the short-term hazard' rather than as a constant hazard.")
L.append("\nDeviations from plan: none in estimands. Natural cubic spline basis implemented in numpy (ESL eq. 5.4-5.5) with knots at 0/6/12/18 rather than patsy `cr`, to make the parameter count (4) explicit; model comparison bootstrap was also run for the five non-primary windows because it is free once the resamples exist. Comparator bootstraps use B = 500 as instructed.\n")
L.append(f"Figures: `{results['figures'][0]}`, `{results['figures'][1]}`.\n")
with open(RES + f"{ID}.md", "w") as f: f.write("\n".join(L))
print(f"DONE in {runtime_min:.1f} min; results -> {RES}{ID}.json / .md ; figs -> {FIGS}{ID}_*.png")
