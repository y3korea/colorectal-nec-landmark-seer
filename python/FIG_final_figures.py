# -*- coding: utf-8 -*-
"""FIG: publication figures 1-6 for the colorectal NEC landmark manuscript.

Reads ONLY the JSON/CSV outputs in results/ (A01, A02, A03, A05, A05_std, A06, A07, A09 grid, A11,
FIG_flow_counts, data/flow_all.csv). No patient-level data are read here. Figures: 180 mm wide, Arial,
300 dpi PNG + PDF, panel labels A/B/C, Okabe-Ito colour-blind-safe palette (validated), risk tables under
landmark curves, open markers where < 100 patients at risk.
Imported as a helper module by FIG_manuscript_final.py (shared palette, JSON loaders, flow rows).
Running it directly re-draws the earlier (pre-IJC) figure set into <NEC_FIG_DIR>/final/: python3 python/FIG_final_figures.py
"""
import json, os, time, csv
import config
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.stats import chi2

T0 = time.time()
BASE = str(config.REPO)
RES = config.RES
OUT = config.FIGS + "final/"

def J(name):
    with open(RES + name) as f: return json.load(f)

A01 = J("A01.json"); A02 = J("A02.json"); A03 = J("A03.json"); A05 = J("A05.json"); A05S = J("A05_std.json")
A06 = J("A06.json"); A07 = J("A07.json"); A11 = J("A11.json"); FLOWC = J("FIG_flow_counts.json")
flow_rows = list(csv.DictReader(open(config.DATA + "flow_all.csv")))

# ---------------------------------------------------------------- style
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7, "axes.labelsize": 7.5, "axes.titlesize": 8, "legend.fontsize": 6.5,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
    "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5, "lines.linewidth": 1.2,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "pdf.fonttype": 42,
    "ps.fonttype": 42, "savefig.dpi": 300, "figure.dpi": 100,
})
MM = 1 / 25.4
W2 = 180 * MM          # 2-column width
# Okabe-Ito (validated with dataviz validator: all checks pass; CVD separation in the 6-8 band for
# PANC/AC -> secondary encoding by marker shape, line style and direct labels)
COL = {"NEC": "#0072B2", "SCLC": "#D55E00", "AC": "#009E73", "PANC": "#CC79A7", "LCNEC": "#E69F00",
       "grey": "#6E6E6E", "light": "#BBBBBB", "black": "#222222"}
MK = {"NEC": "o", "SCLC": "s", "AC": "^", "PANC": "D", "LCNEC": "v"}
LS = {"NEC": "-", "SCLC": "-", "AC": "--", "PANC": "-.", "LCNEC": ":"}
LAB = {"NEC": "Colorectal NEC", "SCLC": "Small-cell lung cancer", "AC": "Colorectal adenocarcinoma",
       "PANC": "Pancreatic adenocarcinoma", "LCNEC": "Lung LCNEC"}
STAGE_COL = {"Distant": "#0072B2", "Regional": "#D55E00", "Localized": "#009E73"}
STAGE_MK = {"Distant": "o", "Regional": "s", "Localized": "^"}
MANIFEST = {"id": "FIG", "figures": {}, "sources": ["A01", "A02", "A03", "A05", "A05_std", "A06", "A07",
            "A09_hazard_grid", "A11", "FIG_flow_counts", "flow_all.csv"]}

def panel_label(ax, s, dx=-0.12, dy=1.04):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=10, fontweight="bold", va="bottom", ha="left")

def save(fig, name):
    fig.savefig(OUT + name + ".png", dpi=300, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(OUT + name + ".pdf", bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)
    print("saved", name)

def curve_with_markers(ax, x, y, n, color, marker, ls="-", label=None, lo=None, hi=None, every=1, ms=3.2,
                       thresh=100, zorder=3, lw=1.2, band_alpha=0.15):
    """line + CI band + markers at every `every` landmark; open marker when n < thresh."""
    x = np.asarray(x, float); y = np.asarray(y, float); n = np.asarray(n, float)
    if lo is not None:
        ax.fill_between(x, lo, hi, color=color, alpha=band_alpha, lw=0, zorder=zorder - 1)
    ax.plot(x, y, color=color, ls=ls, lw=lw, label=label, zorder=zorder)
    sel = np.arange(0, len(x), every)
    full = sel[n[sel] >= thresh]; opn = sel[n[sel] < thresh]
    ax.plot(x[full], y[full], marker=marker, ls="none", ms=ms, mfc=color, mec=color, mew=0.7, zorder=zorder + 1)
    ax.plot(x[opn], y[opn], marker=marker, ls="none", ms=ms, mfc="white", mec=color, mew=0.8, zorder=zorder + 1)

def risk_table(ax_tab, ax_main, rows, xs, title="No. at risk"):
    """rows: list of (label, color, dict landmark->n). ax_tab shares x with ax_main."""
    ax_tab.set_xlim(ax_main.get_xlim()); ax_tab.set_ylim(0, len(rows) + 0.6)
    ax_tab.axis("off")
    ax_tab.text(ax_main.get_xlim()[0], len(rows) + 0.35, title, fontsize=6, ha="left", va="center",
                color=COL["black"], fontweight="bold")
    for i, (lab, col, d) in enumerate(rows):
        yy = len(rows) - 1 - i + 0.5
        ax_tab.text(ax_main.get_xlim()[0] - 0.045 * (ax_main.get_xlim()[1] - ax_main.get_xlim()[0]), yy, lab,
                    fontsize=6, ha="right", va="center", color=col, clip_on=False)
        for x in xs:
            if x in d and d[x] is not None:
                n = d[x]
                ax_tab.text(x, yy, f"{int(n):,}", fontsize=5.6, ha="center", va="center", color=COL["black"],
                            fontweight="normal" if n >= 100 else "light", style="normal" if n >= 100 else "italic")

def fmt(e, lo, hi, d=1):
    return f"{e:+.{d}f} ({lo:+.{d}f} to {hi:+.{d}f})"

# =========================================================================================== Fig 1 flow
def fig1():
    fl = {}
    for r in flow_rows:
        fl.setdefault(r["cohort"], {})[r["step"]] = (int(r["records"]), int(r["patients"]))
    nen = fl["NEN_all"]; lung = fl["LUNG"]; ac = fl["AC"]; panc = fl["PANC"]
    fig = plt.figure(figsize=(W2, 150 * MM))
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 180); ax.set_ylim(0, 150); ax.axis("off")

    def box(x, y, w, h, text, fc="white", ec=COL["black"], lw=0.7, fs=6.5, bold=False, color=COL["black"]):
        p = FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.2", fc=fc, ec=ec, lw=lw)
        ax.add_patch(p)
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=color,
                fontweight="bold" if bold else "normal", linespacing=1.25)

    def arrow(x1, y1, x2, y2, lw=0.8):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=7, lw=lw,
                                     color=COL["black"], shrinkA=0, shrinkB=0))

    # ---- left: colorectal NEN -> NEC
    X0, BW, BH = 4, 62, 9.5
    y = 138
    box(X0, y, BW, BH, f"SEER-17 colorectal neuroendocrine neoplasm records\n(ICD-O-3 8013, 8041-8045, 8240-8246, 8249)\n"
        f"n = {nen['exported'][0]:,} records ({nen['exported'][1]:,} patients)", fc="#EEF3F8", bold=True)
    steps = [("age>=18", "Age < 18 years"), ("site", "Primary site not colon/rectum (C18-C20)"),
             ("diag_confirmation", "No positive histology"), ("intl_primary", "Not a malignant primary (intl. rules)"),
             ("not_autopsy_dco", "Autopsy / death-certificate only"), ("one_per_patient", "Further record of the same patient")]
    prev_key = "malignant"; prev = nen[prev_key][0]
    ytop = y
    for key, lab in steps:
        n_now = nen[key][0]
        excl = prev - n_now
        if key == "one_per_patient":
            excl = nen["survival_known"][0] - nen["one_per_patient"][1]; n_now = nen["one_per_patient"][1]
        ynext = ytop - 15
        arrow(X0 + BW / 2, ytop, X0 + BW / 2, ynext + BH)
        box(X0 + BW + 4, ynext + BH + 1.2, 50, 6.2, f"Excluded: {lab}\nn = {excl:,}", fc="#FAFAFA", ec=COL["grey"], lw=0.5, fs=5.3, color=COL["grey"])
        ax.plot([X0 + BW / 2, X0 + BW + 4], [ynext + BH + 4.3, ynext + BH + 4.3], color=COL["grey"], lw=0.5)
        box(X0, ynext, BW, BH, f"n = {n_now:,}" + (" patients" if key == "one_per_patient" else " records"), fs=6.5)
        ytop = ynext; prev = n_now
    # NEC histology
    ynext = ytop - 15
    arrow(X0 + BW / 2, ytop, X0 + BW / 2, ynext + BH)
    box(X0 + BW + 4, ynext + BH + 1.2, 50, 6.2, f"Excluded: other NEN histology (8240-8246, 8249)\nn = {FLOWC['NEN_all_n'] - FLOWC['NEC_main_n']:,}", fc="#FAFAFA", ec=COL["grey"], lw=0.5, fs=5.3, color=COL["grey"])
    ax.plot([X0 + BW / 2, X0 + BW + 4], [ynext + BH + 4.3, ynext + BH + 4.3], color=COL["grey"], lw=0.5)
    box(X0, ynext - 1, BW, BH + 3, f"Main colorectal NEC cohort (8013, 8041-8045)\nn = {FLOWC['NEC_main_n']:,} patients\nsmall cell 8041-8045: {sum(v for k, v in FLOWC['NEC_hist'].items() if k != '8013'):,}; large cell 8013: {FLOWC['NEC_hist']['8013']:,}",
        fc="#DCE9F5", bold=True, ec=COL["NEC"], lw=1.0, fs=6.2)
    ytop = ynext - 1
    # stage split
    st = FLOWC["NEC_stage"]; yst = ytop - 16
    xs_st = [X0 - 1, X0 + 23, X0 + 47, X0 + 71]
    for xx, (lab, key, fc, ec) in zip(xs_st, [("Distant\n(primary analysis)", "Distant", "#C9DCF0", COL["NEC"]), ("Regional", "Regional", "#F7E3D4", COL["SCLC"]),
                                           ("Localized", "Localized", "#D9F0E8", COL["AC"]), ("Unknown stage", "Unknown", "#F2F2F2", COL["grey"])]):
        ax.plot([X0 + BW / 2, xx + 11], [ytop, yst + 8.5], color=COL["black"], lw=0.6)
        box(xx, yst, 22, 8.5, f"{lab}\nn = {st[key]:,}", fc=fc, ec=ec, lw=0.9 if key == "Distant" else 0.6, bold=(key == "Distant"), fs=6)
    # expanded
    box(X0 - 1, yst - 12.5, 94, 8.5, f"Expanded NEC (sensitivity): + 8246 NOS with registry grade 3/4 (n = {FLOWC['n_8246_g34']:,} of {FLOWC['n_8246_total']:,})\n"
        f"n = {FLOWC['NEC_expanded_n']:,} patients; distant n = {FLOWC['NEC_expanded_distant']:,}", fc="#FAFAFA", ec=COL["grey"], lw=0.6, fs=5.8)

    # ---- right: comparator cohorts (same filters)
    XR, BWR = 124, 52
    ax.text(XR + BWR / 2, 146.5, "Comparator cohorts (same SEER filters, distant stage)", ha="center", va="center", fontsize=7, fontweight="bold")
    yy = 131
    box(XR, yy, BWR, 8.5, f"Lung small-cell family exported\nn = {lung['exported'][0]:,} records -> {lung['one_per_patient'][1]:,} patients", fc="#FAFAFA", ec=COL["grey"], lw=0.6, fs=5.6)
    hw = BWR / 2 - 1
    arrow(XR + hw / 2, yy, XR + hw / 2, yy - 4); arrow(XR + BWR - hw / 2, yy, XR + BWR - hw / 2, yy - 4)
    yy -= 4 + 10
    box(XR, yy, hw, 10, f"SCLC (8041-8045)\nn = {FLOWC['SCLC_n']:,}\nDistant: n = {FLOWC['SCLC_distant']:,}", fc="white", ec=COL["SCLC"], lw=1.0, fs=5.8, bold=True, color=COL["SCLC"])
    box(XR + BWR - hw, yy, hw, 10, f"Lung LCNEC (8013)\nn = {FLOWC['LCNEC_n']:,}\nDistant: n = {FLOWC['LCNEC_distant']:,}", fc="white", ec=COL["LCNEC"], lw=1.0, fs=5.8, bold=True, color=COL["LCNEC"])
    yy -= 14
    for key, top, bot, col, h in [
        ("AC", f"Colorectal adenocarcinoma exported\nn = {ac['exported'][0]:,} records -> {ac['one_per_patient'][1]:,} patients", f"Distant: n = {FLOWC['AC_distant']:,}", COL["AC"], 8.5),
        ("PANC", f"Pancreatic adenocarcinoma exported\nn = {panc['exported'][0]:,} records -> {panc['one_per_patient'][1]:,} patients\nexcl. islet-cell 8150-8157: n = {FLOWC['PANC_excl_islet_n']:,}", f"Distant: n = {FLOWC['PANC_distant']:,}", COL["PANC"], 11)]:
        box(XR, yy, BWR, h, top, fc="#FAFAFA", ec=COL["grey"], lw=0.6, fs=5.6)
        arrow(XR + BWR / 2, yy, XR + BWR / 2, yy - 4)
        yy -= 4 + 8.5
        box(XR, yy, BWR, 8.5, bot, fc="white", ec=col, lw=1.0, fs=6, bold=True, color=col)
        yy -= 14
    ax.text(XR + BWR / 2, 3, "Common filters: age >= 18 y, malignant, positive histology,\nfirst/only malignant primary per international rules,\nnot autopsy/DCO, diagnosed 2000-2023, known survival time,\none record per patient", ha="center", va="bottom", fontsize=5.4, color=COL["grey"], linespacing=1.3)
    save(fig, "Fig1_flow")
    MANIFEST["figures"]["Fig1_flow"] = {"NEN_flow": {k: v for k, v in nen.items()}, "counts": FLOWC}

# ================================================================================= Fig 2 landmark by stage
def fig2():
    fig = plt.figure(figsize=(W2, 95 * MM))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 0.3], hspace=0.42, wspace=0.28, left=0.07, right=0.99, top=0.93, bottom=0.03)
    axA = fig.add_subplot(gs[0, 0]); axAt = fig.add_subplot(gs[1, 0], sharex=axA)
    axB = fig.add_subplot(gs[0, 1]); axBt = fig.add_subplot(gs[1, 1], sharex=axB)
    rows = []
    for key, name in [("monthly_distant", "Distant"), ("monthly_regional", "Regional"), ("monthly_localized", "Localized")]:
        m = A01[key]; x = m["landmarks"]
        curve_with_markers(axA, x, m["cancer_pct"], m["n_at_risk"], STAGE_COL[name], STAGE_MK[name], lo=m["cancer_ci_lo"], hi=m["cancer_ci_hi"],
                           label=f"{name} (n = {m['n']:,})", every=1, ms=2.6)
        rows.append((name, STAGE_COL[name], {x[i]: m["n_at_risk"][i] for i in range(len(x))}))
    axA.set_xlim(-0.5, 24.5); axA.set_ylim(0, 85); axA.set_xticks(range(0, 25, 3))
    axA.set_ylabel("12-month cancer-death probability (%)"); axA.set_xlabel("Landmark (months after diagnosis)")
    axA.legend(loc="upper right", bbox_to_anchor=(1.0, 1.0), handlelength=2.2)
    axA.text(0.02, 0.04, "Filled markers: >= 100 at risk; open: < 100", transform=axA.transAxes, fontsize=5.6, color=COL["grey"])
    risk_table(axAt, axA, rows, list(range(0, 25, 3)))
    panel_label(axA, "A")
    # ---- B
    m = A01["monthly_distant"]; x = np.array(m["landmarks"]); pl = m["plateau_cancer"]
    axB.axvspan(pl["c_ci"][0], pl["c_ci"][1], color=COL["NEC"], alpha=0.08, lw=0, zorder=0)
    axB.axvline(pl["c"], color=COL["NEC"], lw=0.8, ls=":", zorder=1)
    curve_with_markers(axB, x, m["cancer_pct"], m["n_at_risk"], COL["NEC"], "o", lo=m["cancer_ci_lo"], hi=m["cancer_ci_hi"],
                       label="Unadjusted (Aalen-Johansen), 95% CI", ms=2.8)
    s = A05S["original"]
    axB.plot(s["landmarks"], s["std_curve"], color=COL["black"], ls="--", lw=1.1, marker="x", ms=3, mew=0.8,
             label="Directly standardized (age, sex, year,\nlarge-cell histology, rectal site)", zorder=4)
    # broken-stick fit
    xx = np.linspace(0, 18, 200); yy = pl["b0"] + pl["b1_per_month"] * np.maximum(0, xx - pl["c"])
    axB.plot(xx, yy, color=COL["SCLC"], lw=1.0, ls="-", alpha=0.9, label=f"Broken-stick fit (landmarks 0-18)", zorder=2)
    axB.annotate(f"Plateau end {pl['c']:.2f} mo (95% CI {pl['c_ci'][0]:.1f}-{pl['c_ci'][1]:.2f})\nslope after plateau {pl['b1_per_month']:.2f} pp/month\n(95% CI {pl['b1_ci'][0]:.2f} to {pl['b1_ci'][1]:.2f})", xy=(pl["c"], 34), xytext=(10.5, 20),
                 fontsize=5.8, ha="left", arrowprops=dict(arrowstyle="-", lw=0.5, color=COL["grey"]), color=COL["black"], linespacing=1.25)
    axB.set_xlim(-0.5, 24.5); axB.set_ylim(0, 85); axB.set_xticks(range(0, 25, 3))
    axB.set_ylabel("12-month cancer-death probability (%)"); axB.set_xlabel("Landmark (months after diagnosis)")
    axB.legend(loc="lower left", bbox_to_anchor=(0.0, 0.0), handlelength=2.4, fontsize=5.8)
    risk_table(axBt, axB, [("Distant NEC", COL["NEC"], {x[i]: m["n_at_risk"][i] for i in range(len(x))})], list(range(0, 25, 3)))
    panel_label(axB, "B")
    save(fig, "Fig2_landmark_stage")
    MANIFEST["figures"]["Fig2_landmark_stage"] = {
        "A": {k: {"landmarks": A01[k]["landmarks"], "cancer_pct": A01[k]["cancer_pct"], "ci_lo": A01[k]["cancer_ci_lo"], "ci_hi": A01[k]["cancer_ci_hi"],
                  "n_at_risk": A01[k]["n_at_risk"], "n": A01[k]["n"]} for k in ["monthly_distant", "monthly_regional", "monthly_localized"]},
        "B": {"plateau": pl, "standardized": s, "B": A01["B"], "seed": A01["seeds"]["other"]}}

# ======================================================================================== Fig 3 hazard
def poisson_exact(d, pm):
    lo = chi2.ppf(0.025, 2 * d) / 2 / pm if d > 0 else 0.0
    hi = chi2.ppf(0.975, 2 * (d + 1)) / 2 / pm
    return lo, hi

def fig3():
    fig = plt.figure(figsize=(W2, 62 * MM))
    gs = fig.add_gridspec(1, 3, width_ratios=[1, 1.25, 1.05], wspace=0.42, left=0.065, right=0.99, top=0.9, bottom=0.17)
    axA, axB, axC = [fig.add_subplot(gs[0, i]) for i in range(3)]
    # ---- A crude interval rates with exact Poisson CI
    cr = A06["reproduction"]["crude_rates"]; recA = []
    for r in cr:
        a, b = [float(v) for v in r["interval"].split("-")]
        lo, hi = poisson_exact(r["deaths"], r["person_months"])
        axA.plot([a, b], [r["rate"], r["rate"]], color=COL["NEC"], lw=1.8, solid_capstyle="butt", zorder=3)
        axA.fill_between([a, b], [lo, lo], [hi, hi], color=COL["NEC"], alpha=0.15, lw=0, zorder=2)
        recA.append({"interval": r["interval"], "deaths": r["deaths"], "person_months": r["person_months"], "rate": r["rate"], "exact_ci": [lo, hi]})
    axA.set_xlim(0, 24); axA.set_ylim(0, 0.2); axA.set_xticks([0, 3, 6, 9, 12, 15, 18, 24])
    axA.set_xlabel("Months after diagnosis"); axA.set_ylabel("Cancer deaths per person-month")
    rr = A06["reproduction"]["adj_RR"]["r69_36"]
    axA.annotate("", xy=(7.5, cr[2]["rate"] + 0.012), xytext=(4.5, cr[1]["rate"] + 0.012), arrowprops=dict(arrowstyle="->", lw=0.7, color=COL["SCLC"]))
    axA.text(9.0, 0.188, f"Adjusted RR 6-9 vs 3-6 mo:\n{rr['RR']:.2f} ({rr['lo']:.2f}-{rr['hi']:.2f}), P = {rr['p']:.4f}", fontsize=5.6, ha="left", va="top", color=COL["black"])
    axA.text(0.98, 0.03, "Distant NEC, n = 908\nbands: exact Poisson 95% CI", transform=axA.transAxes, fontsize=5.4, ha="right", va="bottom", color=COL["grey"])
    panel_label(axA, "A", dx=-0.3)
    # ---- B smoothed hazards at bw 1.5
    grid = list(csv.DictReader(open(RES + "A09_hazard_grid.csv")))
    recB = {}
    for key in ["AC", "PANC", "SCLC", "NEC"]:
        g = [r for r in grid if r["group"] == key and r["bw_label"] == "bw1.5"]
        mo = np.array([float(r["month"]) for r in g]); h = np.array([float(r["hazard"]) for r in g])
        keep = (mo >= 0.25) & (mo <= 23.75)   # boundary points at 0 and 24 are not interpreted
        if key == "NEC":
            lo = np.array([float(r["lower95"]) for r in g]); hi = np.array([float(r["upper95"]) for r in g])
            axB.fill_between(mo[keep], lo[keep], hi[keep], color=COL["NEC"], alpha=0.15, lw=0)
            recB[key] = {"month": mo.tolist(), "hazard": h.tolist(), "lower95": lo.tolist(), "upper95": hi.tolist()}
        else:
            recB[key] = {"month": mo.tolist(), "hazard": h.tolist()}
        axB.plot(mo[keep], h[keep], color=COL[key], ls=LS[key], lw=1.3 if key == "NEC" else 1.1, label=LAB[key])
        nar = np.array([int(r["n_at_risk"]) for r in g])
    axB.set_xlim(0, 24); axB.set_ylim(0, 0.2); axB.set_xticks([0, 3, 6, 9, 12, 15, 18, 21, 24])
    axB.set_xlabel("Months after diagnosis"); axB.set_ylabel("Smoothed cancer-death hazard\n(per person-month)")
    axB.legend(loc="upper right", fontsize=5.6, handlelength=2.2)
    axB.text(0.02, 0.03, "Epanechnikov kernel, common bandwidth 1.5 mo;\nband: NEC pointwise bootstrap 95% CI (B = 500)", transform=axB.transAxes, fontsize=5.4, color=COL["grey"])
    panel_label(axB, "B", dx=-0.22)
    # ---- C chemo-stratum RR forest
    ci = A06["reproduction"]["chemo_interaction"]
    ivs = ["6-9", "9-15", "15-24"]; recC = []
    ypos = []; y = 0
    for j, iv in enumerate(ivs):
        for k, (grp, col, mk, lab) in enumerate([("chemo_yes", COL["NEC"], "o", f"Chemotherapy recorded (n = {ci['n_chemo']})"),
                                                  ("chemo_no_or_unknown", COL["grey"], "s", f"No/unknown chemotherapy (n = {ci['n_nochemo']})")]):
            e = ci[f"per_{iv}_vs_3-6"][grp]
            yy = -(j * 3 + k)
            axC.plot([e["lo"], e["hi"]], [yy, yy], color=col, lw=0.9)
            axC.plot(e["RR"], yy, marker=mk, color=col, ms=4, mfc=col if grp == "chemo_yes" else "white", mew=0.9, label=lab if j == 0 else None)
            axC.text(4.6, yy, f"{e['RR']:.2f} ({e['lo']:.2f}-{e['hi']:.2f})", fontsize=5.6, va="center", ha="left", color=COL["black"])
            recC.append({"interval": iv, "group": grp, "RR": e["RR"], "lo": e["lo"], "hi": e["hi"], "p": e["p"]})
        ypos.append(-(j * 3 + 0.5))
    axC.axvline(1, color=COL["light"], lw=0.7, ls="--", zorder=0)
    axC.set_xscale("log"); axC.set_xlim(0.1, 4.2); axC.set_xticks([0.125, 0.25, 0.5, 1, 2, 4]); axC.set_xticklabels(["0.125", "0.25", "0.5", "1", "2", "4"])
    axC.set_yticks(ypos); axC.set_yticklabels([f"{iv} mo" for iv in ivs]); axC.set_ylim(-9.6, 1.2)
    axC.set_xlabel("Rate ratio vs 3-6 months (log scale)")
    axC.legend(loc="upper left", bbox_to_anchor=(-0.05, 1.17), fontsize=5.6, ncol=1, handlelength=1.2)
    ig = ci["interaction_global"]
    axC.text(0.99, 0.02, f"Among 3-month survivors (n = {ci['n_survivors_3m']})\nChemotherapy x interval: W = {ig['W']:.2f}, {ig['df']} df, P = {ig['p']:.3f}", transform=axC.transAxes, fontsize=5.4, ha="right", va="bottom", color=COL["grey"])
    for sp in ["left"]: axC.spines[sp].set_visible(False)
    axC.tick_params(axis="y", length=0)
    panel_label(axC, "C", dx=-0.3)
    save(fig, "Fig3_hazard")
    MANIFEST["figures"]["Fig3_hazard"] = {"A": recA, "A_adjRR": rr, "B": recB, "C": recC, "C_interaction": ig}

# ================================================================================ Fig 4 window / TOST
def fig4():
    fig = plt.figure(figsize=(W2, 65 * MM))
    gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.15], wspace=0.45, left=0.07, right=0.99, top=0.9, bottom=0.17)
    axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[0, 1])
    wins = [3, 6, 9, 12, 18, 24]; recA = {}
    off = {"NEC": -0.45, "SCLC": -0.15, "AC": 0.15, "PANC": 0.45}
    for key in ["NEC", "SCLC", "AC", "PANC"]:
        t = A02["tumours"][key]; est = []; lo = []; hi = []
        for w in wins:
            d = t["windows"][str(w)]["differences"]["6-0_cancer_specific"]; est.append(d["estimate"]); lo.append(d["ci"][0]); hi.append(d["ci"][1])
        est, lo, hi = map(np.array, (est, lo, hi)); xw = np.array(wins, float) + off[key]
        axA.errorbar(xw, est, yerr=[est - lo, hi - est], color=COL[key], marker=MK[key], ms=3.5, lw=1.1, ls=LS[key], capsize=1.5, elinewidth=0.7,
                     label=f"{LAB[key]} (n = {t['n_distant']:,})", mfc=COL[key] if key != "AC" else "white", mew=0.8)
        recA[key] = {"windows": wins, "diff_6_0": est.tolist(), "ci_lo": lo.tolist(), "ci_hi": hi.tolist(), "B": t["B"], "n": t["n_distant"]}
    axA.axhline(0, color=COL["light"], lw=0.7, ls="--", zorder=0)
    axA.axvspan(11, 13, color=COL["NEC"], alpha=0.06, lw=0); axA.text(12, -14.6, "primary\nwindow", ha="center", va="bottom", fontsize=5.4, color=COL["grey"])
    axA.set_xticks(wins); axA.set_xlim(1.5, 25.5); axA.set_ylim(-15, 16)
    axA.set_xlabel("Width of the risk window after each landmark (months)"); axA.set_ylabel("Change in conditional cancer-death\nprobability, landmark 6 minus 0 (pp)")
    axA.legend(loc="upper right", fontsize=5.2, handlelength=2.2, bbox_to_anchor=(1.02, 1.03))
    panel_label(axA, "A")
    # ---- B equivalence
    rows = A03["within_NEC"]; recB = []
    order = ["NEC 6 minus 0 months", "NEC 3 minus 0 months", "NEC 6 minus 3 months", "NEC 12 minus 0 months"]
    rows = sorted(rows, key=lambda r: order.index(r["contrast"]))
    for m, col, ls in [(5, COL["SCLC"], "-"), (3, COL["grey"], "--")]:
        axB.axvline(-m, color=col, lw=0.8, ls=ls, zorder=0); axB.axvline(m, color=col, lw=0.8, ls=ls, zorder=0)
    axB.axvspan(-5, 5, color=COL["SCLC"], alpha=0.05, lw=0)
    axB.axvline(0, color=COL["light"], lw=0.6, zorder=0)
    labs = []
    for i, r in enumerate(rows):
        yy = -i
        c90 = r["ci90"]; c95 = r["ci95"]
        axB.plot(c95, [yy, yy], color=COL["NEC"], lw=0.8, solid_capstyle="butt")
        axB.plot(c90, [yy, yy], color=COL["NEC"], lw=3.0, solid_capstyle="butt")
        axB.plot(r["estimate"], yy, marker="o", color="white", mec=COL["NEC"], ms=4.5, mew=1.0, zorder=5)
        p5 = [mm for mm in r["margins"] if mm["margin"] == 5][0]; p3 = [mm for mm in r["margins"] if mm["margin"] == 3][0]
        axB.text(10.5, yy, f"{r['estimate']:+.2f} ({c95[0]:+.1f} to {c95[1]:+.1f})\nTOST P: ±5 {'<0.001' if p5['p_tost'] < 0.0005 else format(p5['p_tost'], '.3f')}; ±3 {'<0.001' if p3['p_tost'] < 0.0005 else format(p3['p_tost'], '.3f')}",
                 fontsize=5.4, va="center", ha="left", color=COL["black"], linespacing=1.2)
        labs.append(r["contrast"].replace("NEC ", "Landmark ").replace(" months", " mo"))
        recB.append({"contrast": r["contrast"], "estimate": r["estimate"], "ci90": c90, "ci95": c95, "p_tost_5": p5["p_tost"], "p_tost_3": p3["p_tost"],
                     "equivalent_5": p5["equivalent_alpha05"], "equivalent_3": p3["equivalent_alpha05"]})
    axB.set_yticks([-i for i in range(len(rows))]); axB.set_yticklabels(labs); axB.set_ylim(-4.3, 0.8)
    axB.set_xlim(-22, 22); axB.set_xticks([-20, -15, -10, -5, -3, 0, 3, 5, 10])
    axB.set_xlabel("Difference in 12-month cancer-death probability (pp)")
    axB.text(5.3, 0.6, "±5 pp margin", fontsize=5.4, color=COL["SCLC"], ha="left"); axB.text(3.2, -4.15, "±3 pp", fontsize=5.4, color=COL["grey"], ha="left", va="bottom")
    axB.text(0.01, 0.98, "Thick: 90% CI (TOST); thin: 95% CI\nDistant NEC n = 908, B = 2000", transform=axB.transAxes, fontsize=5.2, color=COL["grey"], va="top")
    axB.spines["left"].set_visible(False); axB.tick_params(axis="y", length=0)
    panel_label(axB, "B", dx=-0.28)
    save(fig, "Fig4_window_equivalence")
    MANIFEST["figures"]["Fig4_window_equivalence"] = {"A": recA, "B": recB, "B_margins_status": A03["margins_status"]}

# ============================================================================= Fig 5 cross-tumour
def fig5():
    fig = plt.figure(figsize=(W2, 130 * MM))
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 0.34, 1.05], hspace=0.5, wspace=0.34, left=0.07, right=0.99, top=0.96, bottom=0.06)
    axA = fig.add_subplot(gs[0, 0]); axAt = fig.add_subplot(gs[1, 0], sharex=axA)
    axB = fig.add_subplot(gs[0:2, 1]); axC = fig.add_subplot(gs[2, 0]); axD = fig.add_subplot(gs[2, 1])
    # ---- A four-tumour curves
    rows = []; recA = {}
    for key in ["AC", "PANC", "SCLC", "NEC"]:
        m = A11["table3"][key]["cancer"]["monthly"]; x = m["landmarks"]
        curve_with_markers(axA, x, m["pct"], m["n_at_risk"], COL[key], MK[key], ls=LS[key], lo=m["lo"], hi=m["hi"],
                           label=f"{LAB[key]} (n = {A11['table3'][key]['cancer']['n']:,})", every=3, ms=2.8)
        rows.append(({"NEC": "NEC", "SCLC": "SCLC", "AC": "CRC AC", "PANC": "Panc AC"}[key], COL[key], {x[i]: m["n_at_risk"][i] for i in range(len(x))}))
        recA[key] = {"landmarks": x, "pct": m["pct"], "lo": m["lo"], "hi": m["hi"], "n_at_risk": m["n_at_risk"]}
    rows = rows[::-1]
    axA.set_xlim(-1.6, 24.8); axA.set_ylim(8, 85); axA.set_xticks(range(0, 25, 3))
    axA.set_ylabel("12-month cancer-death probability (%)"); axA.set_xlabel("Landmark (months after diagnosis)")
    axA.legend(loc="lower left", fontsize=5.4, handlelength=2.4, bbox_to_anchor=(-0.01, -0.02))
    risk_table(axAt, axA, rows, list(range(0, 25, 3)))
    panel_label(axA, "A", dx=-0.2)
    # ---- B paired-difference forest
    contrasts = ["6-0", "12-6", "12-0"]; comps = ["SCLC", "LCNEC", "AC", "PANC"]; recB = []
    yy = 0; yt = []; yl = []
    for ci_i, c in enumerate(contrasts):
        axB.text(-19.5, yy + 0.55, f"Change {c} months: NEC minus comparator", fontsize=6, fontweight="bold", va="bottom")
        for key in comps:
            src = A11["lcnec"] if key == "LCNEC" else A11["table3"][key]
            d = src["nec_minus_comparator_cancer"][c]
            axB.plot(d["ci"], [yy, yy], color=COL[key], lw=1.0)
            axB.plot(d["diff_pct_points"], yy, marker=MK[key], color=COL[key], ms=4, mfc=COL[key], mew=0.8)
            axB.text(13.3, yy, f"{d['diff_pct_points']:+.1f} ({d['ci'][0]:+.1f} to {d['ci'][1]:+.1f})", fontsize=5.5, va="center", ha="left")
            yt.append(yy); yl.append({"SCLC": "vs SCLC", "LCNEC": "vs lung LCNEC", "AC": "vs colorectal AC", "PANC": "vs pancreatic AC"}[key])
            recB.append({"contrast": c, "comparator": key, "diff": d["diff_pct_points"], "ci": d["ci"], "boot_p": d.get("boot_p_two_sided")})
            yy -= 1
        yy -= 1.4
    axB.axvline(0, color=COL["light"], lw=0.7, ls="--", zorder=0)
    axB.set_yticks(yt); axB.set_yticklabels(yl); axB.set_ylim(yy + 1.0, 1.4); axB.set_xlim(-20, 26); axB.set_xticks([-15, -10, -5, 0, 5, 10])
    axB.set_xlabel("NEC minus comparator, difference in change (pp)")
    axB.spines["left"].set_visible(False); axB.tick_params(axis="y", length=0)
    axB.text(0.99, 0.01, "Positive = NEC declines less than the comparator.\nIndependent bootstraps paired by index, B = 2000", transform=axB.transAxes, fontsize=5.2, ha="right", color=COL["grey"])
    panel_label(axB, "B", dx=-0.26, dy=1.01)
    # ---- C 2x2
    cells = A11["two_by_two"]["cells"]; recC = {}
    pos = {("colorectal", "SC"): 0 - 0.17, ("colorectal", "LC"): 0 + 0.17, ("lung", "SC"): 1 - 0.17, ("lung", "LC"): 1 + 0.17}
    style = {"SC": (COL["NEC"], "o"), "LC": (COL["LCNEC"], "s")}
    for cname, cell in cells.items():
        organ, morph = cname.split("_"); d = cell["cancer"]["change"]["6-0"]; n = cell["cancer"]["n"]
        col, mk = style[morph]; xx = pos[(organ, morph)]
        axC.plot([xx, xx], d["ci"], color=col, lw=1.0)
        axC.plot(xx, d["pct_points"], marker=mk, color=col, ms=5, mfc=col, mew=0.8, label=f"{'Small-cell' if morph == 'SC' else 'Large-cell'} morphology" if organ == "colorectal" else None)
        axC.text(xx + (0.06 if morph == "LC" else -0.06), d["ci"][1] + 0.6, f"{d['pct_points']:+.1f}\n(n = {n:,})", fontsize=5.3, ha="left" if morph == "LC" else "right", va="bottom", color=COL["black"])
        recC[cname] = {"n": n, "diff_6_0": d["pct_points"], "ci": d["ci"], "plateau_end": cell["cancer"]["plateau_end"]["months"]}
    axC.axhline(0, color=COL["light"], lw=0.7, ls="--", zorder=0)
    axC.set_xticks([0, 1]); axC.set_xticklabels(["Colorectal", "Lung"]); axC.set_xlim(-0.6, 1.6); axC.set_ylim(-16, 12)
    axC.set_ylabel("Change 6-0 months in 12-month\ncancer-death probability (pp)")
    axC.legend(loc="upper right", fontsize=5.6)
    tb = A11["two_by_two"]
    axC.text(0.99, 0.03, f"Morphology effect (SC - LC, mean over organs): {tb['morphology_effect_6_0']['estimate']:+.1f} ({tb['morphology_effect_6_0']['ci'][0]:+.1f} to {tb['morphology_effect_6_0']['ci'][1]:+.1f})\n"
             f"Organ effect: {tb['organ_effect_6_0']['estimate']:+.1f} ({tb['organ_effect_6_0']['ci'][0]:+.1f} to {tb['organ_effect_6_0']['ci'][1]:+.1f})\n"
             f"Interaction: {tb['interaction_6_0']['estimate']:+.1f} ({tb['interaction_6_0']['ci'][0]:+.1f} to {tb['interaction_6_0']['ci'][1]:+.1f})",
             transform=axC.transAxes, fontsize=5.3, ha="right", va="bottom", color=COL["black"], linespacing=1.3)
    panel_label(axC, "C", dx=-0.2)
    recC["effects"] = {k: tb[k] for k in ["morphology_effect_6_0", "organ_effect_6_0", "interaction_6_0"]}
    # ---- D RR 6-9 vs 3-6 by tumour
    recD = []
    items = [("NEC", A06["reproduction"]["adj_RR"]["r69_36"], A06["n"]["NEC_distant"]),
             ("SCLC", A06["comparators_histology"]["SCLC"]["r69_36"], A06["comparators_histology"]["SCLC"]["n"]),
             ("AC", A06["comparators_histology"]["AC"]["r69_36"], A06["comparators_histology"]["AC"]["n"]),
             ("PANC", A06["comparators_histology"]["PANC"]["r69_36"], A06["comparators_histology"]["PANC"]["n"])]
    for i, (key, e, n) in enumerate(items):
        yy = -i
        axD.plot([e["lo"], e["hi"]], [yy, yy], color=COL[key], lw=1.0)
        axD.plot(e["RR"], yy, marker=MK[key], color=COL[key], ms=4.5, mfc=COL[key], mew=0.8)
        axD.text(2.15, yy, f"{e['RR']:.2f} ({e['lo']:.2f}-{e['hi']:.2f})", fontsize=5.6, va="center", ha="left")
        recD.append({"tumour": key, "n": n, "RR": e["RR"], "lo": e["lo"], "hi": e["hi"], "p": e["p"]})
    axD.axvline(1, color=COL["light"], lw=0.7, ls="--", zorder=0)
    axD.set_xscale("log"); axD.set_xlim(0.7, 3.3); axD.set_xticks([0.75, 1, 1.5, 2]); axD.set_xticklabels(["0.75", "1", "1.5", "2"]); axD.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
    axD.set_yticks([0, -1, -2, -3]); axD.set_yticklabels([f"{LAB[k]}\n(n = {n:,})" for k, _, n in items]); axD.set_ylim(-5.8, 0.6)
    axD.set_xlabel("Adjusted RR, cancer death 6-9 vs 3-6 months (log scale)")
    axD.spines["left"].set_visible(False); axD.tick_params(axis="y", length=0)
    st = A06["calendar_period"]["shape_tests"]
    def pf(p): return "< 0.001" if p < 0.001 else f"= {p:.3f}"
    sh = {k: st[k]["overall_linear_year"] for k in ["SCLC", "AC", "PANC"]}
    axD.text(0.99, 0.02, "Piecewise-exponential Poisson, cluster-robust SE;\nadjusted for age, age², sex, year (+ histology, site)\n"
             f"Hazard-shape test vs NEC (Wald, 6 df): SCLC P {pf(sh['SCLC']['p'])};\ncolorectal AC P {pf(sh['AC']['p'])}; pancreatic AC P {pf(sh['PANC']['p'])}",
             transform=axD.transAxes, fontsize=5.2, ha="right", va="bottom", color=COL["grey"])
    for k in sh: recD.append({"shape_test_vs_NEC": k, "W": sh[k]["W"], "df": sh[k]["df"], "p": sh[k]["p"]})
    panel_label(axD, "D", dx=-0.3)
    save(fig, "Fig5_cross_tumour")
    MANIFEST["figures"]["Fig5_cross_tumour"] = {"A": recA, "B": recB, "C": recC, "D": recD}

# ================================================================== Fig 6 strata forest + validation
def fig6():
    fig = plt.figure(figsize=(W2, 105 * MM))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.25, 1], height_ratios=[1, 0.2], hspace=0.4, wspace=0.42, left=0.21, right=0.99, top=0.95, bottom=0.05)
    axA = fig.add_subplot(gs[:, 0]); axB = fig.add_subplot(gs[0, 1]); axBt = fig.add_subplot(gs[1, 1], sharex=axB)
    pooled = A05["pooled"]
    axA.axvspan(pooled["diff_6_0"]["ci"][0], pooled["diff_6_0"]["ci"][1], color=COL["NEC"], alpha=0.07, lw=0, zorder=0)
    axA.axvline(pooled["diff_6_0"]["est"], color=COL["NEC"], lw=0.8, ls="--", zorder=0)
    axA.axvline(0, color=COL["light"], lw=0.6, zorder=0)
    recA = []; yy = 0; yt = []; yl = []
    treat_vars = {"Surgery (primary site)", "Chemotherapy recorded"}
    groups = []
    for s in A05["strata"]:
        if not groups or groups[-1][0] != s["variable"]: groups.append((s["variable"], []))
        groups[-1][1].append(s)
    # pooled first
    axA.plot(pooled["diff_6_0"]["ci"], [yy, yy], color=COL["NEC"], lw=1.4)
    axA.plot(pooled["diff_6_0"]["est"], yy, marker="D", color=COL["NEC"], ms=4.5)
    axA.text(17.5, yy, f"{pooled['diff_6_0']['est']:+.1f} ({pooled['diff_6_0']['ci'][0]:+.1f} to {pooled['diff_6_0']['ci'][1]:+.1f})", fontsize=5.4, va="center")
    yt.append(yy); yl.append(f"All distant NEC (n = {pooled['n']})"); yy -= 1.6
    for var, members in groups:
        treat = var in treat_vars
        axA.text(-29.5, yy + 0.45, var + ("  [6 minus 3 mo, 3-month survivors]" if treat else ""), fontsize=5.8, fontweight="bold", va="bottom", clip_on=False)
        for s in members:
            if treat:
                c = s["conditional_on_3mo_survival"]; d = c["diff_6_3"]; n_lbl = c["n"]; nr6 = c["landmark"]["6"]["n_at_risk"]
            else:
                d = s["diff_6_0"]; n_lbl = s["n"]; nr6 = s["landmark"]["6"]["n_at_risk"]
            col = COL["SCLC"] if treat else COL["black"]
            axA.plot(d["ci"], [yy, yy], color=col, lw=0.9)
            axA.plot(d["est"], yy, marker="o", color=col, ms=3.8, mfc=col if nr6 >= 100 else "white", mew=0.9)
            axA.text(17.5, yy, f"{d['est']:+.1f} ({d['ci'][0]:+.1f} to {d['ci'][1]:+.1f})", fontsize=5.4, va="center")
            yt.append(yy); yl.append(f"{s['level']} (n = {n_lbl})")
            recA.append({"variable": var, "level": s["level"], "n": n_lbl, "contrast": "6-3 conditional on 3-mo survival" if treat else "6-0",
                         "est": d["est"], "ci": d["ci"], "n_at_risk_lm6": nr6, "plateau_end": s["plateau_end"]["est"]})
            yy -= 1
        yy -= 0.9
    axA.set_yticks(yt); axA.set_yticklabels(yl, fontsize=5.8); axA.set_ylim(yy + 0.5, 1.0); axA.set_xlim(-30, 32)
    axA.set_xticks([-25, -20, -15, -10, -5, 0, 5, 10, 15]); axA.set_xlabel("Change in 12-month cancer-death probability (pp), 95% bootstrap CI")
    axA.spines["left"].set_visible(False); axA.tick_params(axis="y", length=0)
    axA.text(0.01, 0.005, "Open marker: < 100 at risk at landmark 6. Dashed line and band: pooled estimate and 95% CI.", transform=axA.transAxes, fontsize=5.0, color=COL["grey"])
    panel_label(axA, "A", dx=-0.62, dy=1.0)
    # ---- B derivation vs confirmation
    recB = {}; rows = []
    for key, col, mk, ls in [("derivation_2000_2015", COL["NEC"], "o", "-"), ("confirmation_2016_2023", COL["SCLC"], "s", "--")]:
        p = A07[key]; x = np.arange(25); y = np.array(p["curve_0_24"], float); n = np.array(p["nrisk_0_24"], float)
        keep = n >= 20
        curve_with_markers(axB, x[keep], y[keep], n[keep], col, mk, ls=ls, every=1, ms=2.6,
                           label=f"{p['label'].capitalize()} (n = {p['n']}): plateau end {p['plateau_end']['est']:.2f} mo ({p['plateau_end']['ci'][0]:.2f}-{p['plateau_end']['ci'][1]:.2f})")
        for lm, v in p["landmarks"].items():
            axB.errorbar(float(lm) + (0.15 if key.startswith("conf") else -0.15), v["cif12"], yerr=[[v["cif12"] - v["ci"][0]], [v["ci"][1] - v["cif12"]]], color=col, elinewidth=0.7, capsize=1.6, lw=0, zorder=2)
        rows.append((p["label"].split(" ")[0].capitalize(), col, {int(x[i]): p["nrisk_0_24"][i] for i in range(25)}))
        recB[key] = {"n": p["n"], "curve_0_24": p["curve_0_24"], "nrisk_0_24": p["nrisk_0_24"], "landmarks": p["landmarks"], "change_6_minus_0": p["change_6_minus_0"], "plateau_end": p["plateau_end"]}
    pd_ = A07["period_difference_in_6_minus_0"]
    axB.text(0.02, 0.04, f"Change 6-0: derivation {A07['derivation_2000_2015']['change_6_minus_0']['est']:+.1f} ({A07['derivation_2000_2015']['change_6_minus_0']['ci'][0]:+.1f} to {A07['derivation_2000_2015']['change_6_minus_0']['ci'][1]:+.1f});\n"
             f"confirmation {A07['confirmation_2016_2023']['change_6_minus_0']['est']:+.1f} ({A07['confirmation_2016_2023']['change_6_minus_0']['ci'][0]:+.1f} to {A07['confirmation_2016_2023']['change_6_minus_0']['ci'][1]:+.1f});\n"
             f"difference {pd_['est']:+.1f} ({pd_['ci'][0]:+.1f} to {pd_['ci'][1]:+.1f}), P = {pd_['boot_p']:.2f}", transform=axB.transAxes, fontsize=5.3, color=COL["black"], linespacing=1.25)
    axB.set_xlim(-0.5, 24.5); axB.set_ylim(0, 92); axB.set_xticks(range(0, 25, 3))
    axB.set_ylabel("12-month cancer-death probability (%)"); axB.set_xlabel("Landmark (months after diagnosis)")
    axB.legend(loc="upper right", fontsize=5.4, handlelength=2.4)
    axB.text(0.98, 0.17, "Error bars: 95% bootstrap CI at pre-specified landmarks;\npoints with < 20 at risk suppressed; open markers < 100 at risk", transform=axB.transAxes, fontsize=5.0, ha="right", va="bottom", color=COL["grey"])
    risk_table(axBt, axB, rows, list(range(0, 25, 3)))
    panel_label(axB, "B", dx=-0.2, dy=1.0)
    recB["period_difference_in_6_minus_0"] = pd_
    save(fig, "Fig6_strata_validation")
    MANIFEST["figures"]["Fig6_strata_validation"] = {"A": recA, "A_pooled": pooled, "B": recB}

if __name__ == "__main__":   # guard added 2026-10-08 so that FIG_manuscript_final.py can import the helpers without re-running
    os.makedirs(OUT, exist_ok=True)
    fig1(); fig2(); fig3(); fig4(); fig5(); fig6()
    MANIFEST["runtime_minutes"] = round((time.time() - T0) / 60, 2)
    MANIFEST["output_dir"] = OUT
    with open(RES + "FIG.json", "w") as f: json.dump(MANIFEST, f, indent=1, default=float)
    print("done", MANIFEST["runtime_minutes"], "min")
