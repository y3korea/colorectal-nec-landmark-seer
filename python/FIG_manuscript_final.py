# -*- coding: utf-8 -*-
"""FIG_manuscript_final.py -- final figure files for the International Journal of Cancer submission.

Builds Figures 1-3 and Supplementary Figures S1-S7 from the stored analysis outputs only:
results/{A01,A02,A03,A05,A05_std,A06,A07,A08,A09,A11,A12,FIG_flow_counts}.json, A09_hazard_grid.csv and
data/flow_all.csv (aggregate eligibility-step counts; 78 rows, no patient-level records). No patient-level data are
read, written or copied, and no estimate is recomputed (no bootstrap, no model fit): every plotted or printed number
is read from those files and only formatted here.

Mapping (write/integration_log.md section 12 item 1 and the legends in Manuscript_IJC_v2.md / Supplement_v2.md):
  Figure1  (A) landmark curves by stage, (B) distant curve + standardized curve + broken-stick plateau,
           (C) 6-0 change vs window width, (D) equivalence (TOST) plot
  Figure2  (A) crude interval rates, (B) smoothed hazards at common bandwidth 1.5 months, (C) chemotherapy-stratum RRs
  Figure3  (A) four-tumour monthly curves, (B) NEC-minus-comparator forest incl. lung LCNEC, (C) 2x2 morphology x organ,
           (D) RR 6-9 vs 3-6 months by tumour, incl. lung LCNEC (A12.json, same R code path that reproduces A06 NEC/SCLC)
  SuppFigS1 cohort flow (A NEC, B comparators); S2 strata forest; S3 derivation vs confirmation; S4 six window widths
  (A)-(F); S5 NEC smoothed hazard at three bandwidths; S6 standardized curves (A) covariate sets (B) comparators;
  S7 time-updated chemotherapy (A) landmark curves (B) crude rates split at treatment start.

Format: 180 mm wide; Arial (checked at run time); text >= 6 pt; panel labels "(A)"; Okabe-Ito palette with marker
shape / line style redundancy; filled markers >= 100 at risk, open < 100; PNG 300 dpi + vector PDF (TrueType fonts).
Numbers: round half up (Decimal) and the Unicode minus sign U+2212, as in the manuscript.
Output: <NEC_FIG_DIR>/manuscript_final/  (+ figure_data.json with every plotted value and the source key).
Run: python3 python/FIG_manuscript_final.py   (output: <NEC_FIG_DIR>/manuscript_final/)
"""
import json, os, sys, csv, time, re, textwrap
from decimal import Decimal, ROUND_HALF_UP
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.transforms as mtrans
import matplotlib.font_manager as fm
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, FancyBboxPatch, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import config                                      # noqa: E402
import FIG_final_figures as F                      # noqa: E402  (module is guarded by __main__; import only loads JSON)
from FIG_final_figures import (poisson_exact, COL, MK, LS, STAGE_COL, STAGE_MK,  # noqa: E402
                               A01, A02, A03, A05, A05S, A06, A07, A11, FLOWC, flow_rows, BASE, RES)

T0 = time.time()
OUT = config.MANUSCRIPT_FIG_DIR + os.sep
os.makedirs(OUT, exist_ok=True)
A08 = F.J("A08.json")
A09 = F.J("A09.json")
A12 = F.J("A12.json")
GRID = list(csv.DictReader(open(RES + "A09_hazard_grid.csv")))

# ------------------------------------------------------------------------------------------------ style
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 6.5, "axes.labelsize": 7, "axes.titlesize": 7, "legend.fontsize": 6.2, "legend.title_fontsize": 6.2,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "axes.linewidth": 0.6, "xtick.major.width": 0.6,
    "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5, "lines.linewidth": 1.2,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False, "axes.unicode_minus": True,
    "pdf.fonttype": 42, "ps.fonttype": 42, "savefig.dpi": 300, "figure.dpi": 100, "mathtext.default": "regular",
    "legend.handlelength": 2.4, "legend.borderaxespad": 0.3, "legend.labelspacing": 0.35, "axes.labelpad": 2.5})
FONT_FILE = fm.findfont(fm.FontProperties(family=["sans-serif"]), fallback_to_default=False)
FONT_NAME = fm.FontProperties(fname=FONT_FILE).get_name()
if FONT_NAME != "Arial":   # journal figures use Arial; elsewhere (e.g. Linux/Colab without Arial) fall back with a warning
    if os.environ.get("NEC_REQUIRE_ARIAL", "0") == "1":
        raise SystemExit(f"Arial not resolved (got {FONT_NAME}); see README")
    print(f"WARNING: Arial not available, figures use {FONT_NAME} (layout may differ slightly from the submitted files)")
MANIFEST = {"id": "FIG_manuscript_final", "date": "2026-10-08", "output_dir": config.relpath(OUT), "font": {"name": FONT_NAME, "file": os.path.basename(FONT_FILE)},
            "figures": {}, "panel_text": {}}
MINUS = "−"
LAB = {"NEC": "Colorectal NEC", "SCLC": "Small-cell lung cancer", "AC": "Colorectal adenocarcinoma",
       "PANC": "Pancreatic adenocarcinoma", "LCNEC": "Lung LCNEC"}
STAGE_LS = {"Distant": "-", "Regional": "--", "Localized": ":"}
FS = 6.2        # body text in panels (minimum 6 pt)
GREY = COL["grey"]


# ------------------------------------------------------------------------------------------------ helpers
def r(x, d):
    """round half up to d decimals; Unicode minus; no '-0'"""
    q = Decimal(repr(float(x))).quantize(Decimal(1).scaleb(-d), rounding=ROUND_HALF_UP)
    if q == 0:
        q = abs(q)
    s = f"{q:.{d}f}"
    return s.replace("-", MINUS)


def sg(x, d):
    s = r(x, d)
    return s if s.startswith(MINUS) or Decimal(s) == 0 else "+" + s


def ci(e, lo, hi, d=1, sign=True):
    """manuscript style: +0.53 (−3.55 to 4.67)"""
    return f"{sg(e, d) if sign else r(e, d)} ({r(lo, d)} to {r(hi, d)})"


def rr(e, lo, hi, d=2):
    return f"{r(e, d)} ({r(lo, d)}–{r(hi, d)})"


def phu(p, nd=3):
    return r(p, nd)


def n_(x):
    return f"{int(x):,}"


def mo(c):
    """change-point months: 0 printed as '0', otherwise two decimals (legend style)"""
    return "0" if float(c) == 0 else r(c, 2)


def newfig(w, h):
    fig = plt.figure(figsize=(w / 25.4, h / 25.4))
    fig._wh = (w, h)
    return fig


def axmm(fig, x, y, w, h):
    W, H = fig._wh
    return fig.add_axes([x / W, y / H, w / W, h / H])


def plab(fig, x, y, s, title=None):
    W, H = fig._wh
    fig.text(x / W, y / H, f"({s})", fontsize=9, fontweight="bold", va="top", ha="left")
    if title:
        fig.text((x + 7.2) / W, (y - 0.6) / H, title, fontsize=7, fontweight="bold", va="top", ha="left")


def ftext(fig, x, y, s, **kw):
    W, H = fig._wh
    kw.setdefault("fontsize", FS)
    return fig.text(x / W, y / H, s, **kw)


def save(fig, name):
    texts = []
    for t in fig.findobj(matplotlib.text.Text):
        s = t.get_text()
        if s and s.strip() and t.get_visible():
            texts.append(s)
    MANIFEST["panel_text"][name] = texts
    if config.SYNTHETIC:   # never let figures from synthetic data be mistaken for results (not added for real data)
        fig.text(0.5, 0.5, "SYNTHETIC DATA - NOT STUDY RESULTS", fontsize=22, color="red", alpha=0.35, rotation=30, ha="center", va="center", fontweight="bold")
    fig.savefig(OUT + name + ".png", dpi=300, facecolor="white")
    fig.savefig(OUT + name + ".pdf", dpi=600, facecolor="white")
    plt.close(fig)
    print("saved", name)


def curve(ax, x, y, n, color, marker, ls="-", lo=None, hi=None, every=1, ms=2.8, thresh=100, lw=1.1, z=3, band=0.15):
    x = np.asarray(x, float); y = np.asarray(y, float); n = np.asarray(n, float)
    if lo is not None:
        ax.fill_between(x, lo, hi, color=color, alpha=band, lw=0, zorder=z - 1)
    ax.plot(x, y, color=color, ls=ls, lw=lw, zorder=z)
    sel = np.arange(0, len(x), every)
    full = sel[n[sel] >= thresh]; opn = sel[n[sel] < thresh]
    ax.plot(x[full], y[full], marker=marker, ls="none", ms=ms, mfc=color, mec=color, mew=0.7, zorder=z + 1)
    ax.plot(x[opn], y[opn], marker=marker, ls="none", ms=ms, mfc="white", mec=color, mew=0.8, zorder=z + 1)


def h_line(color, marker=None, ls="-", lw=1.1, mfc=None, ms=3.2, band=False):
    ln = Line2D([], [], color=color, ls=ls, lw=lw, marker=marker, ms=ms, mfc=mfc or color, mec=color, mew=0.7)
    if band:
        return (Patch(facecolor=color, alpha=0.18, lw=0), ln)
    return ln


def h_open_filled():
    return [Line2D([], [], marker="o", ls="none", ms=3.2, mfc=GREY, mec=GREY),
            Line2D([], [], marker="o", ls="none", ms=3.2, mfc="white", mec=GREY, mew=0.8)], ["≥ 100 at risk", "< 100 at risk"]


def risk_table(fig, ax_main, rows, xs, x, y, w, h, nslots=None, title="No. at risk", labgap=4.0, fs=6.0):
    """rows: (label, colour, {landmark: n}); italic = fewer than 100 at risk."""
    ax = axmm(fig, x, y, w, h); ax.set_xlim(ax_main.get_xlim()); ax.axis("off")
    ns = nslots or len(rows)
    ax.set_ylim(0, ns + 1)
    tr = mtrans.blended_transform_factory(ax.transAxes, ax.transData)
    ax.text(0.0, ns + 0.5, title, fontsize=fs, fontweight="bold", ha="left", va="center", transform=tr)
    for i, (lab, col, d) in enumerate(rows):
        yy = ns - 0.5 - i
        ax.text(-labgap / w, yy, lab, fontsize=fs, ha="right", va="center", color=col, transform=tr, clip_on=False)
        for xx in xs:
            if xx in d and d[xx] is not None:
                nn = int(d[xx])
                ax.text(xx, yy, f"{nn:,}", fontsize=fs, ha="center", va="center", color=COL["black"], style="normal" if nn >= 100 else "italic")
    return ax


def boxed(ax, x, y, s, transform=None, **kw):
    kw.setdefault("fontsize", FS)
    kw.setdefault("color", COL["black"])
    return ax.text(x, y, s, transform=transform or ax.transAxes, bbox=dict(boxstyle="square,pad=0.25", fc="white", ec="none", alpha=0.92), **kw)


def grid_curve(group, label):
    g = [rw for rw in GRID if rw["group"] == group and rw["bw_label"] == label]
    m = np.array([float(rw["month"]) for rw in g]); h = np.array([float(rw["hazard"]) for rw in g])
    lo = np.array([float(rw["lower95"]) if rw["lower95"] not in ("", "NA") else np.nan for rw in g])
    hi = np.array([float(rw["upper95"]) if rw["upper95"] not in ("", "NA") else np.nan for rw in g])
    nr = np.array([float(rw["n_at_risk"]) if rw["n_at_risk"] not in ("", "NA") else np.nan for rw in g])
    return m, h, lo, hi, nr


# =============================================================================================== FIGURE 1
def figure1():
    fig = newfig(180, 194)
    axA = axmm(fig, 17, 124, 70, 60); axB = axmm(fig, 106, 124, 70, 60)
    axC = axmm(fig, 16, 26, 57, 60); axD = axmm(fig, 94, 26, 38, 60)
    rec = {}
    # ---------------- (A) landmark curves by stage
    rows, hs, ls_ = [], [], []
    for key, name in [("monthly_distant", "Distant"), ("monthly_regional", "Regional"), ("monthly_localized", "Localized")]:
        m = A01[key]; x = m["landmarks"]
        curve(axA, x, m["cancer_pct"], m["n_at_risk"], STAGE_COL[name], STAGE_MK[name], ls=STAGE_LS[name], lo=m["cancer_ci_lo"], hi=m["cancer_ci_hi"], ms=2.5)
        rows.append((name, STAGE_COL[name], {x[i]: m["n_at_risk"][i] for i in range(len(x))}))
        hs.append(h_line(STAGE_COL[name], STAGE_MK[name], STAGE_LS[name], band=True)); ls_.append(f"{name} (n = {n_(m['n'])})")
        rec[name] = {"source": f"A01.json {key}", "landmarks": x, "cancer_pct": m["cancer_pct"], "ci_lo": m["cancer_ci_lo"], "ci_hi": m["cancer_ci_hi"], "n_at_risk": m["n_at_risk"], "n": m["n"]}
    axA.set_xlim(-0.6, 24.6); axA.set_ylim(0, 100); axA.set_xticks(range(0, 25, 3)); axA.set_yticks(range(0, 101, 20))
    axA.set_ylabel("12-month probability of cancer death (%)"); axA.set_xlabel("Landmark (months after diagnosis)")
    lg = axA.legend(hs, ls_, loc="upper right", title="Stage (95% bootstrap band)", alignment="left"); axA.add_artist(lg)
    hh, ll = h_open_filled(); axA.legend(hh, ll, loc="upper left", handlelength=1.0)
    risk_table(fig, axA, rows, list(range(0, 25, 3)), 17, 99, 70, 15, nslots=3)
    plab(fig, 2, 192, "A")
    # ---------------- (B) distant: unadjusted, standardized, broken stick
    m = A01["monthly_distant"]; x = np.array(m["landmarks"]); pl = m["plateau_cancer"]
    axB.axvspan(pl["c_ci"][0], pl["c_ci"][1], ymin=0.14, ymax=0.76, color=COL["SCLC"], alpha=0.10, lw=0, zorder=0)
    axB.axvline(pl["c"], ymin=0.14, ymax=0.76, color=COL["SCLC"], lw=0.9, ls=":", zorder=1)
    curve(axB, x, m["cancer_pct"], m["n_at_risk"], COL["NEC"], "o", lo=m["cancer_ci_lo"], hi=m["cancer_ci_hi"], ms=2.6)
    s = A05S["original"]
    axB.plot(s["landmarks"], s["std_curve"], color=COL["black"], ls="--", lw=1.0, marker="x", ms=3.0, mew=0.8, zorder=5)
    xx = np.linspace(0, 18, 200); yy = pl["b0"] + pl["b1_per_month"] * np.maximum(0, xx - pl["c"])
    axB.plot(xx, yy, color=COL["SCLC"], lw=1.0, ls="-", zorder=4)
    hB = [h_line(COL["NEC"], "o", "-", band=True), h_line(COL["black"], "x", "--", lw=1.0), h_line(COL["SCLC"], None, "-", lw=1.0),
          (Patch(facecolor=COL["SCLC"], alpha=0.18, lw=0), Line2D([], [], color=COL["SCLC"], ls=":", lw=0.9))]
    lB = ["Unadjusted (Aalen–Johansen), 95% bootstrap band", "Directly standardized (point estimates)", "Broken-stick fit to landmarks 0–18",
          "Change point c (dotted) and its 95% CI (shaded)"]
    axB.legend(hB, lB, loc="upper right")
    boxed(axB, 0.015, 0.03, f"Plateau end c = {r(pl['c'], 2)} months (95% CI {r(pl['c_ci'][0], 2)}–{r(pl['c_ci'][1], 2)})\n"
          f"Slope after c: {r(pl['b1_per_month'], 2)} pp per month ({r(pl['b1_ci'][0], 2)} to {r(pl['b1_ci'][1], 2)})",
          ha="left", va="bottom", linespacing=1.35)
    axB.set_xlim(-0.6, 24.6); axB.set_ylim(0, 100); axB.set_xticks(range(0, 25, 3)); axB.set_yticks(range(0, 101, 20))
    axB.set_ylabel("12-month probability of cancer death (%)"); axB.set_xlabel("Landmark (months after diagnosis)")
    risk_table(fig, axB, [("Distant", COL["NEC"], {x[i]: m["n_at_risk"][i] for i in range(len(x))})], list(range(0, 25, 3)), 106, 99, 70, 15, nslots=3)
    plab(fig, 92, 192, "B")
    rec["B"] = {"source": "A01.json monthly_distant (+plateau_cancer); A05_std.json original.std_curve", "plateau": pl,
                "standardized": {"landmarks": s["landmarks"], "std_curve": s["std_curve"], "covariates": s["covariates"]}}
    # ---------------- (C) 6-0 change by window width
    wins = [3, 6, 9, 12, 18, 24]; off = {"NEC": -0.45, "SCLC": -0.15, "AC": 0.15, "PANC": 0.45}; hC, lC = [], []
    axC.axvspan(11, 13, color=GREY, alpha=0.12, lw=0, zorder=0)
    axC.axhline(0, color=COL["light"], lw=0.7, ls="--", zorder=0)
    for key in ["NEC", "SCLC", "AC", "PANC"]:
        t = A02["tumours"][key]; est, lo, hi = [], [], []
        for w in wins:
            d = t["windows"][str(w)]["differences"]["6-0_cancer_specific"]; est.append(d["estimate"]); lo.append(d["ci"][0]); hi.append(d["ci"][1])
        est, lo, hi = map(np.array, (est, lo, hi)); xw = np.array(wins, float) + off[key]
        mfc = "white" if key in ("AC",) else COL[key]
        axC.errorbar(xw, est, yerr=[est - lo, hi - est], color=COL[key], marker=MK[key], ms=3.2, lw=1.0, ls=LS[key], capsize=1.4, elinewidth=0.7, mfc=mfc, mew=0.8)
        hC.append(h_line(COL[key], MK[key], LS[key], mfc=mfc)); lC.append(f"{LAB[key]} (n = {n_(t['n_distant'])}; B = {n_(t['B'])})")
        rec.setdefault("C", {})[key] = {"source": f"A02.json tumours.{key}.windows.*.differences.6-0_cancer_specific", "windows": wins, "diff_6_0": est.tolist(), "ci_lo": lo.tolist(), "ci_hi": hi.tolist(), "B": t["B"], "n": t["n_distant"]}
    axC.text(12, -16.3, "Pre-specified\n12-month window", ha="center", va="bottom", fontsize=FS, color=GREY, linespacing=1.15)
    axC.set_xticks(wins); axC.set_xlim(1.4, 25.6); axC.set_ylim(-17, 19); axC.set_yticks(range(-15, 16, 5))
    axC.set_xlabel("Width of the risk window after each landmark (months)")
    axC.set_ylabel("Change in conditional probability of cancer\ndeath, landmark 6 minus landmark 0 (pp)")
    axC.legend(hC, lC, loc="upper right", handlelength=2.6)
    plab(fig, 2, 93, "C")
    # ---------------- (D) TOST
    order = ["NEC 6 minus 0 months", "NEC 3 minus 0 months", "NEC 6 minus 3 months", "NEC 12 minus 0 months"]
    rows = sorted(A03["within_NEC"], key=lambda q: order.index(q["contrast"]))
    axD.axvspan(-5, 5, color=COL["SCLC"], alpha=0.07, lw=0, zorder=0)
    for mg, col, ls in [(5, COL["SCLC"], "-"), (3, GREY, "--")]:
        axD.axvline(-mg, color=col, lw=0.8, ls=ls, zorder=1); axD.axvline(mg, color=col, lw=0.8, ls=ls, zorder=1)
    axD.axvline(0, color=COL["light"], lw=0.6, zorder=0)
    tr = mtrans.blended_transform_factory(axD.transAxes, axD.transData)
    xe, x5, x3 = 1.035, (163 - 94) / 38, (171.5 - 94) / 38
    axD.text(xe, 0.62, "Estimate (95% CI), pp", transform=tr, fontsize=FS, fontweight="bold", va="bottom")
    axD.text((x5 + x3) / 2 + 0.07, 0.98, "TOST P", transform=tr, fontsize=FS, fontweight="bold", va="bottom", ha="center")
    axD.text(x5, 0.62, "±5 pp", transform=tr, fontsize=FS, fontweight="bold", va="bottom"); axD.text(x3, 0.62, "±3 pp", transform=tr, fontsize=FS, fontweight="bold", va="bottom")
    recD = []
    for i, q in enumerate(rows):
        y = -i; c90, c95 = q["ci90"], q["ci95"]
        axD.plot(c95, [y, y], color=COL["NEC"], lw=0.8, solid_capstyle="butt", zorder=3)
        axD.plot(c90, [y, y], color=COL["NEC"], lw=3.2, solid_capstyle="butt", zorder=3)
        axD.plot(q["estimate"], y, marker="o", color="white", mec=COL["NEC"], ms=4.2, mew=1.0, zorder=5)
        p5 = [mm for mm in q["margins"] if mm["margin"] == 5][0]["p_tost"]; p3 = [mm for mm in q["margins"] if mm["margin"] == 3][0]["p_tost"]
        axD.text(xe, y, ci(q["estimate"], c95[0], c95[1], d=2), transform=tr, fontsize=FS, va="center")
        axD.text(x5, y, phu(p5), transform=tr, fontsize=FS, va="center"); axD.text(x3, y, phu(p3), transform=tr, fontsize=FS, va="center")
        recD.append({"contrast": q["contrast"], "estimate": q["estimate"], "ci90": c90, "ci95": c95, "p_tost_5": p5, "p_tost_3": p3,
                     "printed": [ci(q["estimate"], c95[0], c95[1], d=2), phu(p5), phu(p3)]})
    axD.set_yticks([0, -1, -2, -3]); axD.set_yticklabels([q["contrast"].replace("NEC ", "").replace(" minus ", " − ") for q in rows])
    axD.set_ylim(-3.55, 0.6); axD.set_xlim(-22, 7.5); axD.set_xticks([-20, -15, -10, -5, 0, 5])
    axD.set_xlabel("Difference between landmarks (pp)")
    axD.spines["left"].set_visible(False); axD.tick_params(axis="y", length=0)
    hD = [Line2D([], [], color=COL["NEC"], lw=3.2), Line2D([], [], color=COL["NEC"], lw=0.8),
          Line2D([], [], marker="o", ls="none", mfc="white", mec=COL["NEC"], ms=4.2, mew=1.0),
          (Patch(facecolor=COL["SCLC"], alpha=0.15, lw=0), Line2D([], [], color=COL["SCLC"], lw=0.8)), Line2D([], [], color=GREY, lw=0.8, ls="--")]
    lD = ["90% CI (TOST criterion)", "95% CI", "Point estimate", "±5 pp margin (primary)", "±3 pp margin"]
    fig.legend(hD, lD, loc="lower left", bbox_to_anchor=(78 / 180, 1.0 / 194), ncol=3, handlelength=2.0, columnspacing=1.2, fontsize=FS)
    plab(fig, 77, 93, "D")
    rec["D"] = {"source": "A03.json within_NEC (seed 20261007, B = 2,000)", "rows": recD, "margins_status": A03["margins_status"]}
    save(fig, "Figure1")
    MANIFEST["figures"]["Figure1"] = {"A": {k: rec[k] for k in ["Distant", "Regional", "Localized"]}, "B": rec["B"], "C": rec["C"], "D": rec["D"]}


# =============================================================================================== FIGURE 2
def figure2():
    fig = newfig(180, 138)
    axA = axmm(fig, 16, 85, 66, 45); axB = axmm(fig, 105, 85, 71, 45); axC = axmm(fig, 30, 12, 64, 40)
    rec = {}
    # ---------------- (A) crude rates
    cr = A06["reproduction"]["crude_rates"]; recA = []
    for q in cr:
        a, b = [float(v) for v in q["interval"].split("-")]
        lo, hi = poisson_exact(q["deaths"], q["person_months"])
        axA.fill_between([a, b], [lo, lo], [hi, hi], color=COL["NEC"], alpha=0.18, lw=0, zorder=2)
        axA.plot([a, b], [q["rate"], q["rate"]], color=COL["NEC"], lw=1.8, solid_capstyle="butt", zorder=3)
        recA.append({"interval": q["interval"], "deaths": q["deaths"], "person_months": q["person_months"], "rate": q["rate"], "exact_ci": [lo, hi]})
    rr69 = A06["reproduction"]["adj_RR"]["r69_36"]
    hol3 = [q for q in A06["multiplicity"]["three_posthoc"] if q["contrast"] == "6-9 vs 3-6"][0]["p_holm"]
    hol21 = [q for q in A06["multiplicity"]["pairwise21"] if q["contrast"] == "6-9 vs 3-6"][0]["p_holm"]
    mono = A06["multiplicity"]["monotone_nonincrease_after_month3"]["p_boot_exact"]
    axA.annotate("", xy=(7.5, cr[2]["rate"] + 0.006), xytext=(4.5, cr[1]["rate"] + 0.006), arrowprops=dict(arrowstyle="-|>", lw=0.7, color=COL["black"], mutation_scale=6))
    axA.text(0.015, 0.985, f"6–9 vs 3–6 months: adjusted RR {rr(rr69['RR'], rr69['lo'], rr69['hi'])}, P = {r(rr69['p'], 4)}\n"
             f"Holm-adjusted P = {r(hol3, 4)} (3 post hoc contrasts), {r(hol21, 3)} (all 21)\n"
             f"Non-increasing rate after month 3 rejected: P = {r(mono, 3)}",
             transform=axA.transAxes, fontsize=FS, va="top", ha="left", linespacing=1.35)
    axA.text(0.985, 0.03, "Bars: rates; boxes: exact Poisson 95% CI", transform=axA.transAxes, fontsize=FS, ha="right", va="bottom", color=GREY)
    axA.set_xlim(0, 24); axA.set_ylim(0, 0.28); axA.set_xticks([0, 3, 6, 9, 12, 15, 18, 24]); axA.set_yticks(np.arange(0, 0.26, 0.05))
    axA.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, p: r(v, 2)))
    axA.set_xlabel("Months after diagnosis"); axA.set_ylabel("Cancer deaths per person-month")
    plab(fig, 2, 136, "A")
    rec["A"] = {"source": "A06.json reproduction.crude_rates, reproduction.adj_RR.r69_36, multiplicity.{three_posthoc,pairwise21,monotone_nonincrease_after_month3}",
                "rates": recA, "adjRR": rr69, "p_holm_3": hol3, "p_holm_21": hol21, "p_monotone_boot": mono}
    # ---------------- (B) smoothed hazards, bandwidth 1.5
    hB, lB, recB = [], [], {}
    for key in ["AC", "PANC", "SCLC", "NEC"]:
        mth, h, lo, hi, nr = grid_curve(key, "bw1.5")
        keep = (mth >= 0.25) & (mth <= 23.75)
        if key == "NEC":
            axB.fill_between(mth[keep], lo[keep], hi[keep], color=COL["NEC"], alpha=0.16, lw=0, zorder=1)
        axB.plot(mth[keep], h[keep], color=COL[key], ls=LS[key], lw=1.4 if key == "NEC" else 1.1, zorder=3)
        mk = np.isin(mth, np.arange(3, 22, 3))
        axB.plot(mth[mk], h[mk], ls="none", marker=MK[key], ms=2.8, color=COL[key], mfc=COL[key] if key != "AC" else "white", mew=0.7, zorder=4)
        c = A09["groups"][key]["curves"]["bw1.5"]
        lab = f"{LAB[key]}: {r(c['ratio_max_min'], 2)}" + (f" ({r(c['boot']['ratio_fixed_ci'][0], 2)}–{r(c['boot']['ratio_fixed_ci'][1], 2)})" if key == "NEC" else "")
        hB.append(h_line(COL[key], MK[key], LS[key], mfc=COL[key] if key != "AC" else "white", band=(key == "NEC"))); lB.append(lab)
        recB[key] = {"month": mth[keep].tolist(), "hazard": h[keep].tolist(), "bw_auto": A09["groups"][key]["bw_auto"], "trough_month": c["min_month"], "trough_haz": c["min_haz"],
                     "peak_month": c["max_month"], "peak_haz": c["max_haz"], "ratio_fixed": c["ratio_max_min"]}
        if key == "NEC":
            recB[key].update({"lower95": lo[keep].tolist(), "upper95": hi[keep].tolist(), "ratio_fixed_ci": c["boot"]["ratio_fixed_ci"], "ratio_window": c["ratio_window"],
                              "ratio_window_ci": c["boot"]["ratio_window_ci"], "trough_month_ci": c["boot"]["min_month_ci"], "peak_month_ci": c["boot"]["max_month_ci"]})
    hB, lB = hB[::-1], lB[::-1]
    cN = A09["groups"]["NEC"]["curves"]["bw1.5"]
    axB.plot([cN["min_month"], cN["max_month"]], [cN["min_haz"], cN["max_haz"]], ls="none", marker="o", ms=5.5, mfc="none", mec=COL["black"], mew=0.9, zorder=6)
    axB.text(cN["min_month"], cN["min_haz"] - 0.013, "trough", ha="center", va="top", fontsize=FS, color=COL["black"], zorder=6)
    axB.text(cN["max_month"], cN["max_haz"] + 0.03, "maximum", ha="center", va="bottom", fontsize=FS, color=COL["black"], zorder=6)
    axB.legend(hB, lB, loc="upper right", title="Peak/trough hazard ratio at fixed months", alignment="left", handlelength=2.6)
    axB.set_xlim(0, 24); axB.set_ylim(0, 0.28); axB.set_xticks(range(0, 25, 3)); axB.set_yticks(np.arange(0, 0.26, 0.05))
    axB.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, p: r(v, 2)))
    axB.set_xlabel("Months after diagnosis"); axB.set_ylabel("Smoothed hazard of cancer death\n(per person-month)")
    ftext(fig, 105, 72, f"Colorectal NEC: trough at {r(cN['min_month'], 0)} months ({r(cN['min_haz'], 3)} per person-month),\n"
          f"maximum at {r(cN['max_month'], 2)} months ({r(cN['max_haz'], 3)}); mean hazard, months 6–9 / 3–6\n"
          f"(trapezoid rule): {rr(cN['ratio_window'], *cN['boot']['ratio_window_ci'])}. Automatic bandwidths: NEC {r(A09['groups']['NEC']['bw_auto'], 2)},\n"
          f"SCLC {r(A09['groups']['SCLC']['bw_auto'], 2)}, colorectal adenocarcinoma {r(A09['groups']['AC']['bw_auto'], 2)}, pancreatic\n"
          f"adenocarcinoma {r(A09['groups']['PANC']['bw_auto'], 2)} months (all curves here: 1.5 months).",
          va="top", ha="left", linespacing=1.35)
    plab(fig, 91, 136, "B")
    rec["B"] = {"source": "A09_hazard_grid.csv bw_label bw1.5 (months 0.25-23.75); A09.json groups.*.curves.bw1.5", "curves": recB}
    # ---------------- (C) chemotherapy-stratum RRs
    chi = A06["reproduction"]["chemo_interaction"]; ivs = ["6-9", "9-15", "15-24"]; recC = []; ypos = []
    tr = mtrans.blended_transform_factory(axC.transAxes, axC.transData)
    styles = [("chemo_yes", COL["NEC"], "o", COL["NEC"], f"Chemotherapy recorded (n = {chi['n_chemo']})"),
              ("chemo_no_or_unknown", GREY, "s", "white", f"No/unknown chemotherapy (n = {chi['n_nochemo']})")]
    for j, iv in enumerate(ivs):
        for k, (grp, col, mk, mfc, lab) in enumerate(styles):
            e = chi[f"per_{iv}_vs_3-6"][grp]; y = -(j * 3 + k)
            axC.plot([e["lo"], e["hi"]], [y, y], color=col, lw=0.9, zorder=3)
            axC.plot(e["RR"], y, marker=mk, color=col, ms=4, mfc=mfc, mew=0.9, zorder=4)
            axC.text(1.03, y, rr(e["RR"], e["lo"], e["hi"]), transform=tr, fontsize=FS, va="center", ha="left")
            recC.append({"interval": iv, "group": grp, "RR": e["RR"], "lo": e["lo"], "hi": e["hi"], "p": e["p"], "printed": rr(e["RR"], e["lo"], e["hi"])})
        ypos.append(-(j * 3 + 0.5))
    axC.text(1.03, 0.75, "RR (95% CI)", transform=tr, fontsize=FS, fontweight="bold", va="center")
    axC.axvline(1, color=COL["light"], lw=0.7, ls="--", zorder=0)
    axC.set_xscale("log"); axC.set_xlim(0.1, 4.2); axC.set_xticks([0.125, 0.25, 0.5, 1, 2, 4]); axC.set_xticklabels(["0.125", "0.25", "0.5", "1", "2", "4"])
    axC.xaxis.set_minor_locator(mticker.NullLocator())
    axC.set_yticks(ypos); axC.set_yticklabels([f"{iv.replace('-', '–')} months" for iv in ivs]); axC.set_ylim(-7.6, 1.3)
    axC.set_xlabel("Rate ratio vs 3–6 months (log scale)")
    axC.spines["left"].set_visible(False); axC.tick_params(axis="y", length=0)
    hC = [Line2D([], [], color=c_, marker=m_, ms=4, mfc=f_, mew=0.9, lw=0.9) for _, c_, m_, f_, _ in styles]
    fig.legend(hC, [s_[4] for s_ in styles], loc="upper left", bbox_to_anchor=(126 / 180, 47 / 138), handlelength=1.8)
    ig = chi["interaction_global"]
    ftext(fig, 127, 35, f"Patients alive at 3 months: n = {chi['n_survivors_3m']}\nReference interval: 3–6 months\n"
          f"Chemotherapy × interval interaction:\nW = {r(ig['W'], 2)} on {ig['df']} df, P = {r(ig['p'], 3)}", va="top", ha="left", linespacing=1.4)
    plab(fig, 2, 58, "C")
    rec["C"] = {"source": "A06.json reproduction.chemo_interaction", "RR": recC, "interaction": ig, "n": [chi["n_survivors_3m"], chi["n_chemo"], chi["n_nochemo"]]}
    save(fig, "Figure2")
    MANIFEST["figures"]["Figure2"] = rec


# =============================================================================================== FIGURE 3
def figure3():
    fig = newfig(180, 178)
    axA = axmm(fig, 33, 116, 56, 54); axB = axmm(fig, 122, 94, 27, 76)
    axC = axmm(fig, 20, 36, 62, 42); axD = axmm(fig, 123, 36, 28, 42)
    rec = {}
    # ---------------- (A) four-tumour monthly curves
    rows, hs, ls_ = [], [], []; recA = {}
    for key in ["NEC", "SCLC", "AC", "PANC"]:
        t = A11["table3"][key]["cancer"]; m = t["monthly"]; x = m["landmarks"]
        curve(axA, x, m["pct"], m["n_at_risk"], COL[key], MK[key], ls=LS[key], lo=m["lo"], hi=m["hi"], every=3, ms=2.8, lw=1.2 if key == "NEC" else 1.0)
        rows.append((LAB[key], COL[key], {x[i]: m["n_at_risk"][i] for i in range(len(x))}))
        hs.append(h_line(COL[key], MK[key], LS[key], band=(key == "NEC"))); ls_.append(f"{LAB[key]} (n = {n_(t['n'])})")
        recA[key] = {"landmarks": x, "pct": m["pct"], "lo": m["lo"], "hi": m["hi"], "n_at_risk": m["n_at_risk"], "plateau_end": t["plateau_end"]}
    axA.set_xlim(-0.6, 24.6); axA.set_ylim(0, 100); axA.set_xticks(range(0, 25, 3)); axA.set_yticks(range(0, 101, 20))
    axA.set_ylabel("12-month probability of cancer death (%)"); axA.set_xlabel("Landmark (months after diagnosis)")
    lg = axA.legend(hs, ls_, loc="lower left", handlelength=2.6); axA.add_artist(lg)
    hh, ll = h_open_filled(); axA.legend(hh, ll, loc="upper left", handlelength=1.0)
    pe = {k: mo(A11["table3"][k]["cancer"]["plateau_end"]["months"]) for k in ["NEC", "SCLC", "AC", "PANC"]}
    axA.text(0.985, 0.985, f"Broken-stick plateau end (months):\ncolorectal NEC {pe['NEC']}; SCLC {pe['SCLC']};\ncolorectal adenocarcinoma {pe['AC']};\npancreatic adenocarcinoma {pe['PANC']}",
             transform=axA.transAxes, ha="right", va="top", fontsize=FS, linespacing=1.3)
    risk_table(fig, axA, rows, [0, 6, 12, 18, 24], 33, 89, 56, 17)
    plab(fig, 2, 176, "A")
    rec["A"] = {"source": "A11.json table3.{NEC,SCLC,AC,PANC}.cancer.monthly / plateau_end", "curves": recA, "plateau_end_printed": pe}
    # ---------------- (B) NEC minus comparator
    contrasts = ["6-0", "12-6", "12-0"]; comps = ["SCLC", "LCNEC", "AC", "PANC"]; recB = []
    y = 0; yt, yl = [], []
    tr = mtrans.blended_transform_factory(axB.transAxes, axB.transData)
    short = {"SCLC": "SCLC", "LCNEC": "Lung LCNEC", "AC": "Colorectal\nadenocarcinoma", "PANC": "Pancreatic\nadenocarcinoma"}
    for c in contrasts:
        axB.text(-1.05, y + 0.45, f"Change {c.replace('-', ' − ')} months", transform=tr, fontsize=FS, fontweight="bold", va="bottom", ha="left", clip_on=False)
        for key in comps:
            src = A11["lcnec"] if key == "LCNEC" else A11["table3"][key]
            d = src["nec_minus_comparator_cancer"][c]
            axB.plot(d["ci"], [y, y], color=COL[key], lw=1.0, zorder=3)
            axB.plot(d["diff_pct_points"], y, marker=MK[key], color=COL[key], ms=3.8, mew=0.7, zorder=4)
            axB.text(1.04, y, ci(d["diff_pct_points"], d["ci"][0], d["ci"][1], d=2), transform=tr, fontsize=FS, va="center")
            yt.append(y); yl.append(short[key])
            recB.append({"contrast": c, "comparator": key, "diff": d["diff_pct_points"], "ci": d["ci"], "boot_p": d.get("boot_p_two_sided"),
                         "printed": ci(d["diff_pct_points"], d["ci"][0], d["ci"][1], d=2),
                         "source": f"A11.json {'lcnec' if key == 'LCNEC' else 'table3.' + key}.nec_minus_comparator_cancer.{c}"})
            y -= 1.1
        y -= 1.0
    axB.text(1.04, 1.05, "Difference, pp (95% CI)", transform=tr, fontsize=FS, fontweight="bold", va="bottom")
    axB.axvline(0, color=COL["light"], lw=0.7, ls="--", zorder=0)
    axB.set_yticks(yt); axB.set_yticklabels(yl, fontsize=6.2, linespacing=1.0); axB.set_ylim(y + 1.4, 1.3); axB.set_xlim(-20, 13.5)
    axB.set_xticks([-10, 0, 10]); axB.set_xticks([-15, -5, 5], minor=True); axB.tick_params(axis="x", which="minor", length=1.8, width=0.5)
    axB.set_xlabel("NEC minus comparator (pp)\npositive = NEC declines less")
    axB.spines["left"].set_visible(False); axB.tick_params(axis="y", length=0)
    plab(fig, 92, 176, "B")
    rec["B"] = {"rows": recB, "B": A11["B"], "seed_rule": A11["seed_rule"]}
    # ---------------- (C) 2 x 2
    cells = A11["two_by_two"]["cells"]; recC = {}
    pos = {("colorectal", "SC"): -0.17, ("colorectal", "LC"): 0.17, ("lung", "SC"): 0.83, ("lung", "LC"): 1.17}
    style = {"SC": (COL["NEC"], "o", "Small-cell morphology"), "LC": (COL["LCNEC"], "s", "Large-cell morphology (8013)")}
    names = {"colorectal_SC": "colorectal small-cell NEC", "colorectal_LC": "colorectal large-cell NEC", "lung_SC": "SCLC", "lung_LC": "lung LCNEC"}
    for cname, cell in cells.items():
        organ, morph = cname.split("_"); d = cell["cancer"]["change"]["6-0"]; n = cell["cancer"]["n"]; col, mk, _ = style[morph]; xx = pos[(organ, morph)]
        axC.plot([xx, xx], d["ci"], color=col, lw=1.1, zorder=3)
        axC.plot(xx, d["pct_points"], marker=mk, color=col, ms=5, mew=0.8, zorder=4)
        boxed(axC, xx + (0.06 if morph == "LC" else -0.06), d["pct_points"], f"{sg(d['pct_points'], 2)}\nn = {n_(n)}", transform=axC.transData,
              ha="left" if morph == "LC" else "right", va="center", linespacing=1.2, zorder=5)
        recC[cname] = {"label": names[cname], "n": n, "diff_6_0": d["pct_points"], "ci": d["ci"], "plateau_end": cell["cancer"]["plateau_end"]["months"]}
    axC.axhline(0, color=COL["light"], lw=0.7, ls="--", zorder=0)
    axC.set_xticks([0, 1]); axC.set_xticklabels(["Colorectal", "Lung"], fontsize=7); axC.set_xlim(-0.65, 1.65); axC.set_ylim(-13, 10); axC.set_yticks(range(-12, 9, 4))
    axC.set_ylabel("Change 6 − 0 months in 12-month\nprobability of cancer death (pp)")
    axC.legend([Line2D([], [], color=style[k][0], marker=style[k][1], ms=5, lw=1.1) for k in ["SC", "LC"]], [style[k][2] for k in ["SC", "LC"]], loc="upper right", handlelength=1.6)
    tb = A11["two_by_two"]; lsl = tb["paired"]["lungSC_minus_lungLC"]["cancer"]["6-0"]
    ftext(fig, 20, 26, "Effects on the 6 − 0 change, pp (95% CI); main effects are\naveraged over the other factor\n"
          f"Morphology (small- minus large-cell): {ci(tb['morphology_effect_6_0']['estimate'], *tb['morphology_effect_6_0']['ci'])}\n"
          f"Organ (colorectal minus lung): {ci(tb['organ_effect_6_0']['estimate'], *tb['organ_effect_6_0']['ci'])}\n"
          f"Morphology × organ interaction: {ci(tb['interaction_6_0']['estimate'], *tb['interaction_6_0']['ci'])}\n"
          f"Lung small-cell minus large-cell: {ci(lsl['diff_pct_points'], *lsl['ci'])}",
          va="top", ha="left", linespacing=1.4)
    plab(fig, 2, 86, "C")
    rec["C"] = {"source": "A11.json two_by_two.cells.*.cancer.change.6-0, two_by_two.{morphology_effect_6_0,organ_effect_6_0,interaction_6_0}, two_by_two.paired.lungSC_minus_lungLC",
                "cells": recC, "effects": {k: tb[k] for k in ["morphology_effect_6_0", "organ_effect_6_0", "interaction_6_0"]}, "lungSC_minus_lungLC": lsl}
    # ---------------- (D) RR 6-9 vs 3-6 by tumour
    # A12 re-runs NEC and SCLC through the same code path as lung LCNEC; guard that it reproduces A06 before plotting LCNEC
    for k_, a6 in [("NEC", A06["reproduction"]["adj_RR"]["r69_36"]), ("SCLC", A06["comparators_histology"]["SCLC"]["r69_36"])]:
        a12 = A12["comparators_same_code_path"][k_]["rr_6_9_vs_3_6"]
        assert all(abs(a12[q] - a6[q]) < 1e-6 for q in ("RR", "lo", "hi")), f"A12 does not reproduce A06 for {k_}"
    st12 = A12["shape_tests"]
    assert abs(st12["NEC_vs_SCLC_all_years"]["p"] - A06["calendar_period"]["shape_tests"]["SCLC"]["overall_linear_year"]["p"]) < 1e-6
    items = [("NEC", A06["reproduction"]["adj_RR"]["r69_36"], A06["n"]["NEC_distant"]), ("SCLC", A06["comparators_histology"]["SCLC"]["r69_36"], A06["comparators_histology"]["SCLC"]["n"]),
             ("LCNEC", A12["lcnec"]["rr_6_9_vs_3_6"], A12["lcnec"]["n"]),
             ("AC", A06["comparators_histology"]["AC"]["r69_36"], A06["comparators_histology"]["AC"]["n"]), ("PANC", A06["comparators_histology"]["PANC"]["r69_36"], A06["comparators_histology"]["PANC"]["n"])]
    trD = mtrans.blended_transform_factory(axD.transAxes, axD.transData); recD = []
    for i, (key, e, n) in enumerate(items):
        y = -i
        axD.plot([e["lo"], e["hi"]], [y, y], color=COL[key], lw=1.1, zorder=3)
        axD.plot(e["RR"], y, marker=MK[key], color=COL[key], ms=4.2, mfc=COL[key], mew=0.8, zorder=4)
        axD.text(1.05, y, rr(e["RR"], e["lo"], e["hi"]), transform=trD, fontsize=FS, va="center")
        recD.append({"tumour": key, "n": n, "RR": e["RR"], "lo": e["lo"], "hi": e["hi"], "p": e["p"], "printed": rr(e["RR"], e["lo"], e["hi"]),
                     "source": "A12.json lcnec.rr_6_9_vs_3_6" if key == "LCNEC" else ("A06.json reproduction.adj_RR.r69_36" if key == "NEC" else f"A06.json comparators_histology.{key}.r69_36")})
    axD.text(1.05, 0.55, "RR (95% CI)", transform=trD, fontsize=FS, fontweight="bold", va="bottom")
    axD.axvline(1, color=COL["light"], lw=0.7, ls="--", zorder=0)
    axD.set_xscale("log"); axD.set_xlim(0.72, 2.2); axD.set_xticks([0.75, 1, 1.5, 2]); axD.set_xticklabels(["0.75", "1", "1.5", "2"]); axD.xaxis.set_minor_locator(mticker.NullLocator())
    axD.set_yticks([-i for i in range(len(items))]); axD.set_yticklabels([f"{LAB[k]}\n(n = {n_(n)})" for k, _, n in items], linespacing=1.1); axD.set_ylim(-len(items) + 0.5, 0.55)
    axD.set_xlabel("Adjusted RR, 6–9 vs 3–6 months\n(log scale)")
    axD.spines["left"].set_visible(False); axD.tick_params(axis="y", length=0)
    st = A06["calendar_period"]["shape_tests"]; sh = {k: st[k]["overall_linear_year"] for k in ["SCLC", "AC", "PANC"]}

    def pp_(p):
        return "P\u00a0<\u00a00.001" if p < 0.001 else f"P\u00a0=\u00a0{r(p, 3)}"
    pa, pp2 = pp_(sh["AC"]["p"]), pp_(sh["PANC"]["p"])
    both = f"colorectal and pancreatic adenocarcinoma {pa}" if pa == pp2 else f"colorectal adenocarcinoma {pa}; pancreatic adenocarcinoma {pp2}"
    ftext(fig, 98, 21.5, textwrap.fill("Piecewise-exponential Poisson models with cluster-robust SE; adjusted for age/10, its square, sex and calendar year "
          "(+ rectal site for colorectal adenocarcinoma; + morphology and rectal site for NEC). Hazard-shape test vs NEC (6-df Wald, pooled model, all years): "
          f"SCLC {pp_(sh['SCLC']['p'])}; {both}; lung LCNEC {pp_(st12['NEC_vs_LCNEC_all_years']['p'])} "
          "(exploratory, not adjusted for multiplicity).", 70), va="top", ha="left", linespacing=1.4, color=COL["black"])
    plab(fig, 92, 86, "D")
    lt = st12["NEC_vs_LCNEC_all_years"]
    rec["D"] = {"source": "A06.json reproduction.adj_RR.r69_36, comparators_histology.*.r69_36, calendar_period.shape_tests.*.overall_linear_year; "
                          "A12.json lcnec.rr_6_9_vs_3_6, shape_tests.NEC_vs_LCNEC_all_years (A12 reproduces A06 NEC/SCLC RRs and the NEC vs SCLC test; asserted)",
                "RR": recD,
                "shape_tests": {**{k: {"W": sh[k]["W"], "df": sh[k]["df"], "p": sh[k]["p"]} for k in sh},
                                "LCNEC": {"W": lt["W"], "df": lt["df"], "p": lt["p"], "n_nec": lt["n_nec"], "n_cmp": lt["n_cmp"], "note": "exploratory; added at revision; not adjusted for multiplicity"}}}
    save(fig, "Figure3")
    MANIFEST["figures"]["Figure3"] = rec


# =============================================================================================== SUPPLEMENTARY FIGURE S1
def suppS1():
    fl = {}
    for q in flow_rows:
        fl.setdefault(q["cohort"], {})[q["step"]] = (int(q["records"]), int(q["patients"]))
    nen, lung, ac, panc = fl["NEN_all"], fl["LUNG"], fl["AC"], fl["PANC"]
    H = 170
    fig = newfig(180, H)
    ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 180); ax.set_ylim(0, H); ax.axis("off")

    def box(x, y, w, h, text, fc="white", ec=COL["black"], lw=0.7, fs=6.2, bold=False, color=COL["black"]):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0,rounding_size=1.0", fc=fc, ec=ec, lw=lw))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs, color=color, fontweight="bold" if bold else "normal", linespacing=1.3)

    def arrow(x1, y1, x2, y2):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>", mutation_scale=6.5, lw=0.7, color=COL["black"], shrinkA=0, shrinkB=0))

    plab(fig, 2, H - 2, "A", "Colorectal NEC cohort")
    X0, BW, BH, XE, WE = 3, 64, 7.0, 71, 52
    yb = H - 26; box(X0, yb, BW, 15, f"SEER 17 colorectal neuroendocrine neoplasm records\nICD-O-3 8013, 8041–8045, 8240–8246, 8249;\n"
                     f"diagnoses 2000–2023\nn = {n_(nen['exported'][0])} records ({n_(nen['exported'][1])} patients)", fc="#EEF3F8", bold=True)
    steps = [("age>=18", "Age < 18 years"), ("site", "Primary site not C18.0, C18.2–C18.9,\nC19.9 or C20.9"), ("diag_confirmation", "No positive histology"),
             ("intl_primary", "Not a primary by international rules"), ("not_autopsy_dco", "Autopsy- or death-certificate-only"), ("one_per_patient", "Further record of the same patient")]
    prev = nen["malignant"][0]; ytop = yb; rec_steps = []
    for key, lab in steps:
        n_now = nen[key][0]; excl = prev - n_now
        if key == "one_per_patient":
            excl = nen["survival_known"][0] - nen["one_per_patient"][1]; n_now = nen["one_per_patient"][1]
        ynext = ytop - 13.5
        arrow(X0 + BW / 2, ytop, X0 + BW / 2, ynext + BH)
        ym = ynext + BH + (ytop - ynext - BH) / 2; hx = 7.2 + 2.6 * lab.count("\n")
        box(XE, ym - hx / 2, WE, hx, f"Excluded: {lab}\nn = {n_(excl)}", fc="#FAFAFA", ec=GREY, lw=0.5, color=COL["black"])
        ax.plot([X0 + BW / 2, XE], [ym, ym], color=GREY, lw=0.5)
        box(X0, ynext, BW, BH, f"n = {n_(n_now)}" + (" patients" if key == "one_per_patient" else " records"), fs=6.5)
        rec_steps.append({"step": key, "excluded": excl, "remaining": n_now}); ytop = ynext; prev = n_now
    ynext = ytop - 17; BHN = 11
    arrow(X0 + BW / 2, ytop, X0 + BW / 2, ynext + BHN)
    ym = ynext + BHN + (ytop - ynext - BHN) / 2
    box(XE, ym - 3.6, WE, 7.2, f"Excluded: other NEN histology (8240–8246, 8249)\nn = {n_(FLOWC['NEN_all_n'] - FLOWC['NEC_main_n'])}", fc="#FAFAFA", ec=GREY, lw=0.5)
    ax.plot([X0 + BW / 2, XE], [ym, ym], color=GREY, lw=0.5)
    n_sc = sum(v for k, v in FLOWC["NEC_hist"].items() if k != "8013")
    box(X0, ynext, BW, BHN, f"Colorectal NEC cohort (8013, 8041–8045)\nn = {n_(FLOWC['NEC_main_n'])} patients\nsmall-cell family 8041–8045: {n_(n_sc)}; large-cell 8013: {n_(FLOWC['NEC_hist']['8013'])}",
        fc="#DCE9F5", bold=True, ec=COL["NEC"], lw=1.0)
    ytop = ynext; st = FLOWC["NEC_stage"]; yst = ytop - 19; wst = 27.5
    for i, (lab, key, fc, ec) in enumerate([("Distant\n(primary analysis)", "Distant", "#C9DCF0", COL["NEC"]), ("Regional", "Regional", "#FFFFFF", COL["black"]),
                                           ("Localized", "Localized", "#FFFFFF", COL["black"]), ("Unknown stage", "Unknown", "#FFFFFF", GREY)]):
        xx = X0 + i * (wst + 2.3)
        ax.plot([X0 + BW / 2, xx + wst / 2], [ytop, yst + 10], color=COL["black"], lw=0.6)
        box(xx, yst, wst, 10, f"{lab}\nn = {n_(st[key])}", fc=fc, ec=ec, lw=1.0 if key == "Distant" else 0.6, bold=(key == "Distant"))
    box(X0, yst - 15, 4 * wst + 3 * 2.3, 10, "Expanded NEC definition (pre-specified sensitivity analysis): main cohort + 8246/3 with registry grade 3 or 4\n"
        f"(n = {n_(FLOWC['n_8246_g34'])} of {n_(FLOWC['n_8246_total'])} such records); expanded cohort n = {n_(FLOWC['NEC_expanded_n'])} patients, distant n = {n_(FLOWC['NEC_expanded_distant'])}",
        fc="#FAFAFA", ec=GREY, lw=0.6)
    ax.text(X0, yst - 19, "Steps that excluded no record (malignant behavior; diagnosis 2000–2023; survival days known) are not drawn.",
            fontsize=FS, color=GREY, ha="left", va="top")
    # ---- (B)
    XR, WR = 127, 51
    plab(fig, 125, H - 2, "B", "Comparator cohorts (same filters)")
    yy = H - 26
    box(XR, yy, WR, 15, f"Lung export\n(ICD-O-3 8041–8045, 8013, 8002, 8040)\nn = {n_(lung['exported'][0])} records ({n_(lung['exported'][1])} patients)\n"
        f"after filters: n = {n_(lung['one_per_patient'][1])} patients", fc="#FAFAFA", ec=GREY, lw=0.6)
    hw = (WR - 2) / 2
    arrow(XR + hw / 2, yy, XR + hw / 2, yy - 5); arrow(XR + WR - hw / 2, yy, XR + WR - hw / 2, yy - 5)
    yy -= 5 + 12
    box(XR, yy, hw, 12, f"SCLC (8041–8045)\nn = {n_(FLOWC['SCLC_n'])}\ndistant: n = {n_(FLOWC['SCLC_distant'])}", ec=COL["SCLC"], lw=1.1, bold=True)
    box(XR + WR - hw, yy, hw, 12, f"Lung LCNEC (8013)\nn = {n_(FLOWC['LCNEC_n'])}\ndistant: n = {n_(FLOWC['LCNEC_distant'])}", ec=COL["LCNEC"], lw=1.1, bold=True)
    n_other = FLOWC["LUNG_hist"].get("8002", 0) + FLOWC["LUNG_hist"].get("8040", 0)
    ax.text(XR, yy - 1.5, f"8002 and 8040 (n = {n_other}): used only in a sensitivity\ndefinition of SCLC", fontsize=FS, color=GREY, ha="left", va="top", linespacing=1.3)
    yy -= 17
    for key, top, h, bot, col in [
        ("AC", f"Colorectal adenocarcinoma export\nn = {n_(ac['exported'][0])} records ({n_(ac['exported'][1])} patients)\nafter filters: n = {n_(ac['one_per_patient'][1])} patients", 12,
         f"Distant: n = {n_(FLOWC['AC_distant'])}", COL["AC"]),
        ("PANC", f"Pancreatic adenocarcinoma export\nn = {n_(panc['exported'][0])} records ({n_(panc['exported'][1])} patients)\nafter filters: n = {n_(panc['one_per_patient'][1])} patients;\n"
                 f"excluding islet-cell codes 8150–8157:\nn = {n_(FLOWC['PANC_excl_islet_n'])}", 17.5, f"Distant: n = {n_(FLOWC['PANC_distant'])}", COL["PANC"])]:
        box(XR, yy - h, WR, h, top, fc="#FAFAFA", ec=GREY, lw=0.6)
        arrow(XR + WR / 2, yy - h, XR + WR / 2, yy - h - 5)
        box(XR, yy - h - 5 - 7, WR, 7, bot, ec=col, lw=1.1, bold=True, fs=6.5)
        yy -= h + 5 + 7 + 7
    ax.text(XR, 36, textwrap.fill("Filters applied in sequence to every cohort: malignant behavior; age ≥ 18 years; primary site; positive histology; "
            "primary by international rules; not autopsy- or death-certificate-only; diagnosis 2000–2023; survival days known; one record per patient. "
            "Patients in both the NEC and a comparator cohort (SCLC 1, colorectal adenocarcinoma 3) were retained.", 46), ha="left", va="top", fontsize=FS, color=GREY, linespacing=1.3)
    save(fig, "SuppFigS1")
    MANIFEST["figures"]["SuppFigS1"] = {"source": "data/flow_all.csv (aggregate step counts), FIG_flow_counts.json", "A_steps": rec_steps, "NEN_flow": nen, "counts": FLOWC,
                                        "comparator_flow": {"LUNG": lung, "AC": ac, "PANC": panc}, "lung_8002_8040": n_other}


# =============================================================================================== SUPPLEMENTARY FIGURE S2
def suppS2():
    fig = newfig(180, 122)
    ax = axmm(fig, 52, 17, 88, 101)
    pooled = A05["pooled"]; pd = pooled["diff_6_0"]
    ax.axvspan(pd["ci"][0], pd["ci"][1], color=COL["NEC"], alpha=0.08, lw=0, zorder=0)
    ax.axvline(pd["est"], color=COL["NEC"], lw=0.8, ls="--", zorder=0); ax.axvline(0, color=COL["light"], lw=0.6, zorder=0)
    tr = mtrans.blended_transform_factory(ax.transAxes, ax.transData)
    rec = []; y = 0; yt, yl = [], []
    treat_vars = {"Surgery (primary site)", "Chemotherapy recorded"}
    groups = []
    for s in A05["strata"]:
        if not groups or groups[-1][0] != s["variable"]:
            groups.append((s["variable"], []))
        groups[-1][1].append(s)
    ax.plot(pd["ci"], [y, y], color=COL["NEC"], lw=1.4); ax.plot(pd["est"], y, marker="D", color=COL["NEC"], ms=4.5)
    ax.text(1.03, y, ci(pd["est"], *pd["ci"], d=2), transform=tr, fontsize=FS, va="center")
    ax.text(1.03, 1.1, "Change, pp (95% CI)", transform=tr, fontsize=FS, va="center", fontweight="bold")
    yt.append(y); yl.append(f"All distant NEC (n = {pooled['n']})"); y -= 1.7
    names = {"Liver metastasis (2010+)": "Liver metastasis (diagnoses 2010 and later)", "Liver-only pattern (2010+)": "Metastatic pattern (diagnoses 2010 and later)",
             "Histology": "Histological subtype", "Race": "Race/ethnicity"}
    for var, members in groups:
        treat = var in treat_vars
        th = ax.text(-0.57, y + 0.5, names.get(var, var), transform=tr, fontsize=FS, fontweight="bold", va="bottom", clip_on=False)
        if treat:
            rend = fig.canvas.get_renderer(); bb = th.get_window_extent(renderer=rend)
            x_end = ax.transAxes.inverted().transform((bb.x1, bb.y0))[0]
            ax.text(x_end + 0.012, y + 0.5, "(landmark 6 − 3 months, patients alive at 3 months)", transform=tr, fontsize=FS, color=COL["SCLC"], va="bottom", clip_on=False)
        for s in members:
            if treat:
                c = s["conditional_on_3mo_survival"]; d = c["diff_6_3"]; n_lbl = c["n"]; nr = c["landmark"]["6"]["n_at_risk"]
            else:
                d = s["diff_6_0"]; n_lbl = s["n"]; nr = s["landmark"]["6"]["n_at_risk"]
            col = COL["SCLC"] if treat else COL["black"]; mk = "s" if treat else "o"
            ax.plot(d["ci"], [y, y], color=col, lw=0.9); ax.plot(d["est"], y, marker=mk, color=col, ms=3.6, mfc=col if nr >= 100 else "white", mew=0.9)
            ax.text(1.03, y, ci(d["est"], *d["ci"], d=2), transform=tr, fontsize=FS, va="center")
            lev = s["level"].replace(">=65", "≥ 65").replace("<65", "< 65").replace("8041-8045", "8041–8045")
            yt.append(y); yl.append(f"{lev} (n = {n_lbl})")
            rec.append({"variable": var, "level": s["level"], "n": n_lbl, "contrast": "6-3 among 3-month survivors" if treat else "6-0", "est": d["est"], "ci": d["ci"],
                        "n_at_risk_later_landmark": nr, "printed": ci(d["est"], *d["ci"], d=2)})
            y -= 1
        y -= 0.95
    ax.set_yticks(yt); ax.set_yticklabels(yl); ax.set_ylim(y + 0.6, 1.5); ax.set_xlim(-26, 16)
    ax.set_xticks([-25, -20, -15, -10, -5, 0, 5, 10, 15]); ax.spines["bottom"].set_bounds(-26, 16)
    ax.set_xlabel("Change in 12-month probability of cancer death, landmark 6 minus landmark 0 months (pp)")
    ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
    ftext(fig, 2, 5.5, "Circles: baseline strata (6 − 0 months); orange squares: treatment strata (6 − 3 months among patients alive at 3 months). Open markers: fewer than 100 patients at risk at the\n"
          "later landmark. Diamond, dashed line and shaded band: whole-cohort estimate and 95% CI. 95% percentile intervals, patient-level bootstrap (B = 2,000 per stratum).",
          va="center", ha="left", color=GREY, linespacing=1.35)
    save(fig, "SuppFigS2")
    MANIFEST["figures"]["SuppFigS2"] = {"source": "A05.json pooled.diff_6_0; strata[*].diff_6_0 / conditional_on_3mo_survival.diff_6_3", "pooled": pooled, "strata": rec}


# =============================================================================================== SUPPLEMENTARY FIGURE S3
def suppS3():
    fig = newfig(180, 112)
    ax = axmm(fig, 27, 41, 149, 66)
    rec = {}; rows = []; hs, ls_ = [], []
    for key, col, mk, ls in [("derivation_2000_2015", COL["NEC"], "o", "-"), ("confirmation_2016_2023", COL["SCLC"], "s", "--")]:
        p = A07[key]; x = np.arange(25); yv = np.array(p["curve_0_24"], float); n = np.array(p["nrisk_0_24"], float); keep = n >= 20
        curve(ax, x[keep], yv[keep], n[keep], col, mk, ls=ls, ms=2.8)
        for lm, v in p["landmarks"].items():
            ax.errorbar(float(lm) + (0.18 if key.startswith("conf") else -0.18), v["cif12"], yerr=[[v["cif12"] - v["ci"][0]], [v["ci"][1] - v["cif12"]]], color=col, elinewidth=0.8, capsize=1.8, lw=0, zorder=2)
        nm = "Derivation" if key.startswith("der") else "Confirmation"; yrs = "2000–2015" if key.startswith("der") else "2016–2023"
        hs.append(h_line(col, mk, ls)); ls_.append(f"{nm}, {yrs} (n = {p['n']}; B = 2,000, seed {p['seed']}): plateau end {r(p['plateau_end']['est'], 2)} months ({r(p['plateau_end']['ci'][0], 2)}–{r(p['plateau_end']['ci'][1], 2)})")
        rows.append((nm, col, {int(x[i]): p["nrisk_0_24"][i] for i in range(25)}))
        rec[key] = {"n": p["n"], "seed": p["seed"], "curve_0_24": p["curve_0_24"], "nrisk_0_24": p["nrisk_0_24"], "landmarks": p["landmarks"], "change_6_minus_0": p["change_6_minus_0"], "plateau_end": p["plateau_end"]}
    hh, ll = h_open_filled()
    ax.legend(hs + hh, ls_ + ll, loc="upper right", handlelength=2.6)
    pdd = A07["period_difference_in_6_minus_0"]; d1 = A07["derivation_2000_2015"]["change_6_minus_0"]; d2 = A07["confirmation_2016_2023"]["change_6_minus_0"]
    ar = A07["adjusted_RR_R"]; r1 = ar["derivation_2000_2015"]["RR_6to9_vs_3to6"]; r2 = ar["confirmation_2016_2023"]["RR_6to9_vs_3to6"]; rq = ar["ratio_of_RR_late_vs_early"]
    boxed(ax, 0.012, 0.03, f"Change 6 − 0 months: derivation {ci(d1['est'], *d1['ci'], d=2)} pp; confirmation {ci(d2['est'], *d2['ci'], d=2)} pp;\n"
          f"difference (confirmation minus derivation) {ci(pdd['est'], *pdd['ci'], d=2)} pp, bootstrap P = {r(pdd['boot_p'], 2)}\n"
          f"Adjusted RR, 6–9 vs 3–6 months: {rr(r1['RR'], r1['lo'], r1['hi'])} in 2000–2015; {rr(r2['RR'], r2['lo'], r2['hi'])} in 2016–2023; ratio {rr(rq['RR'], rq['lo'], rq['hi'])}, P = {r(rq['p'], 3)}",
          ha="left", va="bottom", linespacing=1.4)
    ax.set_xlim(-0.6, 24.6); ax.set_ylim(0, 100); ax.set_xticks(range(0, 25, 3)); ax.set_yticks(range(0, 101, 20))
    ax.set_ylabel("12-month probability of cancer death (%)"); ax.set_xlabel("Landmark (months after diagnosis)")
    risk_table(fig, ax, rows, list(range(0, 25, 3)), 27, 15.5, 149, 14)
    ftext(fig, 27, 6.5, "Lines join point estimates; error bars: 95% bootstrap percentile intervals at the pre-specified landmarks 0, 3, 6 and 12 months;\n"
          "points with fewer than 20 patients at risk are not drawn; numbers at risk in italics: fewer than 100.", va="center", ha="left", color=GREY, linespacing=1.35)
    rec["period_difference_in_6_minus_0"] = pdd; rec["adjusted_RR"] = {"2000_2015": r1, "2016_2023": r2, "ratio": rq}
    save(fig, "SuppFigS3")
    MANIFEST["figures"]["SuppFigS3"] = {"source": "A07.json derivation_2000_2015 / confirmation_2016_2023 / period_difference_in_6_minus_0 / adjusted_RR_R", **rec}


# =============================================================================================== SUPPLEMENTARY FIGURE S4
def suppS4():
    wins = [3, 6, 9, 12, 18, 24]; t = A02["tumours"]["NEC"]; rec = {}
    fig = newfig(180, 130)
    xs = [16, 72.5, 129]; ys = [83, 34]
    for i, w in enumerate(wins):
        ax = axmm(fig, xs[i % 3], ys[i // 3], 47, 36); d = t["windows"][str(w)]
        x = np.array(d["landmark"], float); yv = np.array(d["cancer_specific"], float); n = np.array(d["n_at_risk"], float)
        c = d["plateau_end_cancer_specific"]; cci = d["plateau_end_cancer_specific_ci"]
        y0b, y1b = (0.0, 0.58) if w == 3 else (0.30, 1.0)
        ax.axvspan(cci[0], cci[1], ymin=y0b, ymax=y1b, color=COL["SCLC"], alpha=0.10, lw=0, zorder=0)
        ax.axvline(c, ymin=y0b, ymax=y1b, color=COL["SCLC"], lw=0.9, ls=":", zorder=1)
        curve(ax, x, yv, n, COL["NEC"], "o", lo=d["cancer_specific_ci_lo"], hi=d["cancer_specific_ci_hi"], ms=2.4)
        xx = np.linspace(0, 18, 200); yy = d["broken_stick_intercept"] + d["broken_stick_slope_after_c_per_month"] * np.maximum(0, xx - c)
        ax.plot(xx, yy, color=COL["SCLC"], lw=1.0, ls="-", zorder=2)
        df = d["differences"]["6-0_cancer_specific"]
        txt = f"c = {r(c, 2)} months (95% CI {r(cci[0], 2)}–{r(cci[1], 2)})\nP(c ≥ 3) = {r(d['plateau_prop_boot_ge3'], 2)}\n6 − 0 change: {ci(df['estimate'], *df['ci'], d=2)} pp"
        if w == 3:
            d30 = d["differences"]["3-0_cancer_specific"]; d63 = d["differences"]["6-3_cancer_specific"]
            txt += f"\n3 − 0 change: {ci(d30['estimate'], *d30['ci'])} pp\n6 − 3 change: {ci(d63['estimate'], *d63['ci'])} pp"
            boxed(ax, 0.02, 0.97, txt, ha="left", va="top", linespacing=1.3)
        else:
            boxed(ax, 0.02, 0.03, txt, ha="left", va="bottom", linespacing=1.3)
        ax.set_xlim(-0.5, 18.5); ax.set_ylim(0, 100); ax.set_xticks(range(0, 19, 3)); ax.set_yticks(range(0, 101, 20))
        if i % 3 == 0:
            ax.set_ylabel("Probability of cancer death\nwithin the window (%)")
        if i // 3 == 1:
            ax.set_xlabel("Landmark (months after diagnosis)")
        plab(fig, xs[i % 3] - 13.5, ys[i // 3] + 44, "ABCDEF"[i], f"{w}-month window")
        rec[str(w)] = {"landmarks": d["landmark"], "cancer_specific": d["cancer_specific"], "ci_lo": d["cancer_specific_ci_lo"], "ci_hi": d["cancer_specific_ci_hi"], "n_at_risk": d["n_at_risk"],
                       "c": c, "c_ci": cci, "prop_c_ge_3": d["plateau_prop_boot_ge3"], "b0": d["broken_stick_intercept"], "b1": d["broken_stick_slope_after_c_per_month"],
                       "diff_6_0": df, "diff_3_0": d["differences"]["3-0_cancer_specific"], "diff_6_3": d["differences"]["6-3_cancer_specific"]}
    hh, ll = h_open_filled()
    hS = [h_line(COL["NEC"], "o", "-", band=True), h_line(COL["SCLC"], None, "-", lw=1.0), (Patch(facecolor=COL["SCLC"], alpha=0.18, lw=0), Line2D([], [], color=COL["SCLC"], ls=":", lw=0.9))] + hh
    lS = ["Aalen–Johansen estimate, 95% bootstrap band", "Broken-stick fit (landmarks 0–18)", "Change point c (dotted) and its 95% CI (shaded)"] + ll
    fig.legend(hS, lS, loc="lower left", bbox_to_anchor=(16 / 180, 9 / 130), ncol=3, handlelength=2.4, columnspacing=1.5)
    nr = t["windows"]["12"]["n_at_risk"]; lm = t["windows"]["12"]["landmark"]
    ftext(fig, 16, 4.5, f"Distant colorectal NEC, n = 908; B = 2,000, seed {A02['seed']}, identical resamples across windows. Numbers at risk (identical for every window)\nat landmarks "
          + ", ".join(str(lm[i]) for i in range(0, 19, 3)) + " months: " + ", ".join(n_(nr[i]) for i in range(0, 19, 3)) + ".", va="center", ha="left", color=GREY, linespacing=1.35)
    save(fig, "SuppFigS4")
    MANIFEST["figures"]["SuppFigS4"] = {"source": "A02.json tumours.NEC.windows.{3,6,9,12,18,24}", "B": t["B"], "seed": A02["seed"], "windows": rec}


# =============================================================================================== SUPPLEMENTARY FIGURE S5
def suppS5():
    fig = newfig(180, 98)
    ax = axmm(fig, 16, 22, 161, 72)
    cur = A09["groups"]["NEC"]["curves"]; rec = {}
    ax.axvspan(16, 24, color="#000000", alpha=0.05, lw=0, zorder=0)
    ax.text(20, 0.006, "Fewer than 100 patients at risk", ha="center", va="bottom", fontsize=FS, color=GREY)
    styles = {"half": (COL["SCLC"], "-.", 1.0, f"Half bandwidth ({r(cur['half']['bandwidth'], 3)} months)"),
              "auto": (COL["NEC"], "-", 1.5, f"Automatic bandwidth ({r(cur['auto']['bandwidth'], 3)} months; minimum IMSE)\nwith pointwise 95% bootstrap band (B = {A09['bootstrap']['B']})"),
              "double": (COL["black"], "--", 1.0, f"Double bandwidth ({r(cur['double']['bandwidth'], 3)} months)")}
    hs, ls_ = [], []
    for lab in ["half", "auto", "double"]:
        mth, h, lo, hi, nr = grid_curve("NEC", lab)
        col, ls, lw, name = styles[lab]
        if lab == "auto":
            ax.fill_between(mth, lo, hi, color=col, alpha=0.15, lw=0, zorder=1)
        ax.plot(mth, h, color=col, ls=ls, lw=lw, zorder=3 if lab != "auto" else 4)
        hs.append(h_line(col, None, ls, lw=lw, band=(lab == "auto"))); ls_.append(name)
        c = cur[lab]; b = c["boot"]
        rec[lab] = {"month": mth.tolist(), "hazard": h.tolist(), "bandwidth": c["bandwidth"], "trough_month": c["min_month"], "trough_haz": c["min_haz"], "peak_month": c["max_month"],
                    "peak_haz": c["max_haz"], "ratio_fixed": c["ratio_max_min"], "ratio_fixed_ci": b["ratio_fixed_ci"], "ratio_window": c["ratio_window"], "ratio_window_ci": b["ratio_window_ci"],
                    "trough_month_ci": b["min_month_ci"], "peak_month_ci": b["max_month_ci"]}
        if lab == "auto":
            rec[lab].update({"lower95": lo.tolist(), "upper95": hi.tolist(), "n_at_risk": nr.tolist()})
    a = cur["auto"]
    ax.plot([a["min_month"], a["max_month"]], [a["min_haz"], a["max_haz"]], marker="o", ls="none", ms=5.5, mfc="none", mec=COL["black"], mew=0.9, zorder=6)
    ax.annotate(f"Trough, month {r(a['min_month'], 0)} ({r(a['min_haz'], 3)} per person-month)", xy=(a["min_month"], a["min_haz"]), xytext=(5.0, 0.047), fontsize=FS,
                ha="left", va="center", arrowprops=dict(arrowstyle="-", lw=0.5, color=GREY), color=COL["black"])
    ax.annotate(f"Post-trough maximum, month {r(a['max_month'], 2)}\n({r(a['max_haz'], 3)} per person-month)", xy=(a["max_month"], a["max_haz"]), xytext=(10.3, 0.168), fontsize=FS,
                ha="left", va="center", arrowprops=dict(arrowstyle="-", lw=0.5, color=GREY), color=COL["black"])
    # ratio table (Arial, aligned columns)
    tx = [0.47, 0.585, 0.765]; ty = 0.975; dy = 0.058
    ax.text(tx[0], ty, "Bandwidth", transform=ax.transAxes, fontsize=FS, fontweight="bold", va="top")
    ax.text(tx[1], ty, "Peak/trough ratio,\nfixed months (95% CI)", transform=ax.transAxes, fontsize=FS, fontweight="bold", va="top", linespacing=1.15)
    ax.text(tx[2], ty, "Mean hazard 6–9 / 3–6 months,\ntrapezoid rule (95% CI)", transform=ax.transAxes, fontsize=FS, fontweight="bold", va="top", linespacing=1.15)
    for k, (lab, nm) in enumerate([("half", "Half"), ("auto", "Automatic"), ("double", "Double")]):
        c = cur[lab]; b = c["boot"]; yy = ty - 0.095 - k * dy
        ax.text(tx[0], yy, nm, transform=ax.transAxes, fontsize=FS, va="top")
        ax.text(tx[1], yy, rr(c["ratio_max_min"], *b["ratio_fixed_ci"]), transform=ax.transAxes, fontsize=FS, va="top")
        ax.text(tx[2], yy, rr(c["ratio_window"], *b["ratio_window_ci"]), transform=ax.transAxes, fontsize=FS, va="top")
    ax.set_xlim(0, 24); ax.set_ylim(0, 0.25); ax.set_xticks(range(0, 25, 3)); ax.set_yticks(np.arange(0, 0.26, 0.05))
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, p: r(v, 2)))
    ax.set_xlabel("Months after diagnosis"); ax.set_ylabel("Smoothed hazard of cancer death (per person-month)")
    ax.legend(hs, ls_, loc="upper left", handlelength=2.6, bbox_to_anchor=(0.0, 1.0))
    ftext(fig, 16, 5, f"Distant colorectal NEC, n = {A09['groups']['NEC']['n']} ({A09['groups']['NEC']['cancer_deaths_24m']} cancer deaths within 24 months); Epanechnikov kernel, boundary-corrected (muhaz), "
          f"0.25-month grid;\nmonths 0 and 24 are boundary estimates. Bootstrap B = {A09['bootstrap']['B']}, seed {A09['bootstrap']['seed']}; open circles: trough and post-trough maximum at the automatic bandwidth.",
          va="center", ha="left", color=GREY, linespacing=1.35)
    save(fig, "SuppFigS5")
    MANIFEST["figures"]["SuppFigS5"] = {"source": "A09_hazard_grid.csv group NEC bw_label half/auto/double; A09.json groups.NEC.curves.{half,auto,double}", "curves": rec}


# =============================================================================================== SUPPLEMENTARY FIGURE S6
def suppS6():
    fig = newfig(180, 112)
    axA = axmm(fig, 15, 50, 72, 55); axB = axmm(fig, 104, 50, 72, 55)
    ext = A08["extended"]; lm = np.array(A08["landmarks"], float); rec = {"A": {}, "B": {}}
    aj = ext["aj_908"]; base = ext["base"]
    cc = A08["bootstrap"]["curve_ci"]; lo = np.array([q["lo"] for q in cc]); hi = np.array([q["hi"] for q in cc]); lmc = np.array([q["landmark"] for q in cc], float)
    axA.fill_between(lmc, lo, hi, color="#000000", alpha=0.10, lw=0, zorder=0)
    axA.plot(lm, aj["curve"], color=GREY, ls="--", lw=1.1, zorder=2)
    axA.plot(lm, base["curve"], color=COL["black"], ls="-", lw=1.6, zorder=4)
    hA = [h_line(GREY, None, "--"), (Patch(facecolor="#000000", alpha=0.12, lw=0), Line2D([], [], color=COL["black"], lw=1.6))]
    lA = [f"Unadjusted Aalen–Johansen (6 − 0: {sg(aj['d6_0'], 2)} pp; c = {r(aj['plateau_end'], 2)} months)",
          f"Standardized, original covariates (6 − 0: {sg(base['d6_0'], 2)} pp; c = {r(base['plateau_end'], 2)} months);\nband: 95% bootstrap CI (B = {A08['bootstrap']['B']}; approximate)"]
    sets = [("plus_surg", COL["NEC"], "-", "o", "+ primary-site surgery"), ("plus_race5", COL["AC"], "--", "^", "+ race/ethnicity (5 levels)"),
            ("plus_income3", COL["PANC"], "-.", "D", "+ county-income tertile"), ("plus_chemo", COL["SCLC"], ":", "s", "+ recorded chemotherapy (immortal-time caveat)\n"),
            ("all_combined_908", COL["LCNEC"], "-", "v", "All four added together")]
    for key, col, ls, mk, name in sets:
        e = ext[key]; mfc = "white" if key == "all_combined_908" else col
        axA.plot(lm, e["curve"], color=col, ls=ls, lw=0.9, marker=mk, ms=2.6, markevery=3, mfc=mfc, mew=0.7, zorder=3)
        hA.append(h_line(col, mk, ls, lw=0.9, mfc=mfc)); lA.append(f"{name}{'' if name.endswith(chr(10)) else ' '}(6 − 0: {sg(e['d6_0'], 2)} pp; c = {r(e['plateau_end'], 2)} months)")
        rec["A"][key] = {"curve": e["curve"], "d6_0": e["d6_0"], "plateau_end": e["plateau_end"], "formula": e["formula"], "n_landmarks_divergent_coef": e["n_landmarks_divergent_coef"]}
    rec["A"]["aj_908"] = {"curve": aj["curve"], "d6_0": aj["d6_0"], "plateau_end": aj["plateau_end"]}
    rec["A"]["base"] = {"curve": base["curve"], "d6_0": base["d6_0"], "plateau_end": base["plateau_end"], "curve_ci": cc, "B": A08["bootstrap"]["B"], "seed": A08["bootstrap"]["seed"],
                        "n_replicates_with_divergent_coef": A08["bootstrap"]["n_replicates_with_divergent_coef"]}
    axA.set_xlim(-0.5, 18.5); axA.set_ylim(30, 80); axA.set_xticks(range(0, 19, 3))
    axA.set_xlabel("Landmark (months after diagnosis)"); axA.set_ylabel("12-month probability of cancer death (%)")
    axA.legend(hA, lA, loc="upper left", bbox_to_anchor=(-0.17, -0.2), handlelength=2.6, title=f"Distant colorectal NEC (n = {aj['n']})", alignment="left")
    plab(fig, 2, 110, "A")
    comp = A08["comparators"]; hB, lB = [], []
    for key in ["PANC", "SCLC", "AC"]:
        cd = comp[key]
        axB.plot(lm, cd["aj"], color=COL[key], ls="--", lw=0.9, zorder=2)
        axB.plot(lm, cd["curve"], color=COL[key], ls="-", lw=1.3, marker=MK[key], ms=2.6, markevery=3, mew=0.7, zorder=3)
        d_aj = cd["aj_m6"] - cd["aj_m0"]
        hB.append(h_line(COL[key], MK[key], "-", lw=1.3)); lB.append(f"{LAB[key]} (n = {n_(cd['n_total'])})\n6 − 0: {sg(cd['d6_0'], 2)} pp standardized, {sg(d_aj, 2)} pp unadjusted")
        rec["B"][key] = {"std_curve": cd["curve"], "aj_curve": cd["aj"], "d6_0_std": cd["d6_0"], "d6_0_aj": d_aj, "ref_n": cd["ref_n"], "n": cd["n_total"], "plateau_end_std": cd["plateau_end"]}
    axB.plot(lm, base["curve"], color=COL["black"], ls="-", lw=1.6, zorder=4)
    hB.append(h_line(COL["black"], None, "-", lw=1.6)); lB.append(f"Colorectal NEC, standardized (reference): 6 − 0 {sg(base['d6_0'], 2)} pp")
    hB += [Line2D([], [], color=GREY, lw=1.3), Line2D([], [], color=GREY, lw=0.9, ls="--")]
    lB += ["Solid: directly standardized (original covariate set;\nrandom reference of 5,000 patients)", "Dashed: unadjusted Aalen–Johansen"]
    axB.legend(hB, lB, loc="upper left", bbox_to_anchor=(-0.08, -0.2), handlelength=2.6, title="Comparator tumors (distant stage)", alignment="left")
    axB.set_xlim(-0.5, 18.5); axB.set_ylim(30, 80); axB.set_xticks(range(0, 19, 3)); axB.set_xlabel("Landmark (months after diagnosis)")
    axB.set_ylabel("12-month probability of cancer death (%)")
    plab(fig, 91, 110, "B")
    save(fig, "SuppFigS6")
    MANIFEST["figures"]["SuppFigS6"] = {"source": "A08.json extended.{aj_908,base,plus_surg,plus_race5,plus_income3,plus_chemo,all_combined_908}; bootstrap.curve_ci; comparators.{AC,SCLC,PANC}", **rec}


# =============================================================================================== SUPPLEMENTARY FIGURE S7
def suppS7():
    ct = A05["chemotherapy_time_dependent_supplement"]
    fig = newfig(180, 166)
    rec = {}
    cohorts = [("clean_no_surgery_no_rt", "No primary-site surgery and no radiotherapy"), ("full_cohort", "Full distant cohort")]
    xs = [16, 105]
    for j, (ck, title) in enumerate(cohorts):
        c = ct[ck]; su = c["setup"]; lms = np.array(c["landmarks"], float)
        axA = axmm(fig, xs[j], 108, 71, 46); axB = axmm(fig, xs[j], 42, 71, 38)
        ftext(fig, xs[j] - 1, 158.5, f"{title}\n(n = {su['n']}; chemotherapy recorded in {su['n_chemo_yes']})", fontweight="bold", va="bottom", ha="left", fontsize=6.5, linespacing=1.25)
        rr_ = {"A": {}, "B": {}}
        hA, lA = [], []
        for grp, col, mk, ls, name in [("initiated", COL["NEC"], "o", "-", "Chemotherapy initiated by the landmark"), ("not_initiated", COL["SCLC"], "s", "--", "Not initiated by the landmark (never or not yet)")]:
            L = c["landmark_classification"][grp]["landmark"]
            est = np.array([L[str(int(l))]["est"] for l in lms]); lo = np.array([L[str(int(l))]["ci"][0] for l in lms]); hi = np.array([L[str(int(l))]["ci"][1] for l in lms])
            n = np.array([L[str(int(l))]["n_at_risk"] for l in lms], float); keep = n >= 20
            axA.fill_between(lms[keep], lo[keep], hi[keep], color=col, alpha=0.13, lw=0, zorder=1)
            curve(axA, lms[keep], est[keep], n[keep], col, mk, ls=ls, ms=3.2)
            hA.append(h_line(col, mk, ls, band=True)); lA.append(name + ", 95% bootstrap band")
            rr_["A"][grp] = {"landmarks": lms.tolist(), "est": est.tolist(), "lo": lo.tolist(), "hi": hi.tolist(), "n_at_risk": n.tolist(), "shown": keep.tolist()}
        for grp, col, ls, mk, name in [("chemo_yes", COL["black"], (0, (4, 2)), "x", "From-diagnosis classification: chemotherapy recorded"), ("chemo_no", GREY, ":", "+", "From-diagnosis classification: no/unknown")]:
            nv = c["naive_stratum_landmarks"][grp]; est = np.array(nv["est"], float); n = np.array(nv["n_at_risk"], float); keep = n >= 20
            axA.plot(lms[keep], est[keep], color=col, ls=ls, lw=1.0, marker=mk, ms=3.4, mew=0.9, zorder=2)
            hA.append(h_line(col, mk, ls, lw=1.0)); lA.append(name + " (point estimates)")
            rr_["A"][grp] = {"landmarks": lms.tolist(), "est": est.tolist(), "n_at_risk": n.tolist(), "shown": keep.tolist()}
        axA.set_xlim(-0.4, 12.4); axA.set_ylim(40, 100); axA.set_xticks([0, 1, 2, 3, 6, 9, 12])
        axA.set_xlabel("Landmark (months after diagnosis)")
        if j == 0:
            axA.set_ylabel("12-month probability of cancer death (%)")
        # ---- (B) time-dependent crude rates
        ivs = [q["interval"] for q in c["time_dependent_rates"]]; xi = np.arange(len(ivs))
        hB, lB = [], []
        for k, (key, col, mk, mfc, name) in enumerate([("unexposed", COL["SCLC"], "s", "white", "Person-time before chemotherapy start (never or not yet)"),
                                                     ("exposed", COL["NEC"], "o", COL["NEC"], "Person-time after chemotherapy start")]):
            v = [(q[key]["rate"], q[key]["deaths"], q[key]["person_months"]) for q in c["time_dependent_rates"]]
            cis = [poisson_exact(u[1], u[2]) for u in v]; xp = xi + (-0.14 if k == 0 else 0.14)
            axB.errorbar(xp, [u[0] for u in v], yerr=[[u[0] - q_[0] for u, q_ in zip(v, cis)], [q_[1] - u[0] for u, q_ in zip(v, cis)]], fmt=mk, color=col, mfc=mfc, mec=col, ms=4,
                         mew=0.9, elinewidth=0.8, capsize=1.6, ls="none", zorder=3)
            for xq, u in zip(xp, v):
                axB.text(xq + (-0.07 if k == 0 else 0.07), u[0], f"{u[1]}", fontsize=FS, ha="right" if k == 0 else "left", va="center", color=col)
            hB.append(Line2D([], [], marker=mk, color=col, mfc=mfc, mec=col, ms=4, mew=0.9, ls="none")); lB.append(name)
            rr_["B"][key] = [{"interval": iv, "rate": u[0], "deaths": u[1], "person_months": u[2], "exact_ci": list(q_)} for iv, u, q_ in zip(ivs, v, cis)]
        rp = c["rate_ratios_poisson_ci"]; tdr = c["time_dependent_rates"]
        e36 = [q for q in tdr if q["interval"] == "3-6"][0]["exposed"]["rate"]; e69 = [q for q in tdr if q["interval"] == "6-9"][0]["exposed"]["rate"]
        after3 = sum(q["unexposed"]["deaths"] for q in tdr if q["interval"] != "0-3")
        ftext(fig, xs[j], 31, f"After chemotherapy start: {r(e36, 3)} (3–6 months) and {r(e69, 3)} (6–9 months)\n"
              f"per person-month; RR 6–9 vs 3–6 months {rr(rp['exposed_6_9_vs_3_6']['RR'], *rp['exposed_6_9_vs_3_6']['ci'])}.\n"
              f"Before start: {after3} cancer deaths after month 3.", va="top", ha="left", linespacing=1.35)
        rr_["B"]["rate_ratios"] = rp; rr_["B"]["exposed_rate_3_6"] = e36; rr_["B"]["exposed_rate_6_9"] = e69; rr_["B"]["unexposed_deaths_after_month3"] = after3
        axB.set_xticks(xi); axB.set_xticklabels([iv.replace("-", "–") for iv in ivs]); axB.set_xlim(-0.5, len(ivs) - 0.5); axB.set_ylim(0, 0.75)
        axB.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, p: r(v, 1)))
        axB.set_xlabel("Months after diagnosis")
        if j == 0:
            axB.set_ylabel("Cancer deaths per person-month")
        rec[ck] = rr_
    hh, ll = h_open_filled()
    fig.legend(hA + hh, lA + ll, loc="upper left", bbox_to_anchor=(16 / 180, 99.5 / 166), ncol=2, handlelength=2.6, columnspacing=1.6)
    fig.legend(hB, [lB[0], lB[1]], loc="upper left", bbox_to_anchor=(16 / 180, 19 / 166), ncol=2, handlelength=1.5, columnspacing=1.6,
               title="Numbers next to markers: cancer deaths; error bars: exact Poisson 95% CI", alignment="left")
    ftext(fig, 16, 4.5, "Landmark cells with fewer than 20 patients at risk are not drawn; open markers in (A): fewer than 100 patients at risk.\n"
          "Treatment start: SEER time from diagnosis to first-course treatment.", va="center", ha="left", color=GREY, linespacing=1.35)
    plab(fig, 2, 164, "A"); plab(fig, 2, 86, "B")
    save(fig, "SuppFigS7")
    MANIFEST["figures"]["SuppFigS7"] = {"source": "A05.json chemotherapy_time_dependent_supplement.{clean_no_surgery_no_rt,full_cohort}.{landmark_classification,naive_stratum_landmarks,time_dependent_rates,rate_ratios_poisson_ci}", **rec}


if __name__ == "__main__":
    import warnings
    with warnings.catch_warnings(record=True) as wlist:
        warnings.simplefilter("always")
        figure1(); figure2(); figure3()
        suppS1(); suppS2(); suppS3(); suppS4(); suppS5(); suppS6(); suppS7()
    gl = [str(w.message) for w in wlist if "Glyph" in str(w.message) or "findfont" in str(w.message)]
    MANIFEST["font_warnings"] = gl
    print("font warnings:", gl if gl else "none")
    MANIFEST["runtime_minutes"] = round((time.time() - T0) / 60, 2)
    with open(OUT + "figure_data.json", "w") as f:
        json.dump(MANIFEST, f, indent=1, default=float, ensure_ascii=False)
    print("done", MANIFEST["runtime_minutes"], "min")
