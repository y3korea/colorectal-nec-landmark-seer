# -*- coding: utf-8 -*-
"""graphical_abstract.py -- Graphical Abstract for the International Journal of Cancer submission.

IJC author guidelines 3.5: roughly square, colour, no manuscript title or author names, labels and short phrases only,
no tables or raw primary data; PDF (minimum 600 dpi). Output: vector PDF (TrueType fonts, saved at dpi=600) + PNG 600 dpi.

Reads stored aggregate results only (no patient-level data, nothing re-estimated):
  A11.json  table3.{NEC,SCLC,AC,PANC}.cancer.{monthly,n,plateau_end}   -> left panel curves (same series as Figure 3A)
  A06.json  reproduction.crude_rates (0-3, 3-6, 6-9), reproduction.adj_RR.r69_36, comparators_histology.{SCLC,AC,PANC}.r69_36
  A12.json  lcnec.rr_6_9_vs_3_6 (lung LCNEC; A12 reproduces the A06 NEC/SCLC RRs -- asserted)
Paths come from python/config.py (NEC_RESULTS_DIR, NEC_FIG_DIR), as for every other step of the pipeline.
Output: <NEC_FIG_DIR>/manuscript_final/GraphicalAbstract.{pdf,png} and graphical_abstract_data.json
Run (after run_all.py, or at least after A06, A11 and A12):  python3 python/figures/graphical_abstract.py
"""
import json, os, sys
from decimal import Decimal, ROUND_HALF_UP
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # python/ (config.py)
import config                                      # noqa: E402

RES = config.RES
OUT = config.MANUSCRIPT_FIG_DIR + os.sep
os.makedirs(OUT, exist_ok=True)


def J(name):
    p = RES + name
    if not os.path.exists(p):
        raise SystemExit(f"{p} not found: run the pipeline first (python3 run_all.py, or at least --only A06,A11,A12)")
    with open(p, encoding="utf-8") as f:
        return json.load(f)


A06, A11, A12 = J("A06.json"), J("A11.json"), J("A12.json")

COL = {"NEC": "#0072B2", "SCLC": "#D55E00", "AC": "#009E73", "PANC": "#CC79A7", "LCNEC": "#E69F00",
       "grey": "#6E6E6E", "light": "#BBBBBB", "black": "#222222"}
MK = {"NEC": "o", "SCLC": "s", "AC": "^", "PANC": "D", "LCNEC": "v"}
LS = {"NEC": "-", "SCLC": "-", "AC": "--", "PANC": "-.", "LCNEC": ":"}
LAB = {"NEC": "Colorectal NEC", "SCLC": "Small-cell lung cancer", "AC": "Colorectal adenocarcinoma",
       "PANC": "Pancreatic adenocarcinoma", "LCNEC": "Lung LCNEC"}
MINUS = "−"

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 9,
    "axes.linewidth": 0.8, "xtick.major.width": 0.8, "ytick.major.width": 0.8, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
    "axes.labelsize": 9.5, "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "axes.unicode_minus": True})
FONT_FILE = fm.findfont(fm.FontProperties(family=["sans-serif"]), fallback_to_default=False)
FONT_NAME = fm.FontProperties(fname=FONT_FILE).get_name()
if FONT_NAME != "Arial":   # same policy as python/FIG_manuscript_final.py
    if os.environ.get("NEC_REQUIRE_ARIAL", "0") == "1":
        raise SystemExit(f"Arial not resolved (got {FONT_NAME}); see README")
    print(f"WARNING: Arial not available, the Graphical Abstract uses {FONT_NAME} (layout may differ slightly from the submitted file)")


def r(x, d):
    q = Decimal(repr(float(x))).quantize(Decimal(1).scaleb(-d), rounding=ROUND_HALF_UP)
    return f"{abs(q) if q == 0 else q:.{d}f}".replace("-", MINUS)


def rr(e, d=2):
    return f"{r(e['RR'], d)} ({r(e['lo'], d)}–{r(e['hi'], d)})"


# ---------------------------------------------------------------- guards: A12 must reproduce A06
for k_, a6 in [("NEC", A06["reproduction"]["adj_RR"]["r69_36"]), ("SCLC", A06["comparators_histology"]["SCLC"]["r69_36"])]:
    a12 = A12["comparators_same_code_path"][k_]["rr_6_9_vs_3_6"]
    assert all(abs(a12[q] - a6[q]) < 1e-6 for q in ("RR", "lo", "hi")), k_

W, H = 180, 172                          # mm, roughly square
fig = plt.figure(figsize=(W / 25.4, H / 25.4))
fig.patch.set_facecolor("white")
mm = lambda x, y, w, h: fig.add_axes([x / W, y / H, w / W, h / H])
DATA = {"id": "FIG_graphical_abstract", "date": "2026-10-08", "font": FONT_NAME, "printed": {}, "sources": {}}

# ---------------------------------------------------------------- left: monthly conditional 12-month probability curves
ax = mm(19, 44, 82, 106)
t = A11["table3"]
pe = t["NEC"]["cancer"]["plateau_end"]["months"]
# audit 2026-10-08: shade the pre-specified 0-6-month span (not the descriptive, window-dependent broken-stick hinge)
ax.axvspan(0, 6, color=COL["NEC"], alpha=0.07, lw=0, zorder=0)
ax.text(0.4, 97, "No detectable decline\nover 0–6 months", ha="left", va="top", fontsize=8.5, color=COL["NEC"], fontweight="bold", linespacing=1.15)
curves = {}
for key in ["PANC", "SCLC", "AC", "NEC"]:
    m = t[key]["cancer"]["monthly"]; x = np.array(m["landmarks"], float); y = np.array(m["pct"], float)
    if key == "NEC":
        ax.fill_between(x, m["lo"], m["hi"], color=COL[key], alpha=0.16, lw=0, zorder=2)
    ax.plot(x, y, color=COL[key], ls=LS[key], lw=2.6 if key == "NEC" else 1.6, zorder=4 if key == "NEC" else 3)
    sel = np.arange(0, len(x), 6)
    ax.plot(x[sel], y[sel], ls="none", marker=MK[key], ms=5.5 if key == "NEC" else 4.2, mfc=COL[key], mec="white", mew=0.6, zorder=5)
    curves[key] = {"landmarks": m["landmarks"], "pct": m["pct"], "n": t[key]["cancer"]["n"]}
DATA["sources"]["curves"] = "A11.json table3.{NEC,SCLC,AC,PANC}.cancer.monthly (pct); NEC 95% band lo/hi; plateau_end.months"
mN = t["NEC"]["cancer"]["monthly"]; i0, i6 = mN["landmarks"].index(0), mN["landmarks"].index(6)
p0, p6 = mN["pct"][i0], mN["pct"][i6]
for xx, pp in [(0, p0), (6, p6)]:
    ax.annotate(f"{r(pp, 0)}%", (xx, pp), xytext=(xx + (0.9 if xx == 0 else 0.4), pp + 8.5), fontsize=9, fontweight="bold", color=COL["NEC"],
                ha="left", va="bottom", arrowprops=dict(arrowstyle="-", color=COL["NEC"], lw=0.8, shrinkA=0, shrinkB=3), zorder=6)
ax.set_xlim(-0.5, 24.5); ax.set_ylim(0, 100); ax.set_xticks(range(0, 25, 6)); ax.set_yticks(range(0, 101, 20))
ax.set_xlabel("Months survived since diagnosis")
ax.set_ylabel("Probability of cancer death\nin the next 12 months (%)")
hs = [Line2D([], [], color=COL[k], ls=LS[k], lw=2.6 if k == "NEC" else 1.6, marker=MK[k], ms=5, mfc=COL[k], mec="white") for k in ["NEC", "SCLC", "PANC", "AC"]]
fig.legend(hs, [LAB[k] for k in ["NEC", "SCLC", "PANC", "AC"]], loc="upper left", bbox_to_anchor=(6 / W, 30 / H), ncol=1,
           fontsize=9, handlelength=3.0, labelspacing=0.55)
fig.text(6 / W, 166 / H, "Distant-stage colorectal NEC vs other cancers (SEER)", fontsize=10.5, fontweight="bold", color=COL["black"], va="top")
fig.text(6 / W, 160 / H, "Conditional 12-month risk of cancer death, by time already survived", fontsize=9, color=COL["grey"], va="top")
DATA["printed"]["left"] = {"shaded_span_months": "0-6 (pre-specified landmarks)", "NEC_pct_0": r(p0, 0), "NEC_pct_6": r(p6, 0)}

# ---------------------------------------------------------------- right: three take-home messages
X0, BW = 107, 69
boxes = [(128, 32), (66, 54), (6, 56)]   # (y, h) in mm
bg = fig.add_axes([0, 0, 1, 1], zorder=-10); bg.set_xlim(0, W); bg.set_ylim(0, H); bg.axis("off")
for (y0, h), num in zip(boxes, ["1", "2", "3"]):
    bg.add_patch(FancyBboxPatch((X0, y0), BW, h, boxstyle="round,pad=0,rounding_size=2.2", fc="#F3F7FB", ec="#C9D9E8", lw=0.8, zorder=0))
    bg.add_patch(matplotlib.patches.Circle((X0 + 6, y0 + h - 6), 3.6, fc=COL["NEC"], ec="none", zorder=1))
    fig.text((X0 + 6) / W, (y0 + h - 6) / H, num, ha="center", va="center", color="white", fontsize=10, fontweight="bold")


def head(y0, h, s, fs=10.5):
    fig.text((X0 + 12) / W, (y0 + h - 2.6) / H, s, ha="left", va="top", fontsize=fs, fontweight="bold", color=COL["black"], linespacing=1.12)


# box 1
y0, h = boxes[0]
head(y0, h, "No detectable early decline")
fig.text((X0 + 5) / W, (y0 + h - 13) / H, f"Risk of cancer death in the next\n12 months: {r(p0, 0)}% at diagnosis and\nstill {r(p6, 0)}% after 6 months survived",
         ha="left", va="top", fontsize=9, color=COL["black"], linespacing=1.3)
# box 2: mini bar chart of crude interval rates (dip and rebound)
y0, h = boxes[1]
head(y0, h, "Hazard dips at 3–6 months,\nrebounds at 6–9 months")
cr = {q["interval"]: q["rate"] for q in A06["reproduction"]["crude_rates"]}
ivs = ["0-3", "3-6", "6-9"]
axb = mm(X0 + 8, y0 + 13, 24, 21)
bars = axb.bar(range(3), [cr[i] for i in ivs], color=[COL["NEC"], "#9CC3E0", COL["NEC"]], width=0.68)
axb.set_xticks(range(3)); axb.set_xticklabels([i.replace("-", "–") for i in ivs], fontsize=7.5)
axb.set_yticks([]); axb.spines["left"].set_visible(False); axb.set_ylim(0, 0.17); axb.set_facecolor("none")
axb.tick_params(axis="x", length=0, pad=1.5)
axb.set_xlabel("Months", fontsize=7.5, labelpad=1)
axb.annotate("", xy=(2, cr["6-9"] + 0.012), xytext=(1, cr["3-6"] + 0.012),
             arrowprops=dict(arrowstyle="-|>", color=COL["SCLC"], lw=1.3, connectionstyle="arc3,rad=-0.35", mutation_scale=9))
e = A06["reproduction"]["adj_RR"]["r69_36"]
fig.text((X0 + 36) / W, (y0 + 31) / H, f"Rate ratio\n6–9 vs 3–6 months\n{r(e['RR'], 2)}\n(95% CI {r(e['lo'], 2)}–{r(e['hi'], 2)})",
         ha="left", va="top", fontsize=8.5, color=COL["black"], linespacing=1.25)
fig.text((X0 + 5) / W, (y0 + 3.2) / H, "Bars: cancer deaths per person-month", ha="left", va="bottom", fontsize=7.5, color=COL["grey"])
DATA["printed"]["box2"] = {"crude_rates": {i: cr[i] for i in ivs}, "RR_69_36": rr(e)}
DATA["sources"]["box2"] = "A06.json reproduction.crude_rates, reproduction.adj_RR.r69_36"
# box 3: rebound RR by tumour (dot strip)
y0, h = boxes[2]
head(y0, h, "No difference detected\nfrom small-cell lung cancer;\nabsent in adenocarcinomas\nand lung LCNEC", fs=9.5)
items = [("NEC", A06["reproduction"]["adj_RR"]["r69_36"], "Colorectal NEC"), ("SCLC", A06["comparators_histology"]["SCLC"]["r69_36"], "SCLC"),
         ("LCNEC", A12["lcnec"]["rr_6_9_vs_3_6"], "Lung LCNEC"), ("AC", A06["comparators_histology"]["AC"]["r69_36"], "Colorectal AC"),
         ("PANC", A06["comparators_histology"]["PANC"]["r69_36"], "Pancreatic AC")]
axc = mm(X0 + 28, y0 + 9, 36, 24)
for i, (k, e, lab) in enumerate(items):
    yy = -i
    axc.plot([e["lo"], e["hi"]], [yy, yy], color=COL[k], lw=1.4, solid_capstyle="butt")
    axc.plot(e["RR"], yy, marker=MK[k], color=COL[k], ms=5.5 if k == "NEC" else 4.8, mec="white", mew=0.5)
axc.axvline(1, color=COL["light"], lw=0.8, ls="--", zorder=0)
axc.set_xscale("log"); axc.set_xlim(0.75, 2.0); axc.set_xticks([0.8, 1, 1.5, 2]); axc.set_xticklabels(["0.8", "1", "1.5", "2"], fontsize=7.5)
axc.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
axc.set_yticks([-i for i in range(len(items))]); axc.set_yticklabels([it[2] for it in items], fontsize=7.8)
for tl, (k, _, _) in zip(axc.get_yticklabels(), items):
    tl.set_color(COL[k])
axc.set_ylim(-len(items) + 0.4, 0.6); axc.spines["left"].set_visible(False); axc.tick_params(axis="y", length=0, pad=2)
axc.tick_params(axis="x", pad=1.5); axc.set_facecolor("none")
fig.text((X0 + 45) / W, (y0 + 2.2) / H, "Rate ratio 6–9 vs 3–6 months", ha="center", va="bottom", fontsize=7.5, color=COL["grey"])
DATA["printed"]["box3"] = {k: rr(e) for k, e, _ in items}
DATA["sources"]["box3"] = "A06.json reproduction.adj_RR.r69_36, comparators_histology.{SCLC,AC,PANC}.r69_36; A12.json lcnec.rr_6_9_vs_3_6"
DATA["printed"]["box1"] = {"NEC_pct_0": r(p0, 0), "NEC_pct_6": r(p6, 0)}
DATA["n"] = {k: curves[k]["n"] for k in curves} | {"LCNEC": A12["lcnec"]["n"]}

fig.text(6 / W, 1.2 / H, "AC, adenocarcinoma; LCNEC, large-cell neuroendocrine carcinoma; NEC, neuroendocrine carcinoma; SCLC, small-cell lung cancer",
         ha="left", va="bottom", fontsize=7, color=COL["grey"])
fig.canvas.draw()   # tick labels exist only after a draw
DATA["panel_text"] = [tx.get_text() for tx in fig.findobj(matplotlib.text.Text) if tx.get_text().strip() and tx.get_visible()]
if config.SYNTHETIC:   # never let an image from synthetic data be mistaken for study results (not added for real data)
    fig.text(0.5, 0.5, "SYNTHETIC DATA - NOT STUDY RESULTS", fontsize=22, color="red", alpha=0.35, rotation=30, ha="center", va="center", fontweight="bold")
fig.savefig(OUT + "GraphicalAbstract.pdf", dpi=600, facecolor="white")
fig.savefig(OUT + "GraphicalAbstract.png", dpi=600, facecolor="white")
with open(OUT + "graphical_abstract_data.json", "w", encoding="utf-8") as f:
    json.dump(DATA, f, indent=1, ensure_ascii=False, default=float)
print(f"saved {OUT}GraphicalAbstract.pdf/.png and graphical_abstract_data.json (font {FONT_NAME})")
