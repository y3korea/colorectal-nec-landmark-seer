# -*- coding: utf-8 -*-
"""Aggregate cohort counts used by the cohort-flow figure (Supplementary Figure S1): results/FIG_flow_counts.json.

Reads the derived cohort files in NEC_WORK_DIR and writes counts only (no patient-level data).
Run: python3 python/flow_counts.py   (after Rscript R/00_build_cohorts.R)
"""
import json
import config
import lm


def vc(s):
    return {str(k): int(v) for k, v in s.value_counts().items()}


nen = lm.load("NEN_all"); nec = lm.nec_main(nen); expd = lm.nec_expanded(nen)
g8246 = nen[nen["hist"] == "8246"]
g34 = g8246.grade17.isin(["Poorly differentiated; Grade III", "Undifferentiated; anaplastic; Grade IV"]) | \
    g8246.grade18.isin(["Site-specific grade system category (3)", "Site-specific grade system category (4)"])
out = {"NEN_all_n": int(len(nen)), "NEC_main_n": int(len(nec)), "NEC_stage": vc(nec.stage), "NEC_hist": vc(nec["hist"]),
       "NEC_expanded_n": int(len(expd)), "NEC_expanded_distant": int((expd.stage == "Distant").sum()),
       "n_8246_total": int(len(g8246)), "n_8246_g34": int(g34.sum())}
del nen
lung = lm.load("LUNG"); sclc = lm.sclc_main(lung); lcnec = lung[lung["hist"] == "8013"]
out.update({"LUNG_n": int(len(lung)), "SCLC_n": int(len(sclc)), "SCLC_distant": int((sclc.stage == "Distant").sum()),
            "LCNEC_n": int(len(lcnec)), "LCNEC_distant": int((lcnec.stage == "Distant").sum()), "LUNG_hist": vc(lung["hist"])})
del lung, sclc, lcnec
panc = lm.load("PANC"); pm = lm.panc_main(panc)
out.update({"PANC_n": int(len(panc)), "PANC_excl_islet_n": int(len(pm)), "PANC_distant": int((pm.stage == "Distant").sum())})
del panc, pm
ac = lm.load("AC")
out.update({"AC_n": int(len(ac)), "AC_distant": int((ac.stage == "Distant").sum())})
with open(config.RES + "FIG_flow_counts.json", "w") as f:
    json.dump(out, f, indent=1)
print(json.dumps(out))
