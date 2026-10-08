# A12 lcnec_rate_ratios -- piecewise-exponential (Poisson) analyses adding distant lung large-cell NEC (LCNEC, ICD-O-3 8013)
# as a comparator for distant colorectal NEC. Runs non-interactively: Rscript R/A12_lcnec_rate_ratios.R
# Code path reused verbatim from A06_piecewise_poisson_suite.R (fitpw / rr_per / crude / three / shape_test).
# (0) reproduction anchors: NEC RR 6-9 vs 3-6 = 1.49 (1.17-1.90); NEC vs SCLC shape test W = 12.05, P = 0.061;
#     colorectal small-cell 1.75 (1.29-2.38), large-cell 1.18 (0.79-1.78), SCLC 1.38 (1.34-1.42) (A04/A06)
# (1) LCNEC crude interval rates and adjusted RRs (6-9 vs 3-6; 6-15 vs 0-6; 15-24 vs 6-15), cluster-robust SE
# (2) hazard-shape test NEC vs LCNEC (6-df Wald, pooled model), all years and 2016-2023; ratio of 6-9 vs 3-6 RRs NEC/LCNEC
# (3) 2x2 morphology (small-cell vs large-cell) x organ (colorectum vs lung) for RR 6-9 vs 3-6
# Output is aggregate only (no patient-level rows).
T_START <- Sys.time()
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "config.R"))  # paths, NEC_BOOT, NEC_STRICT
source(file.path(NEC_R_DIR, "pp.R"))
data.table::setDTthreads(1)
# no random numbers are drawn in A12 (no bootstrap), so no seed is set
ID <- "A12"
OUT_JSON <- paste0(RES, ID, ".json"); OUT_MD <- paste0(RES, ID, ".md")

CUTS7 <- c(3, 6, 9, 12, 15, 18); LAB7 <- c("0-3", "3-6", "6-9", "9-12", "12-15", "15-18", "18-24")
CUTS3 <- c(6, 15); LAB3 <- c("0-6", "6-15", "15-24")
fmt <- function(r, d = 2) sprintf(paste0("%.", d, "f (%.", d, "f-%.", d, "f)"), r[["RR"]], r[["lo"]], r[["hi"]])
fp <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 2), sprintf("%.4f", p))
ctr_list <- function(r) list(RR = unname(r[["RR"]]), lo = unname(r[["lo"]]), hi = unname(r[["hi"]]), p = unname(r[["p"]]), logRR = unname(r[["logRR"]]), se = unname(r[["se"]]))

# ---- verbatim from A06 ------------------------------------------------------------------
fitpw <- function(dt, start, cuts, covf, end = 24, extra = NULL) {
  s <- pw(dt, start, cuts, end)
  form <- paste("d ~ 0 + factor(per) +", covf)
  if (!is.null(extra)) form <- paste(form, "+", extra)
  f <- glm(as.formula(form), family = poisson, offset = log(expo), data = s)
  V <- vcovCL(f, cluster = s$pid)
  list(fit = f, V = V, s = s)
}
rr_per <- function(m, a, b) ctr(m$fit, m$V, vec(m$fit, setNames(list(1, -1), c(paste0("factor(per)", a), paste0("factor(per)", b)))))
crude <- function(s) { h <- s[, .(d = sum(d), pm = sum(expo)), by = per][order(per)]; h[, rate := d / pm]; h }
three <- function(dt, covf) {
  m7 <- fitpw(dt, 0, CUTS7, covf); m3 <- fitpw(dt, 0, CUTS3, covf)
  list(m7 = m7, m3 = m3, r69_36 = rr_per(m7, 3, 2), r615_06 = rr_per(m3, 2, 1), r1524_615 = rr_per(m3, 3, 2), crude = crude(m7$s))
}
shape_test <- function(necd, cmpd, yearterm = "year10") {
  X <- rbind(necd[, .(pid, tm, ev, age10, age10sq, male, year10, year, dis = 1L)],
             cmpd[, .(pid, tm, ev, age10, age10sq, male, year10, year, dis = 0L)])
  X[, period3 := cut(year, c(1999, 2007, 2015, 2023), labels = c("2000-2007", "2008-2015", "2016-2023"))]
  s <- pw(X, 0, CUTS7)
  f <- glm(as.formula(paste("d ~ 0 + factor(per) + factor(per):dis + age10 + age10sq + male +", yearterm)), family = poisson, offset = log(expo), data = s)
  V <- vcovCL(f, cluster = s$pid)
  ii <- grep(":dis", names(coef(f))); nmi <- names(coef(f))[ii]
  L <- matrix(0, 6, length(coef(f)))
  for (k in 2:7) { L[k - 1, which(names(coef(f)) == nmi[k])] <- 1; L[k - 1, which(names(coef(f)) == nmi[1])] <- -1 }
  b <- L %*% coef(f); W <- as.numeric(t(b) %*% solve(L %*% V %*% t(L)) %*% b)
  W7 <- wald(f, V, ii)
  rr <- ctr(f, V, vec(f, setNames(list(1, -1), c("factor(per)3:dis", "factor(per)2:dis"))))
  list(W = W, df = 6L, p = pchisq(W, 6, lower.tail = FALSE), W7 = unname(W7["W"]), p7 = unname(W7["p"]),
       ratio_r69_36_NEC_over_cmp = ctr_list(rr), n_nec = nrow(necd), n_cmp = nrow(cmpd))
}
# -----------------------------------------------------------------------------------------
# all seven interval RRs relative to 0-3 and to 3-6 (descriptive), from the 7-interval model
all_rr <- function(m) rbindlist(lapply(2:7, function(k) { r1 <- rr_per(m, k, 1); data.table(interval = LAB7[k], RR_vs_0_3 = r1[["RR"]], lo_vs_0_3 = r1[["lo"]], hi_vs_0_3 = r1[["hi"]]) }))
# difference of two independent log RRs (separate cohorts, no shared patients)
ratio_indep <- function(a, b) { e <- a[["logRR"]] - b[["logRR"]]; se <- sqrt(a[["se"]]^2 + b[["se"]]^2); c(RR = exp(e), lo = exp(e - 1.96 * se), hi = exp(e + 1.96 * se), p = 2 * pnorm(-abs(e / se)), logRR = e, se = se) }

res <- list(id = ID, date = as.character(Sys.Date()), software = list(R = R.version.string, survival = as.character(packageVersion("survival")), sandwich = as.character(packageVersion("sandwich"))),
            methods = list(model = "piecewise-exponential Poisson GLM, offset log(person-months), intervals 0-3/3-6/6-9/9-12/12-15/15-18/18-24 months (and 0-6/6-15/15-24), follow-up truncated at 24 months; outcome cancer-specific death; cluster-robust (vcovCL by patient) SE",
                           covariates = list(lung = COV$SCLC, colorectal_subtype = COV$AC, colorectal_main = COV$NEC, pooled_shape = "age10 + age10sq + male + year10 + factor(per) + factor(per):tumour"),
                           LCNEC_definition = "lung/bronchus, histology 8013, SEER combined summary stage Distant, histologically confirmed (cohort_LUNG)"))
md <- character(0); add <- function(...) md <<- c(md, sprintf(...))

# ------------------------------------------------------------------ load cohorts
cat("Loading cohorts...\n")
NEN <- loadc("NEN_all", "N"); NEC <- nec_main(NEN); rm(NEN)
LU <- loadc("LUNG", "L")
D <- list(NEC = NEC[stage == "Distant"], SCLC = LU[stage == "Distant" & hist %in% SCLC_CODES], LCNEC = LU[stage == "Distant" & hist == "8013"])
rm(LU, NEC)
D$CRC_small <- D$NEC[hist %in% SCLC_CODES]; D$CRC_large <- D$NEC[hist == "8013"]
res$n <- lapply(D, nrow)
cat("n:", unlist(res$n), "\n")
nec_check(res$n$NEC == 908, res$n$SCLC == 77707, res$n$LCNEC == 3185, res$n$CRC_small == 599, res$n$CRC_large == 309)

# ------------------------------------------------------------------ fits
cat("Fitting models...\n")
F <- list(NEC = three(D$NEC, COV$NEC), SCLC = three(D$SCLC, COV$SCLC), LCNEC = three(D$LCNEC, COV$SCLC),
          CRC_small = three(D$CRC_small, COV$AC), CRC_large = three(D$CRC_large, COV$AC))
cat("Shape tests...\n")
sh_sclc <- shape_test(D$NEC, D$SCLC)
sh_lc_all <- shape_test(D$NEC, D$LCNEC)
sh_lc_1623 <- shape_test(D$NEC[year >= 2016], D$LCNEC[year >= 2016])
sh_lc_fac <- shape_test(D$NEC, D$LCNEC, "period3")
# lung large vs small (organ held fixed; 'dis' = LCNEC); colorectal large vs small not needed (A04 covers it)
sh_lung_lc_vs_sc <- shape_test(D$LCNEC, D$SCLC)
# LCNEC restricted to 2016-2023 (same period as the NEC/LCNEC restricted shape test)
lc1623 <- three(D$LCNEC[year >= 2016], COV$SCLC)

# ------------------------------------------------------------------ (0) reproduction
repro <- data.table(
  quantity = c("NEC distant n", "NEC adj RR 6-9 vs 3-6", "  lower 95% CI", "  upper 95% CI", "NEC vs SCLC shape W (6 df)", "NEC vs SCLC shape P",
               "CRC small-cell n", "CRC small-cell RR 6-9 vs 3-6", "  lower", "  upper", "CRC large-cell n", "CRC large-cell RR 6-9 vs 3-6", "  lower", "  upper",
               "SCLC distant n", "SCLC RR 6-9 vs 3-6", "  lower", "  upper", "LCNEC distant n"),
  source = c("A06", "A06", "A06", "A06", "A06", "A06", "A04", "A04", "A04", "A04", "A04", "A04", "A04", "A04", "A06", "A06", "A06", "A06", "task spec"),
  target = c(908, 1.49, 1.17, 1.90, 12.05, 0.061, 599, 1.75, 1.29, 2.38, 309, 1.18, 0.79, 1.78, 77707, 1.38, 1.34, 1.42, 3185),
  ours = c(nrow(D$NEC), F$NEC$r69_36[["RR"]], F$NEC$r69_36[["lo"]], F$NEC$r69_36[["hi"]], sh_sclc$W, sh_sclc$p,
           nrow(D$CRC_small), F$CRC_small$r69_36[["RR"]], F$CRC_small$r69_36[["lo"]], F$CRC_small$r69_36[["hi"]],
           nrow(D$CRC_large), F$CRC_large$r69_36[["RR"]], F$CRC_large$r69_36[["lo"]], F$CRC_large$r69_36[["hi"]],
           nrow(D$SCLC), F$SCLC$r69_36[["RR"]], F$SCLC$r69_36[["lo"]], F$SCLC$r69_36[["hi"]], nrow(D$LCNEC)))
repro[, tol := c(0, rep(0.0051, 3), 0.0051, 0.00051, 0, rep(0.0051, 3), 0, rep(0.0051, 3), 0, rep(0.0051, 3), 0)]
repro[, ok := abs(ours - target) <= tol]
print(repro)
if (!all(repro$ok)) nec_anchor_fail("REPRODUCTION CHECK FAILED: ", paste(repro[ok == FALSE]$quantity, collapse = "; "))
res$reproduction <- list(all_ok = all(repro$ok), table = repro)

# ------------------------------------------------------------------ (1) LCNEC
lc <- F$LCNEC
crude_tab <- function(t3) t3$crude[, .(interval = LAB7[per], deaths = d, person_months = round(pm, 1), rate = rate)]
res$lcnec <- list(n = nrow(D$LCNEC), cancer_deaths_24m = sum(lc$crude$d), crude_rates = crude_tab(lc),
                  rr_6_9_vs_3_6 = ctr_list(lc$r69_36), rr_6_15_vs_0_6 = ctr_list(lc$r615_06), rr_15_24_vs_6_15 = ctr_list(lc$r1524_615),
                  rr_each_interval_vs_0_3 = all_rr(lc$m7),
                  period_2016_2023 = list(n = nrow(D$LCNEC[year >= 2016]), rr_6_9_vs_3_6 = ctr_list(lc1623$r69_36), rr_6_15_vs_0_6 = ctr_list(lc1623$r615_06), rr_15_24_vs_6_15 = ctr_list(lc1623$r1524_615)))
res$comparators_same_code_path <- lapply(c(NEC = "NEC", SCLC = "SCLC", CRC_small = "CRC_small", CRC_large = "CRC_large"), function(g)
  list(n = nrow(D[[g]]), crude_rates = crude_tab(F[[g]]), rr_6_9_vs_3_6 = ctr_list(F[[g]]$r69_36), rr_6_15_vs_0_6 = ctr_list(F[[g]]$r615_06), rr_15_24_vs_6_15 = ctr_list(F[[g]]$r1524_615)))

# ------------------------------------------------------------------ (2) shape tests
res$shape_tests <- list(NEC_vs_SCLC_all_years = sh_sclc, NEC_vs_LCNEC_all_years = sh_lc_all, NEC_vs_LCNEC_year_factor = sh_lc_fac,
                        NEC_vs_LCNEC_2016_2023 = sh_lc_1623, LCNEC_vs_SCLC_all_years_note = "dis = LCNEC; ratio is LCNEC/SCLC", LCNEC_vs_SCLC_all_years = sh_lung_lc_vs_sc)

# ------------------------------------------------------------------ (3) 2x2 morphology x organ
r2 <- list(CRC_small = F$CRC_small$r69_36, CRC_large = F$CRC_large$r69_36, lung_small = F$SCLC$r69_36, lung_large = F$LCNEC$r69_36)
morph_crc <- ratio_indep(r2$CRC_small, r2$CRC_large); morph_lung <- ratio_indep(r2$lung_small, r2$lung_large)
organ_small <- ratio_indep(r2$CRC_small, r2$lung_small); organ_large <- ratio_indep(r2$CRC_large, r2$lung_large)
inter <- ratio_indep(morph_crc, morph_lung)
res$morphology_by_organ <- list(rr_6_9_vs_3_6 = lapply(r2, ctr_list), n = list(CRC_small = nrow(D$CRC_small), CRC_large = nrow(D$CRC_large), lung_small = nrow(D$SCLC), lung_large = nrow(D$LCNEC)),
  ratios = list(small_over_large_colorectal = ctr_list(morph_crc), small_over_large_lung = ctr_list(morph_lung),
                colorectal_over_lung_small = ctr_list(organ_small), colorectal_over_lung_large = ctr_list(organ_large), morphology_x_organ_interaction = ctr_list(inter)),
  note = "Ratios of independent log RRs (disjoint cohorts), SE = sqrt(se1^2 + se2^2); colorectal models adjust for rectum, lung models do not")

# ------------------------------------------------------------------ markdown
add("# A12 -- Distant lung large-cell NEC (LCNEC) as comparator: piecewise-exponential rate ratios")
add("Date %s. %s; survival %s; sandwich %s. Script: R/A12_lcnec_rate_ratios.R (code path of A06 reused verbatim). Aggregate results only.", res$date, R.version.string, res$software$survival, res$software$sandwich)
add("\n## 0. Reproduction check (run first; script stops if any anchor fails)\n")
add("| Quantity | Source | Target | Reproduced | OK |\n|---|---|---|---|---|")
for (i in seq_len(nrow(repro))) add("| %s | %s | %s | %s | %s |", repro$quantity[i], repro$source[i], format(repro$target[i]), ifelse(repro$tol[i] == 0, format(repro$ours[i]), sprintf("%.4f", repro$ours[i])), ifelse(repro$ok[i], "yes", "NO"))
add("\nAll %d anchors reproduced: %s.", nrow(repro), ifelse(all(repro$ok), "yes", "NO"))

add("\n## 1. Distant lung LCNEC (n = %d; %d cancer deaths within 24 months)\n", nrow(D$LCNEC), sum(lc$crude$d))
add("**1a. Crude cancer-death rates per person-month by interval, with colorectal NEC and SCLC for reference**\n")
add("| Interval (months) | LCNEC deaths / PM | LCNEC rate | Colorectal NEC rate | SCLC rate | LCNEC adj RR vs 0-3 |\n|---|---|---|---|---|---|")
ar <- all_rr(lc$m7)
for (k in 1:7) add("| %s | %d / %.0f | %.4f | %.4f | %.4f | %s |", LAB7[k], lc$crude$d[k], lc$crude$pm[k], lc$crude$rate[k], F$NEC$crude$rate[k], F$SCLC$crude$rate[k],
                   ifelse(k == 1, "1 (ref)", sprintf("%.2f (%.2f-%.2f)", ar$RR_vs_0_3[k - 1], ar$lo_vs_0_3[k - 1], ar$hi_vs_0_3[k - 1])))
add("\n**1b. Adjusted rate ratios (age, age^2, sex, year; cluster-robust 95%% CI)**\n")
add("| Tumour (distant) | n | RR 6-9 vs 3-6 | RR 6-15 vs 0-6 | RR 15-24 vs 6-15 |\n|---|---|---|---|---|")
rowf <- function(lab, n, t3) add("| %s | %d | %s, P=%s | %s, P=%s | %s, P=%s |", lab, n, fmt(t3$r69_36), fp(t3$r69_36[["p"]]), fmt(t3$r615_06), fp(t3$r615_06[["p"]]), fmt(t3$r1524_615), fp(t3$r1524_615[["p"]]))
rowf("Lung LCNEC (8013)", nrow(D$LCNEC), lc)
rowf("Lung LCNEC, 2016-2023", nrow(D$LCNEC[year >= 2016]), lc1623)
rowf("Colorectal NEC (main)", nrow(D$NEC), F$NEC)
rowf("SCLC (8041-8045)", nrow(D$SCLC), F$SCLC)
add("\nFootnote: the colorectal NEC row is additionally adjusted for large-cell morphology and rectal site (covariate set of A06, needed to reproduce the A06 anchor); lung rows (LCNEC, SCLC) are adjusted for age, age^2, sex and year only.")

add("\n## 2. Hazard-shape test, colorectal NEC vs lung LCNEC\n")
add("Pooled piecewise model (age10 + age10sq + male + year10 + interval + tumour x interval); 6-df Wald test of the tumour x interval terms relative to 0-3 months; cluster-robust covariance.\n")
add("| Comparison | n (NEC vs comparator) | W (6 df) | P | Ratio of RR 6-9 vs 3-6, NEC / comparator |\n|---|---|---|---|---|")
shrow <- function(lab, s) add("| %s | %d vs %d | %.2f | %s | %s, P=%s |", lab, s$n_nec, s$n_cmp, s$W, fp(s$p), fmt(s$ratio_r69_36_NEC_over_cmp), fp(s$ratio_r69_36_NEC_over_cmp$p))
shrow("NEC vs SCLC, all years (anchor)", sh_sclc)
shrow("NEC vs LCNEC, all years", sh_lc_all)
shrow("NEC vs LCNEC, all years, year as 3-level period", sh_lc_fac)
shrow("NEC vs LCNEC, 2016-2023", sh_lc_1623)
add("| Lung: LCNEC vs SCLC, all years (ratio = LCNEC/SCLC) | %d vs %d | %.2f | %s | %s, P=%s |", sh_lung_lc_vs_sc$n_nec, sh_lung_lc_vs_sc$n_cmp, sh_lung_lc_vs_sc$W, fp(sh_lung_lc_vs_sc$p), fmt(sh_lung_lc_vs_sc$ratio_r69_36_NEC_over_cmp), fp(sh_lung_lc_vs_sc$ratio_r69_36_NEC_over_cmp$p))
add("\nNote on multiplicity: the LCNEC comparisons were added at revision and are exploratory. None of the tests in this table is adjusted for multiplicity; the NEC-vs-LCNEC shape test is reported in three unadjusted versions (all years, year as period factor, 2016-2023), and the 2016-2023 result is one of these three.")

add("\n## 3. Morphology x organ: adjusted RR 6-9 vs 3-6 months (distant disease)\n")
add("| | Small-cell (8041-8045) | Large-cell (8013) | Small / large |\n|---|---|---|---|")
add("| Colorectum | %s, n=%d | %s, n=%d | %s |", fmt(r2$CRC_small), nrow(D$CRC_small), fmt(r2$CRC_large), nrow(D$CRC_large), fmt(morph_crc))
add("| Lung | %s, n=%d | %s, n=%d | %s |", fmt(r2$lung_small), nrow(D$SCLC), fmt(r2$lung_large), nrow(D$LCNEC), fmt(morph_lung))
add("| Colorectum / lung | %s | %s | interaction %s, P=%s |", fmt(organ_small), fmt(organ_large), fmt(inter), fp(inter[["p"]]))
add("\nColorectal models: age10 + age10sq + male + year10 + rectum (as A04); lung models: age10 + age10sq + male + year10. Ratios combine independent log RRs from disjoint cohorts.")

# ---- interpretation (numbers inserted programmatically)
add("\n## 4. Interpretation (for Table 3 and the Discussion)\n")
lcr <- lc$r69_36
reb <- if (lcr[["lo"]] > 1) "a statistically clear rise" else if (lcr[["RR"]] > 1) "a rise whose confidence interval includes 1" else "no rise"
add("1. Distant lung LCNEC (n = %d) shows %s in cancer-death rate from 3-6 to 6-9 months: adjusted RR %s (P = %s); crude rates %.4f, %.4f and %.4f per person-month at 0-3, 3-6 and 6-9 months.",
    nrow(D$LCNEC), reb, fmt(lcr), fp(lcr[["p"]]), lc$crude$rate[1], lc$crude$rate[2], lc$crude$rate[3])
add("2. For comparison the same model gives %s for colorectal NEC and %s for SCLC; beyond 6 months LCNEC rates follow RR 6-15 vs 0-6 = %s and RR 15-24 vs 6-15 = %s (colorectal NEC %s and %s).",
    fmt(F$NEC$r69_36), fmt(F$SCLC$r69_36), fmt(lc$r615_06), fmt(lc$r1524_615), fmt(F$NEC$r615_06), fmt(F$NEC$r1524_615))
add("3. In exploratory comparisons added at revision (not adjusted for multiplicity), the 6-df shape test %s between colorectal NEC and LCNEC over all years (W = %.2f, P = %s) and in 2016-2023 (W = %.2f, P = %s; n = %d vs %d; one of three unadjusted versions, with W = %.2f, P = %s when year enters as a 3-level period factor). The NEC/LCNEC ratio of the 6-9 vs 3-6 RRs is %s (P = %s) over all years and %s in 2016-2023.",
    ifelse(sh_lc_all$p < 0.05, "indicates a difference in hazard shape", "does not clearly indicate a difference in hazard shape"), sh_lc_all$W, fp(sh_lc_all$p), sh_lc_1623$W, fp(sh_lc_1623$p), sh_lc_1623$n_nec, sh_lc_1623$n_cmp, sh_lc_fac$W, fp(sh_lc_fac$p),
    fmt(sh_lc_all$ratio_r69_36_NEC_over_cmp), fp(sh_lc_all$ratio_r69_36_NEC_over_cmp$p), fmt(sh_lc_1623$ratio_r69_36_NEC_over_cmp))
add("4. In the 2 x 2 morphology x organ layout the 6-9 vs 3-6 RR is %s (colorectal small-cell), %s (colorectal large-cell), %s (lung small-cell) and %s (lung large-cell). Small/large ratios are %s in the colorectum and %s in the lung; the morphology x organ interaction is %s (P = %s).",
    fmt(r2$CRC_small), fmt(r2$CRC_large), fmt(r2$lung_small), fmt(r2$lung_large), fmt(morph_crc), fmt(morph_lung), fmt(inter), fp(inter[["p"]]))
add("5. Within the lung, LCNEC differs from SCLC in hazard shape (W = %.2f, P = %s; LCNEC/SCLC ratio of 6-9 vs 3-6 RRs %s).",
    sh_lung_lc_vs_sc$W, fp(sh_lung_lc_vs_sc$p), fmt(sh_lung_lc_vs_sc$ratio_r69_36_NEC_over_cmp))
add("6. Restricted to 2016-2023, LCNEC (n = %d) gives RR 6-9 vs 3-6 = %s.", nrow(D$LCNEC[year >= 2016]), fmt(lc1623$r69_36))
syn7 <- sprintf("7. Synthesis: the 6-9-month rebound %s in lung LCNEC (RR %s). Across the 2 x 2 layout the point estimates are consistent with the rebound being associated with small-cell morphology, but the evidence rests mainly on the lung small- vs large-cell contrast (%s). The colorectal subtype contrast is imprecise (%s), and an organ effect comparable to the morphology effect is not excluded (colorectal/lung %s for small-cell; %s for large-cell). The interaction estimate (%s) is too imprecise to show that the morphology effect is the same in both organs.",
    ifelse(lcr[["lo"]] > 1, "is present", "was not observed"), fmt(lcr), fmt(morph_lung), fmt(morph_crc), fmt(organ_small), fmt(organ_large), fmt(inter))
add("%s", syn7)
add("8. Caveats: SEER does not record immunohistochemical confirmation of large-cell neuroendocrine morphology, so lung 8013 may include some misclassified non-small-cell carcinomas; treatment regimens (SCLC-type vs NSCLC-type chemotherapy) are not distinguished and chemotherapy is recorded only as yes/no-unknown; the colorectal large-cell estimate rests on n = %d. Comparator models do not adjust for treatment.", nrow(D$CRC_large))

add("\n## Revision after audit (2026-10-08)\n")
add("- Wording only; no estimate was recomputed differently. The script was re-run from scratch and all %d reproduction anchors pass.", nrow(repro))
add("- Item 7 (Synthesis) rewritten: removed the claims that the morphology effect is 'of similar size in both organs' and that the rebound follows small-cell morphology 'rather than organ of origin'. A non-significant, wide interaction (%s) cannot show similarity; the colorectal small/large contrast includes 1 and organ effects of comparable size are not excluded. The new wording states that the point estimates are consistent with a small-cell association, resting mainly on the lung contrast.", fmt(inter))
add("- Section 2 and item 3: the LCNEC comparisons are labelled exploratory, added at revision and not adjusted for multiplicity; the 2016-2023 NEC-vs-LCNEC test is identified as one of three unadjusted versions. 'Separates' replaced by 'indicates a difference in hazard shape'.")
add("- Table 1b: footnote added stating that the colorectal NEC row is additionally adjusted for large-cell morphology and rectal site (as A06), while lung rows are adjusted for age, age^2, sex and year only.")
add("- Downstream (for manuscript integration, not changed here): Manuscript v2 currently says that an RR for lung LCNEC 'was not computed' (Figure 3D legend) and that the rebound 'was not analyzed in lung LCNEC' (Discussion). When A12 is integrated, the Figure 3D legend, Table 3 and the Discussion paragraph on small-cell morphology should be updated using the hedged wording of item 7.")
res$interpretation <- list(synthesis = syn7, exploratory_note = "LCNEC comparisons added at revision; exploratory; not adjusted for multiplicity; 2016-2023 NEC-vs-LCNEC test is one of three unadjusted versions.",
  table1b_footnote = "Colorectal NEC row additionally adjusted for large-cell morphology and rectal site (as A06); lung rows adjusted for age, age^2, sex and year only.")
res$revision_after_audit <- list(date = "2026-10-08", changes = c("Synthesis (item 7) reworded to hedged form; removed 'similar size in both organs' and 'rather than organ of origin'",
  "Exploratory / no-multiplicity-adjustment label added to section 2 and item 3", "Footnote on covariate sets added to table 1b",
  "Downstream note: update manuscript Figure 3D legend, Table 3 and Discussion on integration"), numbers_changed = FALSE)
writeLines(md, OUT_MD)
res$runtime_sec <- as.numeric(difftime(Sys.time(), T_START, units = "secs"))
write_json(res, OUT_JSON, auto_unbox = TRUE, pretty = TRUE, digits = NA)
cat("Wrote", OUT_JSON, "and", OUT_MD, "in", round(res$runtime_sec), "s\n")
