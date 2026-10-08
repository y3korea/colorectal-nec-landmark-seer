# A04 (R part): adjusted piecewise-Poisson rate ratios (cluster-robust SE) for case-definition analyses.
# Cohorts: main NEC (distant), expanded NEC (main + 8246 G3/4), 8246-G3/4-only, small-cell (8041-8045) and large-cell (8013) subsets.
# Run: Rscript R/A04_expanded_cohort_codes_subtypes.R  -> writes results/A04_rr.json (merged into A04.json by python/A04_expanded_cohort_codes_subtypes.py)
t0 <- Sys.time()
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "config.R"))  # paths, NEC_BOOT, NEC_STRICT
source(file.path(NEC_R_DIR, "pp.R"))
nen <- loadc("NEN_all", "NEN")
nen[, g3 := grade17 %in% c("Poorly differentiated; Grade III","Undifferentiated; anaplastic; Grade IV") | grade18 %in% c("Site-specific grade system category (3)","Site-specific grade system category (4)")]
nen[, grp := fifelse(hist %in% SCLC_CODES, "small", fifelse(hist == "8013", "large", fifelse(hist == "8246" & g3, "nos8246_g34", "other")))]
nen[, grp := factor(grp, levels = c("small", "large", "nos8246_g34", "other"))]
D <- nen[stage == "Distant"]
cohorts <- list(
  main        = list(dt = nec_main(D),                     cov = "age10 + age10sq + male + year10 + large + rectum"),
  main_grpfac = list(dt = nec_main(D),                     cov = "age10 + age10sq + male + year10 + rectum + grp"),
  expanded    = list(dt = nec_expanded(D),                 cov = "age10 + age10sq + male + year10 + rectum + grp"),
  nos8246_g34 = list(dt = D[grp == "nos8246_g34"],         cov = "age10 + age10sq + male + year10 + rectum"),
  small_cell  = list(dt = D[grp == "small"],               cov = "age10 + age10sq + male + year10 + rectum"),
  large_cell  = list(dt = D[grp == "large"],               cov = "age10 + age10sq + male + year10 + rectum"),
  expanded_2018_2023 = list(dt = nec_expanded(D)[year >= 2018], cov = "age10 + age10sq + male + year10 + rectum + grp"),
  main_2018_2023     = list(dt = nec_main(D)[year >= 2018],     cov = "age10 + age10sq + male + year10 + large + rectum")
)
# Poisson piecewise model on the split data; returns RR for contrast period a vs period b (levels of `per`)
fit_pw <- function(dt, cov, cuts, end, contrasts) {
  s <- pw(dt, 0, cuts, end)
  s[, grp := droplevels(grp)]
  f <- as.formula(paste("d ~ factor(per) +", cov, "+ offset(log(expo))"))
  fit <- glm(f, family = poisson, data = s)
  V <- vcovCL(fit, cluster = s$pid)
  out <- list()
  for (nm in names(contrasts)) {
    ab <- contrasts[[nm]]; v <- vec(fit, setNames(list(1, -1), c(paste0("factor(per)", ab[1]), paste0("factor(per)", ab[2]))))
    # period 1 is the reference (coefficient absent): handle by dropping zero-named terms
    v <- v[names(coef(fit))]; v[is.na(v)] <- 0
    if (ab[2] == 1) { v <- vec(fit, setNames(list(1), paste0("factor(per)", ab[1]))) }
    r <- ctr(fit, V, v)
    out[[nm]] <- list(RR = unname(r["RR"]), lo = unname(r["lo"]), hi = unname(r["hi"]), p = unname(r["p"]),
                      deaths = c(s[per == ab[1], sum(d)], s[per == ab[2], sum(d)]), person_months = c(s[per == ab[1], sum(expo)], s[per == ab[2], sum(expo)]))
  }
  out
}
res <- list()
for (nm in names(cohorts)) {
  dt <- cohorts[[nm]]$dt; cov <- cohorts[[nm]]$cov
  cat(sprintf("%-20s n=%d\n", nm, nrow(dt)))
  # model A: 3-month intervals 0-3,3-6,6-9,9-12,12-15,15-18,18-24 (per 1..7)
  A <- tryCatch(fit_pw(dt, cov, cuts = c(3, 6, 9, 12, 15, 18), end = 24, contrasts = list(rr_6_9_vs_3_6 = c(3, 2), rr_9_12_vs_3_6 = c(4, 2), rr_3_6_vs_0_3 = c(2, 1))), error = function(e) list(error = conditionMessage(e)))
  # model B: 0-6, 6-15, 15-24 (per 1..3)
  B <- tryCatch(fit_pw(dt, cov, cuts = c(6, 15), end = 24, contrasts = list(rr_6_15_vs_0_6 = c(2, 1), rr_15_24_vs_6_15 = c(3, 2))), error = function(e) list(error = conditionMessage(e)))
  res[[nm]] <- list(n = nrow(dt), covariates = cov, model_3mo = A, model_0_6_15_24 = B)
  for (k in names(A)) if (!is.null(A[[k]]$RR)) cat(sprintf("  %-18s RR=%.3f (%.3f-%.3f) p=%.4g\n", k, A[[k]]$RR, A[[k]]$lo, A[[k]]$hi, A[[k]]$p))
  for (k in names(B)) if (!is.null(B[[k]]$RR)) cat(sprintf("  %-18s RR=%.3f (%.3f-%.3f) p=%.4g\n", k, B[[k]]$RR, B[[k]]$lo, B[[k]]$hi, B[[k]]$p))
}
# reproduction check: main cohort RR 6-9 vs 3-6 = 1.49 (1.17-1.90), P = 0.0013
chk <- res$main$model_3mo$rr_6_9_vs_3_6
ok <- abs(chk$RR - 1.49) < 0.01 && abs(chk$lo - 1.17) < 0.01 && abs(chk$hi - 1.90) < 0.01 && abs(chk$p - 0.0013) < 0.0002
cat(sprintf("REPRODUCTION CHECK (R): main RR 6-9 vs 3-6 = %.2f (%.2f-%.2f) P=%.4f  -> %s\n", chk$RR, chk$lo, chk$hi, chk$p, ifelse(ok, "OK", "MISMATCH")))
if (!ok) nec_anchor_fail("Reproduction check failed in R part")
res$reproduction_check <- list(expected = "1.49 (1.17-1.90) P=0.0013", observed = sprintf("%.2f (%.2f-%.2f) P=%.4f", chk$RR, chk$lo, chk$hi, chk$p), ok = ok)
res$runtime_seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
write_json(res, paste0(RES, "A04_rr.json"), auto_unbox = TRUE, digits = 8, pretty = TRUE)
cat("wrote A04_rr.json in", round(res$runtime_seconds, 1), "s\n")
