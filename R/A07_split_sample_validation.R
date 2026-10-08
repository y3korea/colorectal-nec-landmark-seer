# A07 split-sample validation: adjusted piecewise-exponential RR (6-9 vs 3-6 months) by calendar period, distant main NEC.
# Poisson on person-time split at 3,6,9,12,15,18 (follow-up 0-24 m), cluster-robust SE by patient, covariates COV$NEC.
# Writes <results>/A07_rr.json; the Python driver (A07_split_sample_validation.py) merges it into A07.json.
t0 <- Sys.time()
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "config.R"))  # paths, NEC_BOOT, NEC_STRICT
source(file.path(NEC_R_DIR, "pp.R"))
A <- nec_main(loadc("NEN_all", "NEC"))[stage == "Distant"]
CUTS <- c(3, 6, 9, 12, 15, 18)
fitrr <- function(dt, extra = "") {
  s <- pw(dt, 0, CUTS, 24); s[, per := factor(per)]
  f <- glm(as.formula(paste("d ~ per +", COV$NEC, extra, "+ offset(log(expo))")), family = poisson, data = s)
  V <- vcovCL(f, cluster = s$pid)
  list(fit = f, V = V, s = s)
}
summ <- function(dt, label) {
  r <- fitrr(dt)
  c69 <- ctr(r$fit, r$V, vec(r$fit, list(per3 = 1, per2 = -1)))
  c03 <- ctr(r$fit, r$V, vec(r$fit, list(per2 = -1)))          # 0-3 vs 3-6
  c912 <- ctr(r$fit, r$V, vec(r$fit, list(per4 = 1, per2 = -1)))
  deaths <- r$s[, .(deaths = sum(d), pm = sum(expo)), by = per]
  list(label = label, n = nrow(dt), years = range(dt$year),
       RR_6to9_vs_3to6 = as.list(c69[c("RR", "lo", "hi", "p")]),
       RR_0to3_vs_3to6 = as.list(c03[c("RR", "lo", "hi", "p")]),
       RR_9to12_vs_3to6 = as.list(c912[c("RR", "lo", "hi", "p")]),
       interval_deaths = deaths$deaths, interval_person_months = deaths$pm)
}
# reproduction anchor (full cohort)
full <- summ(A, "all 2000-2023")
cat(sprintf("ANCHOR adjusted RR 6-9 vs 3-6 (all): %.2f (%.2f-%.2f) P=%.4f  [manuscript 1.49 (1.17-1.90) P=0.0013]\n",
            full$RR_6to9_vs_3to6$RR, full$RR_6to9_vs_3to6$lo, full$RR_6to9_vs_3to6$hi, full$RR_6to9_vs_3to6$p))
nec_check(abs(full$RR_6to9_vs_3to6$RR - 1.49) < 0.01, abs(full$RR_6to9_vs_3to6$lo - 1.17) < 0.01, abs(full$RR_6to9_vs_3to6$hi - 1.90) < 0.01, nrow(A) == 908)

out <- list(anchor = full,
            derivation_2000_2015 = summ(A[year <= 2015], "2000-2015"),
            confirmation_2016_2023 = summ(A[year >= 2016], "2016-2023"),
            era_2000_2007 = summ(A[year <= 2007], "2000-2007"),
            era_2008_2015 = summ(A[year >= 2008 & year <= 2015], "2008-2015"),
            era_2016_2023 = summ(A[year >= 2016], "2016-2023"))
# period x interval interaction: does the 6-9 vs 3-6 RR differ between periods? (ratio of RRs, cluster-robust Wald)
A[, late := as.integer(year >= 2016)]
r <- fitrr(A, "+ late + per:late")
rat <- ctr(r$fit, r$V, vec(r$fit, list(`per3:late` = 1, `per2:late` = -1)))
out$ratio_of_RR_late_vs_early <- as.list(rat[c("RR", "lo", "hi", "p")])
out$ratio_of_RR_late_vs_early$note <- "RR(6-9 vs 3-6) in 2016-2023 divided by RR in 2000-2015, pooled Poisson with per*late interaction, cluster-robust SE"
out$runtime_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
for (nm in names(out)) if (is.list(out[[nm]]) && !is.null(out[[nm]]$RR_6to9_vs_3to6))
  cat(sprintf("%-22s n=%4d  RR 6-9 vs 3-6 = %.2f (%.2f-%.2f) P=%.4f\n", nm, out[[nm]]$n, out[[nm]]$RR_6to9_vs_3to6$RR, out[[nm]]$RR_6to9_vs_3to6$lo, out[[nm]]$RR_6to9_vs_3to6$hi, out[[nm]]$RR_6to9_vs_3to6$p))
cat(sprintf("ratio of RR (late/early) = %.2f (%.2f-%.2f) P=%.3f\n", rat["RR"], rat["lo"], rat["hi"], rat["p"]))
write_json(out, paste0(RES, "A07_rr.json"), auto_unbox = TRUE, digits = 8, pretty = TRUE)
