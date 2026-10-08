# A05 (part 1/2): standardized distant-NEC landmark curves with original and extended covariate sets.
# Run BEFORE A05_heterogeneity_stratified.py (which merges this output into results/A05.json).
t0 <- Sys.time()
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "config.R"))  # paths, NEC_BOOT, NEC_STRICT
source(file.path(NEC_R_DIR, "pp.R"))
OUT <- paste0(RES, "A05_std.json")
nen <- loadc("NEN_all", "nec"); nec <- nec_main(nen); d <- nec[stage == "Distant"]
d[, nhw := as.integer(race == "Non-Hispanic White")]
d[, liver := as.integer(mets_liver == "Yes")]
cat("distant main NEC n =", nrow(d), " surgery yes/no =", sum(d$surg), "/", sum(d$surg_no), "\n")
nec_check(nrow(d) == 908, sum(d$surg) == 354, sum(d$surg_no) == 554)
LM <- 0:18
run <- function(dt, covf, label) {
  y <- std_curve(dt, covf, ref = dt, LM = LM)
  list(label = label, covariates = covf, n = nrow(dt), landmarks = LM, std_curve = round(y, 2),
       at_0_6_12 = round(y[c(1, 7, 13)], 2), plateau_end = hinge(y, LM))
}
res <- list()
res$original <- run(d, COV$NEC, "original: age10+age10sq+male+year10+large+rectum")
cat(sprintf("REPRODUCTION (std curve, original covariates) 0/6/12 = %.2f/%.2f/%.2f (manuscript 72.11/71.49/58.31); plateau %.2f (manuscript 7.5)\n",
            res$original$at_0_6_12[1], res$original$at_0_6_12[2], res$original$at_0_6_12[3], res$original$plateau_end))
ok <- all(abs(res$original$at_0_6_12 - c(72.11, 71.49, 58.31)) <= 0.05) && res$original$plateau_end == 7.5
if (!ok) nec_anchor_fail("Reproduction check FAILED for standardized curve")
res$plus_surg <- run(d, paste(COV$NEC, "+ surg"), "original + surg (primary-site surgery; unknown=0, none in distant NEC)")
res$plus_nhw  <- run(d, paste(COV$NEC, "+ nhw"),  "original + NHW race indicator")
res$plus_surg_nhw <- run(d, paste(COV$NEC, "+ surg + nhw"), "original + surg + NHW")
d10 <- d[year >= 2010 & mets_liver %in% c("Yes", "No")]
cat("2010+ subset with known liver-metastasis status n =", nrow(d10), " liver yes =", sum(d10$liver), "\n")
res$sub2010_original <- run(d10, COV$NEC, "2010+ subset (known liver status), original covariates")
res$sub2010_plus_liver <- run(d10, paste(COV$NEC, "+ liver"), "2010+ subset, original + liver-metastasis indicator")
res$sub2010_plus_liver_surg_nhw <- run(d10, paste(COV$NEC, "+ liver + surg + nhw"), "2010+ subset, original + liver + surg + NHW")
for (nm in names(res)) cat(sprintf("%-45s n=%4d  0/6/12 = %6.2f/%6.2f/%6.2f  plateau = %5.2f\n", nm, res[[nm]]$n, res[[nm]]$at_0_6_12[1], res[[nm]]$at_0_6_12[2], res[[nm]]$at_0_6_12[3], res[[nm]]$plateau_end))
res$meta <- list(reproduction_ok = ok, runtime_sec = as.numeric(difftime(Sys.time(), t0, units = "secs")),
                 note = "std_curve(): per-landmark multi-state Cox (12-month window, cause-specific hazards for cancer and other death), predicted cancer-death probability averaged over the cohort's own covariate distribution (direct standardization). Plateau end = broken-stick hinge on landmarks 0-18 (grid 0-15 by 0.25).")
write(toJSON(res, auto_unbox = TRUE, digits = 6, pretty = TRUE), OUT)
cat("wrote", OUT, " runtime", round(res$meta$runtime_sec, 1), "s\n")
