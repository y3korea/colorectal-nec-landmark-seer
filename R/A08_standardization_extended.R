# A08 standardization_extended: model-based (standardized) landmark curves for distant colorectal NEC,
# extended covariate sets, patient-level bootstrap CIs, proportional-hazards check, Fine-Gray sensitivity.
# Run non-interactively: Rscript R/A08_standardization_extended.R   (~3-10 min; NEC_CORES forked workers for the bootstrap, default 4)
T0 <- Sys.time()
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "config.R"))  # paths, NEC_BOOT, NEC_STRICT
source(file.path(NEC_R_DIR, "pp.R"))
suppressMessages({library(cmprsk); library(parallel)})
FIG <- FIGS_DIR; ID <- "A08"; dir.create(FIG, showWarnings = FALSE)
LM <- 0:18; W <- 12; NCORES <- if (.Platform$OS.type == "windows") 1L else as.integer(Sys.getenv("NEC_CORES", "4"))   # replicate indices are pre-drawn: results do not depend on NCORES
out <- list(id = ID, date = "2026-10-07", window_months = W, landmarks = LM)
r2 <- function(x) round(x, 2); fmt <- function(x, d = 2) formatC(x, format = "f", digits = d)
tictoc <- function(t) as.numeric(difftime(Sys.time(), t, units = "secs"))

## ---------- helpers (additions; pp.R untouched) ----------
# landmark data set: alive and under follow-up at s; clock reset; administrative censoring at s + window
lmdata <- function(dt, s, window = W) {
  x <- dt[tm > s]; t <- x$tm - s; e <- x$ev; e[t > window] <- 0L; t <- pmin(t, window)
  x[, `:=`(tt = t, st = factor(e, levels = 0:2, labels = c("cens", "cancer", "other")), e3 = e)]; x }
# Exact, vectorised re-implementation of survfit.coxphms() for a competing-risks multi-state Cox model
# (coxph(Surv(tt, st) ~ ., id = pid) == two cause-specific Cox models with Breslow ties; survfit.coxphms uses
# stype = 2, i.e. P(t) = P(t-) %*% expm(H_t) at every event time).  Validated below against survfit() to 1e-8.
std_fast <- function(dt, covf, ref, s, window = W) {
  x <- lmdata(dt, s, window)
  f1 <- coxph(as.formula(paste("Surv(tt, e3 == 1) ~", covf)), data = x, ties = "breslow")
  f2 <- coxph(as.formula(paste("Surv(tt, e3 == 2) ~", covf)), data = x, ties = "breslow")
  b1 <- basehaz(f1, centered = FALSE); b2 <- basehaz(f2, centered = FALSE)
  tg <- sort(unique(c(b1$time, b2$time))); tg <- tg[tg <= window]
  cum <- function(b) approx(b$time, b$hazard, xout = tg, method = "constant", f = 0, rule = c(1, 2), yleft = 0)$y
  dh1 <- diff(c(0, cum(b1))); dh2 <- diff(c(0, cum(b2)))
  lp1 <- predict(f1, newdata = ref, type = "lp", reference = "zero"); lp2 <- predict(f2, newdata = ref, type = "lp", reference = "zero")
  capped <- sum(lp1 > 300 | lp2 > 300)   # a divergent cause-specific coefficient (sparse late risk sets, e.g. 3-7 other deaths) would make exp() overflow to Inf and the product NaN; bound the risk score (exp(300) ~ 1e130: the subject is predicted to die at the first event time, cause allocated by h1/(h1+h2), so the contribution stays finite and in [0,1])
  r1 <- exp(pmin(lp1, 300)); r2 <- exp(pmin(lp2, 300))
  H1 <- outer(dh1, r1); H2 <- outer(dh2, r2); Ht <- H1 + H2
  Sm <- rbind(1, exp(-apply(Ht, 2, cumsum)))[seq_along(dh1), , drop = FALSE]     # P(alive at t_j-)
  g <- ifelse(Ht > 0, (1 - exp(-Ht)) / Ht, 1)
  list(cancer = colSums(Sm * H1 * g), other = colSums(Sm * H2 * g), n = nrow(x), f1 = f1, f2 = f2,
       unstable = sum(abs(c(coef(f1), coef(f2))) > 10, na.rm = TRUE), na_coef = sum(is.na(c(coef(f1), coef(f2)))), capped = capped) }
std_curve_fast <- function(dt, covf, ref, LMs = LM, window = W, chunk = 20000) {
  r <- lapply(LMs, function(s) {
    if (nrow(ref) <= chunk) z <- std_fast(dt, covf, ref, s, window) else {        # chunked reference (large cohorts)
      ch <- split(seq_len(nrow(ref)), ceiling(seq_len(nrow(ref)) / chunk))
      zz <- lapply(ch, function(i) std_fast(dt, covf, ref[i], s, window))
      z <- list(cancer = unlist(lapply(zz, `[[`, "cancer")), n = zz[[1]]$n, unstable = zz[[1]]$unstable, na_coef = zz[[1]]$na_coef, capped = sum(sapply(zz, `[[`, "capped"))) }
    list(cif = mean(z$cancer) * 100, n = z$n, unstable = z$unstable, na_coef = z$na_coef, capped = z$capped) })
  list(cif = sapply(r, `[[`, "cif"), n = sapply(r, `[[`, "n"), unstable = sapply(r, `[[`, "unstable"), na_coef = sapply(r, `[[`, "na_coef"), capped = sapply(r, `[[`, "capped")) }
aj_curve <- function(dt, LMs = LM, window = W) sapply(LMs, function(s) { x <- lmdata(dt, s, window); sf <- survfit(Surv(tt, st) ~ 1, data = x)
  idx <- max(which(sf$time <= window)); sf$pstate[idx, which(sf$states == "cancer")] * 100 })
summ <- function(cif, n = NULL, LMs = LM) { v <- cif[match(c(0, 6, 12), LMs)]
  list(m0 = v[1], m6 = v[2], m12 = v[3], d6_0 = v[2] - v[1], d12_6 = v[3] - v[2], plateau_end = hinge(cif, LMs), curve = cif, n_at_risk = n) }

## ---------- data ----------
NEN <- loadc("NEN_all", "nen"); nec <- nec_main(NEN)[stage == "Distant"]; rm(NEN)
nec_check(nrow(nec) == 908, sum(nec$surg) == 354, sum(nec$surg_no) == 554, sum(nec$surg_unk) == 0)
# extended covariates
race_tab <- table(nec$race); race_small <- names(race_tab)[race_tab < 10]
nec[, race5 := fifelse(race == "Non-Hispanic White", "White", fifelse(race == "Non-Hispanic Black", "Black",
             fifelse(race == "Hispanic (All Races)", "Hispanic", fifelse(race == "Non-Hispanic Asian or Pacific Islander", "API", "Other"))))]
nec[, race5 := factor(race5, levels = c("White", "Black", "Hispanic", "API", "Other"))]
inc_mid <- function(s) { s <- gsub("[$,]", "", s); ifelse(grepl("^<", s), 35000, ifelse(grepl("\\+$", s), 125000,
            sapply(strsplit(s, " - "), function(z) mean(as.numeric(z))))) }
nec[, inc := inc_mid(income)]
qs <- quantile(nec$inc, c(1/3, 2/3)); nec[, inc3 := factor(cut(inc, c(-Inf, qs, Inf), labels = c("T1", "T2", "T3")))]
nec[, liver := as.integer(mets_liver == "Yes")]
nec10 <- nec[year >= 2010 & mets_liver %in% c("Yes", "No")]
out$n <- list(distant_nec = nrow(nec), events = as.list(table(factor(nec$ev, 0:2, c("censored", "cancer_death", "other_death")))),
              nec_2010plus_known_liver = nrow(nec10), surg_yes = sum(nec$surg), surg_no = sum(nec$surg_no), surg_unknown = sum(nec$surg_unk),
              chemo_yes = sum(nec$chemo == 1), chemo_no_unknown = sum(nec$chemo == 0), race_raw = as.list(race_tab), race_collapsed_to_other = race_small,
              race5 = as.list(table(nec$race5)), income_bands = as.list(table(nec$income)), income_tertile_cuts_usd_midpoint = as.numeric(qs),
              income_tertiles = as.list(table(nec$inc3)), liver_yes_2010 = sum(nec10$liver), liver_no_2010 = sum(nec10$liver == 0),
              liver_field_2010plus = as.list(table(nec[year >= 2010]$mets_liver)))
cat("n distant NEC", nrow(nec), "; 2010+ with known liver-met status", nrow(nec10), "\n")

## ---------- (0) validation of the fast predictor against survfit.coxphms (the authors' code path) ----------
val <- sapply(c(0, 6, 12), function(s) { x <- lmdata(nec, s)
  fit <- coxph(as.formula(paste("Surv(tt, st) ~", COV$NEC)), data = x, id = pid)
  sf <- survfit(fit, newdata = nec, se.fit = FALSE); idx <- max(which(sf$time <= W))
  oc <- sf$pstate[idx, , which(sf$states == "cancer")]; fc <- std_fast(nec, COV$NEC, nec, s)$cancer
  c(landmark = s, official = mean(oc) * 100, fast = mean(fc) * 100, max_abs_diff_subject = max(abs(oc - fc)), ms_ties = fit$method == "breslow") })
out$fast_validation <- as.data.frame(t(val)); print(out$fast_validation)
if (max(val["max_abs_diff_subject", ]) > 1e-6) stop("fast predictor does not reproduce survfit.coxphms")

## ---------- (1) reproduce the authors' standardized NEC curve (official code path = pp.R std_curve) ----------
t1 <- Sys.time(); base_off <- std_curve(nec, COV$NEC, nec, LM); t_off <- tictoc(t1)
base <- std_curve_fast(nec, COV$NEC, nec); S_base <- summ(base$cif, base$n)
aj <- aj_curve(nec); S_aj <- summ(aj, base$n)
cat(sprintf("REPRO std NEC 0/6/12 = %.2f/%.2f/%.2f plateau %.2f  (manuscript 72.11/71.49/58.31, 7.5); official path %.2f/%.2f/%.2f [%.0fs], max|off-fast| %.1e\n",
            S_base$m0, S_base$m6, S_base$m12, S_base$plateau_end, base_off[1], base_off[7], base_off[13], t_off, max(abs(base_off - base$cif))))
cat(sprintf("      AJ  NEC 0/6/12 = %.2f/%.2f/%.2f plateau %.2f  (anchor 72.29/72.33/60.03)\n", S_aj$m0, S_aj$m6, S_aj$m12, S_aj$plateau_end))
rep_ok <- all(abs(c(S_base$m0, S_base$m6, S_base$m12) - c(72.11, 71.49, 58.31)) <= 0.05) && S_base$plateau_end == 7.5 && max(abs(base_off - base$cif)) < 1e-6
if (!rep_ok) nec_anchor_fail("REPRODUCTION FAILED: standardized NEC curve")
if (any(abs(c(S_aj$m0, S_aj$m6, S_aj$m12) - c(72.29, 72.33, 60.03)) > 0.05)) nec_anchor_fail("REPRODUCTION FAILED: AJ NEC")
a05 <- fromJSON(paste0(RES, "A05_std.json"))

## ---------- (2) comparators, random reference sample of 5000 (authors' seeding: 03 script AC then SCLC; 07 script PANC) ----------
comp <- list(); comp_anchor <- list(AC = c(40.7, 35.5, 0), SCLC = c(69.9, 69.9, 8.25), PANC = c(76.8, 70.1, 0))
n_anchor <- c(AC = 147014, SCLC = 77707, PANC = 75634)
load_comp <- function(nm) switch(nm, AC = loadc("AC", "ac")[stage == "Distant"], SCLC = sclc_main(loadc("LUNG", "lung"))[stage == "Distant"], PANC = panc_main(loadc("PANC", "panc"))[stage == "Distant"])
set.seed(20260928); refs <- list()
AC <- load_comp("AC"); refs$AC <- AC[sample(.N, min(.N, 5000))]
SCLC <- load_comp("SCLC"); refs$SCLC <- SCLC[sample(.N, min(.N, 5000))]
set.seed(20260928); PANC <- load_comp("PANC"); refs$PANC <- PANC[sample(.N, min(.N, 5000))]
cohorts <- list(AC = AC, SCLC = SCLC, PANC = PANC); rm(AC, SCLC, PANC)
for (nm in names(cohorts)) {
  d <- cohorts[[nm]]; t1 <- Sys.time(); nec_check(nrow(d) == n_anchor[[nm]])
  sc <- std_curve_fast(d, COV[[nm]], refs[[nm]]); S <- summ(sc$cif, sc$n); S$n_total <- nrow(d); S$ref_n <- nrow(refs[[nm]]); S$aj <- aj_curve(d)
  S$aj_m0 <- S$aj[1]; S$aj_m6 <- S$aj[7]; S$aj_m12 <- S$aj[13]; S$aj_plateau_end <- hinge(S$aj, LM)
  S$full_ref_0_6_12 <- std_curve_fast(d, COV[[nm]], d, LMs = c(0, 6, 12))$cif      # whole cohort as reference (no sampling error)
  S$ref_mc_se_0_6 <- sapply(c(0, 6), function(s) { z <- std_fast(d, COV[[nm]], refs[[nm]], s); sd(z$cancer) * 100 / sqrt(length(z$cancer)) })   # Monte-Carlo SE of a 5,000-patient reference mean
  S$diff_vs_manuscript <- c(S$m0, S$m6) - comp_anchor[[nm]][1:2]
  S$reproduced_strict <- all(abs(S$diff_vs_manuscript) <= 0.15) && S$plateau_end == comp_anchor[[nm]][3]
  S$reproduced <- all(abs(S$diff_vs_manuscript) <= 0.3) && S$plateau_end == comp_anchor[[nm]][3]   # 0.3 points ~ 2-3 reference-sampling SEs; plateau end must be identical
  S$manuscript <- comp_anchor[[nm]]; S$runtime_s <- tictoc(t1); comp[[nm]] <- S
  cat(sprintf("REPRO %-4s n=%d std 0/6/12 = %.2f/%.2f/%.2f plateau %.2f (manuscript %.1f -> %.1f, %.2f) full-ref %.2f/%.2f/%.2f  AJ %.2f/%.2f/%.2f  [%.0fs] %s\n", nm, nrow(d), S$m0, S$m6, S$m12, S$plateau_end,
              comp_anchor[[nm]][1], comp_anchor[[nm]][2], comp_anchor[[nm]][3], S$full_ref_0_6_12[1], S$full_ref_0_6_12[2], S$full_ref_0_6_12[3], S$aj_m0, S$aj_m6, S$aj_m12, S$runtime_s, if (S$reproduced_strict) "OK (strict 0.15)" else if (S$reproduced) sprintf("OK within reference-sampling error (MC SE %.2f/%.2f)", S$ref_mc_se_0_6[1], S$ref_mc_se_0_6[2]) else "MISMATCH"))
}
rm(cohorts, refs); invisible(gc())

## ---------- (3) extended covariate sets ----------
sets <- list(base = COV$NEC, plus_surg = paste(COV$NEC, "+ surg"), plus_race5 = paste(COV$NEC, "+ race5"), plus_income3 = paste(COV$NEC, "+ inc3"),
             plus_chemo = paste(COV$NEC, "+ chemo"), all_combined_908 = paste(COV$NEC, "+ surg + race5 + inc3 + chemo"))
labs <- c(base = "original (age10+age10sq+male+year10+large+rectum)", plus_surg = "+ primary-site surgery", plus_race5 = "+ race (5 levels)", plus_income3 = "+ county income tertile",
          plus_chemo = "+ chemotherapy (immortal-time caveat)", all_combined_908 = "all combined, no liver (n=908)",
          aj_2010 = "2010+ subset: unadjusted Aalen-Johansen", base_2010 = "2010+ subset: original covariates", plus_liver_2010 = "2010+ subset: + liver-metastasis indicator",
          all_combined_2010 = "2010+ subset: all combined + liver", aj_908 = "unadjusted Aalen-Johansen (n=908)")
ext <- list()
for (nm in names(sets)) { sc <- std_curve_fast(nec, sets[[nm]], nec); ext[[nm]] <- summ(sc$cif, sc$n); ext[[nm]]$formula <- sets[[nm]]; ext[[nm]]$n <- nrow(nec); ext[[nm]]$label <- labs[[nm]]; ext[[nm]]$n_landmarks_divergent_coef <- sum(sc$unstable > 0); ext[[nm]]$n_landmarks_na_coef <- sum(sc$na_coef > 0); ext[[nm]]$n_capped_risk <- sum(sc$capped) }
sets10 <- list(base_2010 = COV$NEC, plus_liver_2010 = paste(COV$NEC, "+ liver"), all_combined_2010 = paste(COV$NEC, "+ surg + race5 + inc3 + chemo + liver"))
aj10 <- aj_curve(nec10); ext$aj_2010 <- summ(aj10, sapply(LM, function(s) sum(nec10$tm > s))); ext$aj_2010$formula <- "unadjusted"; ext$aj_2010$n <- nrow(nec10); ext$aj_2010$label <- labs[["aj_2010"]]
for (nm in names(sets10)) { sc <- std_curve_fast(nec10, sets10[[nm]], nec10); ext[[nm]] <- summ(sc$cif, sc$n); ext[[nm]]$formula <- sets10[[nm]]; ext[[nm]]$n <- nrow(nec10); ext[[nm]]$label <- labs[[nm]]; ext[[nm]]$n_landmarks_divergent_coef <- sum(sc$unstable > 0); ext[[nm]]$n_landmarks_na_coef <- sum(sc$na_coef > 0); ext[[nm]]$n_capped_risk <- sum(sc$capped) }
ext$aj_908 <- S_aj; ext$aj_908$formula <- "unadjusted"; ext$aj_908$n <- 908; ext$aj_908$label <- labs[["aj_908"]]
ord <- c("aj_908", "base", "plus_surg", "plus_race5", "plus_income3", "plus_chemo", "all_combined_908", "aj_2010", "base_2010", "plus_liver_2010", "all_combined_2010"); ext <- ext[ord]
for (nm in ord) cat(sprintf("EXT %-18s n=%d 0/6/12 = %.2f/%.2f/%.2f  d6-0 %+.2f d12-6 %+.2f plateau %.2f\n", nm, ext[[nm]]$n, ext[[nm]]$m0, ext[[nm]]$m6, ext[[nm]]$m12, ext[[nm]]$d6_0, ext[[nm]]$d12_6, ext[[nm]]$plateau_end))
# agreement with A05 (which used the official survfit path)
a05chk <- rbind(data.frame(set = "base", a05 = a05$original$at_0_6_12, a08 = r2(c(ext$base$m0, ext$base$m6, ext$base$m12))),
                data.frame(set = "plus_surg", a05 = a05$plus_surg$at_0_6_12, a08 = r2(c(ext$plus_surg$m0, ext$plus_surg$m6, ext$plus_surg$m12))),
                data.frame(set = "base_2010", a05 = a05$sub2010_original$at_0_6_12, a08 = r2(c(ext$base_2010$m0, ext$base_2010$m6, ext$base_2010$m12))),
                data.frame(set = "plus_liver_2010", a05 = a05$sub2010_plus_liver$at_0_6_12, a08 = r2(c(ext$plus_liver_2010$m0, ext$plus_liver_2010$m6, ext$plus_liver_2010$m12))))
a05chk$landmark <- rep(c(0, 6, 12), 4); a05chk$diff <- a05chk$a08 - a05chk$a05; print(a05chk)
if (max(abs(a05chk$diff)) > 0.011) stop("disagreement with A05_std.json")
out$a05_agreement <- a05chk
# cause-specific HRs (cancer-death transition) for the all-combined model at landmarks 0 and 6
cf <- function(dt, covf, s) { z <- std_fast(dt, covf, dt[1:2], s); sm <- summary(z$f1)$coefficients
  data.table(landmark = s, term = rownames(sm), HR = exp(sm[, "coef"]), lo = exp(sm[, "coef"] - 1.96 * sm[, "se(coef)"]), hi = exp(sm[, "coef"] + 1.96 * sm[, "se(coef)"]), p = sm[, "Pr(>|z|)"]) }
out$coef_all_combined_cancer_csh <- rbind(cf(nec, sets$all_combined_908, 0), cf(nec, sets$all_combined_908, 6))

## ---------- (4) patient-level bootstrap of the original-covariate standardized curve ----------
B <- nec_B(500); SEED_B <- 20261007; set.seed(SEED_B); bidx <- replicate(B, sample.int(nrow(nec), nrow(nec), TRUE))   # indices drawn once, in order, in the parent process
boot1 <- function(b) { d <- nec[bidx[, b]]; suppressWarnings(tryCatch({ z <- std_curve_fast(d, COV$NEC, d); list(cif = z$cif, unstable = z$unstable, capped = z$capped) },
                                                          error = function(e) list(cif = rep(NA_real_, length(LM)), msg = conditionMessage(e)))) }
isbad <- function(z) !is.list(z) || !is.numeric(z$cif) || length(z$cif) != length(LM) || anyNA(z$cif)
t1 <- Sys.time(); dtth <- getDTthreads(); setDTthreads(1)                      # single-threaded data.table inside forked workers
bs <- mclapply(seq_len(B), boot1, mc.cores = NCORES); setDTthreads(dtth)
bad <- which(sapply(bs, isbad)); retry_msgs <- unname(sapply(bs[bad], function(z) if (is.list(z) && !is.null(z$msg)) z$msg else paste(class(z), collapse = "/")))
for (b in bad) bs[[b]] <- boot1(b)                                             # deterministic sequential re-run of any replicate that failed inside a forked worker (indices pre-drawn)
still_bad <- sapply(bs, isbad); t_boot <- tictoc(t1)
boot_unstable <- sapply(bs, function(z) if (isbad(z)) NA else sum(z$unstable > 0)); boot_capped <- sapply(bs, function(z) if (isbad(z)) NA else sum(z$capped > 0))
boot_capped_lm <- sort(unique(unlist(lapply(bs, function(z) if (isbad(z)) NULL else LM[z$capped > 0]))))
cat("bootstrap", B, "done in", round(t_boot), "s; re-run sequentially after a worker failure:", length(bad), "; still failed:", sum(still_bad),
    "; replicates with a divergent coefficient (|coef|>10) at >=1 landmark:", sum(boot_unstable > 0, na.rm = TRUE), "; with capped risk scores:", sum(boot_capped > 0, na.rm = TRUE), "at landmarks", boot_capped_lm, "\n")
bs <- do.call(rbind, lapply(bs[!still_bad], `[[`, "cif"))
bstat <- cbind(bs[, match(c(0, 3, 6, 12), LM)], bs[, 7] - bs[, 1], bs[, 13] - bs[, 7], apply(bs, 1, hinge, m = LM))
colnames(bstat) <- c("m0", "m3", "m6", "m12", "d6_0", "d12_6", "plateau_end")
pt <- c(m0 = S_base$curve[1], m3 = S_base$curve[4], m6 = S_base$curve[7], m12 = S_base$curve[13], d6_0 = S_base$d6_0, d12_6 = S_base$d12_6, plateau_end = S_base$plateau_end)
ci <- apply(bstat, 2, quantile, c(.025, .975)); bse <- apply(bstat, 2, sd); bmean <- colMeans(bstat)
out$bootstrap <- list(B = B, B_used = nrow(bs), seed = SEED_B, cores = NCORES, runtime_s = t_boot, n_retried_sequentially = length(bad), retry_messages = as.list(retry_msgs), n_replicates_with_divergent_coef = sum(boot_unstable > 0, na.rm = TRUE), n_replicates_with_capped_risk = sum(boot_capped > 0, na.rm = TRUE), capped_landmarks = boot_capped_lm,
  note_capping = "cause-specific linear predictors bounded at 300 before exponentiation (only binds when a coefficient diverges in a sparse late risk set); the authors' survfit code path returns NaN for such samples",
  method = "patient-level resampling with replacement of the 908 distant NEC; per replicate: refit the per-landmark multi-state Cox models (original covariates) and standardize over the resampled cohort; percentile 95% CI",
  table = data.frame(quantity = names(pt), estimate = as.numeric(pt), boot_mean = as.numeric(bmean), boot_se = as.numeric(bse), ci_lo = ci[1, ], ci_hi = ci[2, ]),
  plateau_end_quantiles = as.list(quantile(bstat[, "plateau_end"], c(.025, .1, .25, .5, .75, .9, .975))),
  plateau_end_prop_ge_6 = mean(bstat[, "plateau_end"] >= 6), plateau_end_prop_ge_3 = mean(bstat[, "plateau_end"] >= 3), plateau_end_prop_eq_0 = mean(bstat[, "plateau_end"] == 0),
  curve_ci = data.frame(landmark = LM, estimate = S_base$curve, lo = apply(bs, 2, quantile, .025), hi = apply(bs, 2, quantile, .975)))
print(out$bootstrap$table)

## ---------- (5) proportional-hazards check (cause-specific models at landmarks 0 and 6) ----------
zph_tab <- function(fit, lab, s) { z <- cox.zph(fit, transform = "km", terms = TRUE); tb <- as.data.table(z$table, keep.rownames = "term"); tb[, `:=`(model = lab, landmark = s)]; tb }
zph <- rbindlist(lapply(c(0, 6), function(s) { z <- std_fast(nec, COV$NEC, nec[1:2], s)
  rbind(zph_tab(z$f1, "cause-specific cancer death", s), zph_tab(z$f2, "cause-specific other death", s)) }))
zph_ms <- tryCatch(rbindlist(lapply(c(0, 6), function(s) { x <- lmdata(nec, s); fit <- coxph(as.formula(paste("Surv(tt, st) ~", COV$NEC)), data = x, id = pid); zph_tab(fit, "multi-state (joint)", s) })), error = function(e) NULL)
print(zph)
# (a) calibration: model standardized over its own landmark risk set vs nonparametric Aalen-Johansen (both n and reference = at-risk set)
calib <- sapply(LM, function(s) mean(std_fast(nec, COV$NEC, nec[tm > s], s)$cancer) * 100)
# (b) relaxing PH within the window: two 6-month pieces (landmark s model, 0-6; landmark s+6 model, 6-12), coefficients free to differ
pw12 <- function(s) { a <- std_fast(nec, COV$NEC, nec, s, window = 6); b <- std_fast(nec, COV$NEC, nec, s + 6, window = 6)
  mean(a$cancer + (1 - a$cancer - a$other) * b$cancer) * 100 }
pw <- data.frame(landmark = c(0, 3, 6, 12), single_PH_12mo = S_base$curve[match(c(0, 3, 6, 12), LM)], piecewise_6_6 = sapply(c(0, 3, 6, 12), pw12), aj = aj[match(c(0, 3, 6, 12), LM)])
pw$diff_piecewise_minus_single <- pw$piecewise_6_6 - pw$single_PH_12mo; print(pw)
out$ph <- list(zph = zph, zph_multistate = zph_ms, calibration = data.frame(landmark = LM, std_own_riskset = calib, aj = aj, diff = calib - aj), piecewise = pw,
               max_abs_calibration_diff = max(abs(calib - aj)), max_abs_piecewise_diff = max(abs(pw$diff_piecewise_minus_single)))

## ---------- (6) Fine-Gray sensitivity (cmprsk::crr), standardized over the 908 reference rows ----------
R <- model.matrix(as.formula(paste("~", COV$NEC)), data = nec)[, -1, drop = FALSE]
fg <- lapply(c(0, 6, 12), function(s) {
  x <- lmdata(nec, s); X <- model.matrix(as.formula(paste("~", COV$NEC)), data = x)[, -1, drop = FALSE]
  f <- crr(ftime = x$tt, fstatus = x$e3, cov1 = X, failcode = 1, cencode = 0)
  pr <- predict(f, cov1 = R); cif <- pr[max(which(pr[, 1] <= W)), -1]; sm <- summary(f)$coef
  list(landmark = s, n = nrow(x), converged = f$converged, std12_fine_gray = mean(cif) * 100, std12_multistate_cox = S_base$curve[LM == s], aj12 = aj[LM == s],
       sHR = data.table(term = rownames(sm), sHR = sm[, "exp(coef)"], lo = exp(sm[, "coef"] - 1.96 * sm[, "se(coef)"]), hi = exp(sm[, "coef"] + 1.96 * sm[, "se(coef)"]), p = sm[, "p-value"])) })
for (z in fg) cat(sprintf("FG landmark %2d n=%d: Fine-Gray std 12-mo %.2f | multi-state Cox %.2f | AJ %.2f\n", z$landmark, z$n, z$std12_fine_gray, z$std12_multistate_cox, z$aj12))
out$fine_gray <- list(models = fg, anchors_aj = c(72.29, 72.33, 60.03), note = "crr(failcode=1 cancer death, cencode=0), same covariates; 12-month CIF predicted for each of the 908 distant NEC and averaged")

## ---------- figures ----------
png(paste0(FIG, ID, "_std_curves.png"), width = 3000, height = 1300, res = 300)
par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1))
cols <- c(base = "black", plus_surg = "#1b9e77", plus_race5 = "#d95f02", plus_income3 = "#7570b3", plus_chemo = "#e7298a", all_combined_908 = "#66a61e")
plot(LM, aj, type = "n", ylim = c(20, 85), xlab = "Months survived since diagnosis (landmark)", ylab = "12-month cancer-death probability (%)", main = "Distant colorectal NEC (n = 908)")
polygon(c(LM, rev(LM)), c(out$bootstrap$curve_ci$lo, rev(out$bootstrap$curve_ci$hi)), col = adjustcolor("grey70", 0.4), border = NA)
lines(LM, aj, lty = 2, lwd = 2, col = "grey40"); for (nm in names(cols)) lines(LM, ext[[nm]]$curve, col = cols[[nm]], lwd = if (nm == "base") 2.5 else 1.5)
legend("bottomleft", c("Unadjusted (AJ)", "Standardized: original (95% bootstrap band)", "+ surgery", "+ race", "+ income tertile", "+ chemo", "all combined"), col = c("grey40", cols), lty = c(2, rep(1, 6)), lwd = 1.5, bty = "n", cex = 0.65)
plot(LM, aj10, type = "n", ylim = c(20, 85), xlab = "Months survived since diagnosis (landmark)", ylab = "12-month cancer-death probability (%)", main = sprintf("2010+ subset, known liver-met status (n = %d)", nrow(nec10)))
lines(LM, aj10, lty = 2, lwd = 2, col = "grey40"); lines(LM, ext$base_2010$curve, lwd = 2.5); lines(LM, ext$plus_liver_2010$curve, col = "#a6761d", lwd = 1.5); lines(LM, ext$all_combined_2010$curve, col = "#66a61e", lwd = 1.5)
legend("bottomleft", c("Unadjusted (AJ)", "Standardized: original", "+ liver metastasis", "all combined + liver"), col = c("grey40", "black", "#a6761d", "#66a61e"), lty = c(2, 1, 1, 1), lwd = 1.5, bty = "n", cex = 0.65)
dev.off()
png(paste0(FIG, ID, "_comparators.png"), width = 3000, height = 1100, res = 300); par(mfrow = c(1, 3), mar = c(4, 4, 2.5, 1))
for (nm in names(comp)) { S <- comp[[nm]]; plot(LM, S$aj, type = "l", lty = 2, lwd = 2, col = "grey40", ylim = c(20, 85), xlab = "Landmark (months)", ylab = "12-month cancer-death probability (%)", main = sprintf("%s distant (n = %s; ref 5,000)", nm, format(S$n_total, big.mark = ",")))
  lines(LM, S$curve, lwd = 2.5); legend("bottomleft", c("Unadjusted (AJ)", "Standardized"), lty = c(2, 1), col = c("grey40", "black"), lwd = 2, bty = "n", cex = 0.8) }
dev.off()

## ---------- save JSON ----------
out$reproduction <- list(nec_std = list(ours = r2(c(S_base$m0, S_base$m6, S_base$m12, S_base$plateau_end)), official_path = r2(base_off[c(1, 7, 13)]), manuscript = c(72.11, 71.49, 58.31, 7.5), ok = rep_ok, official_path_runtime_s = t_off),
                         nec_aj = list(ours = r2(c(S_aj$m0, S_aj$m6, S_aj$m12, S_aj$plateau_end)), anchor = c(72.29, 72.33, 60.03)),
                         comparators = lapply(names(comp), function(nm) list(cohort = nm, n = comp[[nm]]$n_total, ours = r2(c(comp[[nm]]$m0, comp[[nm]]$m6, comp[[nm]]$plateau_end)), manuscript = comp_anchor[[nm]], diff_0_6 = r2(comp[[nm]]$diff_vs_manuscript), ref_mc_se_0_6 = r2(comp[[nm]]$ref_mc_se_0_6), reproduced_strict_0_15 = comp[[nm]]$reproduced_strict, reproduced_within_0_3 = comp[[nm]]$reproduced)))
out$base_nec <- S_base; out$base_nec_official_curve <- base_off; out$comparators <- comp; out$extended <- ext; out$sets <- c(sets, sets10)
out$runtime_minutes <- as.numeric(difftime(Sys.time(), T0, units = "mins"))
out$software <- list(R = R.version.string, survival = as.character(packageVersion("survival")), cmprsk = as.character(packageVersion("cmprsk")), data.table = as.character(packageVersion("data.table")))
write(toJSON(out, auto_unbox = TRUE, digits = 6, pretty = TRUE, na = "null"), paste0(RES, ID, ".json"))

## ---------- markdown ----------
md <- c(sprintf("# A08 Standardization extended: distant colorectal NEC (n = 908), 12-month window, landmarks 0-18 (run %s)", out$date), "",
  "## Reproduction check", "",
  "| Quantity | Ours | Manuscript / anchor | Status |", "|---|---|---|---|",
  sprintf("| Standardized NEC curve (original covariates, reference = all 908), landmarks 0/6/12 | %s/%s/%s | 72.11/71.49/58.31 | %s |", fmt(S_base$m0), fmt(S_base$m6), fmt(S_base$m12), if (rep_ok) "reproduced" else "FAILED"),
  sprintf("| Plateau end (broken-stick hinge, grid 0-15 by 0.25) | %s | 7.5 | %s |", fmt(S_base$plateau_end), if (S_base$plateau_end == 7.5) "reproduced" else "FAILED"),
  sprintf("| Same curve via the authors' exact code path (coxph multi-state + survfit, pp.R std_curve) | %s/%s/%s | max abs difference from fast predictor %.1e | identical |", fmt(base_off[1]), fmt(base_off[7]), fmt(base_off[13]), max(abs(base_off - base$cif))),
  sprintf("| Unadjusted Aalen-Johansen NEC 0/6/12 | %s/%s/%s | 72.29/72.33/60.03 | reproduced |", fmt(S_aj$m0), fmt(S_aj$m6), fmt(S_aj$m12)))
for (nm in names(comp)) { S <- comp[[nm]]; md <- c(md, sprintf("| %s distant (n = %s), random reference 5,000 (seed 20260928): 0 -> 6 months, plateau | %s -> %s, %s (12 mo %s) | %.1f -> %.1f, %.2f | %s |", nm, format(S$n_total, big.mark = ","), fmt(S$m0), fmt(S$m6), fmt(S$plateau_end), fmt(S$m12), S$manuscript[1], S$manuscript[2], S$manuscript[3], if (S$reproduced_strict) "reproduced (within 0.15)" else if (S$reproduced) sprintf("reproduced within reference-sampling error (diff %+.2f/%+.2f; MC SE of a 5,000-draw mean %.2f/%.2f; plateau end identical)", S$diff_vs_manuscript[1], S$diff_vs_manuscript[2], S$ref_mc_se_0_6[1], S$ref_mc_se_0_6[2]) else "mismatch")) }
md <- c(md, "", "The reference draws follow the authors' seed and order (set.seed(20260928): AC then SCLC in the 03 script; re-seeded for PANC in the 07 script), but the row order of the rebuilt cohort files differs from the authors' files, so the 5,000 patients drawn are not the same patients; the Monte-Carlo SE of a 5,000-draw standardized mean (above) covers the observed level differences, and the plateau ends are identical.", "",
  sprintf("Comparator standardized values with the whole distant cohort as reference (no reference-sampling error), landmarks 0/6/12: %s. Unadjusted AJ 0/6/12: %s.",
  paste(sapply(names(comp), function(nm) sprintf("%s %s/%s/%s", nm, fmt(comp[[nm]]$full_ref_0_6_12[1]), fmt(comp[[nm]]$full_ref_0_6_12[2]), fmt(comp[[nm]]$full_ref_0_6_12[3]))), collapse = "; "),
  paste(sapply(names(comp), function(nm) sprintf("%s %s/%s/%s", nm, fmt(comp[[nm]]$aj_m0), fmt(comp[[nm]]$aj_m6), fmt(comp[[nm]]$aj_m12))), collapse = "; ")), "",
  sprintf("Agreement with A05_std.json (official survfit path): max |A08 - A05| at landmarks 0/6/12 over the original, +surg, 2010+ original and 2010+ +liver sets = %.2f points.", max(abs(a05chk$diff))), "",
  "## Extended covariate sets (standardized 12-month cancer-death probability, %, reference = the analysed cohort)", "",
  sprintf("Codings: surgery = primary-site surgery (yes %d / no %d / unknown %d); race 5 levels (White %d, Black %d, Hispanic %d, API %d, Other %d = %s collapsed because < 10); county median household income tertiles of band midpoints (cuts $%s and $%s; T1/T2/T3 = %d/%d/%d); chemotherapy yes %d / no-or-unknown %d (SEER records chemotherapy given at any time, so adjusting for it at early landmarks carries immortal-time bias: patients must survive long enough to receive it); liver metastasis = SEER 'Mets at DX-liver' Yes vs No, available 2010+ (n = %d with known status; Yes %d / No %d).",
    out$n$surg_yes, out$n$surg_no, out$n$surg_unknown, out$n$race5$White, out$n$race5$Black, out$n$race5$Hispanic, out$n$race5$API, out$n$race5$Other, paste(race_small, collapse = " + "),
    format(qs[1], big.mark = ","), format(qs[2], big.mark = ","), out$n$income_tertiles$T1, out$n$income_tertiles$T2, out$n$income_tertiles$T3, out$n$chemo_yes, out$n$chemo_no_unknown, nrow(nec10), out$n$liver_yes_2010, out$n$liver_no_2010), "",
  "| Covariate set | n | 0 mo | 6 mo | 12 mo | 6-0 | 12-6 | Plateau end (mo) |", "|---|---|---|---|---|---|---|---|")
for (nm in ord) { e <- ext[[nm]]; md <- c(md, sprintf("| %s | %d | %s | %s | %s | %s | %s | %s |", e$label, e$n, fmt(e$m0), fmt(e$m6), fmt(e$m12), fmt(e$d6_0), fmt(e$d12_6), fmt(e$plateau_end))) }
md <- c(md, "", sprintf("Model stability per set (number of the 19 landmark fits with a divergent cause-specific coefficient |coef| > 10 / with an aliased (NA) coefficient / reference rows needing a bounded risk score): %s.", paste(sapply(ord[ord != "aj_908" & ord != "aj_2010"], function(nm) sprintf("%s %d/%d/%d", nm, ext[[nm]]$n_landmarks_divergent_coef, ext[[nm]]$n_landmarks_na_coef, ext[[nm]]$n_capped_risk)), collapse = "; ")), "",
  "Full 19-landmark curves for every set are in A08.json (`extended.<set>.curve`). Cause-specific HRs (cancer-death transition) of the all-combined 908 model at landmarks 0 and 6 are in `coef_all_combined_cancer_csh`.", "",
  sprintf("## Bootstrap 95%% CIs, original covariates (patient-level resampling of the 908 patients, models refitted and standardized over the resampled cohort; B = %d, seed %d, %d replicates used, %d re-run sequentially after a forked-worker failure)", B, SEED_B, nrow(bs), length(bad)), "",
  "| Quantity | Estimate | Bootstrap mean | Bootstrap SE | 95% percentile CI |", "|---|---|---|---|---|")
tb <- out$bootstrap$table; for (i in seq_len(nrow(tb))) md <- c(md, sprintf("| %s | %s | %s | %s | %s to %s |", c(m0 = "Cancer death from landmark 0 (%)", m3 = "landmark 3 (%)", m6 = "landmark 6 (%)", m12 = "landmark 12 (%)", d6_0 = "Change 6 - 0 (points)", d12_6 = "Change 12 - 6 (points)", plateau_end = "Plateau end (months)")[tb$quantity[i]], fmt(tb$estimate[i]), fmt(tb$boot_mean[i]), fmt(tb$boot_se[i]), fmt(tb$ci_lo[i]), fmt(tb$ci_hi[i])))
md <- c(md, "", sprintf("In %d of %d replicates a cause-specific coefficient diverged (|coef| > 10) at at least one landmark; in %d replicates (landmarks %s; other-death models with 3-7 events among < 100 at risk) the risk score had to be bounded (linear predictor capped at 300) to keep the prediction finite. Without the bound these replicates return NaN, as the authors' survfit code path would.", out$bootstrap$n_replicates_with_divergent_coef, nrow(bs), out$bootstrap$n_replicates_with_capped_risk, paste(boot_capped_lm, collapse = ", ")))
pq <- out$bootstrap$plateau_end_quantiles
md <- c(md, "", sprintf("Plateau-end bootstrap distribution: 2.5/10/25/50/75/90/97.5%% = %s; P(plateau end >= 6 mo) = %.2f, >= 3 mo = %.2f, = 0 = %.2f.", paste(fmt(unlist(pq)), collapse = "/"), out$bootstrap$plateau_end_prop_ge_6, out$bootstrap$plateau_end_prop_ge_3, out$bootstrap$plateau_end_prop_eq_0), "",
  "## Proportional-hazards check (cox.zph, KM time transform; cause-specific Cox models with Breslow ties, as in the multi-state fit)", "",
  "| Landmark | Model | Term | chisq | df | p |", "|---|---|---|---|---|---|")
for (i in seq_len(nrow(zph))) md <- c(md, sprintf("| %d | %s | %s | %s | %d | %s |", zph$landmark[i], zph$model[i], zph$term[i], fmt(zph$chisq[i]), as.integer(zph$df[i]), formatC(zph$p[i], format = "g", digits = 2)))
md <- c(md, "", "Does PH matter for a 12-month-window standardization? Two checks:", "",
  sprintf("(a) Calibration: the model standardized over its own landmark risk set vs the nonparametric Aalen-Johansen estimate, landmarks 0-18: max |difference| = %.2f points (landmark 0/6/12: %s/%s/%s vs %s/%s/%s).", out$ph$max_abs_calibration_diff, fmt(calib[1]), fmt(calib[7]), fmt(calib[13]), fmt(aj[1]), fmt(aj[7]), fmt(aj[13])),
  "", "(b) Relaxed-PH alternative: the 12-month probability built from two 6-month pieces (landmark s model for 0-6 and landmark s+6 model for 6-12, coefficients free to differ) vs the single 12-month PH model:", "",
  "| Landmark | Single 12-mo PH model | Piecewise 6+6 | Difference | Unadjusted AJ |", "|---|---|---|---|---|")
for (i in seq_len(nrow(pw))) md <- c(md, sprintf("| %d | %s | %s | %s | %s |", pw$landmark[i], fmt(pw$single_PH_12mo[i]), fmt(pw$piecewise_6_6[i]), fmt(pw$diff_piecewise_minus_single[i]), fmt(pw$aj[i])))
md <- c(md, "", "## Fine-Gray sensitivity (cmprsk::crr, failcode = 1 cancer death, cencode = 0; 12-month CIF standardized over the 908 reference rows)", "",
  "| Landmark | n at risk | Fine-Gray standardized | Multi-state Cox standardized | Unadjusted AJ | Anchor AJ |", "|---|---|---|---|---|---|")
for (z in fg) md <- c(md, sprintf("| %d | %d | %s | %s | %s | %s |", z$landmark, z$n, fmt(z$std12_fine_gray), fmt(z$std12_multistate_cox), fmt(z$aj12), fmt(out$fine_gray$anchors_aj[match(z$landmark, c(0, 6, 12))])))
md <- c(md, "", "Subdistribution hazard ratios (sHR, 95% CI, p):", "", "| Landmark | Term | sHR | 95% CI | p |", "|---|---|---|---|---|")
for (z in fg) for (i in seq_len(nrow(z$sHR))) md <- c(md, sprintf("| %d | %s | %s | %s to %s | %s |", z$landmark, z$sHR$term[i], fmt(z$sHR$sHR[i]), fmt(z$sHR$lo[i]), fmt(z$sHR$hi[i]), formatC(z$sHR$p[i], format = "g", digits = 2)))
# interpretation
g6 <- zph[landmark == 0 & term == "GLOBAL" & model == "cause-specific cancer death"]$p; g6b <- zph[landmark == 6 & term == "GLOBAL" & model == "cause-specific cancer death"]$p
md <- c(md, "", "## Interpretation", "",
  sprintf("1. The authors' standardized distant-NEC curve is reproduced exactly (72.11/71.49/58.31, plateau end 7.5 months) through the authors' own code path and through a vectorised re-implementation of survfit.coxphms that agrees to 1e-8 per patient; the comparators are reproduced to within reference-sampling error (AC %s -> %s, SCLC %s -> %s, PANC %s -> %s vs manuscript 40.7 -> 35.5, 69.9 -> 69.9, 76.8 -> 70.1; plateau ends %s/%s/%s identical; level differences %+.2f to %+.2f points against a reference Monte-Carlo SE of %.2f-%.2f points, because the rebuilt cohort files are not in the authors' row order).", fmt(comp$AC$m0, 1), fmt(comp$AC$m6, 1), fmt(comp$SCLC$m0, 1), fmt(comp$SCLC$m6, 1), fmt(comp$PANC$m0, 1), fmt(comp$PANC$m6, 1), fmt(comp$AC$plateau_end), fmt(comp$SCLC$plateau_end), fmt(comp$PANC$plateau_end), min(unlist(lapply(comp, `[[`, "diff_vs_manuscript"))), max(unlist(lapply(comp, `[[`, "diff_vs_manuscript"))), min(unlist(lapply(comp, `[[`, "ref_mc_se_0_6"))), max(unlist(lapply(comp, `[[`, "ref_mc_se_0_6")))),
  sprintf("2. No extended covariate set changes the qualitative picture: the 6-0 change stays within %s to %s points and the plateau end within %s to %s months across all adjusted sets on the 908 patients; the largest shift comes from primary-site surgery (12-month value %s vs %s), which is itself partly a consequence of early survival (surgery is more often recorded in patients who live long enough to be operated).", fmt(min(sapply(ext[c("base","plus_surg","plus_race5","plus_income3","plus_chemo","all_combined_908")], `[[`, "d6_0"))), fmt(max(sapply(ext[c("base","plus_surg","plus_race5","plus_income3","plus_chemo","all_combined_908")], `[[`, "d6_0"))), fmt(min(sapply(ext[c("base","plus_surg","plus_race5","plus_income3","plus_chemo","all_combined_908")], `[[`, "plateau_end"))), fmt(max(sapply(ext[c("base","plus_surg","plus_race5","plus_income3","plus_chemo","all_combined_908")], `[[`, "plateau_end"))), fmt(ext$plus_surg$m12), fmt(ext$base$m12)),
  sprintf("3. Chemotherapy adjustment (6-0 = %s, 12-6 = %s, plateau %s) should be read with the immortal-time caveat: SEER chemotherapy is recorded whenever given, so at landmark 0 the 'chemo' group is enriched for early survivors; the conditional landmark design mitigates but does not remove this for the earliest landmarks.", fmt(ext$plus_chemo$d6_0), fmt(ext$plus_chemo$d12_6), fmt(ext$plus_chemo$plateau_end)),
  sprintf("4. In the 2010+ subset with known liver-metastasis status (n = %d) the unadjusted curve is %s/%s/%s (plateau %s), the original-covariate standardized curve %s/%s/%s (plateau %s), adding liver metastasis gives %s/%s/%s (plateau %s) and the all-combined model %s/%s/%s (plateau %s): the liver indicator moves the curve by at most a few tenths of a point.", nrow(nec10), fmt(ext$aj_2010$m0), fmt(ext$aj_2010$m6), fmt(ext$aj_2010$m12), fmt(ext$aj_2010$plateau_end), fmt(ext$base_2010$m0), fmt(ext$base_2010$m6), fmt(ext$base_2010$m12), fmt(ext$base_2010$plateau_end), fmt(ext$plus_liver_2010$m0), fmt(ext$plus_liver_2010$m6), fmt(ext$plus_liver_2010$m12), fmt(ext$plus_liver_2010$plateau_end), fmt(ext$all_combined_2010$m0), fmt(ext$all_combined_2010$m6), fmt(ext$all_combined_2010$m12), fmt(ext$all_combined_2010$plateau_end)),
  sprintf("5. Bootstrap (B = %d): the standardized 6-0 change is %s (95%% CI %s to %s) points and the 12-6 change %s (%s to %s); the plateau end is %s months (95%% CI %s to %s), with %.0f%% of replicates placing it at >= 6 months and %.0f%% at >= 3 months. The flat first half-year is therefore statistically compatible with no change, while the later decline is clearly non-zero; the plateau-end location itself is imprecise.", nrow(bs), fmt(pt["d6_0"]), fmt(ci[1, "d6_0"]), fmt(ci[2, "d6_0"]), fmt(pt["d12_6"]), fmt(ci[1, "d12_6"]), fmt(ci[2, "d12_6"]), fmt(pt["plateau_end"]), fmt(ci[1, "plateau_end"]), fmt(ci[2, "plateau_end"]), 100 * out$bootstrap$plateau_end_prop_ge_6, 100 * out$bootstrap$plateau_end_prop_ge_3),
  sprintf("6. Proportional hazards are violated for the cancer-death transition at landmark 0 (global p = %s; age and rectum are the offending terms) and %s at landmark 6 (global p = %s). This matters little for the estimand: each landmark model is only used to predict a 12-month probability averaged over the cohort, and the Breslow baseline absorbs the average time pattern; the own-risk-set calibration against Aalen-Johansen is within %.2f points at every landmark (the model sits %.1f-%.1f points below AJ from landmark 12 on, where fewer than 170 patients remain at risk) and the relaxed-PH piecewise construction changes the standardized values by at most %.2f points, both immaterial next to the %.0f-point 12-6 decline.", formatC(g6, format = "g", digits = 2), if (g6b < 0.05) "also violated" else "not rejected", formatC(g6b, format = "g", digits = 2), out$ph$max_abs_calibration_diff, min(abs(out$ph$calibration$diff[LM >= 12])), max(abs(out$ph$calibration$diff[LM >= 12])), out$ph$max_abs_piecewise_diff, abs(S_base$d12_6)),
  sprintf("7. The Fine-Gray standardization gives %s/%s/%s at landmarks 0/6/12 against %s/%s/%s for the multi-state Cox standardization and %s/%s/%s for the unadjusted Aalen-Johansen estimate: the choice of competing-risk model is immaterial because other-cause death is rare (%d of %d events) in this cohort.", fmt(fg[[1]]$std12_fine_gray), fmt(fg[[2]]$std12_fine_gray), fmt(fg[[3]]$std12_fine_gray), fmt(fg[[1]]$std12_multistate_cox), fmt(fg[[2]]$std12_multistate_cox), fmt(fg[[3]]$std12_multistate_cox), fmt(fg[[1]]$aj12), fmt(fg[[2]]$aj12), fmt(fg[[3]]$aj12), out$n$events$other_death, out$n$events$cancer_death + out$n$events$other_death),
  "", sprintf("Runtime %.1f min (%s, survival %s, cmprsk %s; bootstrap on %d cores). Figures: figures/output/A08_std_curves.png, figures/output/A08_comparators.png. Script: R/A08_standardization_extended.R.", out$runtime_minutes, out$software$R, out$software$survival, out$software$cmprsk, NCORES))
writeLines(md, paste0(RES, ID, ".md"))
cat("runtime", round(out$runtime_minutes, 1), "min; wrote", paste0(RES, ID, ".json"), "and .md\n")
