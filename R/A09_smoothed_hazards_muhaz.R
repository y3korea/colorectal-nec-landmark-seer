# A09: kernel-smoothed cause-specific hazard of cancer death (muhaz), distant-stage NEC / SCLC / AC / PANC.
# Independent re-implementation of the authors' T02_smooth_hazards.R procedure, plus common-bandwidth curves.
# Run: Rscript R/A09_smoothed_hazards_muhaz.R   (single-threaded, ~1-2 min); then python3 python/A09_write_md.py
# Revision 2026-10-07 (second audit): the 3-6 and 6-9 month interval means (and the 6-9/3-6 window ratio with its bootstrap CI) use the
#   trapezoid rule exactly as the authors' T02_draw_and_report.py (np.trapz(sub.hazard, sub.month)/(hi-lo), sub = month.between(lo, hi));
#   the previous unweighted mean of the same 13 grid points is kept as *_plain for comparison. Seed, resamples and all other quantities unchanged.
T0 <- proc.time()[3]
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "config.R"))  # paths, NEC_BOOT, NEC_STRICT
source(file.path(NEC_R_DIR, "pp.R"))
suppressMessages(library(muhaz))
FIG <- FIGS_DIR
ID <- "A09"; B <- nec_B(500); SEED <- 20260929L

## ---- cohorts (distant stage, same definitions as the authors) ----
coh <- list(
  NEC  = nec_main(loadc("NEN_all", "NEC"))[stage == "Distant"],
  SCLC = sclc_main(loadc("LUNG", "SCLC"))[stage == "Distant"],
  AC   = loadc("AC", "AC")[stage == "Distant"],
  PANC = panc_main(loadc("PANC", "PANC"))[stage == "Distant"])
LAB <- c(NEC = "Colorectal NEC", SCLC = "Small-cell lung cancer", AC = "Colorectal adenocarcinoma", PANC = "Pancreatic adenocarcinoma")
nec_check(nrow(coh$NEC) == 908, nrow(coh$SCLC) == 77707, nrow(coh$AC) == 147014, nrow(coh$PANC) == 75634)

## ---- muhaz wrapper: authors' settings (Epanechnikov, boundary correction both, 97-point grid on 0-24 months) ----
fit_bw <- function(t, e, bw, pilot) muhaz(t, e, min.time = 0, max.time = 24, bw.method = "global", bw.grid = bw, bw.pilot = pilot,
                                          b.cor = "both", n.min.grid = 51, n.est.grid = 97, kern = "epanechnikov")
auto_bw <- function(t, e) {
  pilot <- 24 / (8 * sum(e)^0.2); grid <- seq(0.2 * pilot, 20 * pilot, length.out = 25); rounds <- 0L
  fit <- fit_bw(t, e, grid, pilot)
  while (which.min(fit$glob.imse) == length(grid) && rounds < 3) { grid <- seq(min(grid), max(grid) * 2, length.out = 51); fit <- fit_bw(t, e, grid, pilot); rounds <- rounds + 1L }
  idx <- which.min(fit$glob.imse)
  if (idx %in% c(1, length(grid))) stop("bandwidth optimum at search edge")
  list(fit = fit, bw = fit$bw.glob, pilot = pilot, grid = grid, idx = idx, rounds = rounds, imse = fit$glob.imse)
}
# interval mean of h over [lo, hi] on the est.grid.
#   trapezoid (authors' T02_draw_and_report.py): np.trapz(h, month)/(hi - lo) over the grid points with lo <= month <= hi; on the 0.25-month grid
#   this is sum(w * h) with w = 0.25/(hi-lo) * (1/2, 1, ..., 1, 1/2), i.e. the two endpoints are half-weighted
#   plain: unweighted mean of the same grid points (definition used before the second-audit revision)
imean <- function(g, h, lo, hi, type = c("trapezoid", "plain")) {
  w <- which(g >= lo & g <= hi); gg <- g[w]; hh <- h[w]; k <- length(w)
  if (match.arg(type) == "plain") return(mean(hh))
  sum(diff(gg) * (hh[-1] + hh[-k]) / 2) / (hi - lo)
}
tz_weights <- function(gg, len) { d <- diff(gg); (c(d, 0) + c(0, d)) / 2 / len }   # trapezoid weights: sum(tz_weights(gg, hi-lo) * hh) == imean(..., 'trapezoid')
# features of a curve on the est.grid.
#  primary (authors' reading): trough = minimum of h over months (1, 12]; peak = maximum of h over (trough, 24]; interior flag = trough strictly below both neighbours
#  secondary (strict): first local minimum after month 1 (h[i] < h[i-1] & h[i] <= h[i+1]) and the first local maximum after it
features <- function(g, h) {
  im <- sapply(1:24, function(m) h[which.min(abs(g - m))]); n <- length(h)
  w <- which(g > 1 & g <= 12); imin <- w[which.min(h[w])]; w2 <- which(g > g[imin]); imax <- w2[which.max(h[w2])]
  interior <- imin > 1 && imin < n && h[imin] < h[imin - 1] && h[imin] < h[imin + 1]
  lmin <- NA; lmax <- NA
  for (i in 2:(n - 1)) if (g[i] > 1 && h[i] < h[i - 1] && h[i] <= h[i + 1]) { lmin <- i; break }
  if (!is.na(lmin)) { for (i in (lmin + 1):(n - 1)) if (h[i] > h[i - 1] && h[i] >= h[i + 1]) { lmax <- i; break }
                      if (is.na(lmax) && h[n] > h[n - 1]) lmax <- n }
  list(h_month = as.list(setNames(im, paste0("m", 1:24))),
       min_month = g[imin], min_haz = h[imin], max_month = g[imax], max_haz = h[imax], ratio_max_min = h[imax] / h[imin], trough_interior = interior,
       local_min_month = if (is.na(lmin)) NA else g[lmin], local_min_haz = if (is.na(lmin)) NA else h[lmin],
       local_max_month = if (is.na(lmax)) NA else g[lmax], local_max_haz = if (is.na(lmax)) NA else h[lmax],
       local_ratio = if (is.na(lmin) || is.na(lmax)) NA else h[lmax] / h[lmin],
       mean_3_6 = imean(g, h, 3, 6), mean_6_9 = imean(g, h, 6, 9),                                      # trapezoid rule (authors' definition)
       mean_3_6_plain = imean(g, h, 3, 6, "plain"), mean_6_9_plain = imean(g, h, 6, 9, "plain"))        # unweighted mean of the same 13 grid points
}

## ---- (1) automatic bandwidths + (2) common bandwidths ----
res <- list(); grid_rows <- list(); auto <- list()
for (g in names(coh)) {
  d <- coh[[g]]; t <- pmin(d$tm, 24); e <- as.integer(d$ev == 1 & d$tm <= 24)
  a <- auto_bw(t, e); auto[[g]] <- a
  bws <- c(auto = a$bw, bw1.0 = 1.0, bw1.5 = 1.5); if (g != "NEC") bws <- c(bws, bw0.5 = 0.5)
  if (g == "NEC") bws <- c(bws, half = a$bw / 2, double = a$bw * 2)
  res[[g]] <- list(n = nrow(d), cancer_deaths_24m = sum(e), other_deaths_24m = sum(d$ev == 2 & d$tm <= 24), pilot = a$pilot, bw_auto = a$bw,
                   bw_grid = a$grid, bw_grid_index = a$idx, grid_expansions = a$rounds, glob_imse = a$imse, curves = list())
  for (nm in names(bws)) {
    f <- fit_bw(t, e, bws[[nm]], a$pilot); fe <- features(f$est.grid, f$haz.est)
    res[[g]]$curves[[nm]] <- c(list(bandwidth = bws[[nm]]), fe)
    grid_rows[[paste(g, nm)]] <- data.table(group = g, bw_label = nm, bandwidth = bws[[nm]], month = f$est.grid, hazard = f$haz.est,
                                            n_at_risk = sapply(f$est.grid, function(s) sum(d$tm > s)))  # strictly surviving past s (landmark convention)
  }
  cat(sprintf("%-5s n=%6d deaths24=%6d pilot=%.3f bw_auto=%.4f (grid idx %d/%d, expansions %d)\n", g, nrow(d), sum(e), a$pilot, a$bw, a$idx, length(a$grid), a$rounds))
}

## ---- NEC bootstrap (patient-level, same resample for every bandwidth), pointwise percentile 95% CI ----
d <- coh$NEC; t <- pmin(d$tm, 24); e <- as.integer(d$ev == 1 & d$tm <= 24); n <- length(t); pilot <- auto$NEC$pilot
boot_bws <- c(auto = auto$NEC$bw, bw1.0 = 1.0, bw1.5 = 1.5, half = auto$NEC$bw / 2, double = auto$NEC$bw * 2)
set.seed(SEED); boot <- lapply(boot_bws, function(x) matrix(NA_real_, B, 97))
for (b in seq_len(B)) { idx <- sample.int(n, n, replace = TRUE)
  for (nm in names(boot_bws)) boot[[nm]][b, ] <- fit_bw(t[idx], e[idx], boot_bws[[nm]], pilot)$haz.est }
for (nm in names(boot_bws)) {
  ci <- apply(boot[[nm]], 2, quantile, c(.025, .975)); r <- grid_rows[[paste("NEC", nm)]]
  r[, `:=`(lower95 = ci[1, ], upper95 = ci[2, ])]; grid_rows[[paste("NEC", nm)]] <- r
  g <- r$month; full <- res$NEC$curves[[nm]]
  # (a) re-optimised statistics (trough/peak re-located in every resample): timing + sup/inf ratio, which is selection-biased upward
  fb <- apply(boot[[nm]], 1, function(h) { f <- features(g, h); c(f$min_month, f$max_month, f$ratio_max_min, f$trough_interior) })
  # (b) fixed-location statistics (inferential): h at the full-sample peak month / h at the full-sample trough month, and mean(6-9)/mean(3-6)
  imin <- which.min(abs(g - full$min_month)); imax <- which.min(abs(g - full$max_month)); w36 <- g >= 3 & g <= 6; w69 <- g >= 6 & g <= 9
  r_fixed <- boot[[nm]][, imax] / boot[[nm]][, imin]
  # window ratio 6-9 / 3-6: trapezoid-rule interval means (authors' definition); the weight vectors reproduce features() on the full-sample curve
  tw36 <- tz_weights(g[w36], 6 - 3); tw69 <- tz_weights(g[w69], 9 - 6)
  stopifnot(abs(sum(tw36 * r$hazard[w36]) - full$mean_3_6) < 1e-12, abs(sum(tw69 * r$hazard[w69]) - full$mean_6_9) < 1e-12)
  r_window <- as.vector(boot[[nm]][, w69] %*% tw69) / as.vector(boot[[nm]][, w36] %*% tw36)
  r_window_plain <- rowMeans(boot[[nm]][, w69]) / rowMeans(boot[[nm]][, w36])   # pre-revision definition (unweighted grid-point means), kept for comparison
  q <- function(x, type = 7) unname(quantile(x, c(.025, .975), na.rm = TRUE, type = type))
  res$NEC$curves[[nm]]$ratio_window <- full$mean_6_9 / full$mean_3_6                      # trapezoid (authors' definition)
  res$NEC$curves[[nm]]$ratio_window_plain <- full$mean_6_9_plain / full$mean_3_6_plain    # unweighted grid-point means
  res$NEC$curves[[nm]]$boot <- list(B = B, seed = SEED,
    min_month_ci = q(fb[1, ], type = 1), max_month_ci = q(fb[2, ], type = 1),      # type-1 quantiles: attainable 0.25-month grid values
    prop_trough_interior = mean(fb[4, ] == 1), prop_trough_after_6 = mean(fb[1, ] > 6),
    prop_peak_at_24 = mean(fb[2, ] >= 24), prop_peak_ge_20 = mean(fb[2, ] >= 20),
    ratio_reoptimised_ci = q(fb[3, ]), ratio_reoptimised_boot_mean = mean(fb[3, ]), ratio_reoptimised_boot_median = median(fb[3, ]),
    ratio_fixed_ci = q(r_fixed), ratio_fixed_boot_mean = mean(r_fixed),
    ratio_window_ci = q(r_window), ratio_window_boot_mean = mean(r_window),
    ratio_window_plain_ci = q(r_window_plain), ratio_window_plain_boot_mean = mean(r_window_plain),
    fixed_months = list(trough = full$min_month, peak = full$max_month))
}
GRID <- rbindlist(grid_rows, fill = TRUE); fwrite(GRID, paste0(RES, ID, "_hazard_grid.csv"))

## ---- reproduction check (manuscript / task anchors) ----
nec_auto <- res$NEC$curves$auto
chk <- data.table(item = c("bw NEC", "bw SCLC", "bw AC", "bw PANC", "NEC trough month", "NEC trough hazard", "NEC peak month", "NEC peak hazard",
                           "NEC mean 3-6 (trapezoid)", "NEC mean 6-9 (trapezoid)", "SCLC trough month", "SCLC trough hazard", "SCLC peak month", "SCLC peak hazard",
                           "AC mean 3-6 (trapezoid)", "AC mean 6-9 (trapezoid)", "PANC mean 3-6 (trapezoid)", "PANC mean 6-9 (trapezoid)"),
                  manuscript = c(1.48, 0.34, 0.32, 0.34, 4, 0.084, 8.75, 0.147, 0.095, 0.137, 4.25, 0.069, 10, 0.122, 0.042, 0.035, 0.134, 0.115),
                  ours = c(res$NEC$bw_auto, res$SCLC$bw_auto, res$AC$bw_auto, res$PANC$bw_auto, nec_auto$min_month, nec_auto$min_haz, nec_auto$max_month, nec_auto$max_haz,
                           nec_auto$mean_3_6, nec_auto$mean_6_9, res$SCLC$curves$auto$min_month, res$SCLC$curves$auto$min_haz, res$SCLC$curves$auto$max_month, res$SCLC$curves$auto$max_haz,
                           res$AC$curves$auto$mean_3_6, res$AC$curves$auto$mean_6_9, res$PANC$curves$auto$mean_3_6, res$PANC$curves$auto$mean_6_9),
                  tol = c(.01, .01, .01, .01, .3, .002, .3, .002, .002, .002, .5, .002, .5, .002, .002, .002, .002, .002))
chk[, ok := !is.na(ours) & abs(ours - manuscript) <= tol]
print(chk)
if (!all(chk$ok)) nec_anchor_fail("REPRODUCTION CHECK FAILED")
cat("REPRODUCTION CHECK PASSED (", nrow(chk), "items )\n")
cat(sprintf("NEC auto interval means: 3-6 trapezoid %.4f / plain %.4f; 6-9 trapezoid %.4f / plain %.4f; window ratio trapezoid %.3f / plain %.3f\n",
            nec_auto$mean_3_6, nec_auto$mean_3_6_plain, nec_auto$mean_6_9, nec_auto$mean_6_9_plain, nec_auto$ratio_window, nec_auto$ratio_window_plain))

## ---- figures ----
stagger <- function(y, gap) { o <- order(y); z <- y[o]; for (i in seq_along(z)[-1]) if (z[i] - z[i - 1] < gap) z[i] <- z[i - 1] + gap; y[o] <- z; y }
COL <- c(NEC = "#2a78d6", SCLC = "#eb6834", AC = "#1baf7a", PANC = "#4a3aa7")
png(paste0(FIG, ID, "_panel_bw1.5.png"), width = 7.2, height = 4.6, units = "in", res = 300)
par(mar = c(5.2, 4.2, 2.5, 7.5), xpd = NA, family = "sans", mgp = c(2.3, 0.7, 0))
ymax <- 0.25; plot(NA, xlim = c(0, 24), ylim = c(0, ymax), xlab = "Months since diagnosis", ylab = "Cancer-death hazard (per person-month)", axes = FALSE,
                   main = "Smoothed cause-specific hazard, distant stage, common bandwidth 1.5 months", cex.main = 0.9, font.main = 1)
axis(1, at = seq(0, 24, 3), col = "grey60", col.axis = "grey30"); axis(2, las = 1, col = "grey60", col.axis = "grey30")
abline(h = seq(0.05, ymax, 0.05), col = "grey92", xpd = FALSE)
r <- GRID[group == "NEC" & bw_label == "bw1.5"]; polygon(c(r$month, rev(r$month)), c(r$lower95, rev(r$upper95)), col = adjustcolor(COL["NEC"], 0.15), border = NA)
ends <- sapply(names(coh), function(g) GRID[group == g & bw_label == "bw1.5"]$hazard[97]); ypos <- stagger(ends, 0.014)
for (g in names(coh)) { r <- GRID[group == g & bw_label == "bw1.5"]; lines(r$month, r$hazard, col = COL[g], lwd = ifelse(g == "NEC", 2.6, 2))
  text(24.3, ypos[g], LAB[g], col = COL[g], adj = 0, cex = 0.72) }
legend("topright", inset = c(-0.01, -0.02), legend = c(LAB, "NEC 95% bootstrap CI (B=500)"), col = c(COL, NA), lwd = c(2.6, 2, 2, 2, NA),
       fill = c(NA, NA, NA, NA, adjustcolor(COL["NEC"], 0.15)), border = NA, bty = "n", cex = 0.65, x.intersp = 0.6, seg.len = 1.4)
mtext("Months 0-1 lie in the boundary region and are not interpreted. Other-cause deaths censored; Epanechnikov kernel, boundary correction.", side = 1, line = 4.1, cex = 0.6, col = "grey35", adj = 0)
dev.off()

png(paste0(FIG, ID, "_NEC_bandwidths.png"), width = 7.2, height = 4.6, units = "in", res = 300)
par(mar = c(5.2, 4.2, 2.5, 7.5), xpd = NA, family = "sans", mgp = c(2.3, 0.7, 0))
plot(NA, xlim = c(0, 24), ylim = c(0, 0.3), xlab = "Months since diagnosis", ylab = "Cancer-death hazard (per person-month)", axes = FALSE,
     main = "Distant colorectal NEC (n = 908): smoothed hazard at three bandwidths", cex.main = 0.9, font.main = 1)
axis(1, at = seq(0, 24, 3), col = "grey60", col.axis = "grey30"); axis(2, las = 1, col = "grey60", col.axis = "grey30"); abline(h = seq(0.05, 0.3, 0.05), col = "grey92", xpd = FALSE)
r <- GRID[group == "NEC" & bw_label == "auto"]; polygon(c(r$month, rev(r$month)), c(r$lower95, rev(r$upper95)), col = adjustcolor(COL["NEC"], 0.15), border = NA)
bwset <- c(half = "#eb6834", auto = "#2a78d6", double = "#4a3aa7"); bwlab <- c(half = sprintf("half (%.2f mo)", auto$NEC$bw / 2), auto = sprintf("automatic (%.2f mo)", auto$NEC$bw), double = sprintf("double (%.2f mo)", auto$NEC$bw * 2))
ends <- sapply(names(bwset), function(nm) GRID[group == "NEC" & bw_label == nm]$hazard[97]); ypos <- stagger(ends, 0.017)
for (nm in names(bwset)) { r <- GRID[group == "NEC" & bw_label == nm]; lines(r$month, r$hazard, col = bwset[nm], lwd = ifelse(nm == "auto", 2.6, 1.8))
  text(24.3, ypos[nm], bwlab[nm], col = bwset[nm], adj = 0, cex = 0.7) }
legend("topright", inset = c(-0.01, -0.02), legend = c(bwlab, "95% bootstrap CI (automatic bw, B=500)"), col = c(bwset, NA), lwd = c(1.8, 2.6, 1.8, NA),
       fill = c(NA, NA, NA, adjustcolor(COL["NEC"], 0.15)), border = NA, bty = "n", cex = 0.65, x.intersp = 0.6, seg.len = 1.4)
mtext("Bandwidth chosen by minimum integrated MSE over 25 candidates (pilot 24/(8*deaths^0.2)). Month 0 is a boundary estimate.", side = 1, line = 4.1, cex = 0.6, col = "grey35", adj = 0)
dev.off()

## ---- JSON ----
runtime <- as.numeric(proc.time()[3] - T0) / 60
out <- list(id = ID, date = "2026-10-07", method = "muhaz 1.2.6.5, Epanechnikov, b.cor='both', bw.method='global', n.min.grid=51, n.est.grid=97, 0-24 months; t=min(tm,24), event = cancer death <= 24 months; other-cause deaths censored",
            bootstrap = list(cohort = "NEC", B = B, seed = SEED, resampling = "patient-level, same resample for all bandwidths", ci = "pointwise percentile 95%",
                             ratio_note = "ratio_reoptimised_ci re-locates trough/peak in each resample (sup/inf statistic, selection-biased upward; descriptive only); ratio_fixed_ci holds the full-sample trough/peak months fixed; ratio_window_ci is the trapezoid-rule interval mean h(6-9) / h(3-6) (authors' definition) and ratio_window_plain_ci the unweighted-grid-mean version; month CIs use type-1 quantiles (grid-attainable values)"),
            interval_means = list(definition = "mean_3_6 / mean_6_9 (and ratio_window*): trapezoid rule over the 0.25-month grid points in [lo, hi] divided by (hi - lo), exactly as the authors' T02_draw_and_report.py (np.trapz(sub.hazard, sub.month)/(hi-lo) with sub = month.between(lo, hi)); the two endpoints are half-weighted",
                                  plain = "mean_*_plain / ratio_window_plain*: unweighted mean of the same 13 grid points (definition used before the 2026-10-07 second-audit revision)"),
            revision = "2026-10-07 second audit: interval means and the 6-9/3-6 window ratio (point estimate and bootstrap CI) switched to the trapezoid rule; same seed and resamples; bandwidths, curves, extrema, fixed-month and re-optimised intervals unchanged",
            reproduction_check = chk, groups = res, muhaz_version = as.character(packageVersion("muhaz")), R_version = R.version.string, runtime_minutes = runtime)
write_json(out, paste0(RES, ID, ".json"), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
cat("runtime minutes", round(runtime, 2), "\n")
