# Shared R helpers (independent re-implementation of the authors' piecewise-rate / interaction / standardization script).
# Paths come from config.R (sourced by the calling script before this file).
suppressMessages({library(data.table); library(survival); library(sandwich); library(jsonlite)})
if (!exists("NEC_R_DIR")) stop("source R/config.R before R/pp.R")
NEC_CODES <- c("8013","8041","8042","8043","8044","8045"); SCLC_CODES <- c("8041","8042","8043","8044","8045"); ISLET <- as.character(8150:8157)
loadc <- function(name, tag) {
  A <- fread(paste0(DATA, "cohort_", name, ".csv"), colClasses = list(character = c("pid","hist","s1","s2","seq","grade17","grade18","tt_days","smonths","ntum")))
  A[, pid := paste0(tag, "_", pid)]
  A[, age10 := (age - 65) / 10][, age10sq := age10^2][, year10 := (year - 2012) / 10]
  A[, large := as.integer(hist == "8013")]
  # primary-site surgery, same rule as the authors' table script: no = s1 "00" or s2 "A000"; unknown = "99"/"A990"; not applicable = "98"/"A980"; blank both = unknown; else yes
  A[, surg_no := as.integer(s1 == "00" | s2 == "A000")]
  A[, surg_unk := as.integer(s1 %in% c("99","98") | s2 %in% c("A990","A980") | (s1 %in% c("Blank(s)","") & s2 %in% c("Blank(s)","")))]
  A[, surg := as.integer(surg_no == 0L & surg_unk == 0L)]
  A[, first_primary := as.integer(seq %in% c("One primary only", "1st of 2 or more primaries"))]
  A
}
nec_main <- function(nen) nen[hist %in% NEC_CODES]
nec_expanded <- function(nen) nen[hist %in% NEC_CODES | (hist == "8246" & (grade17 %in% c("Poorly differentiated; Grade III","Undifferentiated; anaplastic; Grade IV") | grade18 %in% c("Site-specific grade system category (3)","Site-specific grade system category (4)")))]
sclc_main <- function(lung) lung[hist %in% SCLC_CODES]
panc_main <- function(panc) panc[!hist %in% ISLET]
COV <- list(NEC = "age10 + age10sq + male + year10 + large + rectum", AC = "age10 + age10sq + male + year10 + rectum", SCLC = "age10 + age10sq + male + year10", PANC = "age10 + age10sq + male + year10")
# piecewise-exponential data: follow-up from `start`, cut at `cuts`, censored at `end`; d = cancer death in interval
pw <- function(dt, start, cuts, end = 24) {
  x <- copy(dt)[tm > start]
  x[, `:=`(t0 = start, t1 = pmin(tm, end), d = as.integer(ev == 1 & tm <= end))]
  s <- as.data.table(survSplit(Surv(t0, t1, d) ~ ., data = as.data.frame(x), cut = cuts, episode = "per"))
  s[, expo := t1 - t0]; s[expo > 0]
}
wald <- function(fit, V, idx) { b <- coef(fit)[idx]; W <- as.numeric(t(b) %*% solve(V[idx, idx]) %*% b); c(W = W, df = length(idx), p = pchisq(W, length(idx), lower.tail = FALSE)) }
ctr <- function(fit, V, v) { e <- sum(v * coef(fit)); se <- sqrt(as.numeric(t(v) %*% V %*% v)); c(RR = exp(e), lo = exp(e - 1.96 * se), hi = exp(e + 1.96 * se), p = 2 * pnorm(-abs(e / se)), logRR = e, se = se) }
vec <- function(fit, nm) { v <- setNames(rep(0, length(coef(fit))), names(coef(fit))); for (k in names(nm)) v[k] <- nm[[k]]; v }
# standardized landmark curve: per-landmark multi-state Cox (12-month window), averaged over reference population
std_curve <- function(dt, covf, ref, LM = 0:18) {
  sapply(LM, function(s) {
    x <- dt[tm > s]; t <- x$tm - s; e <- x$ev; e[t > 12] <- 0L; t <- pmin(t, 12)
    x[, `:=`(tt = t, st = factor(e, levels = 0:2, labels = c("cens", "cancer", "other")))]
    fit <- coxph(as.formula(paste("Surv(tt, st) ~", covf)), data = x, id = pid)
    sf <- survfit(fit, newdata = ref, se.fit = FALSE)
    k <- which(sf$states == "cancer"); idx <- max(which(sf$time <= 12))
    mean(sf$pstate[idx, , k]) * 100 })
}
hinge <- function(y, m = 0:18) { best <- c(Inf, NA); for (c0 in seq(0, 15, 0.25)) { x <- pmax(0, m - c0); sse <- sum(lm.fit(cbind(1, x), y)$residuals^2); if (sse < best[1]) best <- c(sse, c0) }; best[2] }
