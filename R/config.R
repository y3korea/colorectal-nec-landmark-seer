# Shared path / run configuration for every R step of the pipeline (mirror of python/config.py).
# Environment variables (relative values resolved against the working directory; defaults relative to the repo root):
#   NEC_DATA_DIR (./data/raw)  NEC_WORK_DIR (./work)  NEC_RESULTS_DIR (./results)  NEC_FIG_DIR (./figures/output)
#   NEC_BOOT (2000; secondary bootstraps scaled proportionally)  NEC_QUICK (0/1)  NEC_STRICT (1 unless quick/synthetic)
#   NEC_FILE_NEN / NEC_FILE_AC / NEC_FILE_LUNG / NEC_FILE_PANC (export file names inside NEC_DATA_DIR)
# Every R script locates this file with: source(file.path(dirname(<script path>), "config.R")).
local({
  a <- commandArgs(trailingOnly = FALSE); f <- sub("^--file=", "", a[grep("^--file=", a)])
  r_dir <- if (length(f)) dirname(normalizePath(f[1], mustWork = FALSE)) else file.path(getwd(), "R")
  if (nzchar(Sys.getenv("NEC_REPO"))) r_dir <- file.path(Sys.getenv("NEC_REPO"), "R")
  assign("NEC_R_DIR", normalizePath(r_dir, mustWork = FALSE), envir = globalenv())
  assign("NEC_REPO_DIR", normalizePath(file.path(r_dir, ".."), mustWork = FALSE), envir = globalenv())
})
.nec_truthy <- function(v) tolower(trimws(v)) %in% c("1", "true", "yes", "y", "on")
.nec_dir <- function(env, default) {
  v <- trimws(Sys.getenv(env, ""))
  p <- if (nzchar(v)) path.expand(v) else file.path(NEC_REPO_DIR, default)
  normalizePath(p, mustWork = FALSE)
}
NEC_DATA_DIR <- .nec_dir("NEC_DATA_DIR", "data/raw")
NEC_WORK_DIR <- .nec_dir("NEC_WORK_DIR", "work")
NEC_RESULTS_DIR <- .nec_dir("NEC_RESULTS_DIR", "results")
NEC_FIG_DIR <- .nec_dir("NEC_FIG_DIR", "figures/output")
for (.d in c(NEC_WORK_DIR, NEC_RESULTS_DIR, NEC_FIG_DIR)) dir.create(.d, showWarnings = FALSE, recursive = TRUE)
NEC_RAW_FILES <- c(NEN = Sys.getenv("NEC_FILE_NEN", "colorectal_NEN_SEER17_2000_2023_full_20260806.txt"),
                   AC = Sys.getenv("NEC_FILE_AC", "colorectal_adenocarcinoma_SEER17_2000_2023_full_20260928.txt"),
                   LUNG = Sys.getenv("NEC_FILE_LUNG", "SCLC_lung_SEER17_2000_2023_full_20260928.txt"),
                   PANC = Sys.getenv("NEC_FILE_PANC", "pancreatic_adenocarcinoma_SEER17_2000_2023_full_20260928.txt"))
NEC_SYNTHETIC <- file.exists(file.path(NEC_DATA_DIR, "SYNTHETIC_DATA_README.txt"))
NEC_QUICK <- .nec_truthy(Sys.getenv("NEC_QUICK", "0"))
NEC_PAPER_B <- 2000L
NEC_BOOT <- if (nzchar(trimws(Sys.getenv("NEC_BOOT")))) as.integer(Sys.getenv("NEC_BOOT")) else if (NEC_QUICK) 40L else NEC_PAPER_B
NEC_STRICT <- if (nzchar(trimws(Sys.getenv("NEC_STRICT")))) .nec_truthy(Sys.getenv("NEC_STRICT")) else !(NEC_QUICK || NEC_SYNTHETIC)
# replicates for an analysis that used `paper` replicates in the paper (identity at NEC_BOOT = 2000)
nec_B <- function(paper) if (NEC_BOOT == NEC_PAPER_B) as.integer(paper) else max(4L, as.integer(round(paper * NEC_BOOT / NEC_PAPER_B)))
# reproduction anchor not met: stop in strict mode, warn otherwise
nec_anchor_fail <- function(...) { msg <- paste0(...); if (NEC_STRICT) stop(msg, call. = FALSE) else message("[non-strict mode] ANCHOR NOT MET (continuing): ", msg) }
# stopifnot() for reproduction anchors (counts of the authors' cohorts): enforced only in strict mode
nec_check <- function(...) { ok <- tryCatch(all(c(...)), error = function(e) FALSE); if (!isTRUE(ok)) nec_anchor_fail("reproduction count check failed: ", deparse(substitute(list(...)))) ; invisible(ok) }
# trailing-slash variants used by the analysis scripts (paste0 concatenation, as in the original code)
DATA <- paste0(NEC_WORK_DIR, "/"); RES <- paste0(NEC_RESULTS_DIR, "/"); FIGS_DIR <- paste0(NEC_FIG_DIR, "/"); RAW_DIR <- paste0(NEC_DATA_DIR, "/")
