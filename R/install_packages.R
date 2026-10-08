# Install the R packages used by the pipeline (versions of the verified reproduction run in comments).
# Run: Rscript R/install_packages.R
pkgs <- c(data.table = "1.18.4", survival = "3.8.6", sandwich = "3.1.1", jsonlite = "2.0.0", cmprsk = "2.2.12", muhaz = "1.2.6.5")
repo <- "https://cloud.r-project.org"
missing <- names(pkgs)[!vapply(names(pkgs), requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = repo)
# exact versions (optional): set NEC_PIN_VERSIONS=1 to install the verified versions with 'remotes'
if (Sys.getenv("NEC_PIN_VERSIONS") == "1") {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes", repos = repo)
  for (p in names(pkgs)) if (as.character(packageVersion(p)) != pkgs[[p]]) remotes::install_version(p, version = pkgs[[p]], repos = repo, upgrade = "never")
}
for (p in names(pkgs)) cat(sprintf("%-11s installed %-9s (verified run: %s)\n", p, as.character(packageVersion(p)), pkgs[[p]]))
cat(R.version.string, "(verified run: R version 4.5.3)\n")
