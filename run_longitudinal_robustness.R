project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
config_file <- file.path(project_dir, "config", "config.R")
if (!file.exists(config_file)) stop("Copy config/config.example.R to config/config.R and set data_root.")
source(config_file, encoding = "UTF-8")
scripts <- file.path(project_dir, "scripts", "Longitudinal_robustness", c(
  "01_freeze_and_audit.R",
  "02_longitudinal_robustness.R",
  "04_lagged_and_sample_reconciliation.R",
  "05_within_between_synthesis.R"
))
if (any(!file.exists(scripts))) stop("Missing script: ", paste(scripts[!file.exists(scripts)], collapse = ", "))
for (script in scripts) {
  message("Running ", basename(script))
  sys.source(script, envir = new.env(parent = globalenv()), keep.source = FALSE)
}
