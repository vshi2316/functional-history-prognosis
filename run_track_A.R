project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
config_file <- file.path(project_dir, "config", "config.R")
if (!file.exists(config_file)) stop("Copy config/config.example.R to config/config.R and set data_root.")
source(config_file, encoding = "UTF-8")
scripts <- file.path(project_dir, "scripts", "Track_A", c(
  "01_prepare_cohort_inputs.R",
  "02_derive_and_harmonize_domains.R",
  "03_measurement_invariance.R",
  "04_longitudinal_models.R",
  "05_meta_analysis_and_robustness.R"
))
if (any(!file.exists(scripts))) stop("Missing script: ", paste(scripts[!file.exists(scripts)], collapse = ", "))
for (script in scripts) {
  message("Running ", basename(script))
  sys.source(script, envir = new.env(parent = globalenv()), keep.source = FALSE)
}
