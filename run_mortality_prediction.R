project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
config_file <- file.path(project_dir, "config", "config.R")
if (!file.exists(config_file)) stop("Copy config/config.example.R to config/config.R and set data_root.")
source(config_file, encoding = "UTF-8")
scripts <- file.path(project_dir, "scripts", "Mortality_prediction", c(
  "01_build_age50_landmark_dataset.R",
  "02_fit_and_validate_five_year_models.R"
))
if (any(!file.exists(scripts))) stop("Missing script: ", paste(scripts[!file.exists(scripts)], collapse = ", "))
for (script in scripts) {
  message("Running ", basename(script))
  sys.source(script, envir = new.env(parent = globalenv()), keep.source = FALSE)
}
