# Run from the repository root.
stages <- c("track_A", "longitudinal_robustness", "mortality_bridge",
            "mortality_prediction", "clinical_translation", "adl",
            "dementia_bridge")
for (stage in stages) {
  message("Stage: ", stage)
  sys.source(paste0("run_", stage, ".R"), envir = new.env(parent = globalenv()))
}
