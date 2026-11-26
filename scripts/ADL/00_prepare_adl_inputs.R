# Prepare the five-item self-care overlay and the frozen current-score cutpoints.
suppressPackageStartupMessages({library(data.table); library(haven)})
root <- Sys.getenv("DCV_BASE_DIR")
project <- Sys.getenv("DCV_PROJECT_DIR")
out <- Sys.getenv("ADL_OUTPUT_ROOT")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dictionary <- fread(file.path(project, "config", "selfcare_items.csv"))
files <- c(
  HRS = file.path(root, "HRS", "RAND HRS Longitudinal File 2022", "randhrs1992_2022v1_STATA", "randhrs1992_2022v1.dta"),
  ELSA = Sys.getenv("DEMENTIA_ELSA_OUTCOME_INPUT"),
  SHARE = Sys.getenv("DEMENTIA_SHARE_OUTCOME_INPUT"),
  CHARLS = file.path(root, "CHARLS", "02_harmonized", "Regular Waves", "H_CHARLS_D_Data", "H_CHARLS_D_Data.dta"),
  MHAS = file.path(root, "MHAS", "04_gateway_harmonized", "H_MHAS_d.dta")
)
id_fields <- list(HRS = "hhidpn", ELSA = "idauniqc", SHARE = "mergeid",
                  CHARLS = "ID", MHAS = c("rahhidnp", "unhhidnp"))
overlays <- lapply(names(files), function(cohort_name) {
  metadata <- dictionary[cohort == cohort_name]
  header <- names(read_dta(files[[cohort_name]], n_max = 0))
  id_field <- intersect(id_fields[[cohort_name]], header)[1L]
  if (is.na(id_field)) stop("Missing identifier in ", cohort_name)
  variables <- unique(metadata[status == "present", variable])
  if (length(setdiff(variables, header))) stop("Missing self-care items in ", cohort_name)
  raw <- read_dta(files[[cohort_name]], col_select = tidyselect::all_of(c(id_field, variables)))
  records <- lapply(sort(unique(metadata$wave)), function(w) {
    item_names <- metadata[wave == w, variable]
    if (length(item_names) != 5L) stop("Expected five self-care items per wave.")
    items <- vapply(item_names, function(v) {
      if (!v %in% names(raw)) return(rep(NA_real_, nrow(raw)))
      x <- as.numeric(raw[[v]])
      x[!x %in% c(0, 1)] <- NA_real_
      x
    }, numeric(nrow(raw)))
    observed <- rowSums(!is.na(items))
    any_difficulty <- as.integer(rowSums(items == 1, na.rm = TRUE) > 0)
    any_difficulty[observed == 0L] <- NA_integer_
    data.table(cohort = cohort_name, id = as.character(raw[[id_field]]), wave = as.integer(w),
               selfcare5_any = any_difficulty, selfcare5_n_observed = observed)
  })
  rbindlist(records)
})
overlay <- rbindlist(overlays)
if (anyDuplicated(overlay[, .(cohort, id, wave)])) stop("Duplicate self-care participant-wave records.")
saveRDS(overlay, Sys.getenv("ADL_SELFCARE_INPUT"))

# The cutpoint population matches the original clinical-history preparation.
candidate <- as.data.table(readRDS(Sys.getenv("ADL_CANDIDATE_INPUT")))
candidate <- candidate[cohort %in% names(files) & prior_n >= 1L & !is.na(sex_baseline) &
  is.finite(interval_start_age) & is.finite(interval_end_age) &
  is.finite(H_A) & is.finite(D_A) & is.finite(H_F) & is.finite(D_F) &
  is.finite(education_z) & is.finite(calendar_year_c) & interval_end_age > interval_start_age]
setorder(candidate, cohort, subject_id, interval_start_age, interval_end_age)
landmark <- candidate[interval_start_age >= 50, .SD[1L], by = .(cohort, subject_id)]
landmark <- landmark[interval_start_age <= 64]
landmark[, current_F := H_F + D_F]
cutpoints <- landmark[, .(current_p25 = as.numeric(quantile(current_F, 0.25)),
                          current_p75 = as.numeric(quantile(current_F, 0.75)),
                          history_p75 = as.numeric(quantile(H_F, 0.75))), by = cohort]
fwrite(cutpoints, Sys.getenv("ADL_CUTPOINTS_INPUT"))
