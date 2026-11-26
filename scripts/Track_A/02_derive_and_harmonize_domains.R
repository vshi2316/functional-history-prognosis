# derive_domain_variables_v1.R
local({
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

resolve_base_dir <- function() {
  path <- Sys.getenv("DCV_BASE_DIR", unset = "")
  if (!nzchar(path) || !dir.exists(path)) stop("Set data_root in config/config.R.")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

base_dir <- resolve_base_dir()
input_core4_dir <- file.path(Sys.getenv("DCV_ANALYSIS_DIR"), "analysis_ready_core", "ready_v3_core4")
input_hrs_dir <- file.path(Sys.getenv("DCV_ANALYSIS_DIR"), "analysis_ready_core", "ready_v3_hrs")
out_dir <- file.path(Sys.getenv("DCV_ANALYSIS_DIR"), "analysis_ready_core", "domain_ready_v1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

required_files <- c(
  file.path(input_core4_dir, "share_analysis_ready_long_v3_core4.rds"),
  file.path(input_core4_dir, "elsa_analysis_ready_long_v3_core4.rds"),
  file.path(input_core4_dir, "mhas_analysis_ready_long_v3_core4.rds"),
  file.path(input_core4_dir, "charls_analysis_ready_long_v3_core4.rds"),
  file.path(input_hrs_dir, "hrs_analysis_ready_long_v3.rds")
)

missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files) > 0) {
  stop("Missing required input files:\n", paste(missing_files, collapse = "\n"))
}

z_score <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

z_by_cohort <- function(x, cohort) {
  ave(x, cohort, FUN = z_score)
}

mean_na <- function(...) {
  mat <- cbind(...)
  out <- rowMeans(mat, na.rm = TRUE)
  all_missing <- apply(is.na(mat), 1, all)
  out[all_missing] <- NA_real_
  out
}

count_nonmissing <- function(...) {
  mat <- cbind(...)
  rowSums(!is.na(mat))
}

read_inputs <- function() {
  files <- list(
    SHARE = file.path(input_core4_dir, "share_analysis_ready_long_v3_core4.rds"),
    ELSA = file.path(input_core4_dir, "elsa_analysis_ready_long_v3_core4.rds"),
    MHAS = file.path(input_core4_dir, "mhas_analysis_ready_long_v3_core4.rds"),
    CHARLS = file.path(input_core4_dir, "charls_analysis_ready_long_v3_core4.rds"),
    HRS = file.path(input_hrs_dir, "hrs_analysis_ready_long_v3.rds")
  )
  out <- purrr::imap(files, ~ readRDS(.x) %>% mutate(cohort = .y))
  dplyr::bind_rows(out)
}

derive_domains <- function(df) {
  df %>%
    mutate(
      depression_z = z_by_cohort(depression_score, cohort),
      loneliness_z = z_by_cohort(loneliness_score, cohort),
      cognition_z = z_by_cohort(cognition_score, cohort),
      adl_z = z_by_cohort(adl_count, cohort),
      iadl_z = z_by_cohort(iadl_count, cohort),
      srh_z = z_by_cohort(self_rated_health, cohort),

      affective_score_v1 = mean_na(depression_z, loneliness_z),
      affective_component_n = count_nonmissing(depression_score, loneliness_score),

      functional_score_v1 = mean_na(adl_z, iadl_z, srh_z),
      functional_component_n = count_nonmissing(adl_count, iadl_count, self_rated_health),
      functional_any_limitation_v1 = ifelse(
        is.na(adl_any) & is.na(iadl_any),
        NA_real_,
        pmax(adl_any, iadl_any, na.rm = TRUE)
      ),

      basic_cognition_score_v1 = cognition_z,
      basic_cognition_component_n = ifelse(is.na(cognition_score), 0, 1),



      affective_ready_v1 = ifelse(!is.na(depression_score), 1, 0),
      functional_ready_v1 = ifelse(functional_component_n >= 2, 1, 0),
      cognition_ready_v1 = ifelse(!is.na(cognition_score), 1, 0)
    )
}

build_dictionary <- function(df) {
  tibble(
    variable = names(df),
    class = purrr::map_chr(df, ~ paste(class(.x), collapse = ";")),
    non_missing_n = purrr::map_int(df, ~ sum(!is.na(.x))),
    example = purrr::map_chr(df, ~ {
      vals <- unique(.x[!is.na(.x)])
      if (length(vals) == 0) return("")
      paste(head(as.character(vals), 3), collapse = " | ")
    })
  )
}

save_split_outputs <- function(df) {
  split(df, df$cohort) |>
    purrr::iwalk(function(x, nm) {
      saveRDS(x, file.path(out_dir, paste0(tolower(nm), "_domain_ready_v1.rds")))
    })
}

analysis_ready_v3 <- read_inputs()
domain_ready_v1 <- derive_domains(analysis_ready_v3)

save_split_outputs(domain_ready_v1)
saveRDS(domain_ready_v1, file.path(out_dir, "pooled_domain_ready_v1.rds"))
write_csv(build_dictionary(domain_ready_v1), file.path(out_dir, "domain_variable_dictionary_v1.csv"), na = "")

message("Done.")
message("Outputs written to: ", out_dir)
})

# build_domain_coverage_report_v1.R
local({
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

resolve_base_dir <- function() {
  path <- Sys.getenv("DCV_BASE_DIR", unset = "")
  if (!nzchar(path) || !dir.exists(path)) stop("Set data_root in config/config.R.")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

base_dir <- resolve_base_dir()
input_dir <- file.path(Sys.getenv("DCV_ANALYSIS_DIR"), "analysis_ready_core", "domain_ready_v1")
out_dir <- file.path(Sys.getenv("DCV_ANALYSIS_DIR"), "analysis_ready_core", "domain_ready_v1")

pooled_file <- file.path(input_dir, "pooled_domain_ready_v1.rds")
if (!file.exists(pooled_file)) {
  stop("Missing required input file: ", pooled_file)
}

domain_ready <- readRDS(pooled_file)

domain_specs <- tribble(
  ~domain_group, ~ready_var, ~score_var, ~component_var,
  "affective", "affective_ready_v1", "affective_score_v1", "affective_component_n",
  "functional", "functional_ready_v1", "functional_score_v1", "functional_component_n",
  "basic_cognition", "cognition_ready_v1", "basic_cognition_score_v1", "basic_cognition_component_n"
)

domain_coverage_report <- purrr::pmap_dfr(
  domain_specs,
  function(domain_group, ready_var, score_var, component_var) {
    domain_group_value <- domain_group

    domain_ready %>%
      group_by(cohort) %>%
      summarise(
        n_rows = n(),
        ready_n = sum(.data[[ready_var]] == 1, na.rm = TRUE),
        ready_pct = ready_n / n_rows,
        score_non_missing_n = sum(!is.na(.data[[score_var]])),
        score_non_missing_pct = score_non_missing_n / n_rows,
        mean_component_n = mean(.data[[component_var]], na.rm = TRUE),
        median_component_n = stats::median(.data[[component_var]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(domain_group = domain_group_value) %>%
      relocate(domain_group, .before = cohort)
  }
)

cohort_variable_coverage <- domain_ready %>%
  group_by(cohort) %>%
  summarise(
    n_rows = n(),
    across(
      c(
        depression_score, loneliness_score, cognition_score,
        adl_count, iadl_count, self_rated_health
      ),
      ~ sum(!is.na(.x)),
      .names = "{.col}_non_missing_n"
    ),
    .groups = "drop"
  )

write_csv(domain_coverage_report, file.path(out_dir, "domain_coverage_report_v1.csv"), na = "")
write_csv(cohort_variable_coverage, file.path(out_dir, "cohort_variable_coverage_v1.csv"), na = "")

message("Done.")
message("Outputs written to: ", out_dir)
})

# build_harmonization_input_affective_functional_v1.R
local({
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

resolve_base_dir <- function() {
  path <- Sys.getenv("DCV_BASE_DIR", unset = "")
  if (!nzchar(path) || !dir.exists(path)) stop("Set data_root in config/config.R.")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

base_dir <- resolve_base_dir()
input_dir <- file.path(Sys.getenv("DCV_ANALYSIS_DIR"), "analysis_ready_core", "domain_ready_v1")
out_dir <- file.path(Sys.getenv("DCV_ANALYSIS_DIR"), "analysis_ready_core", "harmonization_input_v1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pooled_file <- file.path(input_dir, "pooled_domain_ready_v1.rds")
if (!file.exists(pooled_file)) {
  stop("Missing required input file: ", pooled_file)
}

domain_ready <- readRDS(pooled_file)

count_nonmissing <- function(...) {
  mat <- cbind(...)
  rowSums(!is.na(mat))
}

harmonization_input <- domain_ready %>%
  transmute(
    cohort,
    source_file,
    source_type,
    id,
    wave,
    year,
    country,
    age,
    sex,
    education_years,
    education_level,
    marital_status,
    partnered,
    race_ethnicity,
    weight_person,
    weight_household,
    weight_longitudinal,

    depression_score,
    loneliness_score,
    affective_score_v1,
    affective_component_n,
    affective_ready_v1,

    adl_count,
    iadl_count,
    adl_any,
    iadl_any,
    self_rated_health,
    functional_any_limitation_v1,
    functional_score_v1,
    functional_component_n,
    functional_ready_v1,

    affective_functional_joint_ready_v1 =
      ifelse(affective_ready_v1 == 1 & functional_ready_v1 == 1, 1, 0),
    affective_functional_complete_case_v1 =
      ifelse(
        count_nonmissing(
          affective_score_v1, functional_score_v1,
          age, sex, education_level, marital_status
        ) == 6,
        1, 0
      ),
    affective_source_profile_v1 = case_when(
      !is.na(depression_score) & !is.na(loneliness_score) ~ "depression_plus_loneliness",
      !is.na(depression_score) ~ "depression_only",
      !is.na(loneliness_score) ~ "loneliness_only",
      TRUE ~ NA_character_
    ),
    functional_source_profile_v1 = case_when(
      !is.na(adl_count) & !is.na(iadl_count) & !is.na(self_rated_health) ~ "adl_plus_iadl_plus_srh",
      !is.na(adl_count) & !is.na(iadl_count) ~ "adl_plus_iadl",
      !is.na(adl_count) & !is.na(self_rated_health) ~ "adl_plus_srh",
      !is.na(iadl_count) & !is.na(self_rated_health) ~ "iadl_plus_srh",
      !is.na(adl_count) ~ "adl_only",
      !is.na(iadl_count) ~ "iadl_only",
      !is.na(self_rated_health) ~ "srh_only",
      TRUE ~ NA_character_
    )
  )

split(harmonization_input, harmonization_input$cohort) |>
  purrr::iwalk(function(x, nm) {
    saveRDS(x, file.path(out_dir, paste0(tolower(nm), "_harmonization_input_affective_functional_v1.rds")))
  })

harmonization_summary <- harmonization_input %>%
  group_by(cohort) %>%
  summarise(
    n_rows = n(),
    affective_ready_n = sum(affective_ready_v1 == 1, na.rm = TRUE),
    functional_ready_n = sum(functional_ready_v1 == 1, na.rm = TRUE),
    joint_ready_n = sum(affective_functional_joint_ready_v1 == 1, na.rm = TRUE),
    complete_case_n = sum(affective_functional_complete_case_v1 == 1, na.rm = TRUE),
    affective_score_non_missing_n = sum(!is.na(affective_score_v1)),
    functional_score_non_missing_n = sum(!is.na(functional_score_v1)),
    age_non_missing_n = sum(!is.na(age)),
    sex_non_missing_n = sum(!is.na(sex)),
    education_level_non_missing_n = sum(!is.na(education_level)),
    marital_status_non_missing_n = sum(!is.na(marital_status)),
    .groups = "drop"
  )

harmonization_dictionary <- tibble(
  variable = names(harmonization_input),
  class = purrr::map_chr(harmonization_input, ~ paste(class(.x), collapse = ";")),
  non_missing_n = purrr::map_int(harmonization_input, ~ sum(!is.na(.x))),
  example = purrr::map_chr(harmonization_input, ~ {
    vals <- unique(.x[!is.na(.x)])
    if (length(vals) == 0) return("")
    paste(head(as.character(vals), 3), collapse = " | ")
  })
)

saveRDS(harmonization_input, file.path(out_dir, "pooled_harmonization_input_affective_functional_v1.rds"))
write_csv(harmonization_summary, file.path(out_dir, "harmonization_input_affective_functional_summary_v1.csv"), na = "")
write_csv(harmonization_dictionary, file.path(out_dir, "harmonization_input_affective_functional_dictionary_v1.csv"), na = "")

message("Done.")
message("Outputs written to: ", out_dir)
})

# patch_charls_cognition_v1.R
local({
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(haven)
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

resolve_base_dir <- function() {
  path <- Sys.getenv("DCV_BASE_DIR", unset = "")
  if (!nzchar(path) || !dir.exists(path)) stop("Set data_root in config/config.R.")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

base_dir <- resolve_base_dir()
domain_dir <- file.path(Sys.getenv("DCV_ANALYSIS_DIR"), "analysis_ready_core", "domain_ready_v1")
patch_dir <- file.path(Sys.getenv("DCV_ANALYSIS_DIR"), "analysis_ready_core", "domain_ready_v1_patches")
dir.create(patch_dir, recursive = TRUE, showWarnings = FALSE)

charls_domain_file <- file.path(domain_dir, "charls_domain_ready_v1.rds")
pooled_domain_file <- file.path(domain_dir, "pooled_domain_ready_v1.rds")
charls_raw_file <- file.path(base_dir, "CHARLS", "02_harmonized", "Regular Waves", "H_CHARLS_D_Data", "H_CHARLS_D_Data.dta")

for (f in c(charls_domain_file, pooled_domain_file, charls_raw_file)) {
  if (!file.exists(f)) stop("Missing required input file: ", f)
}

safe_pull <- function(df, nm, type = c("numeric", "character")) {
  type <- match.arg(type)
  if (!nm %in% names(df)) {
    if (type == "character") return(rep(NA_character_, nrow(df)))
    return(rep(NA_real_, nrow(df)))
  }
  x <- df[[nm]]
  if (inherits(x, "labelled")) x <- haven::zap_labels(x)
  if (type == "character") return(as.character(x))
  suppressWarnings(as.numeric(x))
}

na_codes_to_na <- function(x, extra = c(-1, -4, -9, -10, -11, -12, -13, -14, -15, -16, -17, -18, 88, 99, 999, 9999)) {
  x <- suppressWarnings(as.numeric(x))
  x[x %in% extra] <- NA_real_
  x
}

row_sum_or_na <- function(...) {
  mat <- cbind(...)
  out <- rowSums(mat, na.rm = TRUE)
  out[apply(is.na(mat), 1, all)] <- NA_real_
  out
}

z_score <- function(x) {
  x <- na_codes_to_na(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

charls_domain <- readRDS(charls_domain_file)
pooled_domain <- readRDS(pooled_domain_file)
charls_raw <- read_dta(charls_raw_file)

wave_to_year <- c(`1` = 2011, `2` = 2013, `3` = 2015, `4` = 2018)

charls_cognition_patch <- bind_rows(lapply(1:4, function(w) {
  imrc <- na_codes_to_na(safe_pull(charls_raw, paste0("r", w, "imrc")))
  orient <- na_codes_to_na(safe_pull(charls_raw, paste0("r", w, "orient")))
  tr20 <- na_codes_to_na(safe_pull(charls_raw, paste0("r", w, "tr20")))

  patch_raw <- row_sum_or_na(imrc, orient, tr20)

  tibble(
    id = safe_pull(charls_raw, "ID", "character"),
    wave = w,
    year = unname(wave_to_year[as.character(w)]),
    charls_imrc_patch_v1 = imrc,
    charls_orient_patch_v1 = orient,
    charls_tr20_patch_v1 = tr20,
    charls_cognition_patch_raw_v1 = patch_raw,
    charls_cognition_patch_component_n_v1 =
      rowSums(!is.na(cbind(imrc, orient, tr20))),
    charls_cognition_patch_source_v1 = ifelse(
      !is.na(patch_raw),
      "charls_harmonized_imrc_orient_tr20",
      NA_character_
    )
  )
}))

charls_patched <- charls_domain %>%
  left_join(charls_cognition_patch, by = c("id", "wave", "year")) %>%
  mutate(
    cognition_score_before_patch_v1 = cognition_score,
    cognition_score_patched_v1 = dplyr::coalesce(
      cognition_score,
      charls_cognition_patch_raw_v1
    ),
    basic_cognition_score_v1_before_patch = basic_cognition_score_v1,
    basic_cognition_score_v1_patched = z_score(cognition_score_patched_v1),
    cognition_ready_v1_before_patch = cognition_ready_v1,
    cognition_ready_v1_patched = ifelse(!is.na(cognition_score_patched_v1), 1, 0),
    basic_cognition_component_n_patched = ifelse(
      !is.na(cognition_score_patched_v1),
      1,
      0
    ),
    charls_cognition_patch_applied_v1 = ifelse(
      is.na(cognition_score_before_patch_v1) & !is.na(charls_cognition_patch_raw_v1),
      1,
      0
    )
  )

pooled_patched <- pooled_domain %>%
  filter(cohort != "CHARLS") %>%
  bind_rows(charls_patched)

patch_summary <- charls_patched %>%
  group_by(wave, year) %>%
  summarise(
    n_rows = n(),
    cognition_before_non_missing_n = sum(!is.na(cognition_score_before_patch_v1)),
    cognition_patch_non_missing_n = sum(!is.na(charls_cognition_patch_raw_v1)),
    cognition_after_non_missing_n = sum(!is.na(cognition_score_patched_v1)),
    patch_applied_n = sum(charls_cognition_patch_applied_v1 == 1, na.rm = TRUE),
    .groups = "drop"
  )

patch_dictionary <- tibble(
  variable = names(charls_patched),
  class = purrr::map_chr(charls_patched, ~ paste(class(.x), collapse = ";")),
  non_missing_n = purrr::map_int(charls_patched, ~ sum(!is.na(.x))),
  example = purrr::map_chr(charls_patched, ~ {
    vals <- unique(.x[!is.na(.x)])
    if (length(vals) == 0) return("")
    paste(head(as.character(vals), 3), collapse = " | ")
  })
)

saveRDS(charls_patched, file.path(patch_dir, "charls_domain_ready_v1_patched_cognition.rds"))
saveRDS(pooled_patched, file.path(patch_dir, "pooled_domain_ready_v1_patched_charls_cognition.rds"))
write_csv(patch_summary, file.path(patch_dir, "charls_cognition_patch_summary_v1.csv"), na = "")
write_csv(patch_dictionary, file.path(patch_dir, "charls_cognition_patch_dictionary_v1.csv"), na = "")

message("Done.")
message("Outputs written to: ", patch_dir)
})

# harmonize_affective_functional_v1.R
local({
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

resolve_base_dir <- function() {
  path <- Sys.getenv("DCV_BASE_DIR", unset = "")
  if (!nzchar(path) || !dir.exists(path)) stop("Set data_root in config/config.R.")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

base_dir <- resolve_base_dir()
input_dir <- file.path(Sys.getenv("DCV_ANALYSIS_DIR"), "analysis_ready_core", "harmonization_input_v1")
out_dir <- file.path(Sys.getenv("DCV_ANALYSIS_DIR"), "analysis_ready_core", "harmonized_affective_functional_v1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(input_dir, "pooled_harmonization_input_affective_functional_v1.rds")
if (!file.exists(input_file)) {
  stop("Missing required input file: ", input_file)
}

harm_input <- readRDS(input_file)

na_codes_to_na <- function(x, extra = c(-1, -3, -4, -9, -10, -11, -12, -13, -14, -15, -16, -17, -18, 88, 99, 999, 9999)) {
  x <- suppressWarnings(as.numeric(x))
  x[x %in% extra] <- NA_real_
  x
}

z_score <- function(x) {
  x <- na_codes_to_na(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

z_by_group <- function(x, grp) {
  ave(x, grp, FUN = z_score)
}

affective_functional_harmonized <- harm_input %>%
  mutate(
    across(
      c(age, sex, education_years, education_level, marital_status, partnered,
        depression_score, loneliness_score, adl_count, iadl_count, self_rated_health,
        affective_score_v1, functional_score_v1),
      ~ if (cur_column() %in% c("age", "age_baseline")) {
        ifelse(.x %in% c(88, 99), as.numeric(.x), na_codes_to_na(.x))
      } else na_codes_to_na(.x)
    ),
    sex_female_v1 = case_when(
      cohort == "SHARE" & sex %in% c(0, 1) ~ sex,
      sex == 2 ~ 1,
      sex == 1 ~ 0,
      TRUE ~ NA_real_
    ),
    education_level_clean_v1 = education_level,
    marital_status_clean_v1 = marital_status,
    affective_harmonized_v1 = z_by_group(affective_score_v1, cohort),
    functional_harmonized_v1 = z_by_group(functional_score_v1, cohort),
    affective_harmonized_wave_v1 = z_by_group(
      affective_score_v1,
      interaction(cohort, wave, drop = TRUE)
    ),
    functional_harmonized_wave_v1 = z_by_group(
      functional_score_v1,
      interaction(cohort, wave, drop = TRUE)
    ),
    affective_function_joint_harmonized_v1 = rowMeans(
      cbind(affective_harmonized_v1, functional_harmonized_v1),
      na.rm = TRUE
    ),
    affective_function_joint_harmonized_v1 = ifelse(
      is.nan(affective_function_joint_harmonized_v1),
      NA_real_,
      affective_function_joint_harmonized_v1
    ),
    harmonization_eligibility_v1 = ifelse(
      affective_functional_joint_ready_v1 == 1 &
        !is.na(age) &
        !is.na(sex_female_v1) &
        !is.na(education_level_clean_v1),
      1,
      0
    )
  )

split(affective_functional_harmonized, affective_functional_harmonized$cohort) |>
  purrr::iwalk(function(x, nm) {
    saveRDS(x, file.path(out_dir, paste0(tolower(nm), "_affective_functional_harmonized_v1.rds")))
  })

harmonization_summary <- affective_functional_harmonized %>%
  group_by(cohort) %>%
  summarise(
    n_rows = n(),
    affective_non_missing_n = sum(!is.na(affective_harmonized_v1)),
    functional_non_missing_n = sum(!is.na(functional_harmonized_v1)),
    joint_non_missing_n = sum(!is.na(affective_function_joint_harmonized_v1)),
    eligible_n = sum(harmonization_eligibility_v1 == 1, na.rm = TRUE),
    sex_female_non_missing_n = sum(!is.na(sex_female_v1)),
    education_non_missing_n = sum(!is.na(education_level_clean_v1)),
    marital_non_missing_n = sum(!is.na(marital_status_clean_v1)),
    .groups = "drop"
  )

harmonization_dictionary <- tibble(
  variable = names(affective_functional_harmonized),
  class = purrr::map_chr(affective_functional_harmonized, ~ paste(class(.x), collapse = ";")),
  non_missing_n = purrr::map_int(affective_functional_harmonized, ~ sum(!is.na(.x))),
  example = purrr::map_chr(affective_functional_harmonized, ~ {
    vals <- unique(.x[!is.na(.x)])
    if (length(vals) == 0) return("")
    paste(head(as.character(vals), 3), collapse = " | ")
  })
)

saveRDS(affective_functional_harmonized, file.path(out_dir, "pooled_affective_functional_harmonized_v1.rds"))
write_csv(harmonization_summary, file.path(out_dir, "affective_functional_harmonization_summary_v1.csv"), na = "")
write_csv(harmonization_dictionary, file.path(out_dir, "affective_functional_harmonization_dictionary_v1.csv"), na = "")

message("Done.")
message("Outputs written to: ", out_dir)
})

# integrate_charls_cognition_patch_v1.R
local({
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

resolve_base_dir <- function() {
  path <- Sys.getenv("DCV_BASE_DIR", unset = "")
  if (!nzchar(path) || !dir.exists(path)) stop("Set data_root in config/config.R.")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

base_dir <- resolve_base_dir()
domain_dir <- file.path(Sys.getenv("DCV_ANALYSIS_DIR"), "analysis_ready_core", "domain_ready_v1")
patch_dir <- file.path(Sys.getenv("DCV_ANALYSIS_DIR"), "analysis_ready_core", "domain_ready_v1_patches")
out_dir <- file.path(Sys.getenv("DCV_ANALYSIS_DIR"), "analysis_ready_core", "domain_ready_v1_integrated")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pooled_domain_file <- file.path(domain_dir, "pooled_domain_ready_v1.rds")
charls_patch_file <- file.path(patch_dir, "charls_domain_ready_v1_patched_cognition.rds")

for (f in c(pooled_domain_file, charls_patch_file)) {
  if (!file.exists(f)) stop("Missing required input file: ", f)
}

pooled_domain <- readRDS(pooled_domain_file)
charls_patch <- readRDS(charls_patch_file)

integrated <- pooled_domain %>%
  filter(cohort != "CHARLS") %>%
  bind_rows(
    charls_patch %>%
      mutate(
        cognition_score = cognition_score_patched_v1,
        basic_cognition_score_v1 = basic_cognition_score_v1_patched,
        basic_cognition_component_n = basic_cognition_component_n_patched,
        cognition_ready_v1 = cognition_ready_v1_patched
      )
  ) %>%
  arrange(cohort, id, wave, year)

integration_summary <- integrated %>%
  group_by(cohort) %>%
  summarise(
    n_rows = n(),
    cognition_score_non_missing_n = sum(!is.na(cognition_score)),
    basic_cognition_non_missing_n = sum(!is.na(basic_cognition_score_v1)),
    cognition_ready_n = sum(cognition_ready_v1 == 1, na.rm = TRUE),
    .groups = "drop"
  )

charls_delta_summary <- charls_patch %>%
  summarise(
    n_rows = n(),
    cognition_before_non_missing_n = sum(!is.na(cognition_score_before_patch_v1)),
    cognition_after_non_missing_n = sum(!is.na(cognition_score_patched_v1)),
    basic_cognition_before_non_missing_n = sum(!is.na(basic_cognition_score_v1_before_patch)),
    basic_cognition_after_non_missing_n = sum(!is.na(basic_cognition_score_v1_patched)),
    patch_applied_n = sum(charls_cognition_patch_applied_v1 == 1, na.rm = TRUE)
  )

integration_dictionary <- tibble(
  variable = names(integrated),
  class = purrr::map_chr(integrated, ~ paste(class(.x), collapse = ";")),
  non_missing_n = purrr::map_int(integrated, ~ sum(!is.na(.x))),
  example = purrr::map_chr(integrated, ~ {
    vals <- unique(.x[!is.na(.x)])
    if (length(vals) == 0) return("")
    paste(head(as.character(vals), 3), collapse = " | ")
  })
)

saveRDS(integrated, file.path(out_dir, "pooled_domain_ready_v1_integrated_charls_cognition.rds"))
write_csv(integration_summary, file.path(out_dir, "integrated_charls_cognition_summary_v1.csv"), na = "")
write_csv(charls_delta_summary, file.path(out_dir, "charls_cognition_integration_delta_v1.csv"), na = "")
write_csv(integration_dictionary, file.path(out_dir, "integrated_charls_cognition_dictionary_v1.csv"), na = "")

message("Done.")
message("Outputs written to: ", out_dir)
})

# build_basic_cognition_harmonization_input_v1.R
local({
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
  library(readr)
  library(tibble)
})

resolve_base_dir <- function() {
  path <- Sys.getenv("DCV_BASE_DIR", unset = "")
  if (!nzchar(path) || !dir.exists(path)) stop("Set data_root in config/config.R.")
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

base_dir <- resolve_base_dir()
input_dir <- file.path(Sys.getenv("DCV_ANALYSIS_DIR"), "analysis_ready_core", "domain_ready_v1_integrated")
out_dir <- file.path(Sys.getenv("DCV_ANALYSIS_DIR"), "analysis_ready_core", "harmonization_input_v1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(input_dir, "pooled_domain_ready_v1_integrated_charls_cognition.rds")
if (!file.exists(input_file)) {
  stop("Missing required input file: ", input_file)
}

cog_input_base <- readRDS(input_file)

na_codes_to_na <- function(x, extra = c(-1, -3, -4, -9, -10, -11, -12, -13, -14, -15, -16, -17, -18, 88, 99, 999, 9999)) {
  x <- suppressWarnings(as.numeric(x))
  x[x %in% extra] <- NA_real_
  x
}

z_score <- function(x) {
  x <- na_codes_to_na(x)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - mean(x, na.rm = TRUE)) / s
}

z_by_group <- function(x, grp) {
  ave(x, grp, FUN = z_score)
}

count_nonmissing <- function(...) {
  mat <- cbind(...)
  rowSums(!is.na(mat))
}

basic_cognition_harmonization_input <- cog_input_base %>%
  mutate(
    across(
      c(age, sex, education_years, education_level, marital_status,
        cognition_score, basic_cognition_score_v1),
      ~ if (cur_column() %in% c("age", "age_baseline")) {
        ifelse(.x %in% c(88, 99), as.numeric(.x), na_codes_to_na(.x))
      } else na_codes_to_na(.x)
    ),
    sex_female_v1 = case_when(
      cohort == "SHARE" & sex %in% c(0, 1) ~ sex,
      sex == 2 ~ 1,
      sex == 1 ~ 0,
      TRUE ~ NA_real_
    ),
    education_level_clean_v1 = education_level,
    marital_status_clean_v1 = marital_status,
    cognition_harmonized_v1 = z_by_group(cognition_score, cohort),
    cognition_harmonized_wave_v1 = z_by_group(
      cognition_score,
      interaction(cohort, wave, drop = TRUE)
    ),
    cognition_harmonization_eligibility_v1 = ifelse(
      cognition_ready_v1 == 1 &
        !is.na(age) &
        !is.na(sex_female_v1) &
        !is.na(education_level_clean_v1),
      1,
      0
    ),
    cognition_complete_case_v1 = ifelse(
      count_nonmissing(
        cognition_score, age, sex_female_v1,
        education_level_clean_v1, marital_status_clean_v1
      ) == 5,
      1,
      0
    ),
    cognition_source_profile_v1 = case_when(
      cohort == "CHARLS" & !is.na(charls_cognition_patch_source_v1) ~ charls_cognition_patch_source_v1,
      TRUE ~ "cohort_native_cognition_field"
    )
  ) %>%
  transmute(
    cohort,
    source_file,
    source_type,
    id,
    wave,
    year,
    country,
    age,
    sex,
    sex_female_v1,
    education_years,
    education_level,
    education_level_clean_v1,
    marital_status,
    marital_status_clean_v1,
    partnered,
    race_ethnicity,
    weight_person,
    weight_household,
    weight_longitudinal,
    cognition_score,
    basic_cognition_score_v1,
    cognition_harmonized_v1,
    cognition_harmonized_wave_v1,
    basic_cognition_component_n,
    cognition_ready_v1,
    cognition_harmonization_eligibility_v1,
    cognition_complete_case_v1,
    cognition_source_profile_v1
  )

split(basic_cognition_harmonization_input, basic_cognition_harmonization_input$cohort) |>
  purrr::iwalk(function(x, nm) {
    saveRDS(x, file.path(out_dir, paste0(tolower(nm), "_basic_cognition_harmonization_input_v1.rds")))
  })

cognition_summary <- basic_cognition_harmonization_input %>%
  group_by(cohort) %>%
  summarise(
    n_rows = n(),
    cognition_non_missing_n = sum(!is.na(cognition_harmonized_v1)),
    eligibility_n = sum(cognition_harmonization_eligibility_v1 == 1, na.rm = TRUE),
    complete_case_n = sum(cognition_complete_case_v1 == 1, na.rm = TRUE),
    age_non_missing_n = sum(!is.na(age)),
    sex_non_missing_n = sum(!is.na(sex_female_v1)),
    education_non_missing_n = sum(!is.na(education_level_clean_v1)),
    marital_non_missing_n = sum(!is.na(marital_status_clean_v1)),
    .groups = "drop"
  )

cognition_dictionary <- tibble(
  variable = names(basic_cognition_harmonization_input),
  class = purrr::map_chr(basic_cognition_harmonization_input, ~ paste(class(.x), collapse = ";")),
  non_missing_n = purrr::map_int(basic_cognition_harmonization_input, ~ sum(!is.na(.x))),
  example = purrr::map_chr(basic_cognition_harmonization_input, ~ {
    vals <- unique(.x[!is.na(.x)])
    if (length(vals) == 0) return("")
    paste(head(as.character(vals), 3), collapse = " | ")
  })
)

saveRDS(basic_cognition_harmonization_input, file.path(out_dir, "pooled_basic_cognition_harmonization_input_v1.rds"))
write_csv(cognition_summary, file.path(out_dir, "basic_cognition_harmonization_summary_v1.csv"), na = "")
write_csv(cognition_dictionary, file.path(out_dir, "basic_cognition_harmonization_dictionary_v1.csv"), na = "")

message("Done.")
message("Outputs written to: ", out_dir)
})
