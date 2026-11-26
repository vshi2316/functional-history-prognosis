## Within-between effects and random-effects summaries.

options(stringsAsFactors = FALSE)
base_dir <- Sys.getenv("DCV_BASE_DIR", unset = "")
if (!nzchar(base_dir)) stop("Set DCV_BASE_DIR before running this script.")
if (!dir.exists(base_dir)) stop("DCV_BASE_DIR does not exist: ", base_dir)
base_dir <- normalizePath(base_dir, winslash = "/", mustWork = TRUE)
robustness_dir <- file.path(Sys.getenv("DCV_ANALYSIS_DIR"), "robustness_analysis")
out_dir <- file.path(robustness_dir, "04_within_between_manuscript_freeze")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

pkgs <- c("dplyr", "readr", "tibble", "purrr", "nlme", "metafor")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs)) stop("Install before running: ", paste(missing_pkgs, collapse = ", "))
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(purrr)
  library(nlme); library(metafor)
})

input_rds <- file.path(Sys.getenv("DCV_ANALYSIS_DIR"), "analysis_ready_core", "longitudinal_models_v1",
                       "affective_functional_longitudinal_input_v1.rds")
sample_file <- file.path(robustness_dir, "03_lagged_fallback_sample_reconciliation",
                         "sample_branch_reconciliation.csv")
qualification_file <- file.path(robustness_dir, "01_longitudinal_robustness",
                                 "model_qualification_summary.csv")
qualification_effects_file <- file.path(robustness_dir, "01_longitudinal_robustness",
                                         "model_qualification_effects.csv")
lag_file <- file.path(robustness_dir, "03_lagged_fallback_sample_reconciliation",
                      "lagged_target_terms.csv")
required_files <- c(input_rds, sample_file, qualification_file,
                    qualification_effects_file, lag_file)
if (any(!file.exists(required_files))) stop("Missing required files: ", paste(required_files[!file.exists(required_files)], collapse = "; "))

## -------------------------------------------------------------------------
## 1. Freeze and verify sample definitions
## -------------------------------------------------------------------------
sample_branches <- read_csv(sample_file, show_col_types = FALSE)
expected <- tribble(
  ~branch, ~persons_expected, ~observations_expected,
  "affective_repeated_plus_baseline_covariates", 189732, 803080,
  "functional_repeated_plus_baseline_covariates", 208751, 953301,
  "joint_repeated_plus_baseline_covariates", 189718, 802984,
  "harmonization_eligible_repeated", 188187, 791366
)
sample_check <- sample_branches %>% filter(cohort == "TOTAL") %>%
  select(branch, persons_observed = persons, observations_observed = observations) %>%
  right_join(expected, by = "branch") %>%
  mutate(pass = persons_observed == persons_expected & observations_observed == observations_expected)
write_csv(sample_check, file.path(out_dir, "sample_freeze_check.csv"), na = "")
if (!all(sample_check$pass)) stop("Sample counts differ from the frozen input release.")

## -------------------------------------------------------------------------
## 2. Within-between decomposition
## -------------------------------------------------------------------------
dat <- readRDS(input_rds)
numeric_vars <- c("year", "wave", "age_baseline", "sex_baseline", "education_baseline",
                  "affective_harmonized_v1", "functional_harmonized_v1")
dat <- dat %>% mutate(across(all_of(numeric_vars), ~suppressWarnings(as.numeric(.x)))) %>%
  filter(if_all(all_of(c("year", "age_baseline", "sex_baseline", "education_baseline",
                         "affective_harmonized_v1", "functional_harmonized_v1")), ~!is.na(.x))) %>%
  arrange(cohort, subject_id, year, wave)
eligible <- dat %>% count(cohort, subject_id, name = "n_joint") %>% filter(n_joint >= 2)
dat <- dat %>% inner_join(eligible, by = c("cohort", "subject_id")) %>%
  group_by(cohort, subject_id) %>%
  mutate(
    baseline_year = min(year),
    time_since_baseline_years = year - baseline_year,
    affective_between = mean(affective_harmonized_v1),
    affective_within = affective_harmonized_v1 - affective_between,
    functional_between = mean(functional_harmonized_v1),
    functional_within = functional_harmonized_v1 - functional_between
  ) %>% ungroup()

specs <- tribble(
  ~outcome, ~model_id, ~formula_text, ~between_term, ~within_term,
  "functional_harmonized_v1", "functional_on_affective_within_between",
  "functional_harmonized_v1 ~ affective_between + affective_within + time_since_baseline_years + age_baseline + sex_baseline + education_baseline",
  "affective_between", "affective_within",
  "affective_harmonized_v1", "affective_on_functional_within_between",
  "affective_harmonized_v1 ~ functional_between + functional_within + time_since_baseline_years + age_baseline + sex_baseline + education_baseline",
  "functional_between", "functional_within"
)

fit_one <- function(d, cohort_name, outcome, model_id, formula_text, between_term, within_term) {
  x <- d %>% filter(cohort == cohort_name)
  n_obs <- nrow(x); n_subjects <- n_distinct(x$subject_id)
  fit <- tryCatch(
    nlme::lme(fixed = as.formula(formula_text), random = ~1 | subject_id,
      data = x, method = "ML", na.action = na.omit,
      control = nlme::lmeControl(msMaxIter = 200, opt = "optim", returnObject = FALSE)),
    error = function(e) e
  )
  if (inherits(fit, "error")) return(list(
    effects = tibble(cohort = cohort_name, outcome, model_id, term = NA_character_, estimate = NA_real_,
      std_error = NA_real_, p_value = NA_real_, n_obs, n_subjects, status = "failed",
      error_message = conditionMessage(fit)),
    contrast = tibble(cohort = cohort_name, outcome, model_id, contrast = "between_minus_within",
      estimate = NA_real_, std_error = NA_real_, p_value = NA_real_, n_obs, n_subjects,
      status = "failed", error_message = conditionMessage(fit))
  ))
  tt <- as.data.frame(summary(fit)$tTable); tt$term <- rownames(tt); rownames(tt) <- NULL
  effects <- tibble(cohort = cohort_name, outcome, model_id, term = tt$term,
    estimate = tt$Value, std_error = tt$Std.Error, p_value = tt$`p-value`,
    n_obs, n_subjects, status = "ok", error_message = NA_character_)
  b <- nlme::fixef(fit); V <- vcov(fit)
  delta <- unname(b[between_term] - b[within_term])
  delta_se <- sqrt(V[between_term, between_term] + V[within_term, within_term] - 2 * V[between_term, within_term])
  contrast <- tibble(cohort = cohort_name, outcome, model_id, contrast = "between_minus_within",
    estimate = delta, std_error = delta_se, p_value = 2 * pnorm(abs(delta / delta_se), lower.tail = FALSE),
    n_obs, n_subjects, status = "ok", error_message = NA_character_)
  list(effects = effects, contrast = contrast)
}

cohorts <- sort(unique(dat$cohort))
fit_results <- pmap(specs, function(outcome, model_id, formula_text, between_term, within_term) {
  map(cohorts, ~fit_one(dat, .x, outcome, model_id, formula_text, between_term, within_term))
}) %>% flatten()
effects <- map_dfr(fit_results, "effects")
contrasts <- map_dfr(fit_results, "contrast")
target_terms <- c("affective_between", "affective_within", "functional_between", "functional_within")
effects <- effects %>% mutate(p_fdr_sensitivity_family = if_else(term %in% target_terms,
  p.adjust(if_else(term %in% target_terms, p_value, NA_real_), method = "BH"), NA_real_))
contrasts <- contrasts %>% mutate(p_fdr_sensitivity_family = p.adjust(p_value, method = "BH"))
write_csv(effects, file.path(out_dir, "within_between_effects.csv"), na = "")
write_csv(contrasts, file.path(out_dir, "between_minus_within_contrasts.csv"), na = "")

## Random-effects summaries are descriptive because k=5 and heterogeneity is expected.
meta_input <- effects %>% filter(status == "ok", term %in% target_terms, is.finite(estimate),
                                 is.finite(std_error), std_error > 0)
meta_summary <- meta_input %>% group_by(model_id, outcome, term) %>% group_modify(~{
  fit <- metafor::rma.uni(yi = .x$estimate, sei = .x$std_error, method = "REML")
  tibble(k = nrow(.x), estimate = as.numeric(fit$b), std_error = fit$se,
         ci_lb = fit$ci.lb, ci_ub = fit$ci.ub, p_value = fit$pval,
         I2 = fit$I2, tau2 = fit$tau2)
}) %>% ungroup()
write_csv(meta_input, file.path(out_dir, "within_between_meta_input.csv"), na = "")
write_csv(meta_summary, file.path(out_dir, "within_between_meta_summary.csv"), na = "")

## -------------------------------------------------------------------------
## 3. Model qualification and lagged estimates
## -------------------------------------------------------------------------
qualification <- read_csv(qualification_file, show_col_types = FALSE)
qualification_effects <- read_csv(qualification_effects_file, show_col_types = FALSE)
lag_targets <- read_csv(lag_file, show_col_types = FALSE)

model_status <- qualification %>% mutate(family = case_when(
  grepl("_RI$", model_id) ~ "random_intercept",
  grepl("_RS$", model_id) ~ "random_slope",
  grepl("_CAR1$", model_id) ~ "CAR1",
  TRUE ~ "other")) %>% count(family, status, name = "n_models")
write_csv(model_status, file.path(out_dir, "model_status_audit.csv"), na = "")

primary_terms <- qualification_effects %>% filter(
  (model_id == "functional_on_affective_RI" & term == "affective_harmonized_v1") |
  (model_id == "affective_on_functional_RI" & term == "functional_harmonized_v1") |
  (model_id == "functional_on_affective_RS" & term == "affective_harmonized_v1") |
  (model_id == "affective_on_functional_RS" & term == "functional_harmonized_v1"))
write_csv(primary_terms, file.path(out_dir, "RI_RS_target_effects.csv"), na = "")
write_csv(lag_targets, file.path(out_dir, "lag_delta_and_FE_target_effects.csv"), na = "")

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo.txt"))
