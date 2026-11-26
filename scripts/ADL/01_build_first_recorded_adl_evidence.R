#!/usr/bin/env Rscript
# Reproduce the five-cohort first-recorded basic ADL analysis.

options(stringsAsFactors = FALSE, warn = 1)
if (!requireNamespace("data.table", quietly = TRUE)) stop("Install data.table before running this script.")
suppressPackageStartupMessages(library(data.table))

project_dir <- normalizePath(Sys.getenv("DCV_PROJECT_DIR", unset = getwd()),
                             winslash = "/", mustWork = TRUE)
data_root <- Sys.getenv("DCV_BASE_DIR", unset = "")
if (!nzchar(data_root)) stop("Set DCV_BASE_DIR to the private data root.")

frozen <- file.path(data_root, "V17D_frozen_inputs", "v1")
inputs <- c(
  candidate = Sys.getenv("ADL_CANDIDATE_INPUT", unset = file.path(frozen, "V17D_dynamic_candidate_exact_v1.rds")),
  dates = Sys.getenv("ADL_DATES_INPUT", unset = file.path(frozen, "V17D_exact_visit_dates_QC_v1.csv")),
  direct = Sys.getenv("ADL_DIRECT_INPUT", unset = file.path(frozen, "V17D_wave_direct_proxy_status_QC_v1.csv")),
  selfcare = Sys.getenv("ADL_SELFCARE_INPUT", unset = file.path(project_dir, "results", "ADL_first_recorded", "selfcare5_overlay.rds")),
  cutpoints = Sys.getenv("ADL_CUTPOINTS_INPUT"),
  mortality = Sys.getenv("ADL_MORTALITY_INPUT", unset = file.path(frozen, "V17D_mortality_endpoint_main_QC_v1.rds")),
  mortality_charls = Sys.getenv("MORTALITY_ENDPOINT_SUPPORT_INPUT", unset = file.path(frozen, "V17D_mortality_endpoint_CHARLS_QC_v1.rds"))
)
if (!nzchar(inputs["cutpoints"]) || !all(file.exists(inputs))) {
  stop("Set ADL_CUTPOINTS_INPUT and provide the seven ADL input files.")
}
out <- Sys.getenv("ADL_OUTPUT_ROOT", unset = file.path(project_dir, "results", "ADL_first_recorded"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
save_csv <- function(x, name) fwrite(x, file.path(out, name), na = "")

normalise_id <- function(x) {
  x <- toupper(gsub("\\s+", "", trimws(as.character(x))))
  x <- sub("\\.0+$", "", x)
  numeric_id <- !is.na(x) & grepl("^[0-9]+$", x)
  x[numeric_id] <- sub("^0+", "", x[numeric_id])
  x
}

cohorts <- c("HRS", "ELSA", "SHARE", "CHARLS", "MHAS")
dates <- fread(inputs["dates"],
               select = c("cohort", "id_norm", "wave", "interview_year_raw", "interview_month_raw"),
               colClasses = c(id_norm = "character"))
dates <- dates[cohort %in% cohorts & !is.na(interview_year_raw) & interview_month_raw %in% 1:12]
dates[, visit_date := as.Date(sprintf("%04d-%02d-15", as.integer(interview_year_raw),
                                      as.integer(interview_month_raw)))]
calendar <- dates[, .(median_date = as.Date(round(median(as.numeric(visit_date))),
                                     origin = "1970-01-01")), by = .(cohort, wave)]
setorder(calendar, cohort, wave)
calendar[, `:=`(next_wave = shift(wave, type = "lead"),
                next_median_date = shift(median_date, type = "lead")), by = cohort]
calendar[, planned_gap_months :=
  (as.integer(format(next_median_date, "%Y")) - as.integer(format(median_date, "%Y"))) * 12L +
  as.integer(format(next_median_date, "%m")) - as.integer(format(median_date, "%m"))]
calendar[, eligible_pair := !is.na(next_wave) & planned_gap_months >= 12 &
            planned_gap_months <= 36 & !(cohort == "SHARE" & wave %in% c(6L, 7L))]
save_csv(calendar, "wave_calendar.csv")

candidate <- as.data.table(readRDS(inputs["candidate"]))[
  cohort %in% cohorts & prior_n >= 1L & interval_start_age >= 50 &
    interval_start_age < 65 & is.finite(H_F) & is.finite(D_F) &
    is.finite(H_A) & is.finite(D_A) & is.finite(education_z) &
    !is.na(sex_baseline),
  .(cohort, subject_id, id = normalise_id(sub("^[^:]+::", "", subject_id)),
    wave, index_date = as.Date(interval_start_date), age = interval_start_age,
    prior_n, H_F, current_F = H_F + D_F, H_A, D_A,
    sex = factor(sex_baseline), education_z)]
candidate <- candidate[!is.na(index_date) & !is.na(id)]
cutpoints <- fread(inputs["cutpoints"], select = c("cohort", "current_p25", "current_p75"))
candidate <- merge(candidate, cutpoints, by = "cohort")
candidate <- candidate[current_F >= current_p25 - 1e-8 & current_F <= current_p75 + 1e-8]
candidate <- merge(candidate, calendar[eligible_pair == TRUE,
  .(cohort, wave, next_wave, planned_gap_months)], by = c("cohort", "wave"))

ids <- unique(candidate[, .(cohort, id)])
adl <- as.data.table(readRDS(inputs["selfcare"]))
adl[, id := normalise_id(id)]
adl <- adl[ids, on = .(cohort, id), nomatch = 0L,
           .(cohort, id, wave, adl_any = selfcare5_any,
             adl_items = selfcare5_n_observed)]
dates[, `:=`(id = normalise_id(id_norm), id_norm = NULL)]
dates <- dates[ids, on = .(cohort, id), nomatch = 0L, .(cohort, id, wave, visit_date)]
direct <- fread(inputs["direct"],
                select = c("cohort", "id_norm", "wave", "proxy_status"),
                colClasses = c(id_norm = "character"))
direct <- direct[proxy_status == "direct"]
direct[, `:=`(id = normalise_id(id_norm), id_norm = NULL, proxy_status = NULL)]
direct <- direct[ids, on = .(cohort, id), nomatch = 0L]
visits <- merge(merge(dates, direct, by = c("cohort", "id", "wave")),
                adl, by = c("cohort", "id", "wave"), all.x = TRUE)
visits[, adl_complete := adl_items == 5 & adl_any %in% 0:1]

baseline <- visits[, .(cohort, id, wave, baseline_date = visit_date,
                       baseline_adl = adl_any, baseline_complete = adl_complete)]
dat <- merge(candidate, baseline, by = c("cohort", "id", "wave"))
dat <- dat[baseline_complete == TRUE & baseline_adl == 0 & baseline_date == index_date]
setorder(dat, cohort, subject_id, index_date, wave)
dat <- dat[, .SD[1L], by = .(cohort, subject_id)]
next_visit <- visits[, .(cohort, id, next_wave = wave, next_date = visit_date,
                         next_adl = adl_any, next_complete = adl_complete)]
dat <- merge(dat, next_visit, by = c("cohort", "id", "next_wave"), all.x = TRUE)
dat[, observed := !is.na(next_date) & next_complete %in% TRUE]
dat[, outcome := fifelse(observed, as.integer(next_adl == 1), NA_integer_)]

earlier <- merge(dates[, .(cohort, id, wave, visit_date)],
                 adl[, .(cohort, id, wave, adl_any, adl_items)],
                 by = c("cohort", "id", "wave"), all.x = TRUE)
earlier[, complete := adl_items == 5 & adl_any %in% 0:1]
earlier <- merge(earlier, dat[, .(cohort, id, index_date)],
                 by = c("cohort", "id"), allow.cartesian = TRUE)
earlier <- earlier[visit_date < index_date,
                   .(earlier_visits = .N,
                     complete_earlier = sum(complete, na.rm = TRUE),
                     earlier_difficulty = sum(complete & adl_any == 1, na.rm = TRUE)),
                   by = .(cohort, id, index_date)]
dat <- merge(dat, earlier, by = c("cohort", "id", "index_date"), all.x = TRUE)
dat[, eligible_first_record := !is.na(earlier_visits) &
      earlier_visits == complete_earlier & earlier_difficulty == 0]
flow <- dat[, .(current_adl_free = .N,
                eligible_first_record = sum(eligible_first_record),
                next_adl_observed = sum(eligible_first_record & observed),
                first_recorded_difficulty = sum(eligible_first_record & outcome == 1,
                                                na.rm = TRUE)), by = cohort]
flow[, cohort_order := match(cohort, cohorts)]
setorder(flow, cohort_order)
flow[, cohort_order := NULL]
save_csv(flow, "first_recorded_flow.csv")

analysis <- dat[eligible_first_record == TRUE & observed == TRUE]
analysis[, `:=`(year_c = as.integer(format(index_date, "%Y")) - 2010L,
                gap_years = planned_gap_months / 12)]
set.seed(20261110)
hrs_ids <- sort(unique(analysis[cohort == "HRS", subject_id]))
train_ids <- sample(hrs_ids, floor(0.70 * length(hrs_ids)))
analysis[, evaluation_set := fifelse(cohort != "HRS", "cross_cohort_application",
                                     fifelse(subject_id %in% train_ids,
                                             "HRS_development", "HRS_holdout"))]

formula_current <- outcome ~ current_F + H_A + D_A + age + sex +
  education_z + year_c + gap_years
train <- analysis[evaluation_set == "HRS_development"]
model_current <- glm(formula_current, data = train, family = binomial())
model_history <- update(model_current, . ~ . + H_F)
save_csv(rbindlist(list(
  data.table(model = "current", term = names(coef(model_current)),
             estimate = unname(coef(model_current))),
  data.table(model = "current_plus_history", term = names(coef(model_history)),
             estimate = unname(coef(model_history)))
)), "first_recorded_HRS_coefficients.csv")

estimate_association <- function(x) {
  q <- as.numeric(quantile(x$H_F, c(0.25, 0.75)))
  fit <- glm(update(formula_current, . ~ . + H_F), data = x, family = binomial())
  coefficients <- coef(fit)
  estimable <- !is.na(coefficients)
  vc <- vcov(fit)[estimable, estimable, drop = FALSE]
  X0 <- X1 <- model.matrix(fit)
  X0[, "H_F"] <- q[1]
  X1[, "H_F"] <- q[2]
  X0 <- X0[, estimable, drop = FALSE]
  X1 <- X1[, estimable, drop = FALSE]
  p0 <- plogis(drop(X0 %*% coefficients[estimable]))
  p1 <- plogis(drop(X1 %*% coefficients[estimable]))
  grad <- colMeans(X1 * as.numeric(p1 * (1 - p1))) -
    colMeans(X0 * as.numeric(p0 * (1 - p0)))
  rd <- mean(p1 - p0)
  rd_se <- sqrt(drop(t(grad) %*% vc %*% grad))
  or_se <- sqrt(vc["H_F", "H_F"])
  data.table(cohort = x$cohort[1], participants = nrow(x),
             events = sum(x$outcome), H_q25 = q[1], H_q75 = q[2],
             odds_ratio = exp(coefficients["H_F"] * diff(q)),
             or_low = exp((coefficients["H_F"] - 1.96 * or_se) * diff(q)),
             or_high = exp((coefficients["H_F"] + 1.96 * or_se) * diff(q)),
             risk_q25_pct = 100 * mean(p0), risk_q75_pct = 100 * mean(p1),
             risk_difference_per_1000 = 1000 * rd,
             rd_low_per_1000 = 1000 * (rd - 1.96 * rd_se),
             rd_high_per_1000 = 1000 * (rd + 1.96 * rd_se))
}
association <- rbindlist(lapply(cohorts, function(z) {
  estimate_association(analysis[cohort == z])
}))
save_csv(association, "first_recorded_association.csv")

validation <- analysis[evaluation_set != "HRS_development"]
validation[, `:=`(
  risk_current = as.numeric(predict(model_current, newdata = .SD, type = "response")),
  risk_history = as.numeric(predict(model_history, newdata = .SD, type = "response"))
)]
auc <- function(y, p) {
  ranks <- rank(p, ties.method = "average")
  (sum(ranks[y == 1]) - sum(y == 1) * (sum(y == 1) + 1) / 2) /
    (sum(y == 1) * sum(y == 0))
}
brier <- function(y, p) mean((y - p)^2)
assess <- function(x) {
  y <- x$outcome
  p0 <- x$risk_current
  p1 <- x$risk_history
  boot <- replicate(300L, {
    i <- sample.int(length(y), replace = TRUE)
    c(auc(y[i], p1[i]) - auc(y[i], p0[i]),
      brier(y[i], p0[i]) - brier(y[i], p1[i]))
  })
  ai <- quantile(boot[1, ], c(0.025, 0.975))
  bi <- quantile(boot[2, ], c(0.025, 0.975))
  data.table(cohort = x$cohort[1], evaluation_set = x$evaluation_set[1],
             participants = length(y), events = sum(y),
             observed_pct = 100 * mean(y),
             predicted_current_pct = 100 * mean(p0),
             predicted_history_pct = 100 * mean(p1),
             auc_current = auc(y, p0), auc_history = auc(y, p1),
             auc_gain = auc(y, p1) - auc(y, p0),
             auc_gain_low = ai[1], auc_gain_high = ai[2],
             brier_current = brier(y, p0), brier_history = brier(y, p1),
             brier_gain = brier(y, p0) - brier(y, p1),
             brier_gain_low = bi[1], brier_gain_high = bi[2])
}
set.seed(20261111)
performance <- rbindlist(lapply(cohorts, function(z) assess(validation[cohort == z])))
save_csv(performance, "first_recorded_prediction.csv")

mortality_columns <- c("cohort", "subject_id", "primary_event",
                       "death_date_midmonth", "primary_exit_date")
mortality <- rbindlist(list(
  as.data.table(readRDS(inputs["mortality"]))[cohort != "CHARLS", ..mortality_columns],
  as.data.table(readRDS(inputs["mortality_charls"]))[, ..mortality_columns]
))[, .(cohort, subject_id, primary_event,
        death_date = as.Date(death_date_midmonth),
        mortality_exit_date = as.Date(primary_exit_date))]
stopifnot(!anyDuplicated(mortality[, .(cohort, subject_id)]))
audit <- merge(dat[eligible_first_record == TRUE,
                   .(cohort, subject_id, wave, index_date, H_F,
                     next_date, observed, outcome)],
               mortality, by = c("cohort", "subject_id"), all.x = TRUE)
audit <- merge(audit, calendar[eligible_pair == TRUE,
                               .(cohort, wave, scheduled_next = next_median_date)],
               by = c("cohort", "wave"), all.x = TRUE)
audit[, audit_cutoff := fifelse(!is.na(next_date), next_date, scheduled_next)]
stopifnot(!anyNA(audit$audit_cutoff))
audit[, death_before_cutoff := primary_event %in% TRUE &
        !is.na(death_date) & death_date > index_date & death_date < audit_cutoff]
audit[, `:=`(death_before_next = !observed & death_before_cutoff,
             timing_conflict = observed & death_before_cutoff)]
audit[, unresolved := !observed & !death_before_next]
audit[, mortality_covered_to_cutoff := !is.na(mortality_exit_date) &
        mortality_exit_date >= audit_cutoff]
quartiles <- audit[, .(q25 = as.numeric(quantile(H_F, 0.25)),
                       q75 = as.numeric(quantile(H_F, 0.75))), by = cohort]
stopifnot(all(quartiles$q25 < quartiles$q75))
audit <- merge(audit, quartiles, by = "cohort")
audit[, history_group := fifelse(H_F <= q25, "at_or_below_p25",
                                 fifelse(H_F >= q75, "at_or_above_p75", "middle"))]
observation_flow <- audit[, .(eligible = .N,
  next_adl_observed = sum(observed),
  next_adl_events = sum(outcome == 1, na.rm = TRUE),
  death_before_next = sum(death_before_next),
  unresolved = sum(unresolved),
  no_next_direct_visit = sum(is.na(next_date)),
  next_adl_incomplete = sum(!is.na(next_date) & !observed),
  timing_conflicts = sum(timing_conflict),
  missing_mortality_record = sum(is.na(primary_event)),
  unresolved_without_mortality_followup = sum(unresolved & !mortality_covered_to_cutoff)
), by = cohort]
observation_flow[, order := match(cohort, cohorts)]
setorder(observation_flow, order)
observation_flow[, order := NULL]
save_csv(observation_flow, "observation_flow.csv")

strata <- audit[, .(participants = .N,
                    observed = sum(observed),
                    observed_pct = 100 * mean(observed),
                    adl_events_observed = sum(outcome == 1, na.rm = TRUE),
                    death_before_next = sum(death_before_next),
                    death_pct = 100 * mean(death_before_next),
                    unresolved = sum(unresolved),
                    unresolved_pct = 100 * mean(unresolved)),
                by = .(cohort, history_group)]
strata[, cohort_order := match(cohort, cohorts)]
strata[, group_order := match(history_group,
                              c("at_or_below_p25", "middle", "at_or_above_p75"))]
setorder(strata, cohort_order, group_order)
strata[, c("cohort_order", "group_order") := NULL]
save_csv(strata, "observation_by_history.csv")

contrast <- rbindlist(lapply(cohorts, function(z) {
  low <- audit[cohort == z & history_group == "at_or_below_p25"]
  high <- audit[cohort == z & history_group == "at_or_above_p75"]
  difference <- function(field) {
    p0 <- mean(low[[field]])
    p1 <- mean(high[[field]])
    se <- sqrt(p0 * (1 - p0) / nrow(low) + p1 * (1 - p1) / nrow(high))
    c(100 * (p1 - p0), 100 * (p1 - p0 - 1.96 * se),
      100 * (p1 - p0 + 1.96 * se))
  }
  obs <- difference("observed")
  death <- difference("death_before_next")
  data.table(cohort = z, low_n = nrow(low), high_n = nrow(high),
             observation_difference_pp = obs[1],
             observation_low_pp = obs[2], observation_high_pp = obs[3],
             death_difference_pp = death[1],
             death_low_pp = death[2], death_high_pp = death[3])
}))
save_csv(contrast, "observation_history_contrasts.csv")

writeLines(capture.output(sessionInfo()), file.path(out, "adl_sessionInfo.txt"), useBytes = TRUE)
message("ADL results written to ", normalizePath(out, winslash = "/"))
