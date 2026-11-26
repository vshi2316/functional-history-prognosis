# Copy to config/config.R. Run entry points from the repository root.
project_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_root <- Sys.getenv("DCV_BASE_DIR", unset = file.path(project_dir, "data"))
analysis_dir <- file.path(project_dir, "results", "pipeline")
frozen <- file.path(data_root, "V17D_frozen_inputs", "v1")
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)

elsa_file <- file.path(data_root, "ELSA",
  "5050stata_A06F12276F837D82374CA07FCEA89BDC3945439B39E962E7B03B39CEBABE46BC_V1",
  "UKDA-5050-stata", "stata", "stata13_se", "gh_elsa_h.dta")
share_file <- file.path(data_root, "SHARE", "GH_SHARE_g_rel9-0-0_ALL_datasets_stata", "GH_SHARE_g.dta")

Sys.setenv(
  DCV_PROJECT_DIR = project_dir,
  DCV_BASE_DIR = data_root,
  DCV_ANALYSIS_DIR = analysis_dir,
  MORTALITY_LONGITUDINAL_INPUT = file.path(analysis_dir, "analysis_ready_core",
    "longitudinal_models_v1", "affective_functional_longitudinal_input_v1.rds"),
  MORTALITY_EXACT_VISIT_DATE_INPUT = file.path(frozen, "V17D_exact_visit_dates_QC_v1.csv"),
  MORTALITY_WAVE_STATUS_INPUT = file.path(frozen, "V17D_wave_direct_proxy_status_QC_v1.csv"),
  MORTALITY_ENDPOINT_MAIN_INPUT = file.path(frozen, "V17D_mortality_endpoint_main_QC_v1.rds"),
  MORTALITY_ENDPOINT_SUPPORT_INPUT = file.path(frozen, "V17D_mortality_endpoint_CHARLS_QC_v1.rds"),
  MORTALITY_OUTPUT_ROOT = file.path(project_dir, "results", "mortality_bridge"),
  MORTALITY_PREDICTION_CANDIDATE_INPUT = file.path(project_dir, "results", "mortality_bridge", "candidate", "dynamic_candidate_exact.rds"),
  MORTALITY_PREDICTION_OUTPUT_ROOT = file.path(project_dir, "results", "mortality_prediction"),
  CLINICAL_LANDMARK_INPUT = file.path(project_dir, "results", "mortality_prediction", "landmark", "age50_landmark_dataset.rds"),
  CLINICAL_PRIOR1_OUTPUT_ROOT = file.path(project_dir, "results", "clinical_risk_prior1_continuous"),
  CLINICAL_RECALIBRATION_OUTPUT_ROOT = file.path(project_dir, "results", "clinical_risk_recalibration"),
  ADL_CANDIDATE_INPUT = file.path(project_dir, "results", "mortality_bridge", "candidate", "dynamic_candidate_exact.rds"),
  ADL_DATES_INPUT = file.path(frozen, "V17D_exact_visit_dates_QC_v1.csv"),
  ADL_DIRECT_INPUT = file.path(frozen, "V17D_wave_direct_proxy_status_QC_v1.csv"),
  ADL_SELFCARE_INPUT = file.path(project_dir, "results", "ADL_first_recorded", "selfcare5_overlay.rds"),
  ADL_CUTPOINTS_INPUT = file.path(project_dir, "results", "ADL_first_recorded", "cohort_cutpoints.csv"),
  ADL_MORTALITY_INPUT = file.path(frozen, "V17D_mortality_endpoint_main_QC_v1.rds"),
  ADL_OUTPUT_ROOT = file.path(project_dir, "results", "ADL_first_recorded"),
  DEMENTIA_DYNAMIC_CANDIDATE_INPUT = file.path(project_dir, "results", "mortality_bridge", "candidate", "dynamic_candidate_exact.rds"),
  DEMENTIA_ID_CROSSWALK_INPUT = file.path(project_dir, "results", "dementia_bridge", "provider_crosswalk.rds"),
  DEMENTIA_ELSA_OUTCOME_INPUT = elsa_file,
  DEMENTIA_SHARE_OUTCOME_INPUT = share_file,
  DEMENTIA_ELSA_COGNITION_INPUT = file.path(analysis_dir, "analysis_ready_core", "ready_v3_core4", "elsa_analysis_ready_long_v3_core4.rds"),
  DEMENTIA_SHARE_COGNITION_INPUT = file.path(analysis_dir, "analysis_ready_core", "ready_v3_core4", "share_analysis_ready_long_v3_core4.rds"),
  DEMENTIA_OUTPUT_ROOT = file.path(project_dir, "results", "dementia_bridge")
)
