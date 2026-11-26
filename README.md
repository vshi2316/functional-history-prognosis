# Functional history and prognosis across five ageing cohorts

Repository name: `functional-history-prognosis`.

Analysis code for *Functional history distinguishes survival and first-recorded basic activity difficulty at comparable current function*.

The workflow uses CHARLS, ELSA, HRS, MHAS and SHARE. It separates strictly earlier functional burden from the deviation at the current assessment and relates these components to subsequent clinical events. Mortality is the primary outcome. First-recorded basic activity difficulty provides a functional outcome; dementia-related analyses use the ELSA and SHARE records.

## Age cleaning

Age and baseline-age fields retain valid values of 88 and 99 years. The original missing-value codes continue to apply to other numeric fields. The longitudinal preparation and downstream robustness scripts use this variable-specific rule.

The repeated harmonisation-eligible branch includes 188,187 participants and 791,366 observations. The functional branch includes 208,751 participants and 953,301 observations. `Longitudinal_robustness/05` checks the sample sizes for all four longitudinal branches.

## Analysis

For a current score `X`, history `H` is the mean of earlier observed scores and deviation is `D = X - H`. History and deviation use the same cohort-specific scale. The linear predictor can therefore be written as

```text
beta_H * H + beta_D * D
    = beta_D * X + (beta_H - beta_D) * H
```

The coefficient difference estimates the association with earlier burden while holding the current score constant. The longitudinal within-between models are a separate analysis using each person's observed mean and person-centred scores.

The clinical workflow includes cohort harmonisation and longitudinal models, followed by time-updated mortality models. Five-year mortality models are developed in HRS and evaluated with fixed coefficients in other cohorts. A separate analysis examines the central range of current function at ages 50–64. The next-wave activity analysis requires complete five-item assessments with no recorded difficulty through the index visit. Dementia-related cause-specific models retain death as a competing event.

## Files

```text
config/config.example.R            data paths and output paths
config/selfcare_items.csv          five-item definitions by cohort and wave
scripts/Track_A/                   cohort preparation and longitudinal models
scripts/Longitudinal_robustness/    lagged and within-between models
scripts/Mortality/                 dynamic mortality dataset and models
scripts/Mortality_prediction/      five-year models and recalibration
scripts/ADL/                       first-recorded basic activity difficulty
scripts/Dementia/                  dementia-related competing-event models
scripts/SCRIPT_MANIFEST.csv        execution order
run_*.R                           analysis entry points
package-versions.txt               available recorded package versions
LICENSE                            code licence
```

## Software

The mortality and longitudinal runs used R 4.4.2. Recorded package versions are listed in `package-versions.txt`. Install the analysis dependencies once:

```r
install.packages(c(
  "data.table", "haven", "dplyr", "tidyr", "tidyselect", "purrr",
  "readr", "tibble", "nlme", "lavaan", "metafor",
  "fixest", "survival", "sandwich"
))
```

Scripts use `splines` and other packages distributed with R. Each run writes analysis outputs and session information under `results/`.

## Data preparation

Obtain the cohort files from their providers under the applicable access agreements. Provider links are [HRS](https://hrsdata.isr.umich.edu/), [ELSA through the UK Data Service](https://datacatalogue.ukdataservice.ac.uk/studies/study/5050), [SHARE](https://share-eric.eu/data/data-access), [CHARLS](https://charls.pku.edu.cn/en/) and [MHAS](https://www.mhasweb.org/). Harmonised files are available through the [Gateway to Global Aging Data](https://g2aging.org/). Place the files below the private `data_root` directory. Paths below are relative to that directory.

| Cohort | Study files used by the preparation scripts |
| --- | --- |
| HRS | `HRS/RAND HRS Longitudinal File 2022/randhrs1992_2022v1_STATA/randhrs1992_2022v1.dta` and `HRS/Gateway Harmonized HRS/H_HRS_d_stata/H_HRS_d.dta` |
| ELSA | `ELSA/5050stata_A06F12276F837D82374CA07FCEA89BDC3945439B39E962E7B03B39CEBABE46BC_V1/UKDA-5050-stata/stata/stata13_se/gh_elsa_h.dta` |
| SHARE | `SHARE/easySHARE_rel9-0-0_stata/easySHARE_rel9-0-0.dta` and `SHARE/GH_SHARE_g_rel9-0-0_ALL_datasets_stata/GH_SHARE_g.dta` |
| CHARLS | `CHARLS/02_harmonized/Regular Waves/H_CHARLS_D_Data/H_CHARLS_D_Data.dta` |
| MHAS | `MHAS/04_gateway_harmonized/H_MHAS_d.dta` |

The longitudinal preparation covers HRS waves 1–15, ELSA waves 1–10, SHARE release 9, CHARLS waves 1–4 and MHAS waves 1–6. The original provider filenames and variable definitions are retained in the scripts. Functional scores average the available standardised components; the at-least-two-component rule defines the harmonisation readiness flag. HRS input combines RAND with the harmonised HRS file. Participant linkage uses `hhidpn` in HRS, `idauniq` for ELSA longitudinal preparation and `mergeid` in SHARE. ELSA self-care and dementia linkage use `idauniqc`. CHARLS uses `ID`. MHAS uses `rahhidnp`, falling back to `unhhidnp` when required by the release.

### Quality-controlled event inputs

The mortality and observation analyses also require the study's adjudicated event records. The following four files are research inputs from the frozen `V17D_frozen_inputs/v1/` release:

| File | Required fields |
| --- | --- |
| `V17D_exact_visit_dates_QC_v1.csv` | `cohort`, `id_norm`, `wave`, `interview_year_raw`, `interview_month_raw` |
| `V17D_wave_direct_proxy_status_QC_v1.csv` | `cohort`, `id_norm`, `wave`, `proxy_status` |
| `V17D_mortality_endpoint_main_QC_v1.rds` | `cohort`, `subject_id`, `baseline_date_direct`, `baseline_age_direct`, `primary_exit_date`, `primary_event`, `primary_eligible`, `sex_baseline`, `education_baseline`, `death_date_midmonth` |
| `V17D_mortality_endpoint_CHARLS_QC_v1.rds` | the same mortality fields for the repaired CHARLS endpoint records |

Dynamic mortality assembly and the ADL observation audit use the main mortality file for HRS, ELSA, SHARE and MHAS, and replace CHARLS records with the repaired CHARLS file. In the 4,771 CHARLS participants eligible for the ADL observation audit, the two endpoint releases gave identical death dates and observation classifications. The main mortality file comes from `V16_00d_direct_interview_cox_candidate_data.rds`; the repaired CHARLS file comes from `V16_00g_CHARLS_cox_candidate_no_model.rds`. Interview dates come from `V16_00c_extracted_interview_dates_weights.csv`, with interview-source records from the chronology alignment stage. These controlled files are not distributed in this repository.

This code starts event processing from those adjudicated records. The earlier provider-specific mortality adjudication and CHARLS ID repair are not regenerated here. The four QC inputs must be supplied to reproduce this endpoint release. All four local files match the original adjudication outputs and their recorded MD5 checksums. Cohort harmonisation, dynamic exposure construction and all subsequent models are included.

`proxy_status` records `direct`, `proxy_or_exit` or `unknown`. Direct interviews supply repeated exposure measurements. A valid interview year with a missing month is assigned July for dynamic mortality assembly. Observed year-month dates use day 15. The ADL schedule uses records with an observed valid month. Longitudinal keys are participant plus wave; model intervals additionally use their start date.

### ADL inputs

`scripts/ADL/00_prepare_adl_inputs.R` reads the provider's five self-care items and creates `selfcare5_overlay.rds`. Item definitions come from `config/selfcare_items.csv`. Codes 0 and 1 represent no difficulty and difficulty. Other values remain missing under the frozen item rules. A complete assessment contains five observed items.

The same script derives `cohort_cutpoints.csv` from the frozen clinical-history cutpoint population. It selects the first eligible mortality landmark at age 50 or older and retains `start_age <= 64` before calculating the current-function quartiles. The next-wave ADL eligibility window uses `50 <= age < 65`. These distinct age conditions are retained from the original analyses.

Current and next-wave self-care observations require direct interviews. Every dated earlier self-care assessment enters the first-recorded eligibility rule. The scheduled interval is 12–36 months, with SHARE index waves 6 and 7 excluded. The first qualifying index assessment is selected before checking the full earlier self-care record.

### Dementia inputs

The dementia stage builds its provider-ID crosswalk from the dynamic candidate and the harmonised ELSA/SHARE files. It reuses the cohort cognition files generated by `run_track_A.R`.

ELSA uses `r[wave]demene` in waves 1–10. SHARE uses `r[wave]alzdeme` in waves 2 and 4–9. Both sources supply interview status and assessment dates. Dementia status uses explicit 0/1 values. SHARE proxy values 1 and 2 count as proxy interviews; ELSA uses value 1.

## Configuration and execution

Start R in the repository root. Create the local configuration:

```r
file.copy("config/config.example.R", "config/config.R")
```

The default input directory is `data/` under the repository root. To use an existing data directory, set `DCV_BASE_DIR` before running, or edit `data_root` in `config/config.R`. Outputs are written under the repository's `results/` directory. The configuration file and participant-level outputs are excluded by `.gitignore`.

Run the complete clinical workflow:

```r
source("run_all.R", encoding = "UTF-8")
```

From a terminal in the repository root:

```sh
Rscript --vanilla run_all.R
```

Each stage can also run independently after its inputs have been created:

```r
source("run_track_A.R", encoding = "UTF-8")
source("run_longitudinal_robustness.R", encoding = "UTF-8")
source("run_mortality_bridge.R", encoding = "UTF-8")
source("run_mortality_prediction.R", encoding = "UTF-8")
source("run_clinical_translation.R", encoding = "UTF-8")
source("run_adl.R", encoding = "UTF-8")
source("run_dementia_bridge.R", encoding = "UTF-8")
```

`run_track_A.R` creates the longitudinal inputs. Both the robustness and mortality stages use those files. Mortality creates the dynamic candidate required by the prediction, ADL and dementia stages. Run `run_mortality_prediction.R` before `run_clinical_translation.R`. The ADL and dementia entry points include their input-preparation scripts.

The longitudinal and robustness analyses completed on 26 November 2026, with all four branch sample checks passing. The ADL assembly produced 38,254 observed follow-ups and 1,694 first-recorded difficulty events.

Stage entry points load the configuration themselves. Each analysis script runs in a separate R environment. Existing output directories can be reused when rerunning a stage.

## Model settings

Annual-change and reciprocal meta-analyses use REML with Knapp–Hartung intervals. `Track_A/05` writes the overall and cohort-group summaries to `meta_summary_knapp_hartung_v1.csv`. Leave-one-cohort-out and between-within summaries use normal-approximation intervals.

Dynamic mortality uses attained age as the time scale, subject-cluster robust variance and Efron ties. All four affective/functional history and deviation components enter together. Age bands are below 65, 65–74, 75–84 and at least 85 years. HRS provides the four primary global tests with Holm correction. MHAS and SHARE contribute replication estimates; ELSA provides sensitivity estimates. CHARLS contributes event support and five-year evaluation rather than the age-varying mortality model. Three-cohort mortality pooling uses REML with Hartung–Knapp intervals.

The mortality sensitivity models use 12-month delayed entry, one-visit lagging and a restriction to at least three measurements. The specified current-deviation interaction remains in the supplementary analyses.

Five-year mortality prediction compares demographic covariates with current affective/functional scores, then compares current scores with the four history/deviation components. HRS provides a participant-level 70% development split using seed `20260828`. Coefficients and cumulative baseline hazards are applied unchanged during evaluation. Paired participant bootstrap differences use 1,000 draws. The age-restricted continuous analysis retains seed `20260909`; its original age condition is `50 <= start_age <= 64`. Intercept-only recalibration updates the five-year baseline risk within each evaluation cohort while retaining predictor coefficients.

The ADL association model includes current function and earlier function, with affective history/deviation and the specified demographic/time covariates. Standardised risks compare the 75th and 25th history percentiles. Risk-difference intervals use the delta method. HRS model development uses seed `20261110`; held-out and external evaluation uses 300 paired bootstrap draws with seed `20261111`. The observation analysis includes everyone meeting first-recorded eligibility and separates observed reassessment from pre-reassessment death and unresolved status.

Dementia-related cause-specific Cox models compare current scores with history/deviation components. ELSA and SHARE are fitted separately; SHARE models are stratified by country. The entry assessment must explicitly record no dementia. The first subsequent positive report defines detection. Death takes priority when death and detection share a date. The primary analysis retains known deaths within valid dynamic intervals, including deaths after the last dementia assessment. The strict-observation sensitivity additionally caps follow-up at the last valid outcome assessment. Cognition-adjusted and interview-source sensitivities remain included.



## Licence

The existing repository licence is retained in `LICENSE`. Access and redistribution of cohort data follow the data providers' terms.
