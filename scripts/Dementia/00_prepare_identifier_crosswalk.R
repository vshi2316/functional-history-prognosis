# Link the dynamic analysis identifiers to the provider identifiers.
suppressPackageStartupMessages({library(data.table); library(haven)})
normalise_id <- function(x) {
  x <- toupper(gsub("\\s+", "", trimws(as.character(x))))
  x <- sub("\\.0+$", "", x)
  numeric <- !is.na(x) & grepl("^[0-9]+$", x)
  x[numeric] <- sub("^0+", "", x[numeric])
  x
}
candidate <- as.data.table(readRDS(Sys.getenv("DEMENTIA_DYNAMIC_CANDIDATE_INPUT")))
candidate <- unique(candidate[cohort %in% c("ELSA", "SHARE"), .(cohort, candidate_id = as.character(subject_id))])
candidate[, id_key := normalise_id(sub("^[^:]+::", "", candidate_id))]
files <- c(ELSA = Sys.getenv("DEMENTIA_ELSA_OUTCOME_INPUT"),
           SHARE = Sys.getenv("DEMENTIA_SHARE_OUTCOME_INPUT"))
ids <- c(ELSA = "idauniqc", SHARE = "mergeid")
provider <- rbindlist(lapply(names(files), function(z) {
  raw <- read_dta(files[[z]], col_select = tidyselect::all_of(ids[[z]]))
  data.table(cohort = z, provider_id = as.character(raw[[ids[[z]]]]))
}))
provider[, id_key := normalise_id(provider_id)]
if (anyDuplicated(provider[, .(cohort, id_key)])) stop("Provider identifiers are not unique.")
crosswalk <- merge(candidate, provider, by = c("cohort", "id_key"), all.x = TRUE)
if (anyNA(crosswalk$provider_id)) stop("An analysis identifier has no provider match.")
path <- Sys.getenv("DEMENTIA_ID_CROSSWALK_INPUT")
dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
saveRDS(as.data.frame(crosswalk[, .(cohort, candidate_id, provider_id)]), path)
