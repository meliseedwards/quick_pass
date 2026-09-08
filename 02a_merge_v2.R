# =============================================================================
# Script 02 (v2): Merge GP2 Master Key (R12) with NULISA Data
# Date: July 2026
# Description: Loads SE objects from Script 01, joins the GP2 master key
#              using cohort donor ID, adds clinical metadata (diagnosis,
#              age, sex, ancestry) to colData, harmonizes phenotype to
#              PD/Control, and saves updated SE objects.
# =============================================================================



# --- 0. Setup ----------------------------------------------------------------

library(SummarizedExperiment)
library(tidyverse)
library(readxl)

source("~/proteomics_nulisa/scripts/nulisa_pipeline_v2/quick_pass/00_config.R")

dir.create(RESULTS_02_DIR, showWarnings = FALSE, recursive = TRUE)


# --- 1. Load master key ------------------------------------------------------
master_key <- read.csv(file.path(DATA_DIR, MASTER_KEY_FILE), stringsAsFactors = FALSE)

# Filter to study only
study_key <- master_key %>%
  filter(study %in% STUDY_FILTER) %>%
  mutate(clinical_id = as.character(clinical_id))
cat("Study samples in master key:", nrow(study_key), "\n")
stopifnot(nrow(study_key) > 0)

# Select only columns we need for downstream analysis
study_key_clean <- study_key %>%
  select(clinical_id, GP2ID, GP2sampleID, study,
         diagnosis, GP2_PHENO, GP2_phenotype, study_arm, study_type,
         biological_sex_for_qc, age, age_of_onset, age_at_diagnosis,
         age_at_death, age_at_last_follow_up, GP2_phenotype_for_qc,
         family_history_pd, nba_label, nba, wgs_label,
         race_for_qc, biological_sex_for_qc, region_for_qc)%>%
  mutate(join_key = .data[[GP2_JOIN_KEY]])

# Check for duplicate IDs
duplicated_ids <- study_key_clean %>%
  dplyr::count(join_key) %>%
  dplyr::filter(n > 1)
cat("Duplicate join keys in join key:", nrow(duplicated_ids), "\n")

if (nrow(duplicated_ids) > 0) {
  print(head(duplicated_ids))
  stop("Duplicate id(s) in the master key - the join would add rows.")
}


# --- 1b. Load PD treatment status (baseline from PPMI data) ----------------------

# PDTRTMNT (binary: 0=untreated, 1=treated) at baseline visit, joined by ID (PPMI).
if (USE_CURATED) {
  curated_sheets <- excel_sheets(file.path(DATA_DIR, CURATED_FILE))

if (!CURATED_SHEET %in% curated_sheets) {
  stop("Sheet '", CURATED_SHEET, "' is not in ", CURATED_FILE,
       ". Sheets present: ", paste(curated_sheets, collapse = ", "),
       "\nUpdate CURATED_SHEET in the config - it changes with every data cut.")
}

treatment_raw <- read_excel(file.path(DATA_DIR, CURATED_FILE),
                            sheet = CURATED_SHEET, na = c("", "NA")) %>%
  filter(EVENT_ID == "BL") %>%
  mutate(PATNO = as.character(PATNO)) %>%
  select(PATNO, PDTRTMNT)

treatment_bl <- treatment_raw %>% distinct(PATNO, .keep_all = TRUE)

cat("Baseline treatment records loaded:", nrow(treatment_bl),
    "| treated:", sum(treatment_bl$PDTRTMNT == 1, na.rm = TRUE), "\n")
} else {
  treatment_bl <- NULL
  cat("USE_CURATED = FALSE - no curated clinical file for this cohort\n")
}



# --- 2. Helper function: merge metadata into SE -----------------------------
merge_metadata <- function(se_path, panel_name) {
 
  cat("\n=====", panel_name, "=====\n")

  # load SE object from Script 01
  se <- readRDS(se_path)

  # extract current colData; DONOR_ID is the join key
  current_meta <- as.data.frame(colData(se)) %>%
    mutate(DONOR_ID = as.character(DONOR_ID))

  # join with master key 
   merged_meta <- current_meta %>%
    left_join(study_key_clean, by = c("DONOR_ID" = "join_key"))

  if (!is.null(treatment_bl)) {
    merged_meta <- merged_meta %>%
      left_join(treatment_bl, by = c("DONOR_ID" = "PATNO")) %>%
      mutate(PDTRTMNT = factor(PDTRTMNT))
  }

  # v2: neither join may change the number of rows
  stopifnot(nrow(merged_meta) == ncol(se))

  # Report match success
  matched   <- sum(!is.na(merged_meta$GP2ID))
  unmatched <- nrow(merged_meta) - matched
  match_pct <- round(100 * matched / nrow(merged_meta), 1)

  cat("Samples matched to GP2:", matched, "/", nrow(merged_meta),
      "(", match_pct, "%)\n", sep = "")
  cat("Unmatched samples:", unmatched, "\n")
  cat("\nDiagnosis distribution (matched samples):\n")
  print(table(merged_meta$diagnosis, useNA = "always"))
  cat("\nGP2_PHENO distribution:\n")
  print(table(merged_meta$GP2_PHENO, useNA = "always"))
  if (sum(!is.na(merged_meta$age)) > 0) {
    cat("\nAge distribution:\n")
    print(summary(merged_meta$age))
  }


  # --- Harmonize phenotype and sex labels --------------------------------
  # GP2 master key has inconsistent labels. Create clean variables for
  # downstream analysis; original columns are preserved.
  merged_meta <- merged_meta %>%
    mutate(
      sex_clean = case_when(
        trimws(biological_sex_for_qc) %in% c("F", "Female", "1") ~ "Female",
        trimws(biological_sex_for_qc) %in% c("M", "Male", "2")   ~ "Male",
        TRUE ~ NA_character_),
      # PD vs Control from standardized GP2 fields (prodromal/other dropped)
      phenotype_clean = case_when(
        trimws(GP2_phenotype) == "PD" &
          trimws(diagnosis) %in% c("PD", "Idiopathic PD")
        ~ "PD",
        trimws(GP2_phenotype) == "Control" &
          trimws(diagnosis) %in% c("HC", "Control", "Healthy control",
                                   "No PD Nor Other Neurological Disorder")
        ~ "Control",
        TRUE ~ NA_character_
      ),
      # Unified ancestry: prefer NBA label; fall back to WGS label; else NA.
      ancestry = case_when(
      !is.na(nba_label) & nba_label != "" ~ nba_label,
      !is.na(wgs_label) & wgs_label != "" ~ wgs_label,
      TRUE                                ~ NA_character_
    ),
      ancestry_source = case_when(
      !is.na(nba_label) & nba_label != "" ~ "NBA",
      !is.na(wgs_label) & wgs_label != "" ~ "WGS",
      TRUE                                ~ NA_character_
    )
  )
  # Report exclusions (matched to GP2 but no clean phenotype)
  # v2: show WHICH label combinations caused it. A "Parkinson's Disease" vs
  # "Idiopathic PD" mismatch in a new cohort otherwise costs samples silently.
  excluded <- merged_meta %>%
    filter(is.na(phenotype_clean) & !is.na(GP2ID))
  cat("Excluded (matched but no clean PD/Control phenotype):", nrow(excluded), "\n")

  if (nrow(excluded) > 0) {
    print(excluded %>% dplyr::count(GP2_phenotype, diagnosis, sort = TRUE))
  }

  # A GP2_phenotype of PD or Control that did not map is a string mismatch,
  # not a real exclusion. Check.
  unmapped <- merged_meta %>%
    filter(is.na(phenotype_clean), trimws(GP2_phenotype) %in% c("PD", "Control"))
  if (nrow(unmapped) > 0) {
    cat("\nWARNING:", nrow(unmapped),
        "samples have GP2_phenotype PD or Control but did not map.\n",
        " This is almost always an unrecognised diagnosis string:\n")
    print(unmapped %>% dplyr::count(GP2_phenotype, diagnosis, sort = TRUE))
  }

  cat("\nHarmonized phenotype distribution:\n")
  print(table(merged_meta$phenotype_clean, useNA = "always"))
  cat("\nAncestry (unified: NBA then WGS) distribution:\n")
  print(table(merged_meta$ancestry, useNA = "always"))
  cat("\nAncestry source (NBA vs WGS):\n")
  print(table(merged_meta$ancestry_source, useNA = "always"))

  # check the configured ancestry actually exists here
  n_ancestry <- sum(merged_meta$ancestry %in% ANCESTRY_KEEP, na.rm = TRUE)
  cat("Samples matching ANCESTRY_KEEP (", ANCESTRY_KEEP, "): ", n_ancestry, "\n", sep = "")
  if (n_ancestry == 0)
    stop("ANCESTRY_KEEP is '", ANCESTRY_KEEP, "' but no samples have it.")

  # treatment x phenotype, not the marginal table
  if ("PDTRTMNT" %in% colnames(merged_meta)) {
    cat("\nTreatment (PDTRTMNT) x phenotype - DIAGNOSTIC ONLY, not a covariate:\n")
    print(table(merged_meta$PDTRTMNT, merged_meta$phenotype_clean, useNA = "always"))
  }

  # Critical: merged_meta MUST be in same order as colnames(se)
  stopifnot(identical(as.character(merged_meta$SampleName), colnames(se)))
  rownames(merged_meta) <- merged_meta$SampleName
  colData(se) <- DataFrame(merged_meta)

  # Script 01 put the per-target QC into rowData. Confirm it survived, so
  # Script 03 can join the flags from the SE instead of re-reading a CSV.
  stopifnot(all(c("detectability", "low_detect_flag", "zero_prop") %in%
                  colnames(rowData(se))))



  # save updated SE
out_path <- file.path(RESULTS_02_DIR,
                        paste0("se_", COHORT, "_", panel_name, "_with_metadata.rds"))
saveRDS(se, out_path)
  cat("Saved:", basename(out_path), "\n")


# write the unmatched donors so they can be tracked (possibly R13)
  unmatched_ids <- merged_meta %>%
    filter(is.na(GP2ID)) %>%
    select(SampleName, DONOR_ID, PlateID, SampleQC, pct_above_lod)
  if (nrow(unmatched_ids) > 0) {
    write.csv(unmatched_ids,
              file.path(RESULTS_02_DIR,
                        paste0("unmatched_", COHORT, "_", panel_name, ".csv")),
              row.names = FALSE)
    cat("Wrote", nrow(unmatched_ids), "unmatched donor IDs to unmatched_",
        COHORT, "_", panel_name, ".csv\n", sep = "")
  }


  # Return summary stats
  return(list(
    panel      = panel_name,
    n_samples  = nrow(merged_meta),
    n_matched  = matched,
    match_pct  = match_pct,
    n_PD       = sum(merged_meta$phenotype_clean == "PD", na.rm = TRUE),
    n_Control  = sum(merged_meta$phenotype_clean == "Control", na.rm = TRUE),
    n_unmatched = sum(is.na(merged_meta$GP2ID)),
    n_ancestry_keep = n_ancestry,
    # v2: how many samples are actually usable in Script 03 - PD/Control, the
    # right ancestry, and complete covariates.
    n_da_ready = sum(
      merged_meta$phenotype_clean %in% c("PD", "Control") &
      merged_meta$ancestry %in% ANCESTRY_KEEP &
      stats::complete.cases(merged_meta[, DA_COVARIATES, drop = FALSE]),
      na.rm = TRUE)
  ))
}


# --- 3. Merge metadata for each panel ----------------------------------------
SE_IN <- list(Neuro220 = SE_FILTERED)

summaries <- list()
for (p in PANELS) summaries[[p]] <- merge_metadata(SE_IN[[p]], p)



# --- 4. Save merge summary ---------------------------------------------------
summary_df <- bind_rows(summaries)

write.csv(summary_df,
          file.path(RESULTS_02_DIR, paste0(COHORT, "_merge_summary.csv")),
          row.names = FALSE)

print(summary_df)
