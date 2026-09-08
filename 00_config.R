# =============================================================================
# Config for NULISA proteomics pipeline (v2)
# Date: July 2026
# Description: Paths & parameters for NULISA GP2 Pilot proteomics workflow
# =============================================================================



# --- VM Paths -----------------------------------------------------------------
BUCKET_DIR   <- "~/proteomics_bucket/proteomics"
DATA_DIR     <- file.path(BUCKET_DIR, "metadata")
MANIFEST_DIR <- file.path(DATA_DIR, "GP2_nulisa_manifests")
PROJECT_DIR  <- file.path(BUCKET_DIR, "results", "nulisa_v2")



# --- Cohort to analyze --------------------------------------------------------
# change cohort and plate correction, source 00, then run scripts. 
COHORT <- "ntuh" 
PLATE_CORRECTION <- "none" 

# label whether model was adj for plate or not 
PLATE_SUFFIX <- if (PLATE_CORRECTION == "none") "_plate_none" else "_plate_adj"

# "default" = the NA/zero thresholds below; "none" = no NA or zero filtering.
FILTER_MODE <- "none"


if (COHORT == "ntuh") {

  PANELS             <- c("Neuro220")
  MATRIX_TYPE        <- "PLASMA"
  SAMPLE_DET_MIN     <- 0.90            
  ANCESTRY_KEEP      <- "EAS"
  VISIT_KEEP         <- NULL            # no visit column
  ID_METHOD          <- "strip_suffix"  # SampleName minus "_s2" = GP2ID
  ID_COLUMN          <- NULL
  GP2_JOIN_KEY       <- "GP2ID"
  STUDY_FILTER       <- c("NTUH")
  EXCLUDE_FROM_DA <- c("APOE4")   # from TAP report: binary carrier readout
  USE_CURATED        <- FALSE
  SAMPLE_ID_PATTERN  <- "^NTUH_[0-9]+_s[0-9]+$"
  MANIFEST_FILES     <- NULL
  MANIFEST_SAMPLE_COL <- NULL
  MANIFEST_GP2ID_COL  <- NULL
  NPQ_DIR  <- file.path(BUCKET_DIR, "P121 GP2")
  NPQ_FILE <- "P121_BSHRI_NULISAseq_Neuro220_NPQ_06092026.xlsx"

} else if (COHORT == "p118") {
  PANELS             <- c("Neuro220")
  MATRIX_TYPE        <- "PLASMA"
  SAMPLE_DET_MIN     <- 0.90
  ANCESTRY_KEEP      <- "EAS"
  VISIT_KEEP         <- NULL
  ID_METHOD          <- "manifest"
  ID_COLUMN          <- NULL
  GP2_JOIN_KEY       <- "GP2ID"
  STUDY_FILTER       <- c("UMKLM", "KUL")
  EXCLUDE_FROM_DA    <- c("APOE4")
  USE_CURATED        <- FALSE
  SAMPLE_ID_PATTERN  <- NULL
  NPQ_DIR            <- file.path(BUCKET_DIR, "P118 GP2")
  NPQ_FILE           <- "P118_BSHRI_NULISAseq_Neuro220_NPQ_06082026.xlsx"
  MANIFEST_FILES     <- c("UMKLM_selfQCV2_2026-01-23_m11.csv",
                          "KUL_selfQCV2_2026-01-23_m4.csv")
  MANIFEST_SAMPLE_COL <- "sample_id"
  MANIFEST_GP2ID_COL  <- "GP2ID"

} else  stop("Unknown COHORT: ", COHORT)

# Results directories 
COHORT_DIR <- switch(COHORT,
  ntuh     = "results_p121_v2",
  p118     = "results_p118_v2")

RESULTS_BASE <- file.path(PROJECT_DIR, paste0(COHORT_DIR, "_filter_", FILTER_MODE))
RESULTS_DIR  <- paste0(RESULTS_BASE, PLATE_SUFFIX)

RESULTS_01_DIR <- file.path(RESULTS_BASE, "01_load_filter")
RESULTS_02_DIR <- file.path(RESULTS_BASE, "02_merge_gp2")
RESULTS_03_DIR <- file.path(RESULTS_DIR, "03_differential_abundance")
RESULTS_04_DIR <- file.path(RESULTS_DIR, "04_qc_report")

# SC drift script doesn't depend on the DA model, so it lives in RESULTS_BASE
SC_DRIFT_DIR  <- file.path(RESULTS_BASE, "02b_sc_drift")
SC_DRIFT_FILE <- file.path(SC_DRIFT_DIR, paste0("sc_drift_", COHORT, ".csv"))



# --- NULISA data files ----------------------------------------------------------
# release 12
MASTER_KEY_FILE <- "INTERNAL_USE_ONLY_master_key_release12_final_vwb.csv"

# Carrier information in variant report
VARIANT_REPORT_FILE <- "variant_report_r12_final_aug3.csv"


# --- Meta6 GWAS gene list ---------------------------------------------------
META6_FILE <- file.path(DATA_DIR, "meta6_genes.csv")



# --- Pipeline output files (downstream scripts reference these) -------------

# The merged SE for whichever cohort is active. Scripts 04, 06, 08 read this.
SE_FILTERED <- file.path(RESULTS_01_DIR, paste0("se_", COHORT, "_Neuro220_filtered.rds"))
SE_MERGED   <- file.path(RESULTS_02_DIR, paste0("se_", COHORT, "_Neuro220_with_metadata.rds"))


# --- Shared parameters ------------------------------------------------------
FDR_THRESHOLD <- 0.05
OUTLIER_SD_THRESHOLD <- 4



# --- v2 target filters ------------------------------------------------------
ZERO_PROP_MAX <- if (FILTER_MODE == "none") 1 else 0.55
NA_PROP_MAX   <- if (FILTER_MODE == "none") 1 else 0.70

# Does a zero count toward the NA total?
#   FALSE = zeros and NAs are separate states (our reading of the NPQ formula)
ZEROS_COUNT_AS_NA <- FALSE

# Target-level detectability (>=50%), matrix-independent, per Alamar.
# v2 CHANGE: this is now a flag carried into the DA results, not a filter.
DET_THRESH <- 0.50
APPLY_DET_FILTER <- FALSE

# Samples below this detectability (%) ARE excluded from DA and QC. 
# Note Alamar's floor is stricter and matrix-dependent (90% plasma, 70% CSF)
SAMPLE_DETECT_EXCLUDE <- 20

# Minimum non-missing samples before a target is modelled. PLACEHOLDER.
MIN_SAMPLES_PER_TARGET <- 20

# --- v2 sample handling -----------------------------------------------------
# v2 CHANGE: SampleQC == "warning" samples are KEPT and flagged
DROP_SAMPLE_QC_WARN <- FALSE
REMOVE_OUTLIERS <- FALSE

# Some donors contributed more than one sample. Until we know what the visit
# labels mean (0W / 52W) and how to handle donors enrolled under multiple
# sub-study IDs, these donors are excluded rather than picked between.
DROP_REPEAT_DONORS <- TRUE


# --- v2 DA model ------------------------------------------------------------
# Covariates: age and sex only, per Guglielmo. PDTRTMNT is deliberately absent.
DA_COVARIATES <- c("age", "sex_clean")



# --- Special target categories (per NULISAseqR user guide) -------------------
EXEMPT_RARE_CASE <- TRUE

HIGH_ABUNDANCE_TARGETS <- list(
  Inflammation = c("CRP", "KNG1"),
  CNS_Disease  = c("APOE", "CRP"),
  Neuro220     = c("APOA1", "APOA2", "APOE", "APOH", "C1q",
                   "CRP", "PPBP", "PROS1", "SERPINA3")
)

RARE_CASE_TARGETS <- list(
  Inflammation = character(0),
  CNS_Disease  = character(0),
  Neuro220     = c("APOE4", "mHTT-exon1", "pLRRK2-S1292",
                   "pPRKN-S65", "pQ-ATXN3", "pQ-HTT", "pRAB10-T73",
                   "pRAB12-S106", "pRAB29-T71")
)

# --- Vendor QC flags (Alamar Neuro 220 data sheet, Table 3) ------------------
# All of these are flags will be carried into DA results. Not used as filters.

# Cross-reactivity > 1%. Named vector = the reported XR percentage.
CROSS_REACTIVE_TARGETS <- list(
  Neuro220 = c("CCL4" = 1438.1, "CCL3" = 294.3, "Oligo-SNCA" = 99.3,
               "PGK1" = 22.7, "pPRKN-S65" = 10.1, "VSNL1" = 9.9,
               "pQ-ATXN3" = 5.1, "pQ-HTT" = 1.8),
  CNS_Disease  = character(0),
  Inflammation = character(0)
)

# Subset whose cross-reactive partner is ALSO on the panel.
XR_PARTNER_ON_PANEL <- list(
  Neuro220 = c("Oligo-SNCA", "pPRKN-S65", "pQ-ATXN3", "pQ-HTT")
)

# Vendor plasma detectability < 50% in NORMAL controls, excluding rare-case.
VENDOR_LOW_DETECT <- list(
  Neuro220 = c("YWHAZ", "PTN", "SNCB", "SNAP25", "Aβ38", "CGRP",
               "FCGR2B", "UCHL1", "GSDMD", "PARP1", "UBB")
)

# Vendor plasma CV > 20% (intra or inter) per Alamar Data Sheet (Neuro220).
VENDOR_HIGH_CV_PLASMA <- list(
  Neuro220 = c("pPRKN-S65", "pRAB29-T71", "pQ-HTT", "pRAB12-S106",
               "BASP1", "CGRP", "SNCB", "CST5", "pLRRK2-S1292",
               "pRAB10-T73", "pQ-ATXN3", "PARK7")
)


# guards
stopifnot(
  exists("PANELS"), exists("ANCESTRY_KEEP"), exists("MATRIX_TYPE"),
  exists("SAMPLE_ID_PATTERN"), exists("VISIT_KEEP"), exists("ID_METHOD"),
  exists("NPQ_DIR"), exists("NPQ_FILE"),
  ID_METHOD %in% c("column", "strip_suffix", "manifest"),
  (ID_METHOD == "manifest") == !is.null(MANIFEST_FILES),
  is.null(SAMPLE_ID_PATTERN) || nzchar(SAMPLE_ID_PATTERN),
  length(SAMPLE_DET_MIN) == 1, SAMPLE_DET_MIN > 0, SAMPLE_DET_MIN <= 1,
  FILTER_MODE %in% c("default", "none"),
  PLATE_CORRECTION %in% c("none", "covariate", "sva")
)

cat("CONFIG | cohort:", COHORT, "| matrix:", MATRIX_TYPE,
    "| panels:", paste(PANELS, collapse = ", "),
    "| ancestry:", ANCESTRY_KEEP,
    "| plate adj:", PLATE_CORRECTION,
    "| id pattern:", ifelse(is.null(SAMPLE_ID_PATTERN), "none", SAMPLE_ID_PATTERN), "\n")
