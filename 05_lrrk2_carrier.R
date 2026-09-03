# =================================================================================
# Script 05 v2: carrier-stratified analyses
# Date: August 2026
# Description: Asseses whether risk or pathogenic variant alters protein levels
#              Primary: any carrier (e.g., LRRK2: G2385R or R1628P) vs non-carrier
#              Secondary: each variant separately vs non-carrier
# =================================================================================



# --- 0. Setup ---------------------------------------------------------------------

library(SummarizedExperiment)
library(tidyverse)
library(limma)
library(ggrepel)

source("~/proteomics/nulisa_pipeline_v2/quick_pass/00_config.R")
stopifnot(COHORT %in% c("ntuh", "umklm", "kul"))

RESULTS_06_DIR <- file.path(RESULTS_DIR, "06_lrrk2_carrier")
dir.create(RESULTS_06_DIR, showWarnings = FALSE, recursive = TRUE)

RISK_VARIANTS <- c("LRRK2_Gly2385Arg", "LRRK2_Arg1628Pro")

# LRRK2 pathway proteins on the Neuro220 panel - looked at explicitly at the end
LRRK2_PATHWAY <- c("LRRK2", "pLRRK2-S1292",
                   "RAB10", "RAB12", "RAB29",
                   "pRAB10-T73", "pRAB12-S106", "pRAB29-T71")




# --- 1. Base cohort (same filters as Script 03) ------------------------------

se <- readRDS(SE_MERGED_ACTIVE)
meta <- as.data.frame(colData(se))

keep <- meta$phenotype_clean %in% c("PD", "Control") & meta$ancestry %in% ANCESTRY_KEEP
keep[is.na(keep)] <- FALSE
se <- se[, keep]
meta <- as.data.frame(colData(se))

if (!is.null(SAMPLE_DETECT_EXCLUDE)) {
  bad <- !is.na(meta$pct_above_lod) & meta$pct_above_lod < SAMPLE_DETECT_EXCLUDE
  if (any(bad)) { se <- se[, !bad]; meta <- as.data.frame(colData(se)) }
}

cc <- complete.cases(meta[, c("phenotype_clean", DA_COVARIATES), drop = FALSE])
se <- se[, cc]
expr <- assay(se, "npq"); se <- se[rowSums(is.na(expr)) == 0, ]
if (length(EXCLUDE_FROM_DA) > 0) se <- se[!rownames(se) %in% EXCLUDE_FROM_DA, ]

meta <- as.data.frame(colData(se))
cat("Base cohort:", ncol(se), "samples |", nrow(se), "targets\n")



# --- 2. Carrier status -------------------------------------------------------
variants <- read.csv(file.path(DATA_DIR, VARIANT_REPORT_FILE),
                     stringsAsFactors = FALSE) %>%
  dplyr::filter(Gene == "LRRK2", GP2ID %in% meta$DONOR_ID)

cat("\nLRRK2 variants in this cohort:\n")
print(variants %>% dplyr::count(Variant_Name, Pathogenicity, Zygosity))

risk_ids  <- variants %>% filter(Variant_Name %in% RISK_VARIANTS) %>% pull(GP2ID) %>% unique()
other_ids <- variants %>% filter(!Variant_Name %in% RISK_VARIANTS) %>% pull(GP2ID) %>% unique()
# someone carrying a non-risk LRRK2 variant AND a risk variant stays a carrier;
# only those carrying non-risk variants exclusively are dropped
other_only <- setdiff(other_ids, risk_ids)

# who was genotyped at all? non-carrier must mean genotyped-and-negative
stopifnot(identical(as.character(meta$DONOR_ID), as.character(meta$GP2ID)))
cohort_ids <- meta$GP2ID

mk <- read.delim(file.path(DATA_DIR, MASTER_KEY_FILE), stringsAsFactors = FALSE)
genotyped   <- mk$GP2ID[(mk$nba %in% 1 | mk$wgs %in% 1) & mk$GP2ID %in% cohort_ids]
ungenotyped <- setdiff(cohort_ids, genotyped)

cat("\nRisk-variant carriers:", length(risk_ids),
    "| non-risk LRRK2 only (excluded):", length(other_only),
    "| ungenotyped (excluded):", length(ungenotyped), "\n")
if (length(other_only)) cat("  excluded, non-risk variant:", paste(other_only, collapse = ", "), "\n")
if (length(ungenotyped)) cat("  excluded, no NBA or WGS:", paste(ungenotyped, collapse = ", "), "\n")

meta <- meta %>%
  mutate(
    lrrk2_risk  = ifelse(GP2ID %in% risk_ids, "Carrier", "Non-carrier"),
    g2385r      = GP2ID %in% (variants %>% filter(Variant_Name == "LRRK2_Gly2385Arg") %>% pull(GP2ID)),
    r1628p      = GP2ID %in% (variants %>% filter(Variant_Name == "LRRK2_Arg1628Pro") %>% pull(GP2ID)),
    analysable  = !GP2ID %in% c(other_only, ungenotyped)
  )
colData(se) <- DataFrame(meta)

se_c  <- se[, meta$analysable]
metac <- as.data.frame(colData(se_c))
metac$lrrk2_risk <- factor(metac$lrrk2_risk, levels = c("Non-carrier", "Carrier"))
colData(se_c)$lrrk2_risk <- metac$lrrk2_risk

cat("\nAnalysis cohort:", ncol(se_c), "samples\n")
cat("Carrier x phenotype:\n"); print(table(metac$lrrk2_risk, metac$phenotype_clean))
cat("\nCarrier x plate  - if carriers sit on every plate the design is NOT confounded:\n")
print(table(metac$lrrk2_risk, metac$PlateID))
bal <- suppressWarnings(chisq.test(table(metac$lrrk2_risk, metac$PlateID)))
cat("carrier x plate chi-square p =", signif(bal$p.value, 3),
    if (bal$p.value < 0.05) " -> UNBALANCED, interpret with the same caution as PD/Control\n"
    else " -> balanced; carriers are spread across plates\n")
cat("Carriers among PD:", sum(metac$lrrk2_risk == "Carrier" & metac$phenotype_clean == "PD"),
    "| among Control:", sum(metac$lrrk2_risk == "Carrier" & metac$phenotype_clean == "Control"), "\n")



# --- 3. Fit one contrast -----------------------------------------------------
# Phenotype is always adjusted for: these are PD risk alleles, so carriers are
# enriched among cases and the two would otherwise be confounded.
CARRIER_PLATE_ADJUST <- (PLATE_CORRECTION == "covariate")

model_str <- paste("~ carrier + phenotype +", paste(DA_COVARIATES, collapse = " + "),
                   if (CARRIER_PLATE_ADJUST || PLATE_CORRECTION == "covariate") "+ PlateID" else "")
suffix <- if (CARRIER_PLATE_ADJUST) "plateAdjusted" else "plateUnadjusted"

fit_contrast <- function(se_in, group_col, label) {
  m <- as.data.frame(colData(se_in))
  m[[group_col]] <- droplevels(factor(m[[group_col]]))
  if (nlevels(m[[group_col]]) < 2) { cat("SKIP", label, "- only one group\n"); return(NULL) }

  form <- paste("~", group_col, "+ phenotype_clean +", paste(DA_COVARIATES, collapse = " + "))
if (CARRIER_PLATE_ADJUST || PLATE_CORRECTION == "covariate")
    form <- paste(form, "+ PlateID")
  design <- model.matrix(as.formula(form), data = m)
  stopifnot(nrow(design) == ncol(se_in))

  coef_name <- grep(paste0("^", group_col), colnames(design), value = TRUE)[1]
  fit <- eBayes(lmFit(assay(se_in, "npq"), design))
  res <- topTable(fit, coef = coef_name, number = Inf, adjust.method = "BH") %>%
    rownames_to_column("Target") %>%
    mutate(contrast = label,
           n_total = ncol(se_in),
           n_group = sum(m[[group_col]] == levels(m[[group_col]])[2]))

  cat("\n---", label, "---\n")
  cat("n =", ncol(se_in), "|", levels(m[[group_col]])[2], ":",
      sum(m[[group_col]] == levels(m[[group_col]])[2]), "\n")
  cat("Significant at FDR <", FDR_THRESHOLD, ":", sum(res$adj.P.Val < FDR_THRESHOLD), "\n")
  print(res %>% select(Target, logFC, P.Value, adj.P.Val) %>% head(8))
  res
}



# --- 4. Primary + secondary contrasts ----------------------------------------

res_pooled <- fit_contrast(se_c, "lrrk2_risk",
                           "Any LRRK2 risk carrier vs non-carrier")

# single-variant: drop carriers of the OTHER risk variant so the comparison is clean
mk_single <- function(se_in, flag, other_flag, label) {
  m <- as.data.frame(colData(se_in))
  keep <- !m[[other_flag]]                       # remove other-variant carriers
  s <- se_in[, keep]
  mm <- as.data.frame(colData(s))
  colData(s)$variant_group <- factor(ifelse(mm[[flag]], "Carrier", "Non-carrier"),
                                     levels = c("Non-carrier", "Carrier"))
  fit_contrast(s, "variant_group", label)
}
res_g2385r <- mk_single(se_c, "g2385r", "r1628p", "G2385R carrier vs non-carrier")
res_r1628p <- mk_single(se_c, "r1628p", "g2385r", "R1628P carrier vs non-carrier")



# --- 5. LRRK2 pathway proteins, looked at explicitly -------------------------
# The prespecified biological question. Reported regardless of significance,
# with detectability so a p-value on a floor-bound target is not over-read.

det <- as.data.frame(rowData(se_c)) %>% select(Target, detectability)

pathway <- bind_rows(res_pooled, res_g2385r, res_r1628p) %>%
  filter(Target %in% LRRK2_PATHWAY) %>%
  left_join(det, by = "Target") %>%
  select(contrast, Target, logFC, P.Value, adj.P.Val, detectability) %>%
  arrange(contrast, P.Value)

cat("\n\n===== LRRK2 pathway proteins (prespecified) =====\n")
cat("Low detectability means the target sits at the assay floor in most samples;\n")
cat("a p-value there is not the same kind of evidence as one on a well-measured target.\n\n")
print(as.data.frame(pathway), row.names = FALSE)
write.csv(pathway, file.path(RESULTS_06_DIR, "lrrk2_pathway_targets.csv"), row.names = FALSE)


# --- 6. Save + volcano -------------------------------------------------------

all_res <- bind_rows(res_pooled, res_g2385r, res_r1628p) %>% left_join(det, by = "Target")
write.csv(all_res, file.path(RESULTS_06_DIR, "carrier_DA_all_contrasts.csv"), row.names = FALSE)

if (!is.null(res_pooled)) {
  d <- res_pooled %>%
    left_join(det, by = "Target") %>%
    mutate(dir = case_when(adj.P.Val < FDR_THRESHOLD & logFC > 0 ~ "Up in carriers",
                           adj.P.Val < FDR_THRESHOLD & logFC < 0 ~ "Down in carriers",
                           TRUE ~ "Not significant"))
  lab <- d %>% filter(P.Value < 0.05 | Target %in% LRRK2_PATHWAY)
  p <- ggplot(d, aes(logFC, -log10(P.Value), colour = dir)) +
    geom_point(alpha = 0.6, size = 2.5) +
    geom_point(data = filter(d, Target %in% LRRK2_PATHWAY),
               shape = 21, fill = "orange", colour = "black", size = 3.5) +
    geom_text_repel(data = lab, aes(label = Target), colour = "black",
                    size = 3, max.overlaps = 25, box.padding = 0.5, show.legend = FALSE) +
    scale_colour_manual(values = c("Up in carriers" = "tomato",
                                   "Down in carriers" = "steelblue",
                                   "Not significant" = "grey80")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", colour = "grey40") +
    labs(title = "NTUH plasma: LRRK2 risk carriers vs non-carriers",
         subtitle = paste0(trimws(model_str), " | orange = LRRK2 pathway target | ",
                           res_pooled$n_group[1], " carriers vs ",
                           res_pooled$n_total[1] - res_pooled$n_group[1], " non-carriers"),
         x = "log2 fold change (carrier vs non-carrier)", y = "-log10(p)", colour = NULL) +
    theme_bw()
  ggsave(file.path(RESULTS_06_DIR, paste0("volcano_lrrk2_carrier_", suffix, ".png")), p,
         width = 10, height = 8, dpi = 300)
}

cat("\nSaved to", RESULTS_06_DIR, "\n")

# --- 7. Annotated CSV for sharing --------------------------------------------
# Model string and filename are both derived, so the file always states which
# model produced it. Hardcoding these is how the wrong file gets read.

out <- read.csv(file.path(RESULTS_06_DIR, "carrier_DA_all_contrasts.csv")) %>%
  mutate(low_detect_flag = detectability < DET_THRESH,
         model  = trimws(model_str),
         cohort = paste0("NTUH plasma Neuro220, EAS, n=", res_pooled$n_total[1])) %>%
  arrange(contrast, P.Value)

write.csv(out, file.path(RESULTS_06_DIR,
          paste0("carrier_DA_", suffix, "_annotated.csv")), row.names = FALSE)
cat("wrote carrier_DA_", suffix, "_annotated.csv\n", sep = "")

# --- 8. Detectability vs significance ----------------------------------------

d <- out %>%
  filter(contrast == "Any LRRK2 risk carrier vs non-carrier",
         P.Value < 0.05 | Target %in% LRRK2_PATHWAY) %>%
  mutate(grp = ifelse(Target %in% LRRK2_PATHWAY, "LRRK2 pathway", "other nominal hit"))

pdet <- ggplot(d, aes(detectability * 100, -log10(P.Value))) +
  annotate("rect", xmin = -Inf, xmax = 50, ymin = -Inf, ymax = Inf,
           fill = "tomato", alpha = 0.08) +
  geom_point(aes(colour = grp, size = abs(logFC))) +
  geom_text_repel(aes(label = Target), size = 3, max.overlaps = 30) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted") +
  geom_vline(xintercept = 50, linetype = "dashed", colour = "tomato") +
  scale_colour_manual(values = c("LRRK2 pathway" = "darkorange",
                                 "other nominal hit" = "grey40")) +
  labs(title = "NTUH: LRRK2 carrier findings against target detectability",
       subtitle = paste0("Shaded = below Alamar's 50% target-detectability threshold. ",
                         "Model: ", trimws(model_str)),
       x = "Target detectability (%)", y = "-log10(p)",
       colour = NULL, size = "|logFC|") +
  theme_bw()

ggsave(file.path(RESULTS_06_DIR, paste0("detectability_vs_p_", suffix, ".png")),
       pdet, width = 9, height = 7, dpi = 300)
cat("wrote detectability_vs_p_", suffix, ".png\n", sep = "")




# --- 9. Plate-adjusted vs unadjusted comparison ------------------------------
# Only runs once BOTH versions exist. Re-run this script with
# CARRIER_PLATE_ADJUST TRUE and FALSE to populate them.
adj_file <- file.path(RESULTS_06_DIR, "plate_as_covar", "carrier_DA_all_contrasts.csv")
if (!file.exists(adj_file)) {
  cat("\nSkipping model comparison - no plate-adjusted file at\n  ", adj_file, "\n")
} else {
  adj <- read.csv(adj_file)
  una <- read.csv(file.path(RESULTS_06_DIR, "carrier_DA_all_contrasts.csv"))
  cmp <- inner_join(
    adj %>% select(contrast, Target, logFC_adj = logFC, p_adj = P.Value, detectability),
    una %>% select(contrast, Target, logFC_una = logFC, p_una = P.Value),
    by = c("contrast", "Target")) %>%
    filter(contrast == "Any LRRK2 risk carrier vs non-carrier")

  cat("\nPearson r of logFC:", round(cor(cmp$logFC_adj, cmp$logFC_una), 3), "\n")
  cat("nominal hits: adjusted", sum(cmp$p_adj < .05),
      "| unadjusted", sum(cmp$p_una < .05),
      "| both", sum(cmp$p_adj < .05 & cmp$p_una < .05), "\n")

  pcmp <- ggplot(cmp, aes(logFC_una, logFC_adj)) +
    geom_abline(slope = 1, linetype = "dashed", colour = "grey50") +
    geom_hline(yintercept = 0, colour = "grey85") +
    geom_vline(xintercept = 0, colour = "grey85") +
    geom_point(aes(colour = p_adj < .05 | p_una < .05), alpha = 0.7) +
    geom_text_repel(data = filter(cmp, p_adj < .05 | p_una < .05 | Target == "LRRK2"),
                    aes(label = Target), size = 3, max.overlaps = 25) +
    scale_colour_manual(values = c("FALSE" = "grey70", "TRUE" = "tomato"),
                        name = "nominal in either model") +
    labs(title = "NTUH LRRK2 carriers: does plate adjustment change the answer?",
         subtitle = paste0("Each point is one of ", nrow(cmp),
                           " targets. On the dashed line = plate makes no difference."),
         x = "logFC, plate NOT adjusted", y = "logFC, plate adjusted") +
    theme_bw()
  ggsave(file.path(RESULTS_06_DIR, "plate_model_comparison.png"),
         pcmp, width = 9, height = 7, dpi = 300)
  cat("wrote plate_model_comparison.png\n")
}
