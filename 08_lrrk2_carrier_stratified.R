# =============================================================================
# Script 08 (NTUH) v2: LRRK2 carrier analysis WITHIN phenotype (PD Carrier vs
# PD Non-carrier, and Control Carrier vs Control Non-carrier).
# Date: August 2026
# Description: Script 06 pools cases and controls with phenotype as a covariate.
#              This runs the same contrast separately within PD and within
#              Controls.
#
# Run:  source("08_lrrk2_carrier_stratified.R") 
# Note: change config cohort before running new datasets. 
# =============================================================================



# --- 0. Setup ----------------------------------------------------------------

library(SummarizedExperiment)
library(tidyverse)
library(limma)
library(ggrepel)

source("~/proteomics/nulisa_pipeline_v2/quick_pass/00_config.R")
stopifnot(COHORT %in% c("ntuh", "umklm", "kul"))

RESULTS_08_DIR <- file.path(RESULTS_DIR, "08_lrrk2_carrier_stratified")
dir.create(RESULTS_08_DIR, showWarnings = FALSE, recursive = TRUE)

RISK_VARIANTS <- c("LRRK2_Gly2385Arg", "LRRK2_Arg1628Pro")

LRRK2_PATHWAY <- c("LRRK2", "pLRRK2-S1292",
                   "RAB10", "RAB12", "RAB29",
                   "pRAB10-T73", "pRAB12-S106", "pRAB29-T71")

# Set TRUE to include plate. Within PD that reduces the contrast to plate 2 only.
CARRIER_PLATE_ADJUST <- (PLATE_CORRECTION == "covariate")
suffix <- if (CARRIER_PLATE_ADJUST) "plateAdjusted" else "plateUnadjusted"



# --- 1. base cohort ------------------------------------------------------------

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

stopifnot(identical(as.character(meta$DONOR_ID), as.character(meta$GP2ID)))
cohort_ids <- meta$GP2ID

variants <- read.csv(file.path(DATA_DIR, VARIANT_REPORT_FILE),
                     stringsAsFactors = FALSE) %>%
  filter(Gene == "LRRK2", GP2ID %in% meta$DONOR_ID)

risk_ids   <- variants %>% filter(Variant_Name %in% RISK_VARIANTS) %>% pull(GP2ID) %>% unique()
other_ids  <- variants %>% filter(!Variant_Name %in% RISK_VARIANTS) %>% pull(GP2ID) %>% unique()
other_only <- setdiff(other_ids, risk_ids)

mk <- read.delim(file.path(DATA_DIR, MASTER_KEY_FILE), stringsAsFactors = FALSE)
genotyped   <- mk$GP2ID[(mk$nba %in% 1 | mk$wgs %in% 1) & mk$GP2ID %in% cohort_ids]
ungenotyped <- setdiff(cohort_ids, genotyped)

cat("Risk-variant carriers:", length(risk_ids),
    "| non-risk LRRK2 only (excluded):", length(other_only),
    "| ungenotyped (excluded):", length(ungenotyped), "\n")

meta <- meta %>%
  mutate(lrrk2_risk = ifelse(GP2ID %in% risk_ids, "Carrier", "Non-carrier"),
         analysable = !GP2ID %in% c(other_only, ungenotyped))
colData(se) <- DataFrame(meta)

se_c <- se[, meta$analysable]
colData(se_c)$lrrk2_risk <- factor(colData(se_c)$lrrk2_risk,
                                   levels = c("Non-carrier", "Carrier"))
metac <- as.data.frame(colData(se_c))

cat("\nAnalysis cohort:", ncol(se_c), "samples\n")
cat("Carrier x phenotype:\n"); print(table(metac$lrrk2_risk, metac$phenotype_clean))
cat("\nCarrier x plate, WITHIN PD:\n")
print(table(metac$lrrk2_risk[metac$phenotype_clean == "PD"],
            metac$PlateID[metac$phenotype_clean == "PD"]))
cat("\nCarrier x plate, WITHIN Controls:\n")

print(table(metac$lrrk2_risk[metac$phenotype_clean == "Control"],
            metac$PlateID[metac$phenotype_clean == "Control"]))



# --- 3. Fit one contrast -----------------------------------------------------
# phenotype_clean is dropped automatically when it has one level, which is the
# case in every stratified fit here.

fit_contrast <- function(se_in, group_col, label) {
  m <- as.data.frame(colData(se_in))
  m[[group_col]] <- droplevels(factor(m[[group_col]]))
  if (nlevels(m[[group_col]]) < 2) { cat("SKIP", label, "- only one group\n"); return(NULL) }

  form <- paste("~", group_col)
  if (nlevels(droplevels(factor(m$phenotype_clean))) > 1)
    form <- paste(form, "+ phenotype_clean")
  form <- paste(form, "+", paste(DA_COVARIATES, collapse = " + "))
  # plate only enters if it varies AND both groups appear on more than one level
  plate_ok <- nlevels(droplevels(factor(m$PlateID))) > 1
  if ((CARRIER_PLATE_ADJUST || PLATE_CORRECTION == "covariate") && plate_ok)
    form <- paste(form, "+ PlateID")

  design <- model.matrix(as.formula(form), data = m)
  stopifnot(nrow(design) == ncol(se_in))
  if (qr(design)$rank < ncol(design)) {
    cat("SKIP", label, "- design is rank-deficient (a term is confounded)\n"); return(NULL)
  }

  coef_name <- grep(paste0("^", group_col), colnames(design), value = TRUE)[1]
  fit <- eBayes(lmFit(assay(se_in, "npq"), design))
  res <- topTable(fit, coef = coef_name, number = Inf, adjust.method = "BH") %>%
    rownames_to_column("Target") %>%
    mutate(contrast = label,
           model    = trimws(form),
           n_total  = ncol(se_in),
           n_group  = sum(m[[group_col]] == levels(m[[group_col]])[2]))

  cat("\n---", label, "---\n")
  cat("model:", trimws(form), "\n")
  cat("n =", ncol(se_in), "| carriers:",
      sum(m[[group_col]] == levels(m[[group_col]])[2]), "\n")
  cat("Significant at FDR <", FDR_THRESHOLD, ":", sum(res$adj.P.Val < FDR_THRESHOLD), "\n")
  print(res %>% select(Target, logFC, P.Value, adj.P.Val) %>% head(8))
  res
}



# --- 4. Stratified contrasts -------------------------------------------------
se_pd <- se_c[, colData(se_c)$phenotype_clean == "PD"]
colData(se_pd)$lrrk2_risk <- droplevels(colData(se_pd)$lrrk2_risk)
res_pd <- fit_contrast(se_pd, "lrrk2_risk", "PD only: carrier vs non-carrier")

se_ct <- se_c[, colData(se_c)$phenotype_clean == "Control"]
colData(se_ct)$lrrk2_risk <- droplevels(colData(se_ct)$lrrk2_risk)
res_ct <- fit_contrast(se_ct, "lrrk2_risk", "Controls only: carrier vs non-carrier")

stopifnot(!is.null(res_pd), !is.null(res_ct))



# --- 5. LRRK2 pathway targets, side by side ----------------------------------
det <- as.data.frame(rowData(se_c)) %>% select(Target, detectability)

pathway <- bind_rows(res_pd, res_ct) %>%
  filter(Target %in% LRRK2_PATHWAY) %>%
  left_join(det, by = "Target") %>%
  select(contrast, Target, logFC, P.Value, adj.P.Val, detectability, n_total, n_group) %>%
  arrange(Target, contrast)

cat("\n\n===== LRRK2 pathway targets, PD and Controls separately =====\n")
cat("Low detectability = the target sits at the assay floor in most samples.\n\n")
print(as.data.frame(pathway), row.names = FALSE)
write.csv(pathway, file.path(RESULTS_08_DIR,
          paste0("lrrk2_pathway_stratified_", suffix, ".csv")), row.names = FALSE)



# --- 6. Save all results -----------------------------------------------------
all_res <- bind_rows(res_pd, res_ct) %>%
  left_join(det, by = "Target") %>%
  mutate(low_detect_flag = detectability < DET_THRESH,
         cohort = paste0("NTUH plasma Neuro220, EAS, n=", ncol(se_c))) %>%
  arrange(contrast, P.Value)

qc_n <- read.csv(file.path(RESULTS_04_DIR, "qc_target_Neuro220.csv")) %>%
  select(Target, ICC, ICC_residualised, intra_cv, inter_cv)

all_res <- bind_rows(res_pd, res_ct) %>%
  left_join(det, by = "Target") %>%
  left_join(qc_n, by = "Target") %>%
  mutate(low_detect_flag = detectability < DET_THRESH,
         cohort = paste0("NTUH plasma Neuro220, EAS, n=", ncol(se_c))) %>%
  arrange(contrast, P.Value)

write.csv(all_res, file.path(RESULTS_08_DIR,
          paste0("carrier_DA_stratified_", suffix, ".csv")), row.names = FALSE)

cat("\nwrote carrier_DA_stratified_", suffix, ".csv\n", sep = "")


# --- 7. PD vs Control effect sizes, side by side -----------------------------
# Where the two agree, the pooled estimate from Script 06 is a fair summary.
# Where they diverge, the pooled number is averaging two different things.
cmp <- inner_join(
  res_pd %>% select(Target, logFC_pd = logFC, p_pd = P.Value),
  res_ct %>% select(Target, logFC_ct = logFC, p_ct = P.Value),
  by = "Target") %>%
  left_join(det, by = "Target") %>%
  mutate(highlight = Target %in% LRRK2_PATHWAY)

cat("\nPearson r of logFC, PD vs Controls:", round(cor(cmp$logFC_pd, cmp$logFC_ct), 3), "\n")

p_cmp <- ggplot(cmp, aes(logFC_ct, logFC_pd)) +
  geom_abline(slope = 1, linetype = "dashed", colour = "grey50") +
  geom_hline(yintercept = 0, colour = "grey85") +
  geom_vline(xintercept = 0, colour = "grey85") +
  geom_point(aes(colour = highlight), alpha = 0.7, size = 2.2) +
  geom_text_repel(data = filter(cmp, highlight | p_pd < 0.01 | p_ct < 0.01),
                  aes(label = Target), size = 3.2, max.overlaps = 25,
                  box.padding = 0.5, show.legend = FALSE) +
  scale_colour_manual(values = c("FALSE" = "grey70", "TRUE" = "darkorange"),
                      labels = c("other target", "LRRK2 pathway"), name = NULL) +
  labs(title = "LRRK2 carrier effect: within PD vs within Controls",
       subtitle = paste0("Each point is one target. On the dashed line = same effect in ",
                         "both groups.\nPD: ", res_pd$n_group[1], " carriers of ",
                         res_pd$n_total[1], ".  Controls: ", res_ct$n_group[1],
                         " carriers of ", res_ct$n_total[1], ".  Model: ",
                         trimws(res_pd$model[1])),
       x = "logFC, Controls only", y = "logFC, PD only") +
  theme_bw(base_size = 13) + theme(legend.position = "top")

ggsave(file.path(RESULTS_08_DIR, paste0("pd_vs_control_carrier_effect_", suffix, ".png")),
       p_cmp, width = 9, height = 7, dpi = 300)
cat("wrote pd_vs_control_carrier_effect_", suffix, ".png\n", sep = "")

cat("\nSaved to", RESULTS_08_DIR, "\n")



# --- 10. NPQ vs age for the PD-only hits -------------------------------------
# Two regression lines per panel. Parallel lines mean a single linear age term
# is adequate; the vertical gap between them is the carrier effect.
hits_pd <- res_pd %>% filter(adj.P.Val < FDR_THRESHOLD) %>% arrange(P.Value) %>% pull(Target)

mpd <- as.data.frame(colData(se_pd))
age_gap <- mean(mpd$age[mpd$lrrk2_risk == "Carrier"], na.rm = TRUE) -
           mean(mpd$age[mpd$lrrk2_risk == "Non-carrier"], na.rm = TRUE)

dage <- as.data.frame(assay(se_pd, "npq")[hits_pd, , drop = FALSE]) %>%
  rownames_to_column("Target") %>%
  pivot_longer(-Target, names_to = "SampleName", values_to = "NPQ") %>%
  left_join(mpd %>% select(SampleName, age, lrrk2_risk), by = "SampleName") %>%
  left_join(res_pd %>% select(Target, logFC, adj.P.Val), by = "Target") %>%
  mutate(lab = paste0(Target, "  (logFC ", sprintf("%+.2f", logFC),
                      ", FDR ", signif(adj.P.Val, 2), ")"),
         lab = factor(lab, levels = unique(lab[order(match(Target, hits_pd))])))

p_age_npq <- ggplot(dage, aes(age, NPQ, colour = lrrk2_risk, fill = lrrk2_risk)) +
  geom_point(alpha = 0.55, size = 1.9) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 0.9) +
  facet_wrap(~ lab, scales = "free_y") +
  scale_colour_manual(values = c("Non-carrier" = "grey45", "Carrier" = "steelblue")) +
  scale_fill_manual(values = c("Non-carrier" = "grey45", "Carrier" = "steelblue")) +
  labs(title = "Within PD: NPQ vs age, stratified by LRRK2 carrier status",
       subtitle = paste0("NTUH plasma, PD only. ", res_pd$n_group[1], " carriers vs ",
                         res_pd$n_total[1] - res_pd$n_group[1], " non-carriers. ",
                         "Carriers are ", abs(round(age_gap, 1)),
                         " years younger on average.\n",
                        "Model: ", trimws(res_pd$model[1])),
       x = "Age (years)", y = "NPQ (log2)", colour = NULL, fill = NULL) +
  theme_bw(base_size = 13) + theme(legend.position = "top")

ggsave(file.path(RESULTS_08_DIR, paste0("age_npq_PDhits_", suffix, ".png")),
       p_age_npq, width = 11, height = 7, dpi = 300)

cat("wrote age_npq_PDhits_", suffix, ".png\n", sep = "")




# --- test for plate effects within PD ------------------------------------------------

# does plate 1 differ from plates 2-3 among NON-carriers?
# no non-carriers on plate 1, so use PD only and compare plates 2 vs 3
mp <- as.data.frame(colData(se_pd))
for (t in c("ARSA","MAG","GFAP","MOG","NGF","GRN")) {
  y <- assay(se_pd,"npq")[t, ]
  # carriers only, so carrier status is constant - any difference is plate
  cc <- mp$lrrk2_risk == "Carrier"
  p <- summary(lm(y[cc] ~ mp$PlateID[cc]))$coef
  cat(sprintf("%-6s plate effect among carriers: F-test p = %.3f\n", t,
      anova(lm(y[cc] ~ mp$PlateID[cc]))$`Pr(>F)`[1]))
}