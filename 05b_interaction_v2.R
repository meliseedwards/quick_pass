# =============================================================================
# Script 05b (v2): Meta-analysis of the sex x phenotype interaction
# Date: September 2026
# Description: Pools the per-cohort phenotype:sex interaction coefficient
#              (from ~ phenotype * sex + age, no plate) across cohorts with the
#              same random-effects model as 05. Answers: does the PD effect
#              differ between men and women anywhere on the panel?
#              Input: sex_by_phenotype_interaction_all_cohorts.csv 
# =============================================================================

library(tidyverse)
library(metafor)

source("~/proteomics_nulisa/scripts/nulisa_pipeline_v2/quick_pass/00_config_v2.R")

OUT_DIR <- file.path(PROJECT_DIR, "05b_meta_interaction")
dir.create(OUT_DIR, showWarnings = FALSE)

# --- 1. Load ------------------------------------------------------------------
int_all <- read.csv(file.path(PROJECT_DIR, "sex_by_phenotype_interaction_all_cohorts.csv"),
                    stringsAsFactors = FALSE)
needed <- c("cohort_label", "Target", "int_logFC", "int_SE")
stopifnot(all(needed %in% names(int_all)))

int <- int_all %>%
  filter(cohort_label %in% META_COHORTS, !is.na(int_logFC), !is.na(int_SE), int_SE > 0)

dup <- int %>% dplyr::count(Target, cohort_label) %>% filter(n > 1)
if (nrow(dup) > 0) { print(dup); stop("Duplicate Target x cohort rows.") }

cat("Rows in:", nrow(int), "| targets:", n_distinct(int$Target), "\n")
print(int %>% dplyr::count(cohort_label, name = "n_targets"))

# --- 2. One target ------------------------------------------------------------
# Interaction = (PD - Control in men) - (PD - Control in women), log2 scale.
meta_int <- function(d) {
  k <- nrow(d)
  if (k < 2) return(tibble(k = k))
  fit_kh <- rma(yi = d$int_logFC, sei = d$int_SE, method = "REML", test = "knha")
  fit_z  <- rma(yi = d$int_logFC, sei = d$int_SE, method = "REML")
  pooled <- as.numeric(fit_kh$b)
  loo    <- leave1out(fit_kh)
  tibble(k             = k,
         int_logFC     = pooled,
         ci_lb         = fit_kh$ci.lb,
         ci_ub         = fit_kh$ci.ub,
         P_KH          = fit_kh$pval,
         P_z           = fit_z$pval,
         I2            = fit_kh$I2,
         loo_same_sign = all(sign(loo$estimate) == sign(pooled)),
         n_neg         = sum(d$int_logFC < 0),
         n_pos         = sum(d$int_logFC > 0))
}

# --- 3. Run -------------------------------------------------------------------
meta <- int %>%
  group_by(Target) %>%
  group_modify(~ meta_int(.x)) %>%
  ungroup() %>%
  mutate(FDR_KH = p.adjust(P_KH, "BH"),
         FDR_z  = p.adjust(P_z,  "BH")) %>%
  arrange(P_z)

write.csv(meta, file.path(OUT_DIR, "meta_interaction_sex_by_phenotype.csv"), row.names = FALSE)

# --- 4. Report ----------------------------------------------------------------
cat("\nTargets:", nrow(meta),
    "| FDR < 0.05: Knapp-Hartung", sum(meta$FDR_KH < 0.05), "| conventional RE", sum(meta$FDR_z < 0.05),
    "| nominal p < 0.05 (conventional):", sum(meta$P_z < 0.05), "of", nrow(meta), "\n")

cat("\nTop 15:\n")
print(meta %>%
        mutate(across(c(int_logFC, ci_lb, ci_ub), ~ round(.x, 2)), I2 = round(I2)) %>%
        head(15), n = 15, width = Inf)

cat("\nShortlist from the per-cohort runs:\n")
print(meta %>% filter(Target %in% c("FGF21", "GOT1", "CHI3L1")) %>%
        mutate(across(c(int_logFC, ci_lb, ci_ub), ~ round(.x, 2)), I2 = round(I2)), width = Inf)

# --- 5. Forest plots, top 6 --------------------------------------------------
top <- head(meta$Target, 6)
png(file.path(OUT_DIR, "forest_interaction_top6.png"), width = 12, height = 7, units = "in", res = 200)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
for (tg in top) {
  d <- int %>% filter(Target == tg)
  fit <- rma(yi = d$int_logFC, sei = d$int_SE, method = "REML", test = "knha")
  forest(fit, slab = d$cohort_label, xlab = "phenotype x sex interaction (log2)", main = tg, cex = 0.9)
}
dev.off()

cat("\nWrote", OUT_DIR, "\n")