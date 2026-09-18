# =============================================================================
# Script 05 (v2): Random-effects meta-analysis of per-cohort DA results
# Date: September 2026
# Description: For each target, combines the cohort-level lm estimates
#              (logFC, SE) with a random-effects model (metafor::rma, REML
#              random effects; Knapp-Hartung and fixed-effect p reported). 
#              Reports pooled effect, heterogeneity, leave-one-out stability 
#              and cohort-specific QC context.
# =============================================================================



# --- 0. Setup ----------------------------------------------------------------

library(tidyverse)
library(metafor)

source("~/proteomics_nulisa/scripts/nulisa_pipeline_v2/quick_pass/00_config_v2.R")

stopifnot(META_INPUT %in% c("plate_adj", "plate_none", "sva"))
meta_cohorts <- META_COHORTS

META_DIR <- file.path(PROJECT_DIR, paste0("05_meta_", META_INPUT))
dir.create(META_DIR, showWarnings = FALSE, recursive = TRUE)

cat("META | input:", META_INPUT, "| cohorts:", paste(meta_cohorts, collapse = ", "), "\n")



# --- 1. Load the combined DA table -------------------------------------------

da_all <- read.csv(file.path(PROJECT_DIR, paste0("DA_all_cohorts_combined_", META_INPUT, ".csv")),
                   stringsAsFactors = FALSE)

# A misspelled or missing cohort needs to fail loudly here, not silently pool less studies.
missing <- setdiff(meta_cohorts, unique(da_all$cohort_label))
if (length(missing) > 0) stop("Cohort(s) not in combined file: ", paste(missing, collapse = ", "))

# Filter to cohorts of interest and print summary
da <- da_all %>%
  filter(cohort_label %in% meta_cohorts, !is.na(logFC), !is.na(SE), SE > 0)

cat("Rows in:", nrow(da), "| targets:", n_distinct(da$Target), "\n")
print(da %>% dplyr::count(cohort_label, plate_model, name = "n_targets"))



# --- 2. One target -----------------------------------------------------------
# Input: the rows for one target (one per cohort). 
# REML = random effects (Default)
# FE = fixed effects
meta_one <- function(d) {
  k <- nrow(d)
  if (k < 2) {
    return(tibble(k = k, logFC_meta = NA_real_, SE_meta = NA_real_,
                  ci_lb = NA_real_, ci_ub = NA_real_, P_meta = NA_real_,
                  logFC_FE = NA_real_, P_FE = NA_real_,
                  tau2 = NA_real_, I2 = NA_real_, Q_p = NA_real_,
                  loo_max_shift = NA_real_, loo_all_same_sign = NA,
                  n_cohorts_sig = sum(d$adj.P.Val < FDR_THRESHOLD),
                  n_up = sum(d$logFC > 0), n_down = sum(d$logFC < 0)))
  }

  # random effects & fixed effect is kept alongside as seen in the GP2 GRS code.
  fit    <- rma(yi = d$logFC, sei = d$SE, method = "REML") # default               
  fit_kh <- rma(yi = d$logFC, sei = d$SE, method = "REML", test = "knha") # conservative
  fit_fe <- rma(yi = d$logFC, sei = d$SE, method = "FE") # fixed effects

  # leave-one-out: refit without each cohort; how much does the pooled
  # estimate move, and does it keep its sign?
  loo <- leave1out(fit)
  loo_shift <- max(abs(loo$estimate - as.numeric(fit$b)))
  loo_sign  <- all(sign(loo$estimate) == sign(as.numeric(fit$b)))

  tibble(k              = k,
         logFC_meta     = as.numeric(fit$b),
         SE_meta        = fit$se,
         ci_lb          = fit$ci.lb,
         ci_ub          = fit$ci.ub,
         P_meta         = fit$pval,
         logFC_FE       = as.numeric(fit_fe$b),
         P_FE           = fit_fe$pval,
         P_KH = fit_kh$pval,
         tau2           = fit$tau2,
         I2             = fit$I2,
         Q_p            = fit$QEp,
         loo_max_shift  = loo_shift,
         loo_all_same_sign = loo_sign,
         n_cohorts_sig  = sum(d$adj.P.Val < FDR_THRESHOLD),
         n_up           = sum(d$logFC > 0),
         n_down         = sum(d$logFC < 0))
}



# --- 3. Run every target -----------------------------------------------------

meta <- da %>%
  group_by(Target) %>%
  group_modify(~ meta_one(.x)) %>%
  ungroup() %>%
  mutate(adj_P_meta = p.adjust(P_meta, method = "BH"),
         adj_P_KH = p.adjust(P_KH, "BH"), 
         adj_P_FE = p.adjust(P_FE, "BH"),
         bonf_P_meta = p.adjust(P_meta, method = "bonferroni"),
         direction = case_when(adj_P_meta < FDR_THRESHOLD & logFC_meta > 0 ~ "Up in PD",
                               adj_P_meta < FDR_THRESHOLD & logFC_meta < 0 ~ "Down in PD",
                               TRUE ~ "Not significant")) %>%
  arrange(P_meta)

# Per-target QC context from the cohort files: worst detectability across
# cohorts and whether any cohort flagged the target.
qc_ctx <- da %>%
  group_by(Target) %>%
  summarise(min_detectability  = min(detectability, na.rm = TRUE),
            any_low_detect     = any(low_detect_flag, na.rm = TRUE),
            any_vendor_high_cv = any(vendor_high_cv, na.rm = TRUE),
            any_cross_reactive = any(cross_reactive, na.rm = TRUE),
            any_zero_imbalance = any(zero_imbalance_flag, na.rm = TRUE),
            in_meta6           = any(in_meta6),
            .groups = "drop")

meta <- meta %>% left_join(qc_ctx, by = "Target") %>% mutate(meta_input = META_INPUT)
stopifnot(!any(duplicated(meta$Target)))

write.csv(meta, file.path(META_DIR, paste0("meta_results_", META_INPUT, ".csv")), row.names = FALSE)



# --- 4. Report ---------------------------------------------------------------

cat("\nTargets meta-analysed:", nrow(meta),
    "| significant at FDR <", FDR_THRESHOLD, ":", sum(meta$adj_P_meta < FDR_THRESHOLD),
    "| Bonferroni:", sum(meta$bonf_P_meta < 0.05), "\n")
cat("Median I2:", round(median(meta$I2, na.rm = TRUE), 1),
    "| targets with I2 > 50%:", sum(meta$I2 > 50, na.rm = TRUE), "\n")
cat("Significant hits that keep their sign in every leave-one-out refit:",
    sum(meta$adj_P_meta < FDR_THRESHOLD & meta$loo_all_same_sign, na.rm = TRUE), "\n")

cat("\nTop 20:\n")
print(meta %>%
        select(Target, k, logFC_meta, ci_lb, ci_ub, P_meta, adj_P_meta, P_KH, P_FE, I2,
               n_cohorts_sig, n_up, n_down, loo_all_same_sign, any_low_detect) %>%
        mutate(across(c(logFC_meta, ci_lb, ci_ub), ~ round(.x, 2)),
               I2 = round(I2, 0)) %>%
        head(20), n = 20, width = Inf)



# --- 5. Figures --------------------------------------------------------------

# forest plots for the top hits: per-cohort estimates + pooled diamond
top <- meta %>% arrange(P_meta) %>% head(12) %>% pull(Target)
if (length(top) > 0) {
  png(file.path(META_DIR, paste0("forest_top_hits_", META_INPUT, ".png")),
      width = 14, height = 10, units = "in", res = 200)
  par(mfrow = c(3, 4), mar = c(4, 4, 3, 1))
  for (tg in top) {
    d <- da %>% filter(Target == tg)
    fit <- rma(yi = d$logFC, sei = d$SE, method = "REML")
    forest(fit, slab = d$cohort_label, xlab = "log2FC (PD vs Control)", main = tg, cex = 0.9)
  }
  dev.off()
}

# meta volcano
p_vol <- ggplot(meta, aes(logFC_meta, -log10(P_meta), colour = direction)) +
  geom_point(alpha = 0.7, size = 2.2) +
  geom_point(data = filter(meta, I2 > 50), shape = 1, size = 4, colour = "grey30") +
  ggrepel::geom_text_repel(data = filter(meta, adj_P_meta < FDR_THRESHOLD) %>% head(20),
                           aes(label = Target), colour = "black", size = 3, max.overlaps = 30) +
  scale_colour_manual(values = c("Up in PD" = "tomato", "Down in PD" = "steelblue",
                                 "Not significant" = "grey80")) +
  labs(title = paste0("Meta-analysis (", META_INPUT, "): PD vs Control, ",
                      length(meta_cohorts), " cohorts"),
       subtitle = "grey ring = I2 > 50% (heterogeneous across cohorts)",
       x = "pooled log2FC", y = "-log10(p)") +
  theme_bw()
ggsave(file.path(META_DIR, paste0("volcano_meta_", META_INPUT, ".png")), p_vol,
       width = 9, height = 7, dpi = 150)

# per-cohort logFC heatmap for everything significant in the meta-analysis
sig <- meta %>% arrange(P_meta) %>% head(30) %>% pull(Target)
if (length(sig) > 0) {
  hm <- da %>% filter(Target %in% sig) %>%
    mutate(Target = factor(Target, levels = rev(sig)),
           sig_in_cohort = adj.P.Val < FDR_THRESHOLD)
  p_hm <- ggplot(hm, aes(cohort_label, Target, fill = logFC)) +
    geom_tile() +
    geom_point(data = filter(hm, sig_in_cohort), shape = 8, size = 1.5) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "tomato", midpoint = 0,
                         name = "log2FC") +
    labs(title = "Top 30 meta-analysis targets: per-cohort log2FC (* = FDR < 0.05 in that cohort)",
         x = NULL, y = NULL) +
    theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
  ggsave(file.path(META_DIR, paste0("heatmap_meta_sig_", META_INPUT, ".png")), p_hm,
         width = 8, height = max(4, 0.25 * length(sig) + 2), dpi = 150)
}

cat("\nWrote results and figures to", META_DIR, "\n")