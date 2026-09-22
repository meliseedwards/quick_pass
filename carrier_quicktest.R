# =============================================================================
# Quick carrier check: is the encoded protein lower in genetic carriers?
# LRRK2 protein by LRRK2 carrier status, GBA protein by GBA1 carrier status,
# per cohort, R12 variant report. Wilcoxon two-sided; t-test alongside.
# Genotyped status comes from colData (ancestry_source, joined in Script 02).
# ELPD (p111) is not in release 12: expect 0 genotyped.
# Pre-test for the Korea slides; not the full model in 06.
# =============================================================================

library(SummarizedExperiment)
library(tidyverse)
library(patchwork)
library(ggsignif)

source("~/proteomics_nulisa/scripts/nulisa_pipeline_v2/quick_pass/00_config_v2.R")

KOREA_DIR <- file.path(PROJECT_DIR, "korea_talk_2026")
dir.create(KOREA_DIR, showWarnings = FALSE)

GENE_TARGET <- c(LRRK2 = "LRRK2", GBA1 = "GBA")   # variant-report gene -> panel target
COHORTS <- c(p111 = "ELPD", p118a = "UMKLM", p118b = "KUL", p121 = "NTUH", p136 = "TRAPCAF",
             p149 = "CANDAS-SMPD", p143 = "Nigeria-PD", p144_p146 = "LARGEPD+PDGNRTN", 
             p_ppmi = "PPMI")
LRRK2_RISK <- c("LRRK2_Gly2385Arg", "LRRK2_Arg1628Pro")


# change with future runs so titles stop bleeding off page
# base <- theme_bw(base_size = 22) +
#  theme(legend.position = "none", panel.grid.minor = element_blank(),
#        strip.background = element_blank(), strip.text = element_text(size = 22, face = "bold"),
#        plot.title = element_text(size = 17, face = "bold"), plot.subtitle = element_text(size = 13))

# in the variants block: keep the intronic variant, give it its own class, carry zygosity
variants <- read.csv(file.path(DATA_DIR, VARIANT_REPORT_FILE), stringsAsFactors = FALSE) %>%
  filter(Gene %in% names(GENE_TARGET)) %>%
    distinct(GP2ID, Variant_Name, .keep_all = TRUE) %>%
  mutate(class = case_when(Variant_Name == "GBA1_c.1225-34C>A" ~ "GBA1 rs3115534-G (African risk)", 
                           Variant_Name %in% LRRK2_RISK ~ "LRRK2 risk (East Asian)",
                           grepl("athogenic", Pathogenicity) ~ "Pathogenic/likely pathogenic",
                           TRUE ~ "PD risk (common)"))
cat("Variant classes in report:\n"); print(variants %>% dplyr::count(Gene, class))

carrier_test <- function(cid, label) {
  se <- readRDS(file.path(PROJECT_DIR, paste0("results_", cid, "_v2_filter_none"), "02_merge_gp2",
                          paste0("se_", cid, "_Neuro220_with_metadata.rds")))
  meta_all <- as.data.frame(colData(se)) %>% filter(phenotype_clean %in% c("PD", "Control"))
  meta <- meta_all %>% filter(!is.na(GP2ID), !is.na(ancestry_source))   # genotyped only
  n_total <- nrow(meta_all); n_geno <- nrow(meta)
  note <- paste0(n_geno, " of ", n_total, " samples genotyped in R12",
                 if (n_geno == 0) " (cohort not in R12)" else if (n_geno < 0.5 * n_total) " (most not yet in R12)" else "")
  out <- NULL; breakdown <- NULL
  
  for (g in names(GENE_TARGET)) {
    tg  <- GENE_TARGET[[g]]
    ids <- unique(variants$GP2ID[variants$Gene == g])
    d <- meta %>% mutate(carrier = GP2ID %in% ids, NPQ = assay(se, "npq")[tg, SampleName]) %>% filter(!is.na(NPQ))
    n_car <- sum(d$carrier); n_non <- sum(!d$carrier)

    # which variants / classes the carriers in this cohort actually have
    vc <- variants %>% filter(Gene == g, GP2ID %in% d$GP2ID) %>%
      dplyr::count(Variant_Name, Pathogenicity, class, sort = TRUE) %>%
      mutate(cohort = label, gene = g, .before = 1)
    breakdown <- bind_rows(breakdown, vc)
    class_summary <- vc %>% group_by(class) %>% summarise(n = sum(n), .groups = "drop") %>%
      mutate(s = paste0(class, ": ", n)) %>% pull(s) %>% paste(collapse = "; ")

    row <- tibble(cohort = label, gene = g, target = tg, n_carrier = n_car, n_noncarrier = n_non,
                  n_carrier_PD  = sum(d$carrier & d$phenotype_clean == "PD"),
                  n_carrier_CON = sum(d$carrier & d$phenotype_clean == "Control"),   # non-manifesting carriers
                  n_carrier_pathogenic = length(intersect(d$GP2ID[d$carrier], variants$GP2ID[variants$Gene == g & variants$class == "Pathogenic/likely pathogenic"])),
                  n_carrier_risk       = length(intersect(d$GP2ID[d$carrier], variants$GP2ID[variants$Gene == g & variants$class != "Pathogenic/likely pathogenic"])),
                  carrier_classes = class_summary,
                  note = note)
   
    if (n_car >= 3) {
      pd <- d %>% filter(phenotype_clean == "PD")
      row <- row %>% mutate(mean_carrier      = mean(d$NPQ[d$carrier]),   mean_noncarrier   = mean(d$NPQ[!d$carrier]),
                            median_carrier    = median(d$NPQ[d$carrier]), median_noncarrier = median(d$NPQ[!d$carrier]),
                            diff_median = median_carrier - median_noncarrier,
                            p_wilcoxon = wilcox.test(NPQ ~ carrier, data = d)$p.value,
                            p_ttest    = t.test(NPQ ~ carrier, data = d)$p.value,
                            n_carrier_PD_only = sum(pd$carrier),
                            diff_median_PD_only = if (sum(pd$carrier) >= 3 & sum(!pd$carrier) >= 3)
                              median(pd$NPQ[pd$carrier]) - median(pd$NPQ[!pd$carrier]) else NA,
                            p_wilcoxon_PD_only = if (sum(pd$carrier) >= 3 & sum(!pd$carrier) >= 3)
                              wilcox.test(NPQ ~ carrier, data = pd)$p.value else NA)
      cat("  ", label, g, "variants carried:\n"); print(as_tibble(vc %>% select(-cohort, -gene)))
    }
    out <- bind_rows(out, row)
  }
  list(res = out, breakdown = breakdown)
}

runs <- lapply(names(COHORTS), function(x) carrier_test(x, COHORTS[[x]]))
res       <- bind_rows(lapply(runs, `[[`, "res"))
breakdown <- bind_rows(lapply(runs, `[[`, "breakdown"))

print(res %>% mutate(across(where(is.numeric), ~ signif(.x, 3))) %>% as_tibble(), n = 20, width = Inf)
write.csv(res,       file.path(KOREA_DIR, "carrier_quicktest_LRRK2_GBA1.csv"), row.names = FALSE)
write.csv(breakdown, file.path(KOREA_DIR, "carrier_variant_breakdown_by_cohort.csv"), row.names = FALSE)

# optional xlsx, one tab per cohort, if writexl is installed
if (requireNamespace("writexl", quietly = TRUE)) {
  tabs <- split(breakdown, breakdown$cohort)
  writexl::write_xlsx(tabs, file.path(KOREA_DIR, "carrier_variant_breakdown_by_cohort.xlsx"))
  cat("wrote xlsx with", length(tabs), "tabs\n")
} else cat("writexl not installed - CSV only\n")

rs_test <- function(cid, label) {
  se <- readRDS(file.path(PROJECT_DIR, paste0("results_", cid, "_v2_filter_none"), "02_merge_gp2",
                          paste0("se_", cid, "_Neuro220_with_metadata.rds")))
  meta <- as.data.frame(colData(se)) %>%
    filter(phenotype_clean %in% c("PD", "Control"), !is.na(GP2ID), !is.na(ancestry_source))
  rs <- variants %>% filter(Variant_Name == "GBA1_c.1225-34C>A") %>% select(GP2ID, Zygosity)
  d <- meta %>% left_join(rs, by = "GP2ID") %>%
    mutate(G_dose = case_when(Zygosity == "hom" ~ 2L, Zygosity == "het" ~ 1L, TRUE ~ 0L),   # absent = TT
           NPQ = assay(se, "npq")["GBA", SampleName]) %>% filter(!is.na(NPQ))
    stopifnot(nrow(d) <= nrow(meta))   # join must not add samples
  cat("\n", label, "rs3115534 G dosage:\n"); print(table(d$G_dose, d$phenotype_clean))
  fit <- lm(NPQ ~ G_dose + phenotype_clean + age + sex_clean, data = d)
    fit_int <- lm(NPQ ~ G_dose * phenotype_clean + age + sex_clean, data = d)
  cat("  G_dose x phenotype interaction:\n")
  print(round(summary(fit_int)$coefficients[grep("G_dose", rownames(summary(fit_int)$coefficients)), ], 4))
  d %>% group_by(phenotype_clean, G_dose) %>%
    summarise(n = n(), median_GBA = round(median(NPQ), 2), .groups = "drop") %>% print()
  print(round(summary(fit)$coefficients["G_dose", ], 4))
  d %>% group_by(G_dose) %>% summarise(n = n(), median_GBA = median(NPQ), .groups = "drop") %>% print()

  write.csv(d %>% select(GP2ID, phenotype_clean, age, sex_clean, G_dose, NPQ),
            file.path(KOREA_DIR, paste0("rs3115534_dosage_GBA_", cid, ".csv")), row.names = FALSE)
}
for (x in list(c("p136", "TRAPCAF"))) rs_test(x[1], x[2])




# --- Plots for slide deck (Korea) --------------------------------------------


pal  <- c(Control = "#66C5CC", PD = "#F6CF71")     # replace with house Control/PD colours
base <- theme_bw(base_size = 22) +
  theme(legend.position = "none", panel.grid.minor = element_blank(),
        strip.background = element_blank(), strip.text = element_text(size = 22, face = "bold"))
p_lab <- function(p) ifelse(p < 0.001, sprintf("p = %.1e", p), sprintf("p = %.3f", p))

violin_layers <- list(
  geom_violin(aes(fill = phenotype_clean), alpha = 0.35, colour = NA, width = 0.9),
  geom_jitter(aes(colour = phenotype_clean), width = 0.12, alpha = 0.6, size = 1.8),
  stat_summary(fun = median, geom = "crossbar", width = 0.4, linewidth = 0.8),
  scale_fill_manual(values = pal), scale_colour_manual(values = pal),
  facet_wrap(~ phenotype_clean))

# --- (a) UMKLM: LRRK2 by risk-variant carrier status -------------------------
se  <- readRDS(file.path(PROJECT_DIR, "results_p118a_v2_filter_none", "02_merge_gp2", "se_p118a_Neuro220_with_metadata.rds"))
ids <- variants$GP2ID[variants$Gene == "LRRK2"]
d_l <- as.data.frame(colData(se)) %>%
  filter(phenotype_clean %in% c("PD", "Control"), !is.na(GP2ID), !is.na(ancestry_source)) %>%
  mutate(carrier = factor(ifelse(GP2ID %in% ids, "Carrier", "Non-carrier"), levels = c("Non-carrier", "Carrier")),
         NPQ = assay(se, "npq")["LRRK2", SampleName])
n_lab <- d_l %>% dplyr::count(phenotype_clean, carrier) %>%
  mutate(lab = paste0("n = ", n), y = min(d_l$NPQ) - 0.3)

p_a <- ggplot(d_l, aes(carrier, NPQ)) + violin_layers +
  geom_signif(comparisons = list(c("Non-carrier", "Carrier")), test = "wilcox.test",
              map_signif_level = p_lab, textsize = 6, vjust = -0.2, y_position = max(d_l$NPQ) + 0.3) +
  geom_text(data = n_lab, aes(carrier, y, label = lab), size = 6, colour = "grey30") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(title = "LRRK2 protein by LRRK2 risk-variant carrier status",
       subtitle = "UMKLM  |  G2385R / R1628P carriers vs non-carriers, Wilcoxon",
       x = NULL, y = "plasma LRRK2 (NPQ)") + base
ggsave(file.path(KOREA_DIR, "UMKLM_LRRK2_carriers_violin.png"), p_a, width = 11, height = 7, dpi = 300)

# --- (b) TRAPCAF: GBA by rs3115534 G dose ------------------------------------
d_rs <- read.csv(file.path(KOREA_DIR, "rs3115534_dosage_GBA_p136.csv")) %>%
  mutate(genotype = factor(G_dose, levels = 0:2, labels = c("TT", "GT", "GG")))
b <- coef(summary(lm(NPQ ~ G_dose + phenotype_clean + age + sex_clean, data = d_rs)))["G_dose", ]
n_rs <- d_rs %>% dplyr::count(phenotype_clean, genotype) %>%
  mutate(lab = paste0("n = ", n), y = min(d_rs$NPQ) - 0.15)
top <- max(d_rs$NPQ)

p_b <- ggplot(d_rs, aes(genotype, NPQ)) + violin_layers +
  geom_signif(comparisons = list(c("TT", "GT"), c("TT", "GG")), test = "wilcox.test",
              map_signif_level = p_lab, textsize = 6, vjust = -0.2,
              y_position = c(top + 0.15, top + 0.55), step_increase = 0) +
  geom_text(data = n_rs, aes(genotype, y, label = lab), size = 6, colour = "grey30") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(title = "GBA protein by GBA1 rs3115534 genotype (G = risk allele)",
       subtitle = sprintf("TRAPCAF  |  adjusted \u03b2 = %.2f NPQ per G allele, p = %.1e", b["Estimate"], b["Pr(>|t|)"]),
       x = "rs3115534 genotype", y = "plasma GBA (NPQ)") + base
ggsave(file.path(KOREA_DIR, "TRAPCAF_GBA_rs3115534_violin.png"), p_b, width = 12, height = 7, dpi = 300)


# --- (c) PPMI: LRRK2 by pathogenic-variant carrier status -------------------
se  <- readRDS(file.path(PROJECT_DIR, "results_p_ppmi_v2_filter_none", "02_merge_gp2", "se_p_ppmi_Neuro220_with_metadata.rds"))
ids <- variants$GP2ID[variants$Gene == "LRRK2"]
d_pl <- as.data.frame(colData(se)) %>%
  filter(phenotype_clean %in% c("PD", "Control"), !is.na(GP2ID), !is.na(ancestry_source)) %>%
  mutate(carrier = factor(ifelse(GP2ID %in% ids, "Carrier", "Non-carrier"), levels = c("Non-carrier", "Carrier")),
         NPQ = assay(se, "npq")["LRRK2", SampleName])
n_pl <- d_pl %>% dplyr::count(phenotype_clean, carrier) %>%
  mutate(lab = paste0("n = ", n), y = min(d_pl$NPQ) - 0.3)

p_pl <- ggplot(d_pl, aes(carrier, NPQ)) + violin_layers +
  geom_signif(comparisons = list(c("Non-carrier", "Carrier")), test = "wilcox.test",
              map_signif_level = p_lab, textsize = 6, vjust = -0.2, y_position = max(d_pl$NPQ) + 0.3) +
  geom_text(data = n_pl, aes(carrier, y, label = lab), size = 6, colour = "grey30") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(title = "PPMI: LRRK2 protein by LRRK2 pathogenic-variant carrier status",
       subtitle = "G2019S (179) / R1441G (18) carriers vs non-carriers",
       x = NULL, y = "plasma LRRK2 (NPQ)") + base
ggsave(file.path(KOREA_DIR, "PPMI_LRRK2_carriers_violin.png"), p_pl, width = 11, height = 7, dpi = 300)


# --- (d) PPMI: LRRK2 by pathogenic-variant carrier status, pooled ------------
se  <- readRDS(file.path(PROJECT_DIR, "results_p_ppmi_v2_filter_none", "02_merge_gp2", "se_p_ppmi_Neuro220_with_metadata.rds"))
meta_p <- as.data.frame(colData(se)) %>%
  filter(phenotype_clean %in% c("PD", "Control"), !is.na(GP2ID), !is.na(ancestry_source))

d_pl <- meta_p %>%
  mutate(carrier = factor(ifelse(GP2ID %in% variants$GP2ID[variants$Gene == "LRRK2"], "Carrier", "Non-carrier"),
                          levels = c("Non-carrier", "Carrier")),
         NPQ = assay(se, "npq")["LRRK2", SampleName])
n_pl <- d_pl %>% dplyr::count(carrier) %>% mutate(lab = paste0("n = ", n), y = min(d_pl$NPQ) - 0.3)

p_pl <- ggplot(d_pl, aes(carrier, NPQ)) + grey_layers +
  geom_signif(comparisons = list(c("Non-carrier", "Carrier")), test = "wilcox.test",
              map_signif_level = p_lab, textsize = 6, vjust = -0.2, y_position = max(d_pl$NPQ) + 0.3) +
  geom_text(data = n_pl, aes(carrier, y, label = lab), size = 6, colour = "grey30") +
  labs(title = "PPMI: LRRK2 protein by LRRK2 pathogenic-variant carrier status",
       subtitle = "G2019S (179) / R1441G (18) carriers vs non-carriers; PD and controls pooled",
       x = NULL, y = "plasma LRRK2 (NPQ)") + base +
  theme(plot.title = element_text(size = 18), plot.subtitle = element_text(size = 14))
ggsave(file.path(KOREA_DIR, "PPMI_LRRK2_carriers_violin_pooled.png"), p_pl, width = 8, height = 7, dpi = 300)


# --- (e) PPMI: GBA by GBA1 variant carrier status, pooled --------------------
d_pg <- meta_p %>%
  mutate(carrier = factor(ifelse(GP2ID %in% variants$GP2ID[variants$Gene == "GBA1"], "Carrier", "Non-carrier"),
                          levels = c("Non-carrier", "Carrier")),
         NPQ = assay(se, "npq")["GBA", SampleName])
n_pg <- d_pg %>% dplyr::count(carrier) %>% mutate(lab = paste0("n = ", n), y = min(d_pg$NPQ) - 0.3)

p_pg <- ggplot(d_pg, aes(carrier, NPQ)) + grey_layers +
  geom_signif(comparisons = list(c("Non-carrier", "Carrier")), test = "wilcox.test",
              map_signif_level = p_lab, textsize = 6, vjust = -0.2, y_position = max(d_pg$NPQ) + 0.3) +
  geom_text(data = n_pg, aes(carrier, y, label = lab), size = 6, colour = "grey30") +
  labs(title = "PPMI: GBA protein by GBA1 variant carrier status",
       subtitle = "Carriers predominantly N409S (190 of 227); PD and controls pooled",
       x = NULL, y = "plasma GBA (NPQ)") + base +
  theme(plot.title = element_text(size = 18), plot.subtitle = element_text(size = 14))
ggsave(file.path(KOREA_DIR, "PPMI_GBA_carriers_violin_pooled.png"), p_pg, width = 8, height = 7, dpi = 300)


# --- plot without breakdown by phenotype (simple carrier vs noncarrier) -------------------------------------------------

grey_layers <- list(
  geom_violin(fill = "grey80", alpha = 0.6, colour = NA, width = 0.8),
  geom_jitter(width = 0.12, alpha = 0.5, size = 1.8, colour = "grey30"),
  stat_summary(fun = median, geom = "crossbar", width = 0.35, linewidth = 0.8),
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))))

n_l  <- d_l  %>% dplyr::count(carrier)  %>% mutate(lab = paste0("n = ", n), y = min(d_l$NPQ) - 0.3)
p_a2 <- ggplot(d_l, aes(carrier, NPQ)) + grey_layers +
  geom_signif(comparisons = list(c("Non-carrier", "Carrier")), test = "wilcox.test",
              map_signif_level = p_lab, textsize = 6, vjust = -0.2, y_position = max(d_l$NPQ) + 0.3) +
  geom_text(data = n_l, aes(carrier, y, label = lab), size = 6, colour = "grey30") +
  labs(title = "LRRK2 protein by LRRK2 risk-variant carrier status",
       subtitle = "UMKLM: G2385R / R1628P carriers vs non-carriers",
       x = NULL, y = "plasma LRRK2 (NPQ)") + base
ggsave(file.path(KOREA_DIR, "UMKLM_LRRK2_carriers_violin_pooled.png"), p_a2, width = 8, height = 7, dpi = 300)

n_r  <- d_rs %>% dplyr::count(genotype) %>% mutate(lab = paste0("n = ", n), y = min(d_rs$NPQ) - 0.15)
p_b2 <- ggplot(d_rs, aes(genotype, NPQ)) + grey_layers +
  geom_signif(comparisons = list(c("TT", "GT"), c("TT", "GG")), test = "wilcox.test",
              map_signif_level = p_lab, textsize = 6, vjust = -0.2,
              y_position = c(top + 0.15, top + 0.55), step_increase = 0) +
  geom_text(data = n_r, aes(genotype, y, label = lab), size = 6, colour = "grey30") +
  labs(title = "TRAPCAF: GBA protein by GBA1 rs3115534 genotype (G = risk allele)",
       subtitle = sprintf("Adjusted \u03b2 = %.2f NPQ per G allele, p = %.1e (PD and controls pooled)",
                          b["Estimate"], b["Pr(>|t|)"]),
       x = "rs3115534 genotype", y = "plasma GBA (NPQ)") + base
ggsave(file.path(KOREA_DIR, "TRAPCAF_GBA_rs3115534_violin_pooled.png"), p_b2, width = 9, height = 7, dpi = 300)