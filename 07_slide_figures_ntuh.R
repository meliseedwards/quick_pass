# =============================================================================
# Script 07: slide figures for the GP2 All Hands
# Date: August 2026
# Description: Reads CSVs and SEs already written by Scripts 01-06. Runs no
#              analysis. Four figures, each answering one question:
#
#   A. NTUH plate layout by phenotype   - why the case/control contrast is limited
#   B. NTUH DA hits vs QC criteria      - which findings survive QC, and why not
#   C. IL17F Q-Q + ICC (PPMI CSF)       - the ICC statistic tracks what the eye sees
#   D. LRRK2 variant breakdown          - what this cohort can and cannot answer
#
# Run:  source("07_slide_figures.R")
# =============================================================================
library(SummarizedExperiment)
library(tidyverse)
library(ggrepel)
library(read_excel)

PROJECT_DIR <- path.expand("~/proteomics")
NTUH_DIR    <- file.path(PROJECT_DIR, "results_ntuh_v2_plate_none")
CSF_DIR     <- file.path(PROJECT_DIR, "results_csf_ppmi_v2")   # plate-adjusted PPMI
DATA_DIR    <- file.path(PROJECT_DIR, "data")
FIG_DIR     <- file.path(PROJECT_DIR, "slide_figures")

LRRK2_PATHWAY <- c("LRRK2", "pLRRK2-S1292", "RAB10", "RAB12", "RAB29",
                   "pRAB10-T73", "pRAB12-S106", "pRAB29-T71")

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

ICC_MAX <- 10; INTRA_CV_MAX <- 10; INTER_CV_MAX <- 15; DET_THRESH <- 0.50

save_fig <- function(p, name, w = 10, h = 7) {
  ggsave(file.path(FIG_DIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
  cat("wrote ", name, ".png\n", sep = "")
}

theme_set(theme_bw(base_size = 12))

short_plate <- function(x) {
  s <- sub(".*[Pp]late[ _-]?([0-9]+).*", "Plate \\1", x)
  ifelse(grepl("^Plate [0-9]+$", s), s, "Plate 1")   # the odd 2026-06-08 ID
}


# --- A. NTUH plate layout by phenotype ---------------------------------------
plate <- read.csv(file.path(NTUH_DIR, "04_qc_report", "qc_plate_Neuro220.csv"))

pl <- plate %>%
  mutate(plate_label = short_plate(PlateID)) %>%
  select(plate_label, PD = n_PD, Control = n_Control) %>%
  pivot_longer(c(PD, Control), names_to = "phenotype", values_to = "n")

pA <- ggplot(pl, aes(plate_label, n, fill = phenotype)) +
  geom_col(width = 0.65) +
  geom_text(aes(label = ifelse(n > 0, n, "")),
            position = position_stack(vjust = 0.5), colour = "white", size = 5) +
  scale_fill_manual(values = c(Control = "tomato", PD = "steelblue")) +
  labs(title = "NTUH samples were not randomised across plates",
       subtitle = paste0("All controls sit on one plate."),
       x = NULL, y = "Samples", fill = NULL) +
  theme(legend.position = "top")

save_fig(pA, "A_ntuh_plate_layout", 7.5, 6)


# --- B. NTUH DA hits against QC criteria -------------------------------------
# Raw ICC is Alamar's published metric. Residualised ICC (ours, not Alamar's)
# removes phenotype/age/sex first - with this confounded design the gap between
# the two shows which ICC flags are the confound rather than batch.
q <- read.csv(file.path(NTUH_DIR, "04_qc_report", "qc_target_Neuro220.csv"),
              stringsAsFactors = FALSE) %>% filter(sig)

flags <- q %>%
  transmute(Target, adj.P.Val,
    `ICC (raw)`        = !is.na(ICC)              & ICC > ICC_MAX,
    `ICC (adjusted)`   = !is.na(ICC_residualised) & ICC_residualised > ICC_MAX,
    `intra-plate CV`   = is.na(intra_cv)          | intra_cv > INTRA_CV_MAX,
    `inter-plate CV`   = is.na(inter_cv)          | inter_cv > INTER_CV_MAX,
    `detectability`    = is.na(detectability)     | detectability < DET_THRESH) %>%
  pivot_longer(-c(Target, adj.P.Val), names_to = "criterion", values_to = "fails") %>%
  mutate(criterion = factor(criterion, levels = c("ICC (raw)", "ICC (adjusted)",
                              "intra-plate CV", "inter-plate CV", "detectability")),
         Target = factor(Target, levels = q$Target[order(-q$adj.P.Val)]))

n_clean <- flags %>% group_by(Target) %>% summarise(any = any(fails), .groups="drop")

pB <- ggplot(flags, aes(criterion, Target, fill = fails)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  scale_fill_manual(values = c("FALSE" = "grey88", "TRUE" = "tomato"),
                    labels = c("passes", "FAILS"), name = NULL) +
  labs(title = "NTUH DA proteins that survive QC",
       subtitle = paste0(nrow(q), " targets at FDR < 0.05; ", sum(!n_clean$any),
                         " pass every criterion."),
       x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "top", panel.grid = element_blank())
save_fig(pB, "B_ntuh_hits_vs_qc", 8, 8)

cat("  ", sum(!n_clean$any), "of", nrow(q), "significant targets pass all QC criteria\n")


# --- C. IL17F Q-Q, one plate highlighted (PPMI CSF Inflammation) -------------
se_inf <- readRDS(file.path(CSF_DIR, "02_merge_gp2",
                            "se_ppmi_csf_Inflammation_with_metadata.rds"))

qc_inf <- read.csv(file.path(CSF_DIR, "04_qc_report", "qc_target_Inflammation.csv"))
TARGET <- "IL17F"
stopifnot(TARGET %in% rownames(se_inf))

d <- data.frame(NPQ = assay(se_inf, "npq")[TARGET, ],
                PlateID = as.character(colData(se_inf)$PlateID)) %>%
  filter(!is.na(NPQ)) %>%
  mutate(theoretical = qnorm(ppoints(n()))[rank(NPQ, ties.method = "first")])

worst <- d %>% group_by(PlateID) %>%
  summarise(delta = median(NPQ) - median(d$NPQ), n = n(), .groups = "drop") %>%
  arrange(desc(abs(delta))) %>% dplyr::slice(1)

d$highlight <- ifelse(d$PlateID == worst$PlateID,
                      paste0(short_plate(worst$PlateID), " (most deviant)"), "all other plates")

qq <- quantile(d$NPQ, c(.25, .75)); xx <- qnorm(c(.25, .75))
sl <- diff(qq) / diff(xx); ic <- qq[1] - sl * xx[1]
r  <- qc_inf %>% filter(Target == TARGET)

pC <- ggplot(d, aes(theoretical, NPQ)) +
  geom_point(aes(colour = highlight, size = highlight, alpha = highlight)) +
  geom_abline(slope = sl, intercept = ic, colour = "red", linewidth = 0.6) +
  scale_colour_manual(values = setNames(c("grey70", "tomato"),
                        c("all other plates", unique(d$highlight[d$PlateID == worst$PlateID]))),
                      name = NULL) +
  scale_size_manual(values = setNames(c(1, 2.2),
                      c("all other plates", unique(d$highlight[d$PlateID == worst$PlateID]))),
                    guide = "none") +
  scale_alpha_manual(values = setNames(c(0.45, 1),
                       c("all other plates", unique(d$highlight[d$PlateID == worst$PlateID]))),
                     guide = "none") +
  labs(title = paste0(TARGET, " (PPMI CSF, Inflammation panel): the QC statistic finds what the eye sees"),
       subtitle = paste0("Plate ICC = ", round(r$ICC, 1), "%"),
       x = "Theoretical quantiles", y = "NPQ") +
  theme(legend.position = "top")

save_fig(pC, "C_il17f_qq_plate_effect", 9, 7)


# --- D. LRRK2 variant breakdown ----------------------------------------------
# What this cohort can and cannot answer, as a figure rather than a caveat.
se_n  <- readRDS(file.path(NTUH_DIR, "02_merge_gp2", "se_ntuh_Neuro220_with_metadata.rds"))
mn    <- as.data.frame(colData(se_n))
ids   <- mn$GP2ID[mn$phenotype_clean %in% c("PD","Control") & mn$ancestry %in% "EAS"]
n_risk <- length(unique(vr$GP2ID[grepl("risk", vr$Pathogenicity, ignore.case = TRUE)]))

vr <- read.csv(file.path(DATA_DIR,
        "variant_report_files_lara_final_version_precision_med_results_release12_release12_variant_report_updated_final.csv"),
        stringsAsFactors = FALSE) %>%
  filter(GP2ID %in% ids, Gene %in% c("LRRK2", "GBA1"))

vp <- vr %>%
  dplyr::count(Gene, Variant_Name, Pathogenicity, Zygosity) %>%
  mutate(label = sub("^(LRRK2|GBA1)_", "", Variant_Name),
         class = ifelse(grepl("risk", Pathogenicity, ignore.case = TRUE),
                        "PD risk variant", "Pathogenic / likely pathogenic"))

pD <- ggplot(vp, aes(reorder(label, n), n, fill = Zygosity)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            colour = "white", size = 4.5) +
  coord_flip() +
  facet_grid(class ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c(het = "steelblue", hom = "navy")) +
  labs(title = "LRRK2 variants in the NTUH proteomics cohort",
       x = NULL, y = "Carriers", fill = NULL) +
  theme(legend.position = "top", strip.text.y = element_text(angle = 0))

save_fig(pD, "D_lrrk2_variant_breakdown", 8, 5)

cat("\nAll figures in", FIG_DIR, "\n")




# --- E. NTUH volcano coloured by QC status -----------------------------------
# Direction is the x-axis, so colour is used for QC instead. One figure that
# says both "what is significant" and "what can be trusted".
q <- read.csv(file.path(NTUH_DIR, "04_qc_report", "qc_target_Neuro220.csv"),
              stringsAsFactors = FALSE) %>%
  mutate(
    fails = (!is.na(ICC) & ICC > ICC_MAX) |
            (is.na(intra_cv) | intra_cv > INTRA_CV_MAX) |
            (is.na(inter_cv) | inter_cv > INTER_CV_MAX) |
            (is.na(detectability) | detectability < DET_THRESH),
    qc = case_when(adj.P.Val >= FDR_THRESHOLD | is.na(adj.P.Val) ~ "Not significant",
                   fails                                        ~ "Significant, fails >=1 QC",
                   TRUE                                         ~ "Significant, passes all QC"))

# label every FDR-significant target, not just the top 15
lab <- q %>% filter(adj.P.Val < FDR_THRESHOLD)

fdr_line <- if (any(q$adj.P.Val < FDR_THRESHOLD, na.rm = TRUE))
              -log10(max(q$P.Value[q$adj.P.Val < FDR_THRESHOLD], na.rm = TRUE)) else NA
xmin <- min(q$logFC, na.rm = TRUE)

pE <- ggplot(q, aes(logFC, -log10(P.Value), colour = qc)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", colour = "grey55") +
  {if (!is.na(fdr_line)) geom_hline(yintercept = fdr_line, linetype = "dashed",
                                    colour = "grey25", linewidth = 0.7)} +
  geom_point(size = 3, alpha = 0.85) +
  geom_label_repel(data = lab, aes(label = Target), colour = "black",
                   fill = alpha("white", 0.75), label.size = NA,
                   size = 4, max.overlaps = Inf, box.padding = 0.6,
                   min.segment.length = 0, segment.colour = "grey50",
                   show.legend = FALSE) +
  annotate("label", x = xmin, y = -log10(0.05), label = "nominal p = 0.05",
           hjust = 0, size = 3.8, colour = "grey45",
           fill = "white", label.size = NA) +
  {if (!is.na(fdr_line)) annotate("label", x = xmin, y = fdr_line,
           label = "FDR = 0.05", hjust = 0, size = 3.8, colour = "grey25",
           fill = "white", label.size = NA)} +
  scale_colour_manual(values = c(
    "Significant, passes all QC" = "#1b7837",
    "Significant, fails >=1 QC"  = "tomato",
    "Not significant"            = "grey80"), name = NULL) +
  labs(title = "NTUH plasma: PD vs Control",
       subtitle = paste0("Positive logFC = higher in PD. Color = FDR significance and QC status.\n",
                         sum(q$adj.P.Val < FDR_THRESHOLD, na.rm = TRUE),
                         " targets at FDR < 0.05; ",
                         sum(q$qc == "Significant, passes all QC"),
                         " pass all QC (ICC, intra/inter CV, detectability).\n",
                         "Plate not adjusted - adjusting leaves only plate 3",
                         "(22 PD vs 60 Control)."),
       x = "log2 fold change (PD vs Control)", y = "-log10(p)") +
  theme_bw(base_size = 15) + theme(legend.position = "top")

save_fig(pE, "E_ntuh_volcano_by_qc", 11, 9)



# --- F. Violins for prespecified and top-hit targets -------------------------

# Faceted violins with each target's LOD marked. For low-detectability targets
# the LOD line shows visually why a p-value there is not interpretable.
se_n <- readRDS(file.path(NTUH_DIR, "02_merge_gp2", "se_ntuh_Neuro220_with_metadata.rds"))
mn   <- as.data.frame(colData(se_n))
keep <- mn$phenotype_clean %in% c("PD","Control") & mn$ancestry %in% "EAS"
se_n <- se_n[, keep]; mn <- as.data.frame(colData(se_n))

det <- as.data.frame(rowData(se_n)) %>% select(Target, detectability)

violin_panel <- function(targets, fname, title, w = 11, h = 8) {
  tg <- intersect(targets, rownames(se_n))
  d <- as.data.frame(assay(se_n, "npq")[tg, , drop = FALSE]) %>%
    rownames_to_column("Target") %>%
    pivot_longer(-Target, names_to = "SampleName", values_to = "NPQ") %>%
    left_join(mn %>% select(SampleName, phenotype_clean), by = "SampleName") %>%
    left_join(det, by = "Target") %>%
    mutate(lab = paste0(Target, "  (", round(detectability * 100), "% detectable)"))

  ggplot(d, aes(phenotype_clean, NPQ, fill = phenotype_clean)) +
    geom_violin(alpha = 0.45, trim = FALSE) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
    geom_jitter(width = 0.08, alpha = 0.35, size = 1) +
    facet_wrap(~ lab, scales = "free_y") +
    scale_fill_manual(values = c(Control = "tomato", PD = "steelblue")) +
    labs(title = title,
         subtitle = "NTUH plasma, Neuro220. Detectability in each panel header.",
         x = NULL, y = "NPQ (log2)") +
    theme_bw(base_size = 13) + theme(legend.position = "none") -> p
  save_fig(p, fname, w, h)
}

violin_panel(LRRK2_PATHWAY, "F1_violin_lrrk2_pathway",
             "LRRK2 pathway targets: PD vs Control")

top_hits <- read.csv(file.path(NTUH_DIR, "04_qc_report", "qc_target_Neuro220.csv")) %>%
  filter(sig) %>% arrange(P.Value) %>% head(9) %>% pull(Target)

violin_panel(top_hits, "F2_violin_top_hits",
             "Top 9 differentially abundant targets: PD vs Control")

passers <- c("DDC","IL16","CD63","LGALS3","TREM2","IL2")

violin_panel(passers, "F3_violin_qc_passers",
             "The 6 DA targets that pass every QC criterion")




# --- F-alt. Same violins, with the limit of detection marked -----------------
# LOD comes from negative-control wells (assay buffer, no sample), so it marks
# the level of assay background. Values below it may be noise rather than protein.
lod_tab <- read_excel(file.path(DATA_DIR, NTUH_NEURO220_Plasma_FILE), sheet = 1)

names(lod_tab)[names(lod_tab) == "targetLOD_NPQ"] <- "LOD"

lod_tab <- lod_raw %>%
  mutate(LOD = as.numeric(LOD)) %>%
  group_by(Target) %>%
  summarise(lod_min = min(LOD, na.rm = TRUE),
            lod_max = max(LOD, na.rm = TRUE), .groups = "drop")

violin_lod <- function(targets, fname, title, w = 11, h = 8) {
  tg <- intersect(targets, rownames(se_n))
  d <- as.data.frame(assay(se_n, "npq")[tg, , drop = FALSE]) %>%
    rownames_to_column("Target") %>%
    pivot_longer(-Target, names_to = "SampleName", values_to = "NPQ") %>%
    left_join(mn %>% select(SampleName, phenotype_clean), by = "SampleName") %>%
    left_join(det,     by = "Target") %>%
    left_join(lod_tab, by = "Target") %>%
    mutate(lab = paste0(Target, "  (", round(detectability * 100), "% detectable)"))

  p <- ggplot(d, aes(phenotype_clean, NPQ, fill = phenotype_clean)) +
    geom_violin(alpha = 0.45, trim = FALSE) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
    geom_jitter(width = 0.08, alpha = 0.35, size = 1) +
    geom_rect(data = d %>% distinct(lab, lod_min, lod_max),
              aes(xmin = -Inf, xmax = Inf, ymin = lod_min, ymax = lod_max),
              inherit.aes = FALSE, fill = "tomato", alpha = 0.15) +
    facet_wrap(~ lab, scales = "free_y") +
    scale_fill_manual(values = c(Control = "tomato", PD = "steelblue")) +
    labs(title = title,
         subtitle = paste0("NTUH plasma, Neuro220. Shaded red area = LOD range",
                           ", computed from\nnegative-control wells.",
                           "Values below it are indistinguishable from background."),
         x = NULL, y = "NPQ (log2)") +
    theme_bw(base_size = 13) + theme(legend.position = "none")
  save_fig(p, fname, w, h)
}

violin_lod(LRRK2_PATHWAY, "F1b_violin_lrrk2_pathway_LOD",
           "LRRK2 pathway targets: PD vs Control")






# --- D2. Carriers by phenotype -----------------------------------------------
# Same variants, split by PD/Control. Shows what a carrier-stratified analysis
# has to work with in each arm.
risk_ids <- vr %>% filter(Gene == "LRRK2",
                          grepl("risk", Pathogenicity, ignore.case = TRUE)) %>%
  pull(GP2ID) %>% unique()

by_var <- mn %>%
  filter(GP2ID %in% ids) %>%
  mutate(group = case_when(
    GP2ID %in% (vr %>% filter(Variant_Name == "LRRK2_Gly2385Arg") %>% pull(GP2ID)) &
    GP2ID %in% (vr %>% filter(Variant_Name == "LRRK2_Arg1628Pro") %>% pull(GP2ID)) ~ "Both risk variants",
    GP2ID %in% (vr %>% filter(Variant_Name == "LRRK2_Gly2385Arg") %>% pull(GP2ID)) ~ "G2385R only",
    GP2ID %in% (vr %>% filter(Variant_Name == "LRRK2_Arg1628Pro") %>% pull(GP2ID)) ~ "R1628P only",
    GP2ID %in% (vr %>% filter(Variant_Name == "LRRK2_Arg1067Gln") %>% pull(GP2ID)) ~ "Arg1067Gln (pathogenic)",
    TRUE ~ "Non-carrier")) %>%
  dplyr::count(group, phenotype_clean) %>%
  mutate(group = factor(group, levels = c("Non-carrier", "G2385R only", "R1628P only",
                                          "Both risk variants", "Arg1067Gln (pathogenic)"))) 

tot <- by_var %>% group_by(group) %>% summarise(n = sum(n), .groups = "drop")

pD2 <- ggplot(by_var, aes(group, n, fill = phenotype_clean)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = n), position = position_stack(vjust = 0.5),
            colour = "white", size = 4.5) +
  geom_text(data = tot, aes(group, n, label = paste0("n=", n)),
            inherit.aes = FALSE, hjust = -0.25, size = 4, colour = "grey30") +
  coord_flip(clip = "off") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  scale_x_discrete(limits = rev(levels(by_var$group))) +
  scale_fill_manual(values = c(Control = "tomato", PD = "steelblue")) +
  labs(title = "LRRK2 carriers by phenotype, NTUH proteomics cohort",
       subtitle = paste0("n = ", length(ids), " samples. Risk-variant carriers: 53 PD / 21 Control.\n",
                         "Arg1067Gln is pathogenic rather than a risk allele and is excluded\n",
                         "from the carrier analysis (n = 2)."),
       x = NULL, y = "Samples", fill = NULL) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top", plot.margin = margin(10, 30, 10, 10))

save_fig(pD2, "D2_carriers_by_phenotype", 8, 5)




# --- G. NTUH cohort demographics ---------------------------------------------
# Four panels for one slide. Structured so the same code works for any cohort.
library(patchwork)   # if not installed: install.packages("patchwork")

dem <- mn %>% filter(GP2ID %in% ids) %>%
  mutate(carrier = ifelse(GP2ID %in% risk_ids, "LRRK2 risk carrier", "Non-carrier"))

bar_count <- function(d, xvar, title) {
  d %>% dplyr::count(.data[[xvar]], phenotype_clean) %>%
    ggplot(aes(.data[[xvar]], n, fill = phenotype_clean)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = n), position = position_stack(vjust = 0.5),
              colour = "white", size = 4) +
    scale_fill_manual(values = c(Control = "tomato", PD = "steelblue")) +
    labs(title = title, x = NULL, y = "Samples", fill = NULL) +
    theme_bw(base_size = 12)
}

p_pheno <- dem %>% dplyr::count(phenotype_clean) %>%
  ggplot(aes(phenotype_clean, n, fill = phenotype_clean)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = n), vjust = -0.4, size = 4.5) +
  scale_fill_manual(values = c(Control = "tomato", PD = "steelblue")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Phenotype", x = NULL, y = "Samples") +
  theme_bw(base_size = 12) + theme(legend.position = "none")

p_sex     <- bar_count(dem, "sex_clean", "Sex") + theme(legend.position = "none")
p_carrier <- bar_count(dem, "carrier",   "LRRK2 risk carrier status") +
  theme(legend.position = "none", axis.text.x = element_text(size = 9))

p_age <- ggplot(dem, aes(age, fill = phenotype_clean)) +
  geom_histogram(binwidth = 5, colour = "white", position = "stack") +
  scale_fill_manual(values = c(Control = "tomato", PD = "steelblue")) +
  labs(title = "Age at sample collection", x = "Age (years)",
       y = "Samples", fill = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "top")

med <- dem %>% group_by(phenotype_clean) %>%
  summarise(n = n(), med_age = median(age, na.rm = TRUE),
            pct_female = round(100 * mean(sex_clean == "Female", na.rm = TRUE)),
            .groups = "drop")

pG <- (p_pheno | p_sex) / (p_carrier | p_age) +
  plot_annotation(
    title = "NTUH proteomics cohort: composition",
    subtitle = paste0("n = ", nrow(dem), " (EAS, Neuro220 plasma). ",
                      "Median age: PD ", med$med_age[med$phenotype_clean == "PD"],
                      ", Control ", med$med_age[med$phenotype_clean == "Control"], ". ",
                      "Female: PD ", med$pct_female[med$phenotype_clean == "PD"],
                      "%, Control ", med$pct_female[med$phenotype_clean == "Control"], "%."),
    theme = theme(plot.title = element_text(size = 16, face = "bold"),
                  plot.subtitle = element_text(size = 12)))
save_fig(pG, "G_ntuh_demographics", 11, 8)

print(med)




# --- H. Age-associated targets by phenotype ----------------------------------
# Four of the six QC-passing targets rise with age in controls. PD are ~10 years
# older, so the age adjustment is load-bearing for these.
age_tg <- c("TREM2", "LGALS3", "CD63", "IL2", "DDC", "IL16")

dh <- as.data.frame(assay(se_n, "npq")[intersect(age_tg, rownames(se_n)), , drop = FALSE]) %>%
  rownames_to_column("Target") %>%
  pivot_longer(-Target, names_to = "SampleName", values_to = "NPQ") %>%
  left_join(mn %>% select(SampleName, age, phenotype_clean), by = "SampleName")

slopes <- dh %>% filter(phenotype_clean == "Control") %>%
  group_by(Target) %>%
  summarise(cf = list(summary(lm(NPQ ~ age))$coef["age", ]), .groups = "drop") %>%
  mutate(slope = sapply(cf, `[`, 1), p = sapply(cf, `[`, 4),
         lab = paste0(Target, "  (controls: ", sprintf("%+.3f", slope),
                      " NPQ/yr, p = ", signif(p, 2), ")")) %>%
  select(Target, lab)

dh <- dh %>% left_join(slopes, by = "Target")

pH <- ggplot(dh, aes(age, NPQ, colour = phenotype_clean, fill = phenotype_clean)) +
  geom_point(alpha = 0.45, size = 1.6) +
  geom_smooth(method = "lm", se = TRUE, alpha = 0.15, linewidth = 0.8) +
  facet_wrap(~ lab, scales = "free_y") +
  scale_colour_manual(values = c(Control = "tomato", PD = "steelblue")) +
  scale_fill_manual(values = c(Control = "tomato", PD = "steelblue")) +
  labs(title = "Age association in the six QC-passing targets",
       subtitle = paste0("NTUH plasma. Slope and p-value in each header are fitted in ",
                         "CONTROLS only.\nPD are ~10 years older than controls ",
                         "(median 74 vs 64), so age is adjusted for in all DA models."),
       x = "Age (years)", y = "NPQ (log2)", colour = NULL, fill = NULL) +
  theme_bw(base_size = 12) + theme(legend.position = "top")

save_fig(pH, "H_age_association", 11, 7)



# --- 8. Do age or plate explain the PD-only hits? ----------------------------
# For each hit, compare the observed carrier effect to what a 7.6-year age gap
# would produce on its own, using the age slope fitted in PD NON-carriers.
mpd <- as.data.frame(colData(se_pd))
age_gap <- mean(mpd$age[mpd$lrrk2_risk == "Carrier"], na.rm = TRUE) -
           mean(mpd$age[mpd$lrrk2_risk == "Non-carrier"], na.rm = TRUE)

hits <- res_pd %>% filter(adj.P.Val < FDR_THRESHOLD) %>% pull(Target)

chk <- lapply(hits, function(t) {
  y  <- assay(se_pd, "npq")[t, ]
  dd <- data.frame(NPQ = y, age = mpd$age, carrier = mpd$lrrk2_risk)
  cf <- summary(lm(NPQ ~ age, data = subset(dd, carrier == "Non-carrier")))$coef["age", ]
  data.frame(Target = t,
             observed = res_pd$logFC[res_pd$Target == t],
             age_predicted = cf[1] * age_gap,
             age_slope_p = cf[4])
}) %>% bind_rows() %>%
  mutate(pct = round(100 * age_predicted / observed),
         verdict = ifelse(sign(age_predicted) == sign(observed),
                          "age acts in the SAME direction",
                          "age acts in the OPPOSITE direction"))

print(chk, row.names = FALSE)

p_age <- chk %>%
  pivot_longer(c(observed, age_predicted), names_to = "what", values_to = "logFC") %>%
  mutate(what = recode(what,
           observed = "observed carrier effect",
           age_predicted = paste0("predicted from the ", round(age_gap, 1),
                                  "-year age gap alone"))) %>%
  ggplot(aes(reorder(Target, logFC), logFC, fill = what)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_hline(yintercept = 0, colour = "grey40") +
  coord_flip() +
  scale_fill_manual(values = c("grey40", "tomato")) +
  labs(title = "Can age explain the carrier effect within PD?",
       subtitle = paste0("PD carriers are ", abs(round(age_gap, 1)),
                         " years younger than PD non-carriers (p = 4e-05).\n",
                         "Age slopes fitted in PD non-carriers only, so they are ",
                         "independent of the contrast.\nWhere the bars point in ",
                         "OPPOSITE directions, age works against the effect rather ",
                         "than creating it."),
       x = NULL, y = "log2 fold change (carrier vs non-carrier)", fill = NULL) +
  theme_bw(base_size = 13) + theme(legend.position = "top")

ggsave(file.path(RESULTS_08_DIR, paste0("age_vs_carrier_effect_", suffix, ".png")),
       p_age, width = 9, height = 6, dpi = 300)










# lrrk2 by case / control adn carrier -----------------------------------------------------
d <- data.frame(NPQ = assay(se_c,"npq")["LRRK2", ],
                carrier = colData(se_c)$lrrk2_risk,
                pheno   = colData(se_c)$phenotype_clean,
                plate   = colData(se_c)$PlateID)

lrrk2_1 <- ggplot(d, aes(carrier, NPQ, fill = carrier)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 1.8) +
  facet_wrap(~ pheno) +
  labs(title = "LRRK2 by carrier status, split by phenotype",
       y = "NPQ", x = NULL) + theme_bw()

save_fig(lrrk2_1, "Lrrk2_by_carrier_status", 8, 5)

# where does the extra control variance come from?
aggregate(NPQ ~ pheno + carrier, d, function(x) c(n = length(x), sd = round(sd(x), 2)))

d$age <- colData(se_c)$age          # your error — d had no age column

d %>% group_by(pheno, carrier) %>%
  summarise(n = n(), mean = round(mean(NPQ),2), median = round(median(NPQ),2),
            max = round(max(NPQ),2), .groups = "drop")

# rank-based test, insensitive to those extremes
wilcox.test(NPQ ~ carrier, data = subset(d, pheno == "Control"))
wilcox.test(NPQ ~ carrier, data = subset(d, pheno == "PD"))

# who are they, and is it a plate thing?
d %>% filter(pheno == "Control", NPQ > 15) %>% select(carrier, NPQ, plate)






# lrrk2 violin plot for Cornelis ----------------------------------------------------------

d <- data.frame(NPQ = assay(se_c,"npq")["LRRK2", ],
                carrier = colData(se_c)$lrrk2_risk,
                pheno   = colData(se_c)$phenotype_clean)

# medians + n per group, and the rank test per stratum
lab <- d %>% group_by(pheno, carrier) %>%
  summarise(n = n(), med = median(NPQ), .groups = "drop")
wtest <- d %>% group_by(pheno) %>%
  summarise(p = wilcox.test(NPQ ~ carrier)$p.value, .groups = "drop") %>%
  mutate(txt = paste0("Wilcoxon p = ", signif(p, 2)))

yr  <- range(d$NPQ)
pad <- diff(yr) * 0.10

ggplot(d, aes(carrier, NPQ, fill = carrier)) +
  geom_violin(alpha = 0.4, trim = FALSE, width = 0.8) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
  geom_jitter(width = 0.08, alpha = 0.5, size = 1.5) +
  geom_text(data = lab, aes(y = yr[1] - pad, label = paste0("n=", n)),
            size = 3.5, colour = "grey30") +
  geom_text(data = wtest, aes(x = 1.5, y = yr[2] + pad * 0.6, label = txt),
            inherit.aes = FALSE, size = 4) +
  facet_wrap(~ pheno) +
  coord_cartesian(ylim = c(yr[1] - pad * 1.4, yr[2] + pad)) +
  scale_fill_manual(values = c("Non-carrier" = "grey75", "Carrier" = "steelblue")) +
  labs(title = "LRRK2 NPQ lower in LRRK2 risk-variant carriers",
       subtitle = "NTUH plasma, Neuro220. Unadjusted for age, sex or plate.",
       y = "LRRK2 NPQ (log2)", x = NULL, fill = NULL) +
  theme_bw() + theme(legend.position = "none")

ggsave(file.path(RESULTS_06_DIR, "lrrk2_violin_by_carrier_phenotype.png"),
       width = 9, height = 6, dpi = 300)





# --- 9. Volcano, PD-only contrast --------------------------------------------
dv <- res_pd %>% left_join(det, by = "Target") %>%
  mutate(dir = case_when(adj.P.Val < FDR_THRESHOLD & logFC > 0 ~ "Up in carriers",
                         adj.P.Val < FDR_THRESHOLD & logFC < 0 ~ "Down in carriers",
                         TRUE ~ "Not significant"))
labv <- dv %>% filter(adj.P.Val < FDR_THRESHOLD | Target %in% LRRK2_PATHWAY)

p_v <- ggplot(dv, aes(logFC, -log10(P.Value), colour = dir)) +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", colour = "grey55") +
  geom_point(size = 2.8, alpha = 0.8) +
  geom_point(data = filter(dv, Target %in% LRRK2_PATHWAY),
             shape = 21, fill = "orange", colour = "black", size = 3.5) +
  geom_label_repel(data = labv, aes(label = Target), colour = "black",
                   fill = alpha("white", 0.75), label.size = NA, size = 4,
                   max.overlaps = Inf, box.padding = 0.6, min.segment.length = 0,
                   show.legend = FALSE) +
  scale_colour_manual(values = c("Up in carriers" = "tomato",
                                 "Down in carriers" = "steelblue",
                                 "Not significant" = "grey80"), name = NULL) +
  labs(title = "Within PD: LRRK2 risk carriers vs non-carriers",
       subtitle = paste0("NTUH plasma. ", res_pd$n_group[1], " carriers vs ",
                         res_pd$n_total[1] - res_pd$n_group[1], " non-carriers. ",
                         sum(dv$adj.P.Val < FDR_THRESHOLD), " targets at FDR < 0.05.\n",
                         "Model: ", trimws(res_pd$model[1]),
                         " | orange ring = LRRK2 pathway target"),
       x = "log2 fold change (carrier vs non-carrier)", y = "-log10(p)") +
  theme_bw(base_size = 14) + theme(legend.position = "top")

ggsave(file.path(RESULTS_08_DIR, paste0("volcano_PDonly_", suffix, ".png")),
       p_v, width = 10, height = 8, dpi = 300)





# --- investigating cohorts matching to R12 for for cohort fig breakdown ------------------------------

f111 <- "proteomics_P111 GP2_P111_BSHRI_NULISAseq_Neuro220_NPQ_06082026.xlsx"
f118 <- "proteomics_P118 GP2_P118_BSHRI_NULISAseq_Neuro220_NPQ_06082026.xlsx"

for (f in c(f111, f118)) {
  cat("\n=====", f, "=====\n")
  print(excel_sheets(file.path(DATA_DIR, f)))
  d <- read_excel(file.path(DATA_DIR, f), sheet = 1)
  cat("cols:", paste(names(d), collapse = ", "), "\n")
  s <- d %>% filter(SampleType == "Sample") %>% distinct(SampleName)
  cat("study samples:", nrow(s), "\n")
  cat("first 8 names:", paste(head(s$SampleName, 8), collapse = ", "), "\n")
  print(table(d$SampleMatrixType, useNA = "ifany"))
  print(table(d$PlateID))
}

mk <- read.delim(file.path(DATA_DIR, MASTER_KEY_FILE), stringsAsFactors = FALSE)
table(mk$study[grepl("ELPD|UMKLM|KUL", mk$study, ignore.case = TRUE)])

d111 <- read_excel(file.path(DATA_DIR, f111), sheet = 1)
ids111 <- d111 %>% filter(SampleType == "Sample", SampleMatrixType == "PLASMA") %>%
  distinct(SampleName) %>% pull(SampleName)
base111 <- sub("_s[0-9]+$", "", ids111)
cat("plasma samples:", length(ids111), "| match R12:", sum(base111 %in% mk$GP2ID), "\n")


d118 <- read_excel(file.path(DATA_DIR, f118), sheet = 1)
ids118 <- d118 %>% filter(SampleType == "Sample", SampleMatrixType == "PLASMA") %>%
  distinct(SampleName) %>% pull(SampleName)
mk2 <- mk %>% filter(study %in% c("UMKLM","KUL"))

cat("P118 plasma samples:", length(ids118), "\n")
for (col in c("clinical_id","sample_id","alternative_id1","alternative_id2", "FID")) {
  cat(sprintf("  %-16s exact %3d | trimmed %3d\n", col,
      sum(ids118 %in% mk2[[col]]),
      sum(trimws(gsub("\\s+"," ",ids118)) %in% trimws(mk2[[col]]))))
}
head(mk2 %>% select(GP2ID, clinical_id, alternative_id1, alternative_id2), 5)


si118 <- read_excel(file.path(DATA_DIR, f118), sheet = "Sample Information")
names(si118); head(si118, 5)




# --- K. Carrier x plate layout, NTUH ------------------------------------------
# Companion to the LRRK2 carrier slide: shows why the unadjusted estimate
# cannot be separated from plate.
pl_carrier <- tribble(
  ~plate,     ~group,         ~n,
  "Plate 1",  "Carrier",      15,
  "Plate 1",  "Non-carrier",   0,
  "Plate 2",  "Carrier",      38,
  "Plate 2",  "Non-carrier",  37,
  "Plate 3",  "Carrier",      21,
  "Plate 3",  "Non-carrier",  59
) %>% mutate(group = factor(group, levels = c("Non-carrier", "Carrier")))

pK <- ggplot(pl_carrier, aes(plate, n, fill = group)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = n), position = position_dodge(width = 0.7),
            vjust = -0.35, size = 4.5) +
  scale_fill_manual(values = c("Non-carrier" = "grey60", "Carrier" = "steelblue")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(title = "Carriers are not evenly spread across plates",
       subtitle = paste0("Plate 1 holds 15 carriers and no non-carriers."),
       x = NULL, y = "Samples", fill = NULL) +
  theme_bw(base_size = 13) +
  theme(legend.position = "top", legend.justification = "left")

save_fig(pK, "K_carrier_plate_layout", 6.5, 5)





# aug 4 - umklm lrrk2 violin plot --------
library(ggplot2); library(dplyr)

se <- readRDS(SE_MERGED_ACTIVE)
m  <- as.data.frame(colData(se))
keep <- m$phenotype_clean %in% c("PD","Control") & m$ancestry %in% "EAS"
se <- se[, keep]; m <- as.data.frame(colData(se))

v <- read.csv(file.path(DATA_DIR, VARIANT_REPORT_FILE), stringsAsFactors = FALSE)
g2385 <- v %>% filter(Variant_Name == "LRRK2_Gly2385Arg") %>% pull(GP2ID)
r1628 <- v %>% filter(Variant_Name == "LRRK2_Arg1628Pro") %>% pull(GP2ID)

d <- data.frame(NPQ = assay(se, "npq")["LRRK2", ], DONOR_ID = m$DONOR_ID) %>%
  mutate(group = case_when(DONOR_ID %in% g2385 ~ "G2385R",
                           DONOR_ID %in% r1628 ~ "R1628P",
                           TRUE ~ "Non-carrier"),
         group = factor(group, levels = c("Non-carrier", "G2385R", "R1628P")))

lab <- d %>% dplyr::count(group) %>% mutate(lab = paste0("n=", n))
yr <- range(d$NPQ); pad <- diff(yr) * 0.10

ggplot(d, aes(group, NPQ, fill = group)) +
  geom_violin(alpha = 0.45, trim = FALSE) +
  geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
  geom_jitter(width = 0.09, alpha = 0.45, size = 1.6) +
  geom_text(data = lab, aes(y = yr[1] - pad, label = lab),
            size = 3.8, colour = "grey30") +
  coord_cartesian(ylim = c(yr[1] - pad * 1.4, yr[2] + pad * 0.4)) +
  scale_fill_manual(values = c("Non-carrier" = "grey75",
                               "G2385R" = "steelblue", "R1628P" = "#7fb3d5")) +
  labs(title = "LRRK2 protein is lower in G2385R carriers",
       subtitle = paste0("UMKLM plasma, Neuro220, EAS. LRRK2 is 99.6% detectable.\n",
                         "Plate-adjusted DA model: G2385R logFC -0.95 (p = 0.029); ",
                         "R1628P -0.28 (p = 0.33).\nUnadjusted for age, sex or plate ",
                         "in this fig."),
       x = NULL, y = "LRRK2 NPQ (log2)") +
  theme_bw(base_size = 13) + theme(legend.position = "none")

ggsave("~/proteomics/slide_figures/umklm_lrrk2_by_variant.png", width = 8, height = 6, dpi = 300)
