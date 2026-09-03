# =============================================================================
# Bridging sample reproducibility across plates (P118)
# The same 30 CSF samples were run on Plate 1 and Plate 5. Each point is one
# target in one sample, measured twice. Points on the line agree perfectly.
# =============================================================================

library(tidyverse)
library(readxl)

source("~/proteomics/nulisa_pipeline_v2/quick_pass/00_config.R")

# P118 file is the same for umklm and kul
npq <- read_excel(file.path(DATA_DIR, UMKLM_NEURO220_Plasma_FILE),
                  sheet = 1, na = "NA")

# CSF wells are named CSF11_P0118_Plate1_N220.xml / CSF11_P118_Plate5_N220.xml
csf <- npq %>%
  dplyr::filter(SampleMatrixType == "CSF") %>%
  mutate(sid   = str_extract(SampleName, "^CSF[0-9]+"),
         plate = case_when(
           PlateID == "P0118_Plate1_N220.xml" ~ "p1",
           PlateID == "P118_Plate5_N220.xml"  ~ "p5",
           TRUE ~ NA_character_)) %>%
  dplyr::filter(!is.na(plate))

stopifnot(!any(duplicated(csf[c("sid", "Target", "plate")])))

cat("CSF samples:", dplyr::n_distinct(csf$sid),
    "| plates:", paste(unique(csf$PlateID), collapse = ", "), "\n")

pairs <- csf %>%
  mutate(NPQ = as.numeric(NPQ)) %>%
  select(sid, Target, plate, NPQ) %>%
  pivot_wider(names_from = plate, values_from = NPQ) %>%
  drop_na(p1, p5)

cat("Paired measurements:", nrow(pairs),
    "|", dplyr::n_distinct(pairs$sid), "samples x",
    dplyr::n_distinct(pairs$Target), "targets\n")

# --- summary stats -----------------------------------------------------------
r_all <- cor(pairs$p1, pairs$p5)

per_sample <- pairs %>%
  group_by(sid) %>%
  summarise(r = cor(p1, p5), .groups = "drop")

diffs <- abs(pairs$p1 - pairs$p5)

cat("\nPearson r (all):", round(r_all, 4), "\n")
cat("Per-sample r: min", round(min(per_sample$r), 4),
    "median", round(median(per_sample$r), 4),
    "max", round(max(per_sample$r), 4), "\n")
cat("|difference| NPQ: median", round(median(diffs), 3),
    "| 75th", round(quantile(diffs, .75), 3),
    "| 95th", round(quantile(diffs, .95), 3), "\n")
cat("Mean SIGNED difference:", round(mean(pairs$p1 - pairs$p5), 4),
    "(near zero = no systematic plate shift)\n")

# worst-reproducing targets, so the tail is not hidden by the correlation
worst <- pairs %>%
  group_by(Target) %>%
  summarise(median_diff = median(abs(p1 - p5)), .groups = "drop") %>%
  arrange(desc(median_diff))

cat("\nWorst 5 targets by median |difference|:\n")
print(head(worst, 5))

# --- plot --------------------------------------------------------------------
lims <- range(c(pairs$p1, pairs$p5)) + c(-0.5, 0.5)

caption <- sprintf(
  "%s CSF bridging samples x %s targets = %s paired measurements\nPearson r = %.3f   |   per-sample r: %.3f-%.3f\nmedian |difference| = %.2f NPQ   |   dashed line = perfect agreement",
  dplyr::n_distinct(pairs$sid), dplyr::n_distinct(pairs$Target),
  format(nrow(pairs), big.mark = ","),
  r_all, min(per_sample$r), max(per_sample$r), median(diffs))

p <- ggplot(pairs, aes(p1, p5)) +
  geom_point(alpha = 0.28, size = 1, colour = "#526A83") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "#AF6458", linewidth = 0.6) +
  annotate("label", x = lims[1], y = lims[2], label = caption,
           hjust = 0, vjust = 1, size = 3.2, label.size = 0.2,
           fill = "white", colour = "grey20") +
  coord_equal(xlim = lims, ylim = lims) +
  labs(x = "NPQ on Plate 1", y = "NPQ on Plate 5",
       title = "Bridging samples reproduce across plates (P118)") +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(file.path(SC_DRIFT_DIR, "bridging_reproducibility_P118.png"),
       p, width = 7.5, height = 7, dpi = 300)

cat("\nSaved bridging_reproducibility_P118.png\n")
p