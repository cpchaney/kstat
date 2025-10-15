#!/usr/bin/env Rscript

# ----------------------------
# CellPhoneDB Postprocessing
# ----------------------------

library(tidyverse)

# ---- Config ----
time_stamp <- "09_21_2025_093144"
output_dir <- "output/cellphonedb"

# ---- Load Data ----
interaction_scores <- read.table(
  file.path(output_dir, paste0("statistical_analysis_interaction_scores_", time_stamp, ".txt")),
  sep = "\t", header = TRUE
) %>%
  filter(directionality == "Ligand-Receptor") %>%
  select(id_cp_interaction, sender.receiver) %>%
  rename(interaction_score = sender.receiver)

means <- read.table(
  file.path(output_dir, paste0("statistical_analysis_means_", time_stamp, ".txt")),
  sep = "\t", header = TRUE
) %>%
  filter(directionality == "Ligand-Receptor") %>%
  select(id_cp_interaction, sender.receiver) %>%
  rename(mean = sender.receiver)

pvalues <- read.table(
  file.path(output_dir, paste0("statistical_analysis_pvalues_", time_stamp, ".txt")),
  sep = "\t", header = TRUE
) %>%
  filter(directionality == "Ligand-Receptor") %>%
  select(id_cp_interaction, sender.receiver) %>%
  rename(pvalue = sender.receiver)

# ---- Merge Results ----
summary_table <- list(interaction_scores, means, pvalues) %>%
  reduce(full_join, by = "id_cp_interaction") %>%
  arrange(desc(interaction_score), desc(mean), pvalue)

# ---- Add Ligand-Receptor Components ----
component_table <- read.csv("data/ligand_receptor_components.csv") %>%
  group_by(id_cp_interaction, component_type) %>%
  summarise(
    symbol = str_c(mgi_symbol, collapse = "+"),
    .groups = "drop"
  ) %>%
  pivot_wider(
    names_from = component_type,
    values_from = symbol
  ) %>%
  rename(
    ligand = ligand,
    receptor = receptor
  )

# ---- Final Merge ----
summary_table <- full_join(summary_table, component_table, by = "id_cp_interaction") %>%
  mutate(interaction = paste0(ligand, "-", receptor))

# ---- Save Output ----
output_path <- file.path(output_dir, paste0("summary_table_", time_stamp, ".csv"))

write.csv(summary_table, file = output_path, row.names = FALSE, quote = FALSE)

cat("? Summary table written to:", output_path, "\n")
