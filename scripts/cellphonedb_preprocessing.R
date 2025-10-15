#!/usr/bin/env Rscript

# --------------------------- #
#  CellPhoneDB Gene Mapping   #
# --------------------------- #

# Load required utilities
source("~/research/tools/r/sc-rna-seq_util.R")

# Define input/output paths
input_file  <- "~/LabNGS/resources/cellphonedb/v5.0.0/gene_input.csv"
output_file <- "~/LabNGS/resources/cellphonedb/v5.0.0/gene_orthologue.csv"

# --------------------------- #
#  Load Input                 #
# --------------------------- #

gene_input <- read.csv(input_file, stringsAsFactors = FALSE)

# --------------------------- #
#  Map HGNC ? MGI Symbols     #
# --------------------------- #

# Ensure function mapHumanToMouse() is available from sourced utility script
mapped_genes <- mapHumanToMouse(gene_input$hgnc_symbol)

# --------------------------- #
#  Write Output               #
# --------------------------- #

# Create output DataFrame
gene_orthologue <- data.frame(
  hgnc_symbol = gene_input$hgnc_symbol,
  mgi_symbol  = mapped_genes,
  stringsAsFactors = FALSE
)

# Save to CSV
write.csv(
  gene_orthologue,
  file = output_file,
  quote = FALSE,
  row.names = FALSE
)

cat("? Saved orthologue mapping to:", output_file, "\n")
