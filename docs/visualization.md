<!-- Part of the KSTAT Tutorial. Return to [Main Tutorial](tutorial.md) -->
# 3. Visualization

This section shows how to generate spatial plots of measured and inferred gene, regulon, and pathway activity.

---

## 3.1. Visualizing Measured Landmark Expression 
To visualize the spatial distribution of a landmark gene, we use the spatial count matrix and overlay the signal on the tissue background. 

--- 

### 3.1.1 Script: `mouse_e18.5d_kidney_visualization.R` 
This script loads the precomputed sparse landmark matrix and plots the expression of a selected gene using `ggplot2`. 
```r 
# Load configuration and libraries 
stage <- "e18.5d" 
config <- fromJSON(file = "config/mouse_e18.5d.json") 
landmark_counts <- Matrix::readMM("data/mouse_e18.5d_landmarks_matrix.mtx") 
rownames(landmark_counts) <- config$landmarks  

# Select gene to plot 
gene_of_interest <- "Mgp" 
signal_vector <- landmark_counts[gene_of_interest, ] 

# Convert to spatial coordinates 
signal_table <- getSignalCoordinates( channel_name = gene_of_interest, v = signal_vector, height = config$scaled_height, width = config$scaled_width ) 

# Generate spatial plot 
getSpatialPlot( signal_table, width = config$scaled_width, height = config$scaled_height, point_size = 0.5, channel_color_scale = setNames("red", gene_of_interest), background_file = paste0("image/", config$experiment_name, "_background_scaled.jpg"), plot_legend = FALSE ) 
``` 
--- 

### 3.1.2 Output 
This will generate a 2D spatial plot of the landmark **Mgp**, highlighting all spatial bins where it is detected. 
- The gene is visualized as red dots 
- A tissue background is overlaid from the corresponding scaled image 
- Plot aspect ratio and coloring are optimized for clarity

### 3.1.3 Saving the Plot 
To export the spatial plot to a high-resolution TIFF file for use in figures or publications, use `ggsave()`: 
```r 
output_file <- here::here(
  "image",
  paste0(
    config$experiment_name,
    "_",
    stringr::str_to_lower(gene_of_interest),
    "_measured_expression.tiff"
  )
) 
ggsave( 
	filename = output_file, 
	plot = combined_landmark_plot, 
	device = "tiff", 
	height = 11, 
	width = 8.5, 
	units = "in", 
	dpi = 300, 
	compression = "lzw" 
) 
``` 
This will save the plot to the `image/` directory with a resolution of 300 DPI, suitable for print or publication.

## 3.2. Visualizing Inferred Spatial Expression
KSTAT infers the spatial expression of genes across the tissue by projecting single-cell profiles onto spatial bins using the model output. This section shows how to visualize the inferred expression of a gene using the precomputed probabilities. 

--- 
### 3.2.1 Loading Inferred Data 
We begin by loading the cell-bin probability matrix and the unique bin index key used to map dense spatial vectors to valid coordinates. 
```r 
cells_bins_probabilities <- Matrix::readMM( 
	"data/mouse_e18.5d_kidney_cells_bins_probabilities.mtx" 
) 
unique_bins_key <- read.csv( 
	"data/mouse_e18.5d_kidney_unique_bins_key.csv" 
)$index + 1 
``` 

--- 

### 3.2.2 Projecting a Gene Onto the Tissue 
To visualize a gene (e.g., `Cspg4`), we extract its imputed expression vector and use the `getSpatialSignal()` function to compute its inferred spatial distribution. 
```r 
gene_of_interest <- "Cspg4" 
spatial_signal <- getSpatialSignal( 
	as.numeric(expression_matrix[gene_of_interest, , drop = FALSE]), 
	cells_bins_probabilities, 
	mapped.bins = NULL, 
	log = FALSE 
) 
``` 

--- 

### 3.2.3 Plotting the Inferred Signal 
We then use `plotThresholdedSignal()` to plot a denoised and thresholded version of the inferred signal over the tissue background. 
```r 
inferred_plot <- plotThresholdedSignal( 
	spatial_signal = spatial_signal[unique_bins_key], 
	signal_identifier = gene_of_interest, 
	height = config$scaled_height, 
	width = config$scaled_width, 
	point_size = 0.5, 
	denoise = TRUE, 
	percentile = 96, 
	background_file = "image/mouse_e18.5d_kidney_background_scaled.jpg", 
	plot_legend = FALSE 
) 
print(inferred_plot) 
``` 

--- 

### 3.2.4 Saving the Plot 
To export the inferred expression map for this gene: 
```r 
output_file <- paste0( 
	"image/mouse_e18.5d_kidney_", 
	tolower(gene_of_interest), 
	"_inferred_expression.tiff" 
)
ggsave( 
	filename = output_file, 
	plot = inferred_plot + theme(text = element_text(family = "Arial")), 
	device = "tiff", 
	height = 11, 
	width = 8.5, 
	units = "in", 
	dpi = 300, 
	compression = "lzw" 
) 
``` 
The result is a publication-ready TIFF file saved to the `image/` directory.

## 3.3. Visualizing Inferred Regulon Activity 
In addition to gene expression, KSTAT supports the spatial projection of **regulon activity**, such as transcription factor target gene modules. This section describes how to visualize the inferred spatial distribution of a regulon's activity. 

--- 

### 3.3.1 Load the Regulon Activity Matrix 
The `regulons_auc_matrix` is typically stored in the metadata of a `SingleCellExperiment` object, and represents the AUC scores for regulon activation per cell. 
```r 
regulons_auc_matrix <- t(metadata(sce)[["regulons_auc_matrix"]])[, colnames(sce)] 
``` 

--- 

### 3.3.2 Projecting a Regulon (e.g. `Lef1`) 
We select a regulon of interest and compute its spatial distribution by projecting it onto spatial bins using the inferred cell-bin probabilities. 
```r 
regulon_of_interest <- "Lef1" 
spatial_signal <- getSpatialSignal( 
	as.numeric(regulons_auc_matrix[regulon_of_interest, , drop = FALSE]), 
	cells_bins_probabilities, 
	mapped.bins = NULL, 
	log = FALSE ) 
``` 

--- 

### 3.3.3 Visualizing the Activity 
We generate a denoised, thresholded spatial plot of the projected regulon signal: 
```r 
p <- plotThresholdedSignal( 
	spatial_signal = spatial_signal[unique_bins_key], 
	signal_identifier = regulon_of_interest, 
	height = config$scaled_height, 
	width = config$scaled_width, 
	point_size = 0.5, 
	denoise = TRUE, 
	percentile = 96, 
	background_file = paste0("image/", config$experiment_name, "_background_scaled.jpg"), 
	plot_legend = FALSE 
) 
print(p) 
``` 

--- 

### 3.3.4 Saving the Plot 
The plot is saved to the `image/` directory using a descriptive filename: 
```r 
output_file <- paste0( 
	"image/", 
	config$experiment_name, 
	"_", 
	tolower(regulon_of_interest), 
	"_inferred_activity.tiff" 
) 
ggsave( 
	filename = output_file, 
	plot = p + theme(text = element_text(family = "Arial")), 
	device = "tiff", 
	height = 11, 
	width = 8.5, 
	units = "in", 
	dpi = 300, 
	compression = "lzw" 
) 
```
## 3.4. Visualizing Inferred Gene Set Enrichment 
KSTAT enables the visualization of spatial enrichment for gene sets, such as those from MSigDB (e.g. HALLMARK pathways). These are inferred from single-cell profiles and mapped to the tissue using the same probability-based projection used for genes and regulons. 

--- 

### 3.4.1 Load Gene Set AUC Matrix 
The gene set enrichment scores (AUC values) are typically stored in the metadata of the `SingleCellExperiment` object: 
```r 
gene_sets_auc_matrix <- metadata(sce)[["gene_sets_auc_matrix"]][, colnames(sce)] 
``` 

--- 

### 3.4.2 Select a Gene Set 
Here we visualize the **G2M Checkpoint** signature from the MSigDB Hallmark collection: 
```r 
gene_set_of_interest <- "HALLMARK_G2M_CHECKPOINT" 
``` 

---

### 3.4.3 Compute the Spatial Projection 
```r 
spatial_signal <- getSpatialSignal( 
	as.numeric(gene_sets_auc_matrix[gene_set_of_interest, , drop = FALSE]), 
	cells_bins_probabilities, 
	mapped.bins = NULL, 
	log = FALSE 
) 
``` 

--- 

### 3.4.4 Generate the Plot 
```r 
p <- plotThresholdedSignal( 
	spatial_signal = spatial_signal[unique_bins_key], 
	signal_identifier = gene_set_of_interest, 
	height = config$scaled_height, w
	idth = config$scaled_width, 
	point_size = 0.5, 
	denoise = TRUE, 
	percentile = 98, 
	background_file = paste0("image/", config$experiment_name, "_background_scaled.jpg"), 
	plot_legend = FALSE 
) 
print(p) 
``` 

--- 

### 3.4.5 Save the Plot 
We generate a simplified filename to avoid long or invalid file names, and export the image at native resolution
```r 
gene_set_label <- tolower( 
	paste(substr(strsplit(gene_set_of_interest, "_")[[1]], 1, 3), collapse = "_") 
) 
output_file <- paste0( 
	"image/", 
	config$experiment_name, 
	"_", 
	gene_set_label, 
	"_inferred_activity.tiff" 
) 
ggsave( 
	filename = output_file, 
	plot = p + theme(text = element_text(family = "Arial")), 
	device = "tiff", 
	height = config$scaled_height, 
	width = config$scaled_width, 
	units = "px", 
	dpi = 300, 
	compression = "lzw"
) 
```
