# kstat
The Kidney Spatial Transcriptome Analysis Tool

We developed the Kidney Spatial Transcriptome Analysis Tool (KSTAT), enabling researchers to identify cell locations, predict cell-cell communication, and map gene pathway activity.

## Overview

KSTAT (Kidney Spatial Transcriptome Analysis Tool) enables integrated analysis of single-nucleus RNA sequencing (snRNA-seq) and spatial transcriptomics data. It supports:

- Cell type localization
- Prediction of cell-cell communication
- Spatial mapping of gene pathway activity

KSTAT combines R and Python components for efficient, modular analysis.

## Project Structure
```text
kstat/
|-- config/ # Configuration files (JSON/YAML)
|-- data/ # Raw and preprocessed datasets
|-- image/ # Output images from analyses
|-- output/ # Intermediate results
|-- python/
| |-- kstat/ # Core Python modules (e.g. preprocessing.py)
| |-- scripts/ # Execution scripts (e.g. run_preprocessing.py)
|-- R/ # R scripts (Seurat workflows, plotting)
|-- sql/ # SQL queries or schemas (if applicable)
|-- docs/ # Markdown tutorials and documentation
|-- README.md # Project overview and quick start
