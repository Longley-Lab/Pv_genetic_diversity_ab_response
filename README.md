# Antibody Diversity and Diagnostic Performance Analysis

## Overview
This project contains two R scripts used to analyze and visualize antibody response diversity against *Plasmodium vivax* antigens. The analyses compare responses across clinical cohorts from Thailand and Brazil, along with multiple negative control sources. These scripts are designed to support figures and supplementary results for a peer-reviewed publication evaluating antigen utility and immune diversity.

## Scripts and Purpose

### 1. `ROC_analysis_pv_diversity.R`
Performs ROC curve analysis to evaluate the **diagnostic performance** of antigen candidates across cohorts.

- **Generates:**
  - ROC curves (`.svg`) for each antigen family
  - Area Under the Curve (AUC) tables with confidence intervals and p-values (`.xlsx`)
- **Key Features:**
  - Standardized, high-resolution ROC plots with clean group-wise legends
  - Group-wise antigen comparison via bootstrap ROC testing
  - Custom color palettes for each antigen family

### 2. `Boxplot_pv_diversity.R`
Generates boxplots to visualize the **distribution and variability** of antibody responses across infection categories.

- **Generates:**
  - One boxplot per antigen (`.png`)
  - A separate, publication-ready color legend (`legend_two_column_layout.png`)
- **Key Features:**
  - Log-transformed data visualized across 12 well-defined infection/control categories
  - Threshold lines, region annotations, and spacing for clear comparison
  - Color scheme matched to ROC plots for consistency

## Input Files (Shared by Both Scripts)

Place the following CSV files in your working directory:

- `Thai_IgG_geneticdiversity.csv` – *Antibody data from Thai cohort*
- `Brazil_IgG_geneticdiversity.csv` – *Antibody data from Brazilian cohort*
- `NC_IgG_geneticdiversity.csv` – *Negative control data (Thai, Brazilian, Australian VBDR)*

These files must include:
- **Column 10**: Time since last *P. vivax* infection (used for categorization)
- **Columns 17–39**: Antibody measurements

## Dependencies

Ensure these R packages are installed before running the scripts:
```r
install.packages(c("ggplot2", "dplyr", "pROC", "ROCR", "MASS", "gridExtra", "writexl", "svglite"))
```

## Usage

1. Open **RStudio** or your R environment.
2. Set your working directory to where the CSVs and scripts are located:
   ```r
   setwd("path/to/your/data")
   ```
3. Source either script:
   ```r
   source("ROC_analysis_pv_diversity.R")
   source("Boxplot_pv_diversity.R")
   ```

## Output Summary

| Script                    | Output Types                          | File Format        |
|--------------------------|----------------------------------------|--------------------|
| `ROC_analysis_pv_diversity.R` | ROC curves per antigen group            | `.svg`             |
|                          | AUC summary tables                     | `.xlsx`            |
| `Boxplot_pv_diversity.R`     | Individual boxplots for each antigen    | `.png`             |
|                          | Combined two-column color legend       | `.png`             |

All outputs are generated in high resolution, suitable for direct use in publication figures.

## Visual Consistency
- Antigen names and color palettes are harmonized across both scripts.
- Font sizes, legend spacing, and annotation placement have been optimized for clarity in manuscripts.