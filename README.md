# antibody response analysis against different *Plasmodium vivax* antigen haplotypes

## Overview
This project contains two R scripts used to analyze and visualize antibody response against different *Plasmodium vivax* antigen haplotypes.

## Scripts and Purpose

### 1. `Boxplot_pv_diversity.R`
Generates boxplots to visualize the **distribution and variability** of antibody responses across infection categories.

- **Generates:**
  - One boxplot per antigen (`.png`)
  - A separate color legend (`legend_two_column_layout.png`)

### 2. `ROC_analysis_pv_diversity.R`
Performs ROC curve analysis to evaluate the **diagnostic performance** of antigen candidates across cohorts.

- **Generates:**
  - ROC curves (`.svg`) for each antigen family
  - Area Under the Curve (AUC) tables with confidence intervals and p-values (`.xlsx`)


## Input Files (Shared by Both Scripts)

Place the following CSV files in your working directory:

- `Thai_IgG_geneticdiversity.csv` – *Antibody data from Thai cohort*
- `Brazil_IgG_geneticdiversity.csv` – *Antibody data from Brazilian cohort*
- `NC_IgG_geneticdiversity.csv` – *Negative control antibody data (Thai Red Cross, Brazilian Red Cross, Australian Red Cross, Australian VBDR)*

These files include:

- Epidemiological information previously collected from two year-long observational cohort studies in Thailand [https://doi.org/10.1016/j.ijpara.2019.01.004.] 
  and Brazil [https://doi.org/10.1590/0074-02760210330]. (column 1-16)
- Epidemiological information of negative control sample (column 1-16)
- Antibody response against different _P. vivax_ antigen haplotypes (column 17-39)


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
| `Boxplot_pv_diversity.R`     | Individual boxplots for each antigen groups   | `.png`             |
|                          | Combined two-column color legend       | `.png`             |
| `ROC_analysis_pv_diversity.R` | ROC curves per antigen group            | `.svg`             |
|                          | AUC stat summary tables                     | `.xlsx`            |




