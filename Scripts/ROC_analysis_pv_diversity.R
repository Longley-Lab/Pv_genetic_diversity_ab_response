# Load required libraries
library(MASS)
library(ROCR)
library(pROC)
library(ggplot2)
library(svglite)
library(writexl)

# ========== CONFIGURATION ========== #

# Set working directories
setwd('')

data_sources <- list(
  Thai = list(
    data_path = "Input_Files/Thai_IgG_geneticdiversity.csv",
    label = "Thai"
  ),
  Brazil = list(
    data_path = "Input_Files/Brazil_IgG_geneticdiversity.csv",
    label = "Brazilian"
  )
)

control_path <- "Input_Files/NC_IgG_geneticdiversity.csv"

# Create 'output' folder if it doesn't exist
if (!dir.exists("output")) {
  dir.create("output")
  cat("Created 'output' folder\n")
} else {
  cat("'output' folder already exists\n")
}


# Antigen group definitions
# Safe names for internal use
antigen_groups <- list(
  DBP    = c(8, 5, 4, 10),
  MSP5   = c(23, 15, 14),
  MSP8   = c(1, 3),
  RBP2A  = c(12, 9, 13),
  RBP2B  = c(21, 18, 19, 20),
  RIPR   = c(22, 16, 17),
  PTEX   = c(6, 11),
  fam_a  = c(7, 2)
)

# Clean display titles (e.g., for plot titles and file names)
antigen_titles <- list(
  DBP    = "DBPII",
  MSP5   = "MSP5",
  MSP8   = "MSP8",
  RBP2A  = "RBP2a",
  RBP2B  = "RBP2b",
  RIPR   = "RIPR",
  PTEX   = "PTEX150",
  fam_a  = "Pv-fam-a"
)

colors_list <- list(
  DBP    = c("#1f77b4", "#2a9df4", "#6baed6", "#a6c8e0"),
  MSP5   = c("#ff7f0e", "#ffa64d", "#cc6600"),
  MSP8   = c("#2ca02c", "#66c266"),
  RBP2A  = c("#d62728", "#e65555", "#a31d1d"),
  RBP2B  = c("#9467bd", "#af5fbf", "#684b96", "#7e57c2"),
  RIPR   = c("#8c564b", "#c97a6d", "#5f3a32"),
  PTEX   = c("#e377c2", "#b3509e"),
  fam_a  = c("#7f7f7f", "#4d4d4d")
)

ant_names_short <- c("MSP8 L34 Sal1","Pv-fam-a SEM-8 Hap1", "MSP8 SEM-10 Hap1", 
                     "DBPII SEM-35 Hap1", "DBPII AH", "PTEX150 L18 Sal1", 
                     "Pv-fam-a L02 Sal1", "DBPII Sal1", "RBP2a SEM-28 Hap1", 
                     "DBPII SEM-37 Hap2", "PTEX150 SEM-39 Hap1", "RBP2a Sal1", 
                     "RBP2a SEM-29 Hap9","MSP5 SEM-33", "MSP5 SEM-31", "RIPR SEM-41",
                     "RIPR SEM-42", "RBP2b SEM-2", "RBP2b SEM-4", "RBP2b SEM-6",
                     "RBP2b P25", "RIPR", "MSP5 L19 Sal1")  # same as before
ant_names_hap   <- c("Sal-1","Hap1", "Hap1", 
                     "Hap1", "AH", "Sal-1", 
                     "Sal-1", "Sal-1", "Hap1", 
                     "Hap2", "Hap1", "Sal-1", 
                     "Hap9","Hap2", "Hap1", "Hap1",
                     "Hap2", "Hap1", "Hap2", "Hap3",
                     "Sal-1", "Sal-1", "Sal-1")  # same as before


# ========== FUNCTION DEFINITIONS ========== #

categorize_data <- function(data, is_control = FALSE) {
  if (!is_control) {
    cat <- rep("brazil_never", nrow(data))
    cat[data[, 10] == 0] <- "brazil_current"
    cat[data[, 10] > 0 & data[, 10] <= 9 * 30] <- "brazil_recent"
    cat[data[, 10] > 9 * 30] <- "brazil_old"
  } else {
    cat <- rep("VBDR", nrow(data))
    prefixes <- substr(data[, 1], 1, 3)
    cat[prefixes == "TRC"] <- "TRC"
    cat[prefixes == "BRC"] <- "BRC"
    cat[prefixes == "ARC"] <- "ARC"
  }
  return(cat)
}

prepare_AB_matrix <- function(dataset, control, ab_cols = 17:39) {
  rbind(dataset[, ab_cols], control[, ab_cols]) |> 
    log()
}

generate_bin_cat <- function(categorized_vector) {
  cat_bin <- rep("old", length(categorized_vector))
  cat_bin[categorized_vector %in% c("brazil_current", "brazil_recent")] <- "new"
  return(ifelse(cat_bin == "old", 0, 1))
}

run_ROC_analysis <- function(AB_matrix, bin_vector) {
  AB_matrix$bin_cat_num <- bin_vector
  roc_results <- list()
  for (i in 1:(ncol(AB_matrix) - 1)) {
    df <- data.frame(sero = AB_matrix[, i], pcr = AB_matrix$bin_cat_num) |> na.omit()
    roc_results[[paste0("ROC_Column_", i)]] <- roc(df$pcr, df$sero)
  }
  return(roc_results)
}

calculate_AUC <- function(indices, roc_results) {
  base_roc <- roc_results[[paste0("ROC_Column_", indices[1])]]
  auc_df <- data.frame(
    antigen_name = ant_names_short[indices[1]],
    AUC = auc(base_roc),
    CI_lower = ci.auc(base_roc)[1],
    CI_upper = ci.auc(base_roc)[3],
    p_value = NA
  )
  
  for (i in indices[-1]) {
    roc_i <- roc_results[[paste0("ROC_Column_", i)]]
    test <- roc.test(base_roc, roc_i, method = "bootstrap", alternative = "two.sided")
    ci_i <- ci.auc(roc_i)
    auc_df <- rbind(auc_df, data.frame(
      antigen_name = ant_names_short[i],
      AUC = auc(roc_i),
      CI_lower = ci_i[1],
      CI_upper = ci_i[3],
      p_value = test$p.value
    ))
  }
  return(auc_df)
}

plot_ROC <- function(indices, colors, title, filename, roc_results) {
  all_roc_data_list <- list()
  
  for (i in indices) {
    roc_result <- roc_results[[paste("ROC_Column", i, sep = "_")]]
    roc_plot_df <- data.frame(
      Specificity = 1 - roc_result$specificities,
      Sensitivity = roc_result$sensitivities,
      legend = ant_names_hap[i]
    )
    all_roc_data_list[[length(all_roc_data_list) + 1]] <- roc_plot_df
  }
  
  all_roc_data <- do.call(rbind, all_roc_data_list)
  all_roc_data$legend <- factor(all_roc_data$legend, levels = unique(all_roc_data$legend))
  
  plot <- ggplot(all_roc_data, aes(x = Specificity, y = Sensitivity, color = legend)) +
    geom_line(size = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    labs(title = title, x = "1 - Specificity", y = "Sensitivity", color = "Legend") +
    scale_color_manual(values = colors) +
    theme_minimal() +
    coord_cartesian(clip = "off") +
    theme(
      legend.position = c(1.0, 0.25),
      legend.justification = "left",
      legend.title = element_text(size = 30, face = "bold"),
      legend.text = element_text(size = 28),
      legend.key.size = unit(0.5, "cm"),
      legend.key.height = unit(1, "cm"),
      legend.spacing.y = unit(2, "cm"),
      plot.title = element_text(hjust = 0.5, size = 34, face = "bold"),
      panel.grid.major = element_line(color = "gray90"),
      plot.margin = margin(1, 5, 1, 1, "cm"),
      axis.title.x = element_text(size = 28),
      axis.title.y = element_text(size = 28),
      axis.text.x = element_text(size = 24),    
      axis.text.y = element_text(size = 24) 
    ) +
    guides(color = guide_legend(byrow = TRUE))
  
  print(plot)
  ggsave(filename, plot = plot, width = 8, height = 6, dpi = 600)
}


# ========== RUN ANALYSIS ========== #

control_data <- read.csv(control_path)
control_cat <- categorize_data(control_data, TRUE)

for (region in names(data_sources)) {
  message("Processing ", region)
  
  dataset <- read.csv(data_sources[[region]]$data_path)
  cat_vector <- categorize_data(dataset)
  inf_cat <- c(cat_vector, control_cat)
  
  AB <- prepare_AB_matrix(dataset, control_data)
  bin_cat_num <- generate_bin_cat(inf_cat)
  
  if (length(bin_cat_num) != nrow(AB)) stop("Mismatch between categories and data matrix")
  
  roc_results <- run_ROC_analysis(AB, bin_cat_num)
  
  for (group in names(antigen_groups)) {
    indices <- antigen_groups[[group]]
    colors <- colors_list[[group]]
    title <- paste(data_sources[[region]]$label, antigen_titles[[group]])
    filename <- file.path("output", paste0("ROC_", region, "_", group, ".svg"))
    AUC_table <- calculate_AUC(indices, roc_results)
    write_xlsx(AUC_table, file.path("output", paste0("AUC_", region, "_", group, ".xlsx")))
    plot_ROC(indices, colors, title, filename, roc_results)
  }
}
