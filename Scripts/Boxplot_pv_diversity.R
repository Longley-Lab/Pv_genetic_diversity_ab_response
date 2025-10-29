library(ggplot2)
library(dplyr)
library(MASS)
library(ROCR)
library(gridExtra)
library(lme4)
library(lmerTest)
library(emmeans)
library(tidyr)
library(stringr)
library(purrr)



###################################
###################################
##                               ##
##  ####    ####  ######  ####   ##
##  ## ##  ##  ##   ##   ##  ##  ##
##  ##  ## ######   ##   ######  ##
##  ## ##  ##  ##   ##   ##  ##  ##
##  ####   ##  ##   ##   ##  ##  ##
##                               ##
###################################
###################################

setwd('C:/Users/wenco/Desktop/Pv_genetic_diversity_ab_response-main')

Thai_data = read.csv("Input_Files/Thai_IgG_geneticdiversity.csv")

brazil_data = read.csv("Input_Files/Brazil_IgG_geneticdiversity.csv")

control_data = read.csv("Input_Files/NC_IgG_geneticdiversity.csv")


# Create 'output' folder if it doesn't exist
if (!dir.exists("output")) {
  dir.create("output")
  cat("Created 'output' folder\n")
} else {
  cat("'output' folder already exists\n")
}



###################################
###################################
##                               ##
##  CATEGORISATION               ##
##                               ##
###################################
###################################

###################################
## Brazil
# 10 as column 10 has the time since last Pv infection
# can change here to 6*30 instead of 9*30, or 3*30 for 3 months, or 1*30 for 1 month

Thai_cat <- rep("Thai_never", nrow(Thai_data) )

Thai_cat[which( Thai_data[,10] == 0 )] <- "Thai_current"

Thai_cat[intersect( which(Thai_data[,10]>0), which(Thai_data[,10]<=9*30) )] <- "Thai_recent"

Thai_cat[which(Thai_data[,10]>9*30)] <- "Thai_old"




###################################
## Brazil
# 10 as column 10 has the time since last Pv infection
# can change here to 6*30 instead of 9*30, or 3*30 for 3 months, or 1*30 for 1 month

brazil_cat <- rep("brazil_never", nrow(brazil_data) )

brazil_cat[which( brazil_data[,10] == 0 )] <- "brazil_current"

brazil_cat[intersect( which(brazil_data[,10]>0), which(brazil_data[,10]<=9*30) )] <- "brazil_recent"

brazil_cat[which(brazil_data[,10]>9*30)] <- "brazil_old"




###################################
## Controls

control_cat <- rep("VBDR", nrow(control_data) )

for(i in 1:length(control_cat))
{
  if( substr( control_data[i,1], 1, 3 ) == "TRC" )
  {
    control_cat[i] <- "TRC"
  }
  if( substr( control_data[i,1], 1, 3 ) == "BRC" )
  {
    control_cat[i] <- "BRC"
  }
  if( substr( control_data[i,1], 1, 3 ) == "ARC" )
  {
    control_cat[i] <- "ARC"
  }
}

####################################
## Put together Brazil and control data
##in this bit the numbers refer to the columns where the antibody data is...so if you add to 81 (starting at 82) it will be [,17:95(example)]
#control data and brazil data need to have the same number of columns - each antibody column needs to line up in the columns
# 150 is it will drop ab measurments when more than 150 samples missing for that protein

AB <- rbind( Thai_data[,17:39], brazil_data[,17:39], control_data[,17:39] )

AB <- log(AB)

#ant_drop <- c()

#for(j in 1:ncol(AB))
#{
#  if( length(which(is.na(AB[,j]))) > 150 )
#  {
#    ant_drop = c(ant_drop, j)
#  }
#}

#AB = AB[,-ant_drop]

N_ant <- ncol(AB)
ant_names <- colnames(AB)

View(N_ant)
View(ant_names)

############################################
## Create shortened antibody names
#either add new shortened names, or replace shortened names...

ant_names_short = c("MSP8 L34 sal1","Pv-fam-a SEM-8 Hap1", "MSP8 SEM-10 Hap1", 
                    "DBPII SEM-35 Hap1", "DBPII AH", "PTEX150 L18 sal1", 
                    "Pv-fam-a L02 sal1", "Pv DBPII sal1", "RBP2a SEM-28 Hap1", 
                    "DBPII SEM-37 Hap2", "PTEX150 SEM-39 Hap1", "RBP2a sal1", 
                    "RBP2a SEM-29 Hap9","MSP5 SEM-33", "MSP5 SEM-31", "RIPR SEM-41",
                    "RIPR SEM-42", "RBP2b SEM-2", "RBP2b SEM-4", "RBP2b SEM-6",
                    "RBP2b P25", "Pv RIPR", "MSP5 L19 sal1")

ant_names_hap = c("MSP8 Sal-1","Pv-fam-a Hap1", "MSP8 Hap1", 
                  "DBPII Hap1", "DBPII AH", "PTEX150 Sal-1", 
                  "Pv-fam-a Sal-1", "DBPII Sal-1", "RBP2a Hap1", 
                  "DBPII Hap2", "PTEX150 Hap1", "RBP2a Sal-1", 
                  "RBP2a Hap9","MSP5 Hap2", "MSP5 Hap1", "RIPR Hap1",
                  "RIPR Hap2", "RBP2b Hap1", "RBP2b Hap2", "RBP2b Hap3",
                  "RBP2b Sal-1", "RIPR Sal-1", "MSP5 Sal-1")

inf_cat <- c( Thai_cat, brazil_cat, control_cat )

N_part <- length(inf_cat)



View(ant_names_short)
matrix(ant_names_short)


# Prepare the data (Keep it in wide format)
selected_cols <- c(1:23)  # Antibodies to plot
antibody_data <- AB[, selected_cols]  # Extract antibody data
antibody_data$Category <- inf_cat  # Add infection category as a column

# Set the order of the x-axis categories
antibody_data$Category <- factor(antibody_data$Category, levels = c(
  "Thai_current", "Thai_recent", "Thai_old", "Thai_never",
  "brazil_current", "brazil_recent", "brazil_old", "brazil_never",
  "TRC", "BRC", "ARC",  "VBDR"
))


# Define color mapping for categories
color_map <- c(
  "brazil_current" = "#CC79A7",     # Reddish Purple
  "brazil_recent" = "#D55E00",      # Vermillion (red-orange)
  "brazil_old" = "#E69F00",         # Orange
  "brazil_never" = "#F0E442",       # Yellow
  
  "Thai_current" = "#CC79A7",       # Reddish Purple
  "Thai_recent" = "#D55E00",        # Vermillion
  "Thai_old" = "#E69F00",           # Orange
  "Thai_never" = "#F0E442",         # Yellow
  
  # Similar blue shades for control categories
  "TRC" = "#0072B2",                # Blue
  "ARC" = "#56B4E9",                # Sky Blue
  "BRC" = "#4682B4",                # Steel Blue
  "VBDR" = "#5F9EA0"                # Cadet Blue
)


# Define explicit y-axis breaks and labels
y_breaks <- log(c(1e-5, 1e-4, 1e-3, 1e-2))
y_labels <- c("1.E-5", "1.E-4", "1.E-3", "1.E-2")

for (i in selected_cols ) {
  
  p <- ggplot(antibody_data, aes(x = Category, y = AB[, i], fill = Category)) +
    geom_boxplot(na.rm = TRUE) +
    
  
    annotate("text", x = 2.5, y = log(5e-6), label = "Thailand", size = 12, vjust = 1) +
    annotate("text", x = 6.5, y = log(5e-6), label = "Brazil", size = 12, vjust = 1) +
    annotate("text", x = 10.5, y = log(5e-6), label = "Controls", size = 12, vjust = 1) +
    

    geom_segment(x = 0.38, xend = 12.6, y = log(7e-6), yend = log(7e-6), 
                 color = "black", linewidth = 0.8) +
    geom_segment(x = 2.5, xend = 2.5, y = log(7e-6), yend = log(6e-6), 
                 color = "black", linewidth = 0.8) +
    geom_segment(x = 6.5, xend = 6.5, y = log(7e-6), yend = log(6e-6), 
                 color = "black", linewidth = 0.8) +
    geom_segment(x = 10.5, xend = 10.5, y = log(7e-6), yend = log(6e-6), 
                 color = "black", linewidth = 0.8) +
    geom_segment(x = 0.40, xend = 0.40, y = log(7e-6), yend = log(3e-2), 
                 color = "black", linewidth = 0.8) +
    geom_segment(x = 0.38, xend = 12.6, y = log(3e-2), yend = log(3e-2), 
                 color = "black", linewidth = 0.8) +
    geom_segment(x = 12.58, xend = 12.58, y = log(7e-6), yend = log(3e-2), 
                 color = "black", linewidth = 0.8) +
    
    scale_fill_manual(values = color_map) +
    scale_x_discrete(labels = rep("", 12)) +
    scale_y_continuous(
      breaks = y_breaks, labels = y_labels,
      limits = log(c(5e-6, 3e-2)), expand = c(0, 0)
    ) +
    geom_hline(yintercept = log(1.95e-5), linetype = "dashed") +
    geom_hline(yintercept = log(0.02), linetype = "dashed") +
    coord_cartesian(clip = "off") +  
    
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),  
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 32),
      legend.position = "none", 
      plot.title = element_text(face = "bold", size = 42, hjust = 0.5),
      plot.margin = margin(t = 10, r = 10, b = 70, l = 10),
      panel.grid = element_blank(),
      axis.ticks.y = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.25, "cm"),
      axis.text.y = element_text(size = 26)
    ) +
    ylab("RAU") +
    ggtitle(paste(ant_names_hap[i]))
   
  ggsave(
    filename = paste0("output/", ant_names_short[i], ".png"),
    plot = p,
    width = 8, height = 6, units = "in",
    dpi = 300
  )
  
  print(p)  # Print each plot separately
}



# Create dummy data for the legend
legend_data <- data.frame(
  Category = factor(c(
    "infected", 
    "infected 1-9 months", 
    "infected 9-12 months",
    "no detected infection",
    "negative controls: Thai RC",
    "negative controls: Brazilian RC",
    "negative controls: Aus RC",
    "negative controls: VBDR"
  ), levels = c(
    "infected", 
    "infected 1-9 months", 
    "infected 9-12 months",
    "no detected infection",
    "negative controls: Thai RC",
    "negative controls: Brazilian RC",
    "negative controls: Aus RC",
    "negative controls: VBDR"
  )),
  x = 1:8,  # Dummy values for x-axis
  y = 1:8   # Dummy values for y-axis
)

#  Define color mapping as per specification
legend_colors <- c(
  "infected" = "#CC79A7",                     # Reddish Purple
  "infected 1-9 months" = "#D55E00",           # Vermillion
  "infected 9-12 months" = "#E69F00",          # Orange
  "no detected infection" = "#F0E442",         # Yellow
  "negative controls: Thai RC" = "#0072B2",    # Blue
  "negative controls: Brazilian RC" = "#4682B4", # Steel Blue
  "negative controls: Aus RC" = "#56B4E9",     # Sky Blue
  "negative controls: VBDR" = "#5F9EA0"                          # Cadet Blue
)

# Create dummy plot to extract legend with two-column layout
legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = Category)) +
  geom_point(size = 5) +
  scale_color_manual(values = legend_colors) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 16),
    legend.key.size = unit(1, "cm"),
    legend.spacing.x = unit(1, "cm")
  ) +
  guides(
    color = guide_legend(ncol = 2) 
  )

# Save the legend plot as a separate high-resolution PNG
ggsave(
  filename = "output/legend_two_column_layout.png",
  plot = legend_plot,
  width = 10, height = 3, units = "in",
  dpi = 300
)

# 🎨 View the plot in RStudio
print(legend_plot)



## =========================
## LMM: Sal-1 vs Haplotypes
## =========================

# 1) Build subject, cohort, and simplified category vectors --------------------

# Make subject IDs unique across datasets
subject_ids <- c(paste0("T_", Thai_data$sampleid),
                 paste0("B_", brazil_data$sampleid),
                 paste0("C_", control_data$sampleid))

stopifnot(length(subject_ids) == nrow(AB))
stopifnot(length(inf_cat)    == nrow(AB))

# Derive cohort from your inf_cat codes
cohort_vec <- dplyr::case_when(
  grepl("^Thai_",   inf_cat) ~ "Thailand",
  grepl("^brazil_", inf_cat) ~ "Brazil",
  TRUE                       ~ "Control"
)

# Derive category (remove cohort prefix for Thai/Brazil; keep control codes)
category_vec <- inf_cat
category_vec <- gsub("^Thai_",   "", category_vec)
category_vec <- gsub("^brazil_", "", category_vec)

# 2) Map each selected antibody column to antigen & variant --------------------

# Use the names you already provided in ant_names_hap (one per selected column)
stopifnot(length(selected_cols) == length(ant_names_hap))

hap_map <- tibble::tibble(
  col_index = selected_cols,
  hap_label = ant_names_hap
) |>
  dplyr::mutate(
    antigen = sub(" [^ ]+$", "", hap_label),   # text before the last space
    variant = sub("^.* ",     "", hap_label)   # last token: e.g., Sal-1, Hap1, AH
  )

# 3) Build a long data frame ---------------------------------------------------

# Keep only the selected columns from AB and add id/cohort/category
wide_df <- as.data.frame(AB[, selected_cols])
colnames(wide_df) <- paste0("ab_", seq_along(selected_cols))  # temp names

wide_df$subject  <- subject_ids
wide_df$cohort   <- factor(cohort_vec, levels = c("Thailand", "Brazil", "Control"))
wide_df$category <- factor(category_vec, levels = c(
  "current","recent","old","never",  # from Thai_/brazil_
  "TRC","BRC","ARC","VBDR"           # controls
))

# Pivot to long and attach antigen/variant labels from hap_map
long_df <- wide_df |>
  tidyr::pivot_longer(
    cols = dplyr::starts_with("ab_"),
    names_to = "ab_col",
    values_to = "log_value"  # already natural log scale per your AB <- log(AB)
  ) |>
  dplyr::mutate(col_index = as.integer(sub("ab_", "", ab_col))) |>
  dplyr::left_join(hap_map |>
                     dplyr::mutate(ab_col = paste0("ab_", dplyr::row_number())),
                   by = "ab_col") |>
  dplyr::select(subject, cohort, category, antigen, variant, log_value)

# Ensure factors with desired references
long_df$variant <- factor(long_df$variant)
if (!("Sal-1" %in% levels(long_df$variant))) {
  stop("No 'Sal-1' level found in variant factor. Check ant_names_hap.")
}
long_df$variant <- relevel(long_df$variant, ref = "Sal-1")



# ----------------------------------------------------------
# 4) Simple per-antigen LMM + Sal-1 vs others via emmeans
# ----------------------------------------------------------
# make sure reference is Sal-1
# Ensure Sal-1 is the reference
# make sure reference is Sal-1

# ----------------------------------------------------------
# 4) Unpaired linear models + emmeans (no helpers)
# ----------------------------------------------------------

# Ensure reference levels
# ----------------------------------------------------------
# 4) Unpaired linear models + emmeans (t-based CIs; no CI normalizer)
# ----------------------------------------------------------

# Ensure reference levels
long_df$variant  <- stats::relevel(factor(long_df$variant), ref = "Sal-1")
long_df$category <- droplevels(factor(long_df$category))
long_df$cohort   <- droplevels(factor(long_df$cohort))

# minimal helper to coerce numeric
numify <- function(x) suppressWarnings(as.numeric(as.character(x)))

run_lm_unpaired <- function(dat, antigen_name, overall_weights = "proportional") {
  df1 <- dplyr::filter(dat, antigen == antigen_name, !is.na(log_value))
  if (nrow(df1) == 0L) return(NULL)
  
  # ---------- A) Sal-1 vs haplotypes within category × cohort ----------
  mA <- lm(log_value ~ variant * category * cohort, data = df1)
  
  emm_s <- emmeans::emmeans(mA, ~ variant | category * cohort)
  ct_s  <- emmeans::contrast(emm_s, "trt.vs.ctrl", ref = "Sal-1")
  out_s <- as.data.frame(summary(ct_s, infer = c(TRUE, TRUE), adjust = "BH"),
                         stringsAsFactors = FALSE)
  
  # coerce numeric and keep only rows with CIs
  for (cn in c("estimate","SE","df","t.ratio","p.value","lower.CL","upper.CL")) {
    if (cn %in% names(out_s)) out_s[[cn]] <- numify(out_s[[cn]])
  }
  out_s <- dplyr::filter(out_s, !is.na(estimate), !is.na(lower.CL), !is.na(upper.CL))
  
  # back-transform (natural log → ratio)
  out_s$GMR     <- exp(out_s$estimate)
  out_s$GMR_lo  <- exp(out_s$lower.CL)
  out_s$GMR_hi  <- exp(out_s$upper.CL)
  out_s$antigen <- antigen_name
  out_s$haplotype_vs <- sub(" - Sal-1$", "", out_s$contrast)
  
  out_s <- dplyr::relocate(out_s, antigen, category, cohort, contrast, haplotype_vs,
                           estimate, SE, df, lower.CL, upper.CL, t.ratio, p.value,
                           GMR, GMR_lo, GMR_hi)
  
  # ---------- B) Overall by cohort (averaging over categories) ----------
  mB <- lm(log_value ~ variant * cohort, data = df1)
  
  emm_o <- emmeans::emmeans(mB, ~ variant | cohort, weights = overall_weights)
  ct_o  <- emmeans::contrast(emm_o, "trt.vs.ctrl", ref = "Sal-1")
  out_o <- as.data.frame(summary(ct_o, infer = c(TRUE, TRUE), adjust = "BH"),
                         stringsAsFactors = FALSE)
  
  for (cn in c("estimate","SE","df","t.ratio","p.value","lower.CL","upper.CL")) {
    if (cn %in% names(out_o)) out_o[[cn]] <- numify(out_o[[cn]])
  }
  out_o <- dplyr::filter(out_o, !is.na(estimate), !is.na(lower.CL), !is.na(upper.CL))
  
  out_o$GMR     <- exp(out_o$estimate)
  out_o$GMR_lo  <- exp(out_o$lower.CL)
  out_o$GMR_hi  <- exp(out_o$upper.CL)
  out_o$antigen <- antigen_name
  out_o$haplotype_vs <- sub(" - Sal-1$", "", out_o$contrast)
  
  out_o <- dplyr::relocate(out_o, antigen, cohort, contrast, haplotype_vs,
                           estimate, SE, df, lower.CL, upper.CL, t.ratio, p.value,
                           GMR, GMR_lo, GMR_hi)
  
  # ---------- C) Country pairwise (Thailand vs Brazil) within variant × category ----------
  # (kept simple; you can ignore this CSV if not needed)
  emm_c <- emmeans::emmeans(mA, ~ cohort | variant * category)
  ct_c  <- emmeans::contrast(emm_c, "pairwise")
  out_c <- as.data.frame(summary(ct_c, infer = c(TRUE, TRUE), adjust = "BH"),
                         stringsAsFactors = FALSE)
  
  for (cn in c("estimate","SE","df","t.ratio","p.value","lower.CL","upper.CL")) {
    if (cn %in% names(out_c)) out_c[[cn]] <- numify(out_c[[cn]])
  }
  out_c <- dplyr::filter(out_c, !is.na(estimate), !is.na(lower.CL), !is.na(upper.CL))
  
  out_c$GM_ratio <- exp(out_c$estimate)  # TH/BR ratio if that’s the contrast order
  out_c$GM_lo    <- exp(out_c$lower.CL)
  out_c$GM_hi    <- exp(out_c$upper.CL)
  out_c$antigen  <- antigen_name
  
  out_c <- dplyr::relocate(out_c, antigen, variant, category, contrast,
                           estimate, SE, df, lower.CL, upper.CL, t.ratio, p.value,
                           GM_ratio, GM_lo, GM_hi)
  
  list(detail = out_s, overall_by_cohort = out_o, country_comp = out_c)
}

# Run across antigens & bind
all_antigens <- sort(unique(hap_map$antigen))
res <- purrr::map(all_antigens, ~ run_lm_unpaired(long_df, .x))

detail_tbl   <- dplyr::bind_rows(lapply(res, `[[`, "detail"))
overall_tbl  <- dplyr::bind_rows(lapply(res, `[[`, "overall_by_cohort"))
country_tbl  <- dplyr::bind_rows(lapply(res, `[[`, "country_comp"))

# BH-FDR per antigen
if (nrow(detail_tbl))  detail_tbl  <- detail_tbl  %>% dplyr::group_by(antigen) %>% dplyr::mutate(q_BH = p.adjust(p.value, "BH")) %>% dplyr::ungroup()
if (nrow(overall_tbl)) overall_tbl <- overall_tbl %>% dplyr::group_by(antigen) %>% dplyr::mutate(q_BH = p.adjust(p.value, "BH")) %>% dplyr::ungroup()
if (nrow(country_tbl)) country_tbl <- country_tbl %>% dplyr::group_by(antigen, variant, category) %>% dplyr::mutate(q_BH = p.adjust(p.value, "BH")) %>% dplyr::ungroup()

# Save
dir.create("output/lmm", showWarnings = FALSE, recursive = TRUE)
if (nrow(detail_tbl))   write.csv(detail_tbl,  "output/lmm/UNPAIRED_LM_BY_STRATUM.csv", row.names = FALSE)
if (nrow(overall_tbl))  write.csv(overall_tbl, "output/lmm/UNPAIRED_LM_OVERALL_BY_COHORT.csv", row.names = FALSE)
if (nrow(country_tbl))  write.csv(country_tbl, "output/lmm/UNPAIRED_LM_COUNTRY_PAIRWISE.csv", row.names = FALSE)

cat("Done. Files written (if non-empty):\n",
    "- output/lmm/UNPAIRED_LM_BY_STRATUM.csv\n",
    "- output/lmm/UNPAIRED_LM_OVERALL_BY_COHORT.csv\n",
    "- output/lmm/UNPAIRED_LM_COUNTRY_PAIRWISE.csv\n")
