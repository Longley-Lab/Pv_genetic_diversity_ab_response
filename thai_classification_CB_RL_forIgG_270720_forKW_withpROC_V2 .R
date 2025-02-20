library(MASS)
library(ROCR)
library(pROC)
library(ggplot2)
library(svglite)
library(writexl)



#Function

calculate_AUC <- function(...) {
  
  AUC <- data.frame(antigen_name = character(), AUC = numeric(), p_value = numeric(), stringsAsFactors = FALSE)
  
  indices <- c(...)
  
  roc_result_first <- roc_results[[paste("ROC_Column", indices[1], sep = "_")]]
  first_antigen_name <- ant_names_short[indices[1]]
  
  AUC <- rbind(AUC, data.frame(antigen_name = first_antigen_name, AUC = auc(roc_result_first), p_value = NA))  # No comparison for the first
  
  for (i in indices[-1]) {
    
    roc_result <- roc_results[[paste("ROC_Column", i, sep = "_")]]
    
    current_auc <- auc(roc_result)
    
    # Perform the roc.test comparing current ROC result to the first ROC result
    test_result <- roc.test(roc_result_first, roc_result, method = "bootstrap")
    
    # Extract the p-value from the test
    p_value <- test_result$p.value
    
    # Add the AUC and p-value to the data frame
    AUC <- rbind(AUC, data.frame(antigen_name = ant_names_short[i], AUC = current_auc, p_value = p_value))
  }
  return(AUC)
}


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

setwd()

brazil_data = read.csv("20240629_Thai_IgG_geneticdiversity.csv")

control_data = read.csv("20240628_NC_IgG_geneticdiversity.csv")



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

# so for Brazil IgM going to try and running without dropping any - normally it drops 5 ehime proteins from IgG but this doesn't matter for this analyis
# as not using IgG only for RBP2b IgG. Then it should keep my IgM data for L19, L41 and L34 where I had to exclude data from 5 plates (so more than 150 missing)

AB <- rbind( brazil_data[,17:39], control_data[,17:39] )

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
#I just added the names of the ones I was using because there were only 9
ant_names_short = c("MSP8 L34 Sal1","Pv-fam-a SEM-8 Hap1", "MSP8 SEM-10 Hap1", 
                    "DBPII SEM-35 Hap1", "DBPII AH", "PTEX150 L18 Sal1", 
                    "Pv-fam-a L02 Sal1", "DBPII Sal1", "RBP2a SEM-28 Hap1", 
                    "DBPII SEM-37 Hap2", "PTEX150 SEM-39 Hap1", "RBP2a Sal1", 
                    "RBP2a SEM-29 Hap9","MSP5 SEM-33", "MSP5 SEM-31", "RIPR SEM-41",
                    "RIPR SEM-42", "RBP2b SEM-2", "RBP2b SEM-4", "RBP2b SEM-6",
                    "RBP2b P25", "RIPR", "MSP5 L19 Sal1")
                    
inf_cat <- c( brazil_cat, control_cat )

N_part <- length(inf_cat)



View(ant_names_short)
matrix(ant_names_short)

###################################
## Binary category

bin_cat = rep("old", N_part)
bin_cat[which(inf_cat=="brazil_current")] = "new"
bin_cat[which(inf_cat=="brazil_recent")]  = "new"



###################################
## Trim down to top 8 antigens - Note Rhea 18Jul2022 dont think we need this section
#the numbers correspond to which columns the antibodies are in in the spreadsheet
#N_ant is the total number you are looking at so if there are eg 25, replace the 8 with 25
# first just going to look at my 20 IgM responses - still calling in "top_8" just so i don't mess with code too much 

#top_8 <- c(L01_IgM, L02_IgM, L12_IgM, L19_IgM, L28_IgM, L32_IgM, L34_IgM, L36_IgM, L41_IgM, L54_IgM, L55_IgM, RBP2b_IgM, DBPII_AH_IgM, DBP_Sal1_IgM, L23_IgM, L31_IgM, PVX_092995_IgM, PVX_087885_IgM, PvEBP_IgM, PvRIPR_IgM)
#top_8 <- c(1,2,8,12,16,20,27,28,30,31,34,39,40,50,51,55,60,62,63,65)

#N_ant <- 20
#AB <- AB[,top_8]
#ant_names_short <- ant_names_short[top_8]

#View(top_8)
#View(ant_names_short)

#This part is to create ROC using pROC package

# Convert old or new into 0 or 1
bin_cat_num <- ifelse(bin_cat == "old", 0, 1)
print(bin_cat_num)  

# Check if lengths match
if (length(bin_cat_num) == nrow(AB)) {
  # Add bin_cat_num as a new column
  AB$bin_cat_num <- bin_cat_num
  
  # Rename the dataframe
  AB_pcr <- AB
  
  # Print the updated dataframe
  print(AB_pcr)
} else {
  stop("Length of bin_cat_num does not match the number of rows in AB")
}





#########################################################


# Initialize a list to store ROC results
roc_results <- list()

# For loop to iterate through columns 1 to 10
for (i in 1:N_ant) {
  # Extract the sero column (predictor)
  sero_column <- AB_pcr[, i]
  
  # Extract the pcr column (response)
  pcr_column <- AB_pcr[, ncol(AB_pcr)]
  
  # Merge sero and pcr results together
  roc_df <- data.frame(sero_column, pcr_column)
  
  # Clean the NA rows
  roc_df_cleaned <- na.omit(roc_df)
  
  # Extract cleaned predictor and response vectors
  sero_predictor <- roc_df_cleaned$sero_column
  pcr_response <- roc_df_cleaned$pcr_column
  
  # Compute the ROC curve
  roc_result <- roc(pcr_response, sero_predictor)
  
  # Store the ROC result in the list
  roc_results[[paste("ROC_Column", i, sep = "_")]] <- roc_result
  
}

######DBP#############################
# Store the colors you want to use
my_colors_DBP <- c("#1f77b4", "#2a9df4", "#6baed6", "#a6c8e0")

# Initialize an empty data frame to store all ROC data
all_roc_data_list <- list()

counter <- 1

for (i in c(8, 5, 4, 10)) {
  
  
  roc_result <- roc_results[[paste("ROC_Column", i, sep = "_")]]
  
  roc_plot_df <- data.frame(
    Specificity = 1 - roc_result$specificities,
    Sensitivity = roc_result$sensitivities,
    legend = ant_names_short[i]  # Label for the ROC curve
  )
  
  
  all_roc_data_list[[counter]] <- roc_plot_df
  

  counter <- counter + 1
}


all_roc_data <- do.call(rbind, all_roc_data_list)




all_roc_data$legend <- factor(all_roc_data$legend, levels = unique(all_roc_data$legend))

# Plot all ROC curves on the same axes and ensure the legend appears
plot <- ggplot(all_roc_data, aes(x = Specificity, y = Sensitivity, color = legend)) +
  geom_line(size = 1.2) +  # Line thickness
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "ROC Curves for DBP",
       x = "1 - Specificity", y = "Sensitivity", color = "Legend") +  # Add legend title
  scale_color_manual(values = my_colors_DBP) +  # Assign custom colors
  theme_minimal() +
  coord_cartesian(clip = "off") +
  theme(legend.position = c(1.0,0.12),
        legend.justification = "left",
        legend.title = element_text(size = 10),  # Reduce legend title font size
        legend.text = element_text(size = 8),   # Reduce legend text font size
        legend.key.size = unit(0.5, "cm"),     # Reduce size of the legend keys (color boxes)
        plot.title = element_text(hjust = 0.5), # Center the title
        panel.grid.major = element_line(color = "gray90"),
        plot.margin = margin(1, 4, 1, 1, "cm")) # Lighten grid lines
        


# Print the plot
print(plot)


#export the image
ggsave("ROC_Curves_DBP_Thai.svg", plot = plot, width = 8, height = 6, dpi = 600)

#export auc

AUC_DBP <- calculate_AUC(8, 5, 4, 10)

write_xlsx(AUC_DBP, "AUC_DBP.xlsx")



###########MSP5###############################
my_colors_MSP5 <- c("#ff7f0e", "#ffa64d", "#cc6600")

# Initialize an empty data frame to store all ROC data
all_roc_data_list <- list()

counter <- 1

for (i in c(23, 15, 14)) {
  
  
  roc_result <- roc_results[[paste("ROC_Column", i, sep = "_")]]
  
  roc_plot_df <- data.frame(
    Specificity = 1 - roc_result$specificities,
    Sensitivity = roc_result$sensitivities,
    legend = ant_names_short[i]  # Label for the ROC curve
  )
  
  
  all_roc_data_list[[counter]] <- roc_plot_df
  
  
  counter <- counter + 1
}


all_roc_data <- do.call(rbind, all_roc_data_list)




all_roc_data$legend <- factor(all_roc_data$legend, levels = unique(all_roc_data$legend))

# Plot all ROC curves on the same axes and ensure the legend appears
plot <- ggplot(all_roc_data, aes(x = Specificity, y = Sensitivity, color = legend)) +
  geom_line(size = 1.2) +  # Line thickness
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "ROC Curves for MSP5",
       x = "1 - Specificity", y = "Sensitivity", color = "Legend") +  # Add legend title
  scale_color_manual(values = my_colors_MSP5) +  # Assign custom colors
  theme_minimal() +
  coord_cartesian(clip = "off") +
  theme(legend.position = c(1.0,0.12),
        legend.justification = "left",
        legend.title = element_text(size = 10),  # Reduce legend title font size
        legend.text = element_text(size = 8),   # Reduce legend text font size
        legend.key.size = unit(0.5, "cm"),     # Reduce size of the legend keys (color boxes)
        plot.title = element_text(hjust = 0.5), # Center the title
        panel.grid.major = element_line(color = "gray90"),
        plot.margin = margin(1, 4, 1, 1, "cm")) # Lighten grid lines



# Print the plot
print(plot)

#export the image
ggsave("ROC_Curves_MSP5_Thai.svg", plot = plot, width = 8, height = 6, dpi = 600)

#export auc

AUC_MSP5 <- calculate_AUC(23, 15, 14)

write_xlsx(AUC_MSP5, "AUC_MSP5.xlsx")



############MSP8#################################
my_colors_MSP8 <- c("#2ca02c", "#66c266")

# Initialize an empty data frame to store all ROC data
all_roc_data_list <- list()

counter <- 1

for (i in c(1, 3)) {
  
  
  roc_result <- roc_results[[paste("ROC_Column", i, sep = "_")]]
  
  roc_plot_df <- data.frame(
    Specificity = 1 - roc_result$specificities,
    Sensitivity = roc_result$sensitivities,
    legend = ant_names_short[i]  # Label for the ROC curve
  )
  
  
  all_roc_data_list[[counter]] <- roc_plot_df
  
  
  counter <- counter + 1
}


all_roc_data <- do.call(rbind, all_roc_data_list)




all_roc_data$legend <- factor(all_roc_data$legend, levels = unique(all_roc_data$legend))

# Plot all ROC curves on the same axes and ensure the legend appears
plot <- ggplot(all_roc_data, aes(x = Specificity, y = Sensitivity, color = legend)) +
  geom_line(size = 1.2) +  # Line thickness
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "ROC Curves for MSP8",
       x = "1 - Specificity", y = "Sensitivity", color = "Legend") +  # Add legend title
  scale_color_manual(values = my_colors_MSP8) +  # Assign custom colors
  theme_minimal() +
  coord_cartesian(clip = "off") +
  theme(legend.position = c(1.0,0.12),
        legend.justification = "left",
        legend.title = element_text(size = 10),  # Reduce legend title font size
        legend.text = element_text(size = 8),   # Reduce legend text font size
        legend.key.size = unit(0.5, "cm"),     # Reduce size of the legend keys (color boxes)
        plot.title = element_text(hjust = 0.5), # Center the title
        panel.grid.major = element_line(color = "gray90"),
        plot.margin = margin(1, 4, 1, 1, "cm")) # Lighten grid lines



# Print the plot
print(plot)

#export the image
ggsave("ROC_Curves_MSP8_Thai.svg", plot = plot, width = 8, height = 6, dpi = 600)

#Export auc

AUC_MSP8 <- calculate_AUC(1, 3)

write_xlsx(AUC_MSP8, "AUC_MSP8.xlsx")


######################RBP2A########################
my_colors_RBP2A <- c("#d62728", "#e65555", "#a31d1d")

# Initialize an empty data frame to store all ROC data
all_roc_data_list <- list()

counter <- 1

for (i in c(12, 9, 13)) {
  
  
  roc_result <- roc_results[[paste("ROC_Column", i, sep = "_")]]
  
  roc_plot_df <- data.frame(
    Specificity = 1 - roc_result$specificities,
    Sensitivity = roc_result$sensitivities,
    legend = ant_names_short[i]  # Label for the ROC curve
  )
  
  
  all_roc_data_list[[counter]] <- roc_plot_df
  
  
  counter <- counter + 1
}


all_roc_data <- do.call(rbind, all_roc_data_list)




all_roc_data$legend <- factor(all_roc_data$legend, levels = unique(all_roc_data$legend))

# Plot all ROC curves on the same axes and ensure the legend appears
plot <- ggplot(all_roc_data, aes(x = Specificity, y = Sensitivity, color = legend)) +
  geom_line(size = 1.2) +  # Line thickness
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "ROC Curves for RBP2A",
       x = "1 - Specificity", y = "Sensitivity", color = "Legend") +  # Add legend title
  scale_color_manual(values = my_colors_RBP2A) +  # Assign custom colors
  theme_minimal() +
  coord_cartesian(clip = "off") +
  theme(legend.position = c(1.0,0.12),
        legend.justification = "left",
        legend.title = element_text(size = 10),  # Reduce legend title font size
        legend.text = element_text(size = 8),   # Reduce legend text font size
        legend.key.size = unit(0.5, "cm"),     # Reduce size of the legend keys (color boxes)
        plot.title = element_text(hjust = 0.5), # Center the title
        panel.grid.major = element_line(color = "gray90"),
        plot.margin = margin(1, 4, 1, 1, "cm")) # Lighten grid lines



# Print the plot
print(plot)

#export the image
ggsave("ROC_Curves_RBP2A_Thai.svg", plot = plot, width = 8, height = 6, dpi = 600)


#Export auc

AUC_RBP2A <- calculate_AUC(12, 9, 13)

write_xlsx(AUC_RBP2A, "AUC_RBP2A.xlsx")


###############RBP2B#########################################
my_colors_RBP2B <- c("#9467bd", "#af5fbf", "#684b96", "#7e57c2")

# Initialize an empty data frame to store all ROC data
all_roc_data_list <- list()

counter <- 1

for (i in c(21, 18, 19, 20)) {
  
  
  roc_result <- roc_results[[paste("ROC_Column", i, sep = "_")]]
  
  roc_plot_df <- data.frame(
    Specificity = 1 - roc_result$specificities,
    Sensitivity = roc_result$sensitivities,
    legend = ant_names_short[i]  # Label for the ROC curve
  )
  
  
  all_roc_data_list[[counter]] <- roc_plot_df
  
  
  counter <- counter + 1
}


all_roc_data <- do.call(rbind, all_roc_data_list)




all_roc_data$legend <- factor(all_roc_data$legend, levels = unique(all_roc_data$legend))

# Plot all ROC curves on the same axes and ensure the legend appears
plot <- ggplot(all_roc_data, aes(x = Specificity, y = Sensitivity, color = legend)) +
  geom_line(size = 1.2) +  # Line thickness
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "ROC Curves for RBP2B",
       x = "1 - Specificity", y = "Sensitivity", color = "Legend") +  # Add legend title
  scale_color_manual(values = my_colors_RBP2B) +  # Assign custom colors
  theme_minimal() +
  coord_cartesian(clip = "off") +
  theme(legend.position = c(1.0,0.12),
        legend.justification = "left",
        legend.title = element_text(size = 10),  # Reduce legend title font size
        legend.text = element_text(size = 8),   # Reduce legend text font size
        legend.key.size = unit(0.5, "cm"),     # Reduce size of the legend keys (color boxes)
        plot.title = element_text(hjust = 0.5), # Center the title
        panel.grid.major = element_line(color = "gray90"),
        plot.margin = margin(1, 4, 1, 1, "cm")) # Lighten grid lines



# Print the plot
print(plot)

#export the image
ggsave("ROC_Curves_RBP2B_Thai.svg", plot = plot, width = 8, height = 6, dpi = 600)

#Export auc

AUC_RBP2B <- calculate_AUC(21, 18, 19, 20)

write_xlsx(AUC_RBP2B, "AUC_RBP2B.xlsx")


#########################RIPR#############################
my_colors_RIPR <- c("#8c564b", "#c97a6d", "#5f3a32")

# Initialize an empty data frame to store all ROC data
all_roc_data_list <- list()

counter <- 1

for (i in c(22, 16, 17)) {
  
  
  roc_result <- roc_results[[paste("ROC_Column", i, sep = "_")]]
  
  roc_plot_df <- data.frame(
    Specificity = 1 - roc_result$specificities,
    Sensitivity = roc_result$sensitivities,
    legend = ant_names_short[i]  # Label for the ROC curve
  )
  
  
  all_roc_data_list[[counter]] <- roc_plot_df
  
  
  counter <- counter + 1
}


all_roc_data <- do.call(rbind, all_roc_data_list)




all_roc_data$legend <- factor(all_roc_data$legend, levels = unique(all_roc_data$legend))

# Plot all ROC curves on the same axes and ensure the legend appears
plot <- ggplot(all_roc_data, aes(x = Specificity, y = Sensitivity, color = legend)) +
  geom_line(size = 1.2) +  # Line thickness
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "ROC Curves for RIPR",
       x = "1 - Specificity", y = "Sensitivity", color = "Legend") +  # Add legend title
  scale_color_manual(values = my_colors_RIPR) +  # Assign custom colors
  theme_minimal() +
  coord_cartesian(clip = "off") +
  theme(legend.position = c(1.0,0.12),
        legend.justification = "left",
        legend.title = element_text(size = 10),  # Reduce legend title font size
        legend.text = element_text(size = 8),   # Reduce legend text font size
        legend.key.size = unit(0.5, "cm"),     # Reduce size of the legend keys (color boxes)
        plot.title = element_text(hjust = 0.5), # Center the title
        panel.grid.major = element_line(color = "gray90"),
        plot.margin = margin(1, 4, 1, 1, "cm")) # Lighten grid lines



# Print the plot
print(plot)

#export the image
ggsave("ROC_Curves_RIPR_Thai.svg", plot = plot, width = 8, height = 6, dpi = 600)

#Export auc

AUC_RIPR <- calculate_AUC(22, 16, 17)

write_xlsx(AUC_RIPR, "AUC_RIPR.xlsx")


#####################PTEX150##################################
my_colors_PTEX <- c("#e377c2", "#b3509e")

# Initialize an empty data frame to store all ROC data
all_roc_data_list <- list()

counter <- 1

for (i in c(6, 11)) {
  
  
  roc_result <- roc_results[[paste("ROC_Column", i, sep = "_")]]
  
  roc_plot_df <- data.frame(
    Specificity = 1 - roc_result$specificities,
    Sensitivity = roc_result$sensitivities,
    legend = ant_names_short[i]  # Label for the ROC curve
  )
  
  
  all_roc_data_list[[counter]] <- roc_plot_df
  
  
  counter <- counter + 1
}


all_roc_data <- do.call(rbind, all_roc_data_list)




all_roc_data$legend <- factor(all_roc_data$legend, levels = unique(all_roc_data$legend))

# Plot all ROC curves on the same axes and ensure the legend appears
plot <- ggplot(all_roc_data, aes(x = Specificity, y = Sensitivity, color = legend)) +
  geom_line(size = 1.2) +  # Line thickness
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "ROC Curves for PTEX150",
       x = "1 - Specificity", y = "Sensitivity", color = "Legend") +  # Add legend title
  scale_color_manual(values = my_colors_PTEX) +  # Assign custom colors
  theme_minimal() +
  coord_cartesian(clip = "off") +
  theme(legend.position = c(1.0,0.12),
        legend.justification = "left",
        legend.title = element_text(size = 10),  # Reduce legend title font size
        legend.text = element_text(size = 8),   # Reduce legend text font size
        legend.key.size = unit(0.5, "cm"),     # Reduce size of the legend keys (color boxes)
        plot.title = element_text(hjust = 0.5), # Center the title
        panel.grid.major = element_line(color = "gray90"),
        plot.margin = margin(1, 4, 1, 1, "cm")) # Lighten grid lines



# Print the plot
print(plot)

#export the image
ggsave("ROC_Curves_PTEX_Thai.svg", plot = plot, width = 8, height = 6, dpi = 600)

#Export auc
AUC_PTEX <- calculate_AUC(6, 11)

write_xlsx(AUC_PTEX, "AUC_PTEX.xlsx")

#######################fam-a####################
my_colors_fam_a <- c("#7f7f7f", "#4d4d4d")

# Initialize an empty data frame to store all ROC data
all_roc_data_list <- list()

counter <- 1

for (i in c(7, 2)) {
  
  
  roc_result <- roc_results[[paste("ROC_Column", i, sep = "_")]]
  
  roc_plot_df <- data.frame(
    Specificity = 1 - roc_result$specificities,
    Sensitivity = roc_result$sensitivities,
    legend = ant_names_short[i]  # Label for the ROC curve
  )
  
  
  all_roc_data_list[[counter]] <- roc_plot_df
  
  
  counter <- counter + 1
}


all_roc_data <- do.call(rbind, all_roc_data_list)




all_roc_data$legend <- factor(all_roc_data$legend, levels = unique(all_roc_data$legend))

# Plot all ROC curves on the same axes and ensure the legend appears
plot <- ggplot(all_roc_data, aes(x = Specificity, y = Sensitivity, color = legend)) +
  geom_line(size = 1.2) +  # Line thickness
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "ROC Curves for fam-a",
       x = "1 - Specificity", y = "Sensitivity", color = "Legend") +  # Add legend title
  scale_color_manual(values = my_colors_fam_a) +  # Assign custom colors
  theme_minimal() +
  coord_cartesian(clip = "off") +
  theme(legend.position = c(1.0,0.12),
        legend.justification = "left",
        legend.title = element_text(size = 10),  # Reduce legend title font size
        legend.text = element_text(size = 8),   # Reduce legend text font size
        legend.key.size = unit(0.5, "cm"),     # Reduce size of the legend keys (color boxes)
        plot.title = element_text(hjust = 0.5), # Center the title
        panel.grid.major = element_line(color = "gray90"),
        plot.margin = margin(1, 4, 1, 1, "cm")) # Lighten grid lines



# Print the plot
print(plot)

#export the image
ggsave("ROC_Curves_fam_a_Thai.svg", plot = plot, width = 8, height = 6, dpi = 600)

#Export AUC

AUC_fam_a <- calculate_AUC(7, 2)
# Write the AUC and p-value data to an Excel file
write_xlsx(AUC_fam_a, "AUC_fam_a.xlsx")

############################################




