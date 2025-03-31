library(ggplot2)
library(dplyr)
library(MASS)
library(ROCR)
library(gridExtra)


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

Thai_data = read.csv("Thai_IgG_geneticdiversity.csv")

brazil_data = read.csv("Brazil_IgG_geneticdiversity.csv")

control_data = read.csv("NC_IgG_geneticdiversity.csv")



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
                  "fam-a Sal-1", "DBPII Sal-1", "RBP2A Hap1", 
                  "DBPII Hap2", "PTEX150 Hap1", "RBP2A Sal-1", 
                  "RBP2A Hap9","MSP5 Hap2", "MSP5 Hap1", "RIPR Hap1",
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
    
  
    annotate("text", x = 2.5, y = log(5e-6), label = "Thailand", size = 6, vjust = 1) +
    annotate("text", x = 6.5, y = log(5e-6), label = "Brazil", size = 6, vjust = 1) +
    annotate("text", x = 10.5, y = log(5e-6), label = "Controls", size = 6, vjust = 1) +
    

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
      axis.title.y = element_text(size = 14),
      legend.position = "none", 
      plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
      plot.margin = margin(t = 10, r = 10, b = 70, l = 10),
      panel.grid = element_blank(),
      axis.ticks.y = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.25, "cm")
    ) +
    ylab("RAU") +
    ggtitle(paste(ant_names_hap[i]))
   
  ggsave(
    filename = paste(ant_names_short[i], ".png"),
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
    "VBDR"
  ), levels = c(
    "infected", 
    "infected 1-9 months", 
    "infected 9-12 months",
    "no detected infection",
    "negative controls: Thai RC",
    "negative controls: Brazilian RC",
    "negative controls: Aus RC",
    "VBDR"
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
  "VBDR" = "#5F9EA0"                          # Cadet Blue
)

# Create dummy plot to extract legend with two-column layout
legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = Category)) +
  geom_point(size = 5) +
  scale_color_manual(values = legend_colors) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.key.size = unit(1, "cm"),
    legend.spacing.x = unit(1, "cm")
  ) +
  guides(
    color = guide_legend(ncol = 2) 
  )

# Save the legend plot as a separate high-resolution PNG
ggsave(
  filename = "legend_two_column_layout.png",
  plot = legend_plot,
  width = 10, height = 3, units = "in",
  dpi = 300
)

# 🎨 View the plot in RStudio
print(legend_plot)


