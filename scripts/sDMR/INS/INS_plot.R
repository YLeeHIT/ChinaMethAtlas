# Load required libraries
library(ggplot2)   # For plotting
library(ggforce)   # For geometric functions like geom_circle
library(tidyr)     # For data manipulation
library(dplyr)     # For data manipulation
library(ggpubr)    # For publication-ready plots
library(ggtern)    # For ternary plots (if necessary)

# Function: plot_InsDensity
# Description: This function reads methylation difference data from a file, processes the data, 
# and generates a density plot to visualize the methylation difference between body and upstream/downstream regions.
# Args:
#   file_path: Path to the data file
#   cutoff: The threshold for filtering out small methylation differences

plot_InsDensity <- function(file_path, cutoff) {
  # Read in the data
  all_data <- read.table(file_path, header = FALSE)
  colnames(all_data) <- c("Chr", "Start", "End", "Body", "Up", "Down", 
                          "UpNum", "DownNum", "SBody", "SUp", "SDown", "Number")
  
  # Calculate methylation differences (non-standardized data)
  all_data <- all_data %>% 
    mutate(mp = (Body - Up), 
           md = (Body - Down)) %>%
    mutate(category = case_when(
      mp > 0 & md > 0 ~ "High",
      mp < 0 & md < 0 ~ "Low",
      mp >= 0 & md <= 0 ~ "Other1",
      mp <= 0 & md >= 0 ~ "Other2",
      TRUE ~ "NA"
    ))
  
  # Calculate standardized methylation differences
  all_data <- all_data %>% 
    mutate(Smp = (SBody - SUp), 
           Smd = (SBody - SDown)) %>%
    mutate(Scategory = case_when(
      Smp > 0 & Smd > 0 ~ "High",
      Smp < 0 & Smd < 0 ~ "Low",
      Smp >= 0 & Smd <= 0 ~ "Other1",
      Smp <= 0 & Smd >= 0 ~ "Other2",
      TRUE ~ "NA"
    ))
  
  # Count categories before filtering
  category_counts <- all_data %>%
    group_by(Scategory) %>%
    summarise(count = n())
  
  # Filter out data points with small methylation differences
  all_filtered <- all_data %>%
    filter(!(sqrt(Smp^2 + Smd^2) <= cutoff))
  
  # Count categories after filtering
  category_counts_filter <- all_filtered %>%
    group_by(Scategory) %>%
    summarise(count = n())
  
  # Calculate the proportion of Number column
  all_data <- all_data %>% mutate(proportion = Number / 106)
  
  # Generate the density plot
  ggplot(all_data, aes(x = Smd, y = Smp)) +
    # Add density plot filled with different levels
    geom_density_2d_filled(show.legend = FALSE, bins = 50) +
    scale_fill_manual(values = grey_col3) +  # Use a custom grey color palette
    geom_hline(yintercept = 0, linetype = "solid", color = "#303030") +
    geom_vline(xintercept = 0, linetype = "solid", color = "#303030") +
    theme_bw() +
    scale_x_continuous(labels = scales::percent) +
    scale_y_continuous(labels = scales::percent) +
    theme(axis.title = element_text(color = 'black', size = 16),
          axis.text = element_text(size = 14, color = "black"),
          panel.grid = element_blank(),
          legend.position = "none") +
    labs(x = "Methylation Difference", y = "Methylation Difference")
}

# Function: plot_InsDensity2
# Description: Similar to plot_InsDensity but uses a different approach for visualizing methylation differences.
# Args:
#   file_path: Path to the data file
#   cutoff: The threshold for filtering out small methylation differences

plot_InsDensity2 <- function(file_path, cutoff) {
  # Read in the data
  all_data <- read.table(file_path, header = FALSE)
  colnames(all_data) <- c("Chr", "Start", "End", "Body", "Up", "Down", 
                          "UpNum", "DownNum", "SBody", "SUp", "SDown", "Number")
  
  # Calculate standardized methylation differences
  all_data <- all_data %>% 
    mutate(Smp = (SBody - SUp), 
           Smd = (SBody - SDown)) %>%
    mutate(Scategory = case_when(
      Smp > 0 & Smd > 0 ~ "High",
      Smp < 0 & Smd < 0 ~ "Low",
      Smp >= 0 & Smd <= 0 ~ "Other1",
      Smp <= 0 & Smd >= 0 ~ "Other2",
      TRUE ~ "NA"
    ))
  
  # Calculate the proportion of Number column
  all_data <- all_data %>% mutate(proportion = Number / 106)
  
  # Generate the density plot
  grey_col3 <- colorRampPalette(c("white", "#808080", "#606060", "#404040", "#303030", 
                                  "#202020", "#101010", "black"))(50)
  
  ggplot(all_data, aes(x = Smd, y = Smp)) +
    geom_density_2d_filled(show.legend = FALSE, bins = 50) +  # Density filled plot
    scale_fill_manual(values = grey_col3) +  # Custom grey color scale
    geom_hline(yintercept = 0, linetype = "solid", color = "#303030") +
    geom_vline(xintercept = 0, linetype = "solid", color = "#303030") +
    theme_bw() +
    scale_x_continuous(labels = scales::percent) +
    scale_y_continuous(labels = scales::percent) +
    theme(axis.title = element_text(color = 'black', size = 16),
          axis.text = element_text(size = 14, color = "black"),
          panel.grid = element_blank(),
          legend.position = "none") +
    labs(x = "Methylation Difference", y = "Methylation Difference")
}

# Function: plot_InsAll
# Description: This function processes methylation data, filters and categorizes based on methylation difference,
# and generates a scatter plot visualizing methylation differences between "Body", "Up", and "Down" regions.
# It also writes the results to output files.
# Args:
#   file_path: Path to the data file.
#   cutoff: The cutoff for filtering methylation differences (below which points are excluded from the plot).

plot_InsAll <- function(file_path, cutoff) {
  all_data <- read.table(file_path, header = FALSE)
  colnames(all_data) <- c("Chr", "Start", "End", "Body", "Up", "Down", 
                          "UpNum", "DownNum", "SBody", "SUp", "SDown", "Number")
  
  # Process standardized data (SBody vs SUp, SBody vs SDown)
  all_data <- all_data %>%
    mutate(Smp = (SBody - SUp), 
           Smd = (SBody - SDown)) %>%
    mutate(Scategory = case_when(
      Smp > 0 & Smd > 0 ~ "High",
      Smp < 0 & Smd < 0 ~ "Low",
      Smp >= 0 & Smd <= 0 ~ "Other1",
      Smp <= 0 & Smd >= 0 ~ "Other2",
      TRUE ~ "NA"
    ))
  
  # Convert "Number" to proportion
  all_data <- all_data %>% mutate(proportion = Number)
  
  # Generate scatter plot
  ggplot(all_data, aes(x = Smp, y = Smd, color = Scategory, size = proportion)) +
    geom_point(alpha = 0.4) +  # Scatter plot with transparency
    scale_color_manual(values = c("High" = "#BF1D2D", "Low" = "#293890", 
                                 "Other1" = "#6f6f6f", "Other2" = "#6f6f6f")) +
    scale_size_continuous(breaks = seq(0, 1, by = 0.2), range = c(1, 5)) +  # Customize point size
    theme_bw() +  # Black-and-white theme for clarity
    geom_hline(yintercept = 0, linetype = "solid") +  # Add horizontal line at y = 0
    geom_vline(xintercept = 0, linetype = "solid") +  # Add vertical line at x = 0
    scale_x_continuous(labels = scales::percent, limits = c(-1.2, 1.2)) +  # Customize x-axis
    scale_y_continuous(labels = scales::percent, limits = c(-1.2, 1.2)) +  # Customize y-axis
    theme(axis.text = element_text(color = 'black', size = 14), 
          axis.title = element_text(color = 'black', size = 16), 
          plot.title = element_text(size = 16, color = "black"), 
          panel.grid = element_blank(), 
          legend.position = "none") +  # Customize plot appearance
    labs(x = "Methylation Difference", y = "Methylation Difference") +  # Labels
    geom_circle(aes(x0 = 0, y0 = 0, r = 0.5), color = "black", linewidth = 0.5, linetype = "dashed")  # Add dashed circle
}
