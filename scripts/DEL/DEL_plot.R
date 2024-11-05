#!/usr/bin/env Rscript

# DEL Methylation Analysis and Visualization Script
# This script performs DEL methylation level analysis and visualization for a given dataset.
# It includes calculating standardized methylation differences, categorizing methylation levels,
# and generating scatter and density plots based on methylation differences.

# Required Libraries
library(dplyr)
library(ggplot2)
library(ggforce)

# Function: plot_DEL_scatter
# Arguments:
# - file_path: Path to the DEL methylation data file.
# - cutoff: Threshold for filtering DELs based on methylation differences (default is 0.5).
# Description:
# - This function reads DEL methylation data, calculates methylation differences, categorizes them, and generates a scatter plot.

plot_DEL_scatter <- function(file_path, cutoff, circle_radius = 0.5) {
    # Load data
    all_data <- read.table(file_path, header = FALSE)
    colnames(all_data) <- c("Chr", "Start", "End", "Body", "Up", "Down", "SBody", "SUp", "SDown", "Number", "Len")
    
    # Calculate non-standardized and standardized methylation differences
    all_data <- all_data %>% 
        mutate(mp = (Body - Up) / 100, 
               md = (Body - Down) / 100) %>%
        mutate(category = case_when(
                mp > 0 & md > 0 ~ "High",
                mp < 0 & md < 0 ~ "Low",
                mp >= 0 & md <= 0 ~ "Other1",
                mp <= 0 & md >= 0 ~ "Other2",
                TRUE ~ "NA")) %>%
        mutate(Smp = (SBody - SUp), 
               Smd = (SBody - SDown)) %>%
        mutate(Scategory = case_when(
                Smp > 0 & Smd > 0 ~ "High",
                Smp < 0 & Smd < 0 ~ "Low",
                Smp >= 0 & Smd <= 0 ~ "Other1",
                Smp <= 0 & Smd >= 0 ~ "Other2",
                TRUE ~ "NA")) %>%
        mutate(proportion = Number / 106)  # Normalized by the total sample count (106)
                      
        # Scatter Plot
        ggplot(all_data, aes(x = Smp, y = Smd, color = Scategory, size = proportion)) +
            geom_point(alpha = 0.4) +
            scale_color_manual(values = c("High" = "#BF1D2D", "Low" = "#293890", "Other1" = "#6f6f6f", "Other2" = "#6f6f6f")) +
            scale_size_continuous(breaks = seq(0, 1, by = 0.2), range = c(1, 5)) +
            theme_bw() +
            geom_hline(yintercept = 0, linetype = "solid") +
            geom_vline(xintercept = 0, linetype = "solid") +
            scale_x_continuous(labels = scales::percent, limits = c(-1.2, 1.2)) +
            scale_y_continuous(labels = scales::percent, limits = c(-1.2, 1.2)) +
            theme(
                  axis.text = element_text(color = 'black', size = 14),
                  axis.title = element_text(color = 'black', size = 16),
                  plot.title = element_text(size = 16, color = "black"),
                  panel.grid = element_blank(),
                  legend.position = "none"
                  ) +
        labs(x = "Methylation Difference", y = "Methylation Difference") +
        geom_circle(aes(x0 = 0, y0 = 0, r = cutoff), color = "black", linewidth = 0.6, linetype = "dashed")
}

# Function: plot_DEL_density
# Arguments:
# - file_path: Path to the DEL methylation data file.
# - cutoff: Threshold for filtering DELs based on methylation differences (default is 0.5).
# Description:
# - This function reads DEL methylation data, calculates standardized methylation differences,
#   and generates a density plot with KDE-based contours.

plot_DEL_density <- function(file_path, cutoff = 0.5) {
    # Load data
    all_data <- read.table(file_path, header = FALSE)
    colnames(all_data) <- c("Chr", "Start", "End", "Body", "Up", "Down", "SBody", "SUp", "SDown", "Number", "Len")
    
    # Calculate standardized methylation differences and proportion
    all_data <- all_data %>% 
        mutate(Smp = (SBody - SUp), 
               Smd = (SBody - SDown)) %>%
        mutate(proportion = Number / 106)  # Normalized by the total sample count (106)
        
        # Density Plot
        grey_col <- colorRampPalette(c("white", "#808080", "#606060", "#404040", "#303030", "#202020", "black"))(50)
        ggplot(all_data, aes(x = Smd, y = Smp)) +
            geom_density_2d_filled(show.legend = FALSE, bins = 50) +
            scale_fill_manual(values = grey_col) +
            geom_hline(yintercept = 0, linetype = "solid") +
            geom_vline(xintercept = 0, linetype = "solid") +
            theme_bw() +
            scale_x_continuous(labels = scales::percent) +
            scale_y_continuous(labels = scales::percent) +
            theme(
                  axis.title = element_text(color = 'black', size = 16),
                  axis.text = element_text(size = 14, color = "black"),
                  panel.grid = element_blank(),
                  legend.position = "none"
                  ) +
        geom_circle(aes(x0 = 0, y0 = 0, r = cutoff), color = "black", linewidth = 0.3, linetype = "dashed") +
        labs(x = "Methylation Difference", y = "Methylation Difference")
}

# Example usage
# plot_DEL_scatter("path/to/datafile.cpg", cutoff = 0.5)
# plot_DEL_density("path/to/datafile.cpg", cutoff = 0.5)
