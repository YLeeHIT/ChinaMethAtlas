#!/usr/bin/env Rscript

# PCA and Anosim Analysis Script for Population Data
# This script performs PCA analysis and Anosim analysis on population data.
# - Step 1: Load data and metadata (classification labels).
# - Step 2: Perform PCA and visualize the results with custom biplot.
# - Step 3: Conduct Anosim analysis to assess the similarity between groups.

# Load required libraries
library(RColorBrewer)  # For color palettes
library(PCAtools)      # For PCA analysis and plotting
library(dplyr)         # For data manipulation
library(vegan)         # For Anosim analysis

# Set file paths
infile <- "rawdata.txt"        # Input raw data file
inlab <- "pca_lab.txt"    # Input classification labels file

# Step 1: Load Data and Metadata
# Read classification labels and convert to a data frame
labdata <- read.table(inlab, sep = "\t", header = TRUE) %>% as.data.frame()

# Read raw data and convert to data frame
rawdata <- read.table(infile, sep = "\t", header = FALSE) %>% as.data.frame()
newdata <- t(rawdata)  # Transpose data for PCA and Anosim analysis

# Step 2: PCA Analysis and Custom Plotting
# Perform PCA with PCAtools
p_div <- pca(newdata, metadata = labdata, removeVar = 0.1, scale = TRUE)

# Customize and plot PCA biplot
pca_plot <- biplot(
    p_div, pointSize = 2, lab = NULL, x = "PC1", y = "PC2",
    colLegendTitle = "Region", colby = 'Division',
    colkey = c(North = '#0072B5FF', South = '#BC3C29FF', Xizang = '#E18727FF'),
    xlab = "PC1", ylab = "PC2",
    gridlines.major = FALSE, gridlines.minor = FALSE, borderWidth = 0.4,
    axisLabSize = 8, ellipse = TRUE, ellipseType = "t", ellipseLevel = 0.9,
    ellipseFill = FALSE, ellipseLineSize = 0.2, legendPosition = 'none',
    legendLabSize = 7, legendIconSize = 2, legendTitleSize = 8) +
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Save PCA plot to a file
ggsave("PCA_plot.png", plot = pca_plot, width = 8, height = 6)

# Step 3: Anosim Analysis
# Calculate similarity matrix using Euclidean distance
similarity_matrix <- vegdist(newdata)

# Perform Anosim analysis based on group labels
group <- labdata$Division
anosim_result <- anosim(similarity_matrix, group)

# Output Anosim results
cat("ANOSIM Analysis Results:\n")
print(summary(anosim_result))
