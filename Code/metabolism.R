# ==============================================================================
# Preprocessing, Scaling, and PCA of Metabolomics Data
# ==============================================================================

# Clear workspace
rm(list = ls())

# Load required libraries
library(openxlsx)
library(factoextra)
library(ggplot2)
library(dplyr)

# ------------------------------------------------------------------------------
# 1. Load Data and Align Metadata
# ------------------------------------------------------------------------------
metainfo <- openxlsx::read.xlsx("../Data/000.clinical.process.addCell.xlsx")
dat <- read.csv("./org_mix/all_metab_data.txt", sep = "\t")

# Extract expression matrix (assuming columns 23 to end contain sample data)
dat_meta <- dat[, 23:ncol(dat)]
rownames(dat_meta) <- dat$metab_id

# Align metabolomics data columns with metadata sample IDs (uID)
dat_meta <- dat_meta[, as.character(metainfo$uID)]

# ------------------------------------------------------------------------------
# 2. Normalization and Filtering
# ------------------------------------------------------------------------------
# Convert intensities to a relative abundance scale (CPM equivalents)
dat_meta <- apply(dat_meta, 2, function(x) { 10^6 * x / sum(x) })

# Filter 1: Retain metabolites present (value > 0) in more than 5 samples
dat_meta <- dat_meta[apply(dat_meta, 1, function(x) { sum(x > 0) > 5 }), ]

# Filter 2: Retain the top 75% most variable metabolites
# Calculate Coefficient of Variation (CV) = Standard Deviation / Mean
dat_cv <- apply(dat_meta, 1, function(x) { sd(x) / mean(x) })
dat_cv <- sort(dat_cv, decreasing = TRUE)

# Subset data based on the highest CV scores
dat_meta <- dat_meta[names(dat_cv) %>% head(length(dat_cv) * 0.75), ]

# ------------------------------------------------------------------------------
# 3. Timepoint Selection and Re-filtering
# ------------------------------------------------------------------------------
# Define timepoints of interest (Note: Verify "013_M12" vs "0013_M12")
select_time <- c(
  "001_BL", "002_D1", "003_D14", "004_M1", 
  "005_M2", "006_M3", "009_M6", "013_M12"
)

# Subset metadata and metabolomics matrix for selected timepoints
metainfo <- metainfo[metainfo$timeIDcontiue %in% select_time, ]
dat_meta <- dat_meta[, metainfo$uID]

# Re-apply Filter 1 after dropping samples to ensure robustness
dat_meta <- dat_meta[apply(dat_meta, 1, function(x) { sum(x > 0) > 5 }), ]

# ------------------------------------------------------------------------------
# 4. Data Transformation and Scaling Options
# ------------------------------------------------------------------------------
# Log1p transformation (log(x + 1)) to stabilize variance for right-skewed data
dat_meta_log <- log1p(dat_meta)

# Transpose matrix for PCA (prcomp expects rows = samples, columns = features)
X <- t(dat_meta_log) 

# Define custom scaling functions
scale_ctr <- function(mat) {
  # Mean-centering only
  scale(mat, center = TRUE, scale = FALSE)
}

scale_uv <- function(mat) {
  # Unit Variance / Autoscaling (Divide by Standard Deviation)
  scale(mat, center = TRUE, scale = TRUE)
}

scale_par <- function(mat) {
  # Pareto scaling (Divide by square root of Standard Deviation)
  sds <- apply(mat, 2, sd, na.rm = TRUE)
  scale(mat, center = TRUE, scale = sqrt(sds))
}

# ------------------------------------------------------------------------------
# 5. Principal Component Analysis (PCA) and Visualization
# ------------------------------------------------------------------------------
# Ensure the grouping factor is explicitly defined and ordered for plotting
group_factor <- factor(metainfo$timeIDcontiue, levels = select_time)

# Define a function to run PCA and generate a standardized plot
plot_pca <- function(Xscaled, title_txt) {
  # Run PCA (scaling is set to FALSE because inputs are pre-scaled)
  pcs <- prcomp(Xscaled, center = FALSE, scale. = FALSE) 
  
  p <- fviz_pca_ind(
    pcs,
    geom = "point", 
    col.ind = group_factor, # Color by clinical group
    palette = c("grey20", "#feb24c", "red", "#7FC97F", "#BEAED4", "#FDC086", "#c51b8a", "#386CB0"),
    addEllipses = TRUE,
    ellipse.level = 0.95,
    ellipse.type = "confidence", 
    legend.title = "Groups", 
    repel = TRUE
  ) +
    theme_bw(base_size = 12) +
    theme(
      aspect.ratio = 1,
      legend.title = element_blank(),
      plot.title = element_text(hjust = 0.5),
      axis.text = element_text(family = "Arial"), 
      axis.title = element_text(family = "Arial"),
      panel.grid = element_blank(),
      panel.border = element_rect(fill = NA, colour = "grey87")
    ) + 
    # Standardize shapes across the 8 groups if shape mapping is added later
    scale_shape_manual(values = 1:8) + 
    ggtitle(title_txt)
  
  return(p)
}

# Generate PCA plots for unscaled (log only) and the three scaling methods
p_log <- plot_pca(t(dat_meta_log), "PCA-log")
p_ctr <- plot_pca(scale_ctr(X), "PCA-CTR")
p_uv  <- plot_pca(scale_uv(X), "PCA-UV")
p_par <- plot_pca(scale_par(X), "PCA-PAR")

# Save the Unit Variance (UV) scaled PCA plot
ggsave("./metabolism/001.uv_pva.pdf", plot = p_uv, width = 5.79, height = 4.21)