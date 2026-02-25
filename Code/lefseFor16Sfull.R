=# ==============================================================================
# 16S rRNA Microbiome Analysis: LEfSe (Biomarker Discovery)
# ==============================================================================

# Clear workspace
rm(list = ls())

# Load required libraries
library(openxlsx)
library(dplyr)
library(ggplot2)
library(SummarizedExperiment)
library(lefser)

# ------------------------------------------------------------------------------
# 1. Load Data
# ------------------------------------------------------------------------------
# Load relative abundance data
loadtmp <- readRDS("../Data/002.16Sfull/001.species.dat.Rds")
dat_rela <- loadtmp$dat_rela

# [FIX]: Define dat_rela_select before filtering it
dat_rela_select <- dat_rela

# Filter out features with zero abundance across all samples
dat_rela_select <- dat_rela_select[rowSums(dat_rela_select) > 0, ]

# Load and align clinical metadata
metainfo <- openxlsx::read.xlsx("../Data/000.clinical.process.xlsx", sheet = "S16")
metainfo_select <- metainfo
rownames(metainfo_select) <- metainfo_select$uID

# Subset metadata to match the columns of the abundance matrix
metainfo_select <- metainfo_select[colnames(dat_rela_select), ]

# Note: all_patient is defined but currently unused. Kept for potential downstream use.
all_patient <- unique(metainfo_select$patID)

# ------------------------------------------------------------------------------
# 2. Process Taxonomic Lineage
# ------------------------------------------------------------------------------
dat_name <- read.csv("../Data/002.16Sfull/ASV_Taxon_Origin_asv.full.xls", sep = "\t")

# Keep relevant taxonomy columns (1-8) and drop the redundant 2nd column
dat_name <- dat_name[, 1:8]
dat_name <- dat_name[, -2]

# Create a MetaPhlAn-style lineage string (separated by '|')
dat_name$merger <- apply(dat_name, 1, paste, sep = "", collapse = "|")
dat_name <- unique(dat_name)

# Align taxonomy rows with abundance matrix rows
dat_name <- dat_name[match(rownames(dat_rela_select), dat_name$species), ]
rownames(dat_name) <- dat_name$species

# [FIX]: Robust sanity check to ensure perfect row alignment. 
# Script will halt here with an error if rows do not match perfectly.
stopifnot(all(rownames(dat_name) == rownames(dat_rela_select)))

# ------------------------------------------------------------------------------
# 3. Construct SummarizedExperiment Object
# ------------------------------------------------------------------------------
# lefser requires a SummarizedExperiment object as input
counts_df <- as.data.frame(dat_rela_select)
coldata_df <- as(metainfo_select, "data.frame")

objectD <- SummarizedExperiment(
  assays = list(counts = counts_df),
  colData = coldata_df,
  rowData = dat_name[, -8] # Exclude the 'merger' column from rowData
)

# ------------------------------------------------------------------------------
# 4. Iterative LEfSe Analysis (All Timepoints vs Baseline)
# ------------------------------------------------------------------------------
all_time <- c(
  "002_D1", "003_D14", "004_M1", "005_M2", "006_M3", 
  "007_M4", "008_M5", "009_M6", "010_M8", "011_M9", 
  "012_M10", "013_M12"
)

result <- list()

# Loop through each timepoint to compare against Baseline (001_BL)
for (each in all_time) {
  
  # Subset SummarizedExperiment object for Baseline and current timepoint
  objectD_select <- objectD[, objectD$timeIDcontiue %in% c("001_BL", each)]
  
  # Run LEfSe (set seed for reproducibility of the internal Kruskal-Wallis/Wilcoxon tests)
  set.seed(1234)
  res <- lefser(objectD_select, groupCol = "timeIDcontiue")
  
  # Skip to the next iteration if no significant biomarkers are found
  if (nrow(res) == 0) { next }
  
  # Assign group directionality based on LDA score
  # Positive score = enriched in the later timepoint; Negative = enriched in Baseline
  res$group <- ifelse(res$scores > 0, each, "001_BL")
  
  # Lock factor levels to maintain sorting in the plot
  res$Names <- factor(res$Names, levels = as.character(res$Names))
  
  # Generate custom LDA bar plot
  p <- ggplot(as.data.frame(res), aes(x = Names, y = scores)) + 
    geom_bar(stat = "identity", aes(fill = group), width = 0.618) + 
    scale_fill_manual(values = c("#7da1cc", "#e57a77")) + 
    coord_flip() +
    ylab("LDA SCORE (log 10)") + 
    theme(
      axis.title.y = element_blank(), 
      axis.title.x = element_text(size = 11, face = "bold"), 
      axis.text.y  = element_text(vjust = 0.7, size = 9, face = "bold"), 
      axis.text.x  = element_text(vjust = 0.7, size = 9, face = "bold"), 
      plot.title   = element_text(hjust = 0.5, size = 13, face = "bold")
    )
  
  # Save plot (dynamically named by timepoint)
  ggsave(
    sprintf("./Code/001.16S/012.lefse/%s.pdf", each), 
    plot = p, 
    width = 5.82, 
    height = 6.32
  )                 
  
  # Store results
  result[[each]] <- res
}

# ------------------------------------------------------------------------------
# 5. Export Final Results
# ------------------------------------------------------------------------------
# Append supporting data matrices to the output list
result[["dat_name"]] <- dat_name
result[["counts"]]   <- data.frame(row = rownames(counts_df), counts_df, check.names = FALSE)
result[["coldata"]]  <- coldata_df

# Save compiled results to Excel
openxlsx::write.xlsx(result, "./001.16S/012.lefse/total.xlsx")