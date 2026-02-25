# ==============================================================================
# Differential Abundance and GSEA Analysis of Metagenomic KO Data
# ==============================================================================

# Clear workspace (Note: consider removing this line if the script is meant 
# to be sourced by other pipeline scripts, but it is fine for standalone runs)
rm(list = ls())

# Load required libraries
library(openxlsx)
library(limma)
library(clusterProfiler)
library(enrichplot)
library(colorspace)
library(dplyr) # Required for the pipe (%>%) operator

# ------------------------------------------------------------------------------
# 1. Load and Preprocess Clinical Metadata
# ------------------------------------------------------------------------------
metainfo <- openxlsx::read.xlsx("./Data/000.clinical.process.addCell.xlsx")

# Define valid timepoints and subset metadata
valid_timepoints <- c(
  "001_BL", "002_D1", "003_D14", "004_M1", 
  "005_M2", "006_M3", "009_M6", "0013_M12"
)
metainfo <- metainfo[metainfo$timeIDcontiue %in% valid_timepoints, ]
rownames(metainfo) <- metainfo$metaID

# ------------------------------------------------------------------------------
# 2. Load and Preprocess Metagenomic KO Data
# ------------------------------------------------------------------------------
dat <- read.delim(
  gzfile("./Data/003.metagenome/001.KO_unstratified.tsv.gz"), 
  row.names = 1, 
  check.names = FALSE
)

# Retain non-functional rows (e.g., UNMAPPED, UNINTEGRATED) but clean column names
colnames(dat) <- sub("_genefamilies", "", colnames(dat))

# Remove taxa-level stratifications (rows containing "g__")
dat <- dat[!grepl("g__", rownames(dat)), ]

# Align abundance data columns with metadata rows
dat <- dat[, metainfo$metaID]

# Convert counts to Counts Per Million (CPM)
dat <- apply(dat, 2, function(x) { 10^6 * x / sum(x) })

# Filter low-abundance features: keep features present (value > 0) in > 20% of samples
dat <- dat[apply(dat, 1, function(x) { sum(x > 0) > length(x) * 0.2 }), ]

# Apply Log2 transformation with pseudo-count to handle zeros
pseudo_count <- 1e-6
y <- log2(dat + pseudo_count)

# Verify sample alignment (will throw an error if mismatched)
stopifnot(all(colnames(y) == rownames(metainfo)))

# ------------------------------------------------------------------------------
# 3. Linear Modeling (limma)
# ------------------------------------------------------------------------------
# Define grouping factor based on timepoints
group <- factor(
  metainfo$timeIDcontiue,
  levels = valid_timepoints 
)
# Ensure R-compatible names for columns (e.g., adds 'X' before numbers)
levels(group) <- make.names(levels(group))

# Create design matrix (without intercept)
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# Fit linear model and apply empirical Bayes smoothing with trend 
# (Trend = TRUE is recommended for log-intensity data)
fit <- lmFit(y, design)                 
fit <- eBayes(fit, trend = TRUE)

# ------------------------------------------------------------------------------
# 4. Define and Compute Contrasts
# ------------------------------------------------------------------------------
cont <- makeContrasts(
  D1_vs_BL   = X002_D1 - X001_BL,
  D14_vs_BL  = X003_D14 - X001_BL,
  M1_vs_BL   = X004_M1 - X001_BL,
  M2_vs_BL   = X005_M2 - X001_BL,
  M3_vs_BL   = X006_M3 - X001_BL,
  M6_vs_BL   = X009_M6 - X001_BL,
  M12_vs_BL  = X0013_M12 - X001_BL, 
  M1_vs_D14  = X004_M1 - X003_D14,
  levels = design
)

fit2 <- contrasts.fit(fit, cont)
fit2 <- eBayes(fit2, trend = TRUE)

# ------------------------------------------------------------------------------
# 5. Extract Differential Analysis Results
# ------------------------------------------------------------------------------
result <- list()
contrast_names <- colnames(cont)

# Loop through all contrasts to extract topTable results dynamically
for (cmp in contrast_names) {
  tt <- topTable(fit2, coef = cmp, n = Inf, adjust.method = "BH")
  tt$pathway <- rownames(tt)   # Append feature names
  tt$type <- cmp               # Record contrast type
  result[[cmp]] <- tt
}

# Append processed CPM data and metadata to the output list
result[["cpm_dat"]] <- data.frame(pathway = rownames(dat), dat, check.names = FALSE)
result[["metainfo"]] <- metainfo

# Save complete results
saveRDS(result, "./002.metagenome/001.KO.pathway.RDS")
openxlsx::write.xlsx(result, "./002.metagenome/001.DF.KO.xlsx")

# ------------------------------------------------------------------------------
# 6. Gene Set Enrichment Analysis (GSEA)
# ------------------------------------------------------------------------------
# Load KEGG database mappings
kegg_dat <- readRDS("./Data/KEGG/ko_long.RDS")

# Filter for relevant functional categories (Metabolism only)
kegg_dat <- kegg_dat[!kegg_dat$A %in% c("A09160 Human Diseases"), ]
kegg_dat <- kegg_dat[kegg_dat$A %in% c("A09100 Metabolism"), ]

# Extract unique pathway terms and names
kegg_dat_name <- kegg_dat[, c("C_code", "C_name")] %>% unique()
colnames(kegg_dat_name) <- c("term", "name")

# Extract gene-to-pathway mapping
kegg_dat_select <- kegg_dat[, c("C_code", "K")] 
colnames(kegg_dat_select) <- c("term", "gene")

# Example: Run GSEA for the D14_vs_BL contrast
select_dat <- result[["D14_vs_BL"]]

# Create a ranked named vector of logFC values
select_dat_gene <- select_dat$logFC
names(select_dat_gene) <- select_dat$pathway 
select_dat_gene <- sort(select_dat_gene, decreasing = TRUE)

# Filter KEGG annotations to include only genes present in the dataset
kegg_dat_select <- kegg_dat_select[kegg_dat_select$gene %in% rownames(select_dat), ]

# Execute GSEA
result_GSEA <- GSEA(
  geneList = select_dat_gene,
  TERM2GENE = kegg_dat_select,
  TERM2NAME = kegg_dat_name
)

# Inspect top GSEA results
head(result_GSEA@result)

# ------------------------------------------------------------------------------
# 7. Visualization
# ------------------------------------------------------------------------------
# Specify pathway IDs of interest for plotting
paths <- c("ko00360", "ko00121", "ko00640") 

# Generate standard GSEA plot
enrichplot::gseaplot2(result_GSEA, pathwayId = paths, pvalue_table = TRUE)

# Generate customized multi-pathway GSEA plot with distinct colors
enrichplot::gseaplot2(
  result_GSEA, 
  pathwayId = paths, 
  color = colorspace::rainbow_hcl(length(paths)), 
  subplots = c(1, 2), 
  pvalue_table = TRUE
)