library(SingleCellExperiment)
library(scran)
library(scater)
library(dplyr)
library(readr)

setwd("C:/Users/49152/Downloads/Multi-omics/")

# -------------------------------
# Step 1: Read counts
# -------------------------------
raw_counts <- read_tsv("Github_data/Transcriptomics/C10SVEC_singlecells_Counts.txt", show_col_types = FALSE)
colnames(raw_counts)[1] <- "Gene"

count_matrix <- as.matrix(raw_counts[, -1])
rownames(count_matrix) <- raw_counts$Gene

sce <- SingleCellExperiment(assays = list(counts = count_matrix))

# -------------------------------
# Step 2: Metadata parsing
# -------------------------------
cell_names <- colnames(sce)
cell_info <- strsplit(cell_names, "_")

chip <- sapply(cell_info, function(x) x[1])       # first part = chip
cell_line <- sapply(cell_info, function(x) x[2])  # second part = cell line

# Special case: "07J_SVEC_C10"
cell_line[cell_names == "07J_SVEC_C10"] <- "SVEC"

colData(sce)$Chip <- chip
colData(sce)$CellLine <- cell_line

# -------------------------------
# Step 3: Filtering
# -------------------------------
keep_genes <- rowSums(counts(sce) > 0) > 3  # expressed in at least 4 cells
sce <- sce[keep_genes, ]

# Count number of detected genes (≥ 5 counts) per cell
detected_genes <- apply(counts(sce), 2, function(x) sum(x >= 5))

# Add to colData
colData(sce)$DetectedGenes5 <- detected_genes

# Apply thresholds by CellLine
keep_cells <- (colData(sce)$CellLine == "C10"  & colData(sce)$DetectedGenes5 >= 3000) |
  (colData(sce)$CellLine == "SVEC" & colData(sce)$DetectedGenes5 >= 2000)

sce <- sce[, keep_cells]

# -------------------------------
# Step 4: Library size normalization + log-transform
# -------------------------------
sce <- computeSumFactors(sce)
sce <- logNormCounts(sce)   # stores normalized log counts in logcounts(sce)

# -------------------------------
# Step 5: HVG selection
# -------------------------------
dec <- modelGeneVar(sce)
top_hvgs <- getTopHVGs(dec, n = 2000)
sce <- sce[top_hvgs, ]

# -------------------------------
# Step 6: Export
# -------------------------------
rna_matrix <- as.matrix(logcounts(sce))

# Convert to data frame with Gene column
rna_df <- data.frame(Gene = rownames(rna_matrix), rna_matrix, check.names = FALSE)

write.csv(rna_df,
          file = "MOFA/input/1_Single_cell/transcriptomics_for_MOFA.csv",
          row.names = FALSE, quote = FALSE)

metadata <- as.data.frame(colData(sce))

# add rownames (sample IDs) as a proper column
metadata$SampleID <- rownames(metadata)

# move SampleID to the first column
metadata <- metadata[, c("SampleID", setdiff(names(metadata), "SampleID"))]

# write out with SampleID as column name
write.csv(metadata,
          file = "MOFA/input/1_Single_cell/transcriptomics_metadata.csv",
          row.names = FALSE, quote = FALSE)

# -------------------------------
# Step 7: PCA diagnostic
# -------------------------------
pca <- prcomp(t(rna_matrix), scale. = FALSE)  # PCA on cells
plot(pca$x[,1], pca$x[,2],
     col = ifelse(metadata$CellLine == "C10", "red", "blue"),
     pch = 16,
     xlab = "PC1", ylab = "PC2",
     main = "RNA PCA: C10 (red) vs SVEC (blue)")