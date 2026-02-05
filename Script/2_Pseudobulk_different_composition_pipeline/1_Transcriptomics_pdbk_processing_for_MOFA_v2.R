library(SingleCellExperiment)
library(scran)
library(scater)
library(dplyr)
library(tibble)
library(readr)
library(scuttle)   # for aggregateAcrossCells
library(edgeR)
library(limma)
library(ggplot2)
library(factoextra)  # for explained variance visualization

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
# Step 3.5: Create pseudobulks with redistribution of small groups
# -------------------------------

set.seed(42)  # for reproducibility

cells_per_bulk <- 3  # desired number of cells per pseudobulk

meta <- as.data.frame(colData(sce))

meta <- meta %>%
  dplyr::group_by(Chip, CellLine) %>%
  dplyr::group_modify(~{
    df <- .x   # the cells in the group
    group_vals <- .y   # tibble of the current group's Chip and CellLine
    
    n_cells <- nrow(df)
    n_full <- n_cells %/% cells_per_bulk
    n_leftover <- n_cells %% cells_per_bulk
    
    # Assign cells to bulk indices (1:n_full)
    if (n_full > 0) {
      bulk_ids <- rep(1:n_full, each = cells_per_bulk)[1:n_cells]
    } else {
      bulk_ids <- rep(1, n_cells)
    }
    
    # Handle leftovers: redistribute evenly to existing bulks
    if (n_leftover > 0 && n_full > 0) {
      leftover_idx <- (n_full * cells_per_bulk + 1):n_cells
      redistribute <- rep(1:n_full, length.out = n_leftover)
      bulk_ids[leftover_idx] <- redistribute
    }
    
    df$BulkIndex <- bulk_ids
    df$PseudoBulkID <- paste(group_vals$CellLine, group_vals$Chip,
                             sprintf("bulk%02d", df$BulkIndex),
                             sep = "_")
    df
  }) %>%
  dplyr::ungroup()

# Attach pseudobulk IDs back to colData
colData(sce)$PseudoBulkID <- meta$PseudoBulkID

# 4. Aggregate counts within each pseudobulk group
# This sums raw counts per gene within each pseudobulk
# Aggregate counts per pseudobulk
sce_bulk <- aggregateAcrossCells(
  sce,
  ids = colData(sce)$PseudoBulkID,
  statistics = "sum"
)

# -------------------------------
# Step 5: Update metadata and store single-cell IDs per pseudobulk
# -------------------------------

# Move rownames to a column
coldata_df <- as.data.frame(colData(sce)) %>%
  rownames_to_column(var = "SingleCellID")

# Aggregate metadata per pseudobulk, keeping actual single-cell IDs
bulk_meta <- coldata_df %>%
  group_by(PseudoBulkID) %>%
  summarise(
    CellLine = dplyr::first(CellLine),
    Chip = dplyr::first(Chip),
    nCells = n(),
    SingleCellIDs = paste(SingleCellID, collapse = ";"),
    .groups = "drop"
  )

# rownames = PseudoBulkID
rownames(bulk_meta) <- bulk_meta$PseudoBulkID

# Replace sce with the pseudobulked object for downstream processing
sce <- sce_bulk

cat("Number of pseudobulks generated:", ncol(sce), "\n")
table(bulk_meta$CellLine)

# Sum of counts per pseudobulk (pseudobulked object)
bulk_totals <- colSums(counts(sce))

# Sum of counts per pseudobulk computed manually from original single-cell counts
manual_totals <- tapply(
  colSums(counts(sce_bulk)),  # you could also use original 'sce'
  colData(sce)$PseudoBulkID,
  sum
)

# Compare
all.equal(as.numeric(bulk_totals), as.numeric(manual_totals))

# 1) TMM normalization + voom (recommended for pseudobulk)
dge <- DGEList(counts = counts(sce))
dge <- calcNormFactors(dge, method = "TMM")

# design only intercept because we do NOT adjust for Chip here (we want factors that may capture cell-line differences)
design <- model.matrix(~1, data = as.data.frame(colData(sce)))

# voom: converts to log2-CPM with observation-level weights
v <- voom(dge, design = design, plot = FALSE)
logexpr <- v$E   # log2-CPM matrix (genes x samples)

# Put voom logcounts into sce for downstream steps
assay(sce, "logcounts_voom") <- logexpr

# 2) HVG selection: use variance modelling on the voom matrix
# scran has helpers — we compute per-gene variance statistics from the matrix
# If your scran version lacks modelGeneVarRow, you can wrap the matrix into SCE and call modelGeneVar.
tmp_sce <- SingleCellExperiment(assays = list(logcounts = assay(sce, "logcounts_voom")))
dec <- modelGeneVar(tmp_sce, assay.type = "logcounts")  # variance modeling on the voom logcounts
# Pick top HVGs (adjust n to your data; with few samples be conservative)
n_hvgs <- 2000
n_hvgs <- min(nrow(sce), n_hvgs)  # safety
top_hvgs <- getTopHVGs(dec, n = n_hvgs)

# Subset sce to HVGs and use the voom-corrected matrix as final expression
sce <- sce[top_hvgs, ]
final_rna <- assay(sce, "logcounts_voom")  # rows = genes (HVGs), cols = pseudobulks

# 3) Centering and scaling: GLOBAL (do NOT center per group)
# Center each gene (row) to zero mean across all pseudobulks, then optionally scale to unit variance.
row_means <- rowMeans(final_rna, na.rm = TRUE)
row_sds <- matrixStats::rowSds(final_rna, na.rm = TRUE)
# Avoid division by zero
row_sds[row_sds == 0] <- 1

final_rna_z <- sweep(final_rna, 1, row_means, FUN = "-")
final_rna_z <- sweep(final_rna_z, 1, row_sds, FUN = "/")

# Transpose: PCA expects samples as rows and genes as columns
pca_res <- prcomp(t(final_rna_z), center = FALSE, scale. = FALSE)

# Variance explained by each PC
var_explained <- pca_res$sdev^2 / sum(pca_res$sdev^2)
cat("Variance explained by first 5 PCs:\n")
print(round(100 * var_explained[1:5], 2))

# Create a dataframe for plotting
pca_df <- data.frame(
  Sample = colnames(final_rna_z),
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  CellLine = bulk_meta[colnames(final_rna_z), "CellLine"],
  Chip = bulk_meta[colnames(final_rna_z), "Chip"]
)

# PC1 vs PC2 colored by CellLine
ggplot(pca_df, aes(x = PC1, y = PC2, color = CellLine, shape = Chip)) +
  geom_point(size = 3) +
  theme_classic(base_size = 14) +
  labs(title = "PCA of pseudobulk transcriptomics (HVGs)",
       subtitle = "Colored by CellLine, shaped by Chip",
       x = paste0("PC1 (", round(100 * var_explained[1], 1), "%)"),
       y = paste0("PC2 (", round(100 * var_explained[2], 1), "%)")) +
  scale_color_brewer(palette = "Set1")

# Optionally: visualize explained variance curve
fviz_eig(pca_res, addlabels = TRUE, ylim = c(0, 50)) +
  ggtitle("PCA scree plot (variance explained per component)")

# Save final matrix for MOFA (genes x samples)
rna_matrix_for_mofa <- final_rna_z

# Export to CSV for MOFA input (as before)
rna_df <- data.frame(Gene = rownames(rna_matrix_for_mofa), rna_matrix_for_mofa, check.names = FALSE)
write.csv(rna_df,
          file = "MOFA/input/2_Pseudobulk_different_composition/transcriptomics_pdbk_different_for_MOFA.csv",
          row.names = FALSE, quote = FALSE)

# Prepare metadata (SampleID etc.)
# bulk_meta must exist from your earlier steps; ensure it matches sce column order
metadata <- bulk_meta[ colnames(sce), ]
metadata$SampleID <- rownames(metadata)
metadata <- metadata[, c("SampleID", setdiff(names(metadata), "SampleID"))]
write.csv(metadata,
          file = "MOFA/input/2_Pseudobulk_different_composition/transcriptomics_pdbk_different_metadata.csv",
          row.names = FALSE, quote = FALSE)
