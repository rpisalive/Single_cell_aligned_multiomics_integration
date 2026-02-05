library(dplyr)
library(tibble)
library(pheatmap)
library(MOFA2)

# Load models
m_single <- load_model("C:/Users/49152/Downloads/Multi-omics/MOFA/output/single_cell_trained_model.hdf5")
m_pbulk <- load_model("C:/Users/49152/Downloads/Multi-omics/MOFA/output/pdbk_different_trained_model.hdf5")

Z_single <- get_factors(m_single, factors = "all", as.data.frame = FALSE)[[1]]  

meta_pb <- read.csv(
  "C:/Users/49152/Downloads/Multi-omics/MOFA/input/2_Pseudobulk_different_composition/transcriptomics_pdbk_different_metadata.csv",
  stringsAsFactors = FALSE
)

# Check
head(meta_pb)
# Columns: SampleID, PseudoBulkID, CellLine, Chip, nCells, SingleCellIDs

# Split SingleCellIDs into vectors
meta_pb <- meta_pb %>%
  rowwise() %>%
  mutate(SingleCellIDs_list = strsplit(SingleCellIDs, ";")) %>%
  ungroup()

# Initialize an empty list for aggregated factors
Z_single_agg <- lapply(meta_pb$PseudoBulkID, function(pb_id) {
  cells <- meta_pb$SingleCellIDs_list[[which(meta_pb$PseudoBulkID == pb_id)]]
  # Average the factors across the cells in this pseudobulk
  colMeans(Z_single[cells, , drop = FALSE], na.rm = TRUE)
})

# Convert to a data frame
Z_single_agg <- do.call(rbind, Z_single_agg)
rownames(Z_single_agg) <- meta_pb$PseudoBulkID

Z_pbulk <- get_factors(m_pbulk, factors = "all", as.data.frame = FALSE)[[1]]  

# Align pseudobulk IDs
common_pseudobulks <- intersect(rownames(Z_single_agg), rownames(Z_pbulk))
Z_single_agg <- Z_single_agg[common_pseudobulks, , drop = FALSE]
Z_pbulk <- Z_pbulk[common_pseudobulks, , drop = FALSE]

cor_mat <- cor(Z_single_agg, Z_pbulk, use = "pairwise.complete.obs")

# Heatmap
pheatmap(cor_mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue","white","red"))(50),
         main = "Factor correlation: pseudobulk vs aggregated single-cell")

# Overall similarity
mean_abs_corr <- mean(abs(cor_mat))
cat("Mean absolute correlation:", round(mean_abs_corr, 3), "\n")
