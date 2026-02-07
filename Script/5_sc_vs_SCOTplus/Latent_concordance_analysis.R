library(MOFA2)
library(pheatmap)
library(RColorBrewer)
library(clue) # for Hungarian matching

# ----------------------------
# Load models
# ----------------------------
m_single <- load_model("C:/Users/49152/Downloads/Multi-omics/MOFA/output/single_cell_trained_model.hdf5")
m_scot   <- load_model("C:/Users/49152/Downloads/Multi-omics/SCOT_plus/MOFA_output/SCOTplus_trained_model.hdf5")

# ----------------------------
# Extract factor matrices
# ----------------------------
Z_single <- get_factors(m_single, factors = "all", as.data.frame = FALSE)[[1]]
Z_scot   <- get_factors(m_scot, factors = "all", as.data.frame = FALSE)[[1]]

# ----------------------------
# Align shared cells
# ----------------------------
common_cells <- intersect(rownames(Z_single), rownames(Z_scot))
Z_single_aligned <- Z_single[common_cells, , drop = FALSE]
Z_scot_aligned   <- Z_scot[common_cells, , drop = FALSE]

# Optional: rename factors
colnames(Z_single_aligned) <- paste0("Single_F", seq_len(ncol(Z_single_aligned)))
colnames(Z_scot_aligned)   <- paste0("SCOT_F", seq_len(ncol(Z_scot_aligned)))

# ----------------------------
# Compute correlation matrix
# ----------------------------
cor_mat <- cor(Z_single_aligned, Z_scot_aligned, use = "pairwise.complete.obs", method = "pearson")

# ----------------------------
# Optional: reorder SCOT+ factors for best match
# ----------------------------
# Hungarian algorithm maximizes sum of correlations
cost_mat <- 1 - abs(cor_mat) # Hungarian minimizes cost
assignment <- solve_LSAP(cost_mat)
cor_mat_ordered <- cor_mat[, assignment]

# ----------------------------
# Heatmap colors
# ----------------------------
breaks <- seq(-1, 1, length.out = 101)
colors <- colorRampPalette(c("blue", "white", "red"))(100)

# ----------------------------
# Heatmap
# ----------------------------
pheatmap(
  cor_mat_ordered,
  cluster_rows = TRUE,       # cluster reference factors
  cluster_cols = FALSE,      # keep SCOT+ factors reordered by Hungarian matching
  color = colors,
  breaks = breaks,
  main = "Concordance of MOFA+ latent factors between reference and SCOT+ models",
  fontsize_row = 10,
  fontsize_col = 10,
  display_numbers = TRUE
)

# ----------------------------
# Summary statistics
# ----------------------------
mean_abs_corr <- mean(abs(cor_mat), na.rm = TRUE)
cat("Mean absolute factor correlation:", round(mean_abs_corr, 3), "\n")

# Max correlation per reference factor
max_corr <- apply(cor_mat, 1, function(x) max(abs(x), na.rm = TRUE))
cat("Mean of max correlations per reference factor:", round(mean(max_corr), 3), "\n")

# Optional: count how many factors are strongly matched (>0.5)
strong_matches <- sum(max_corr > 0.5)
cat("Number of reference factors with strong match (>0.5):", strong_matches, "\n")

# =========================================================
# FEATURE-LEVEL (LOADING) CORRELATION ANALYSIS
# =========================================================

# ----------------------------
# Extract loadings (weights)
# ----------------------------
W_single <- get_weights(m_single, views = "all", factors = "all")
W_scot   <- get_weights(m_scot, views = "all", factors = "all")

# Views assumed named "RNA" and "Protein"
# Check names(W_single) if needed

# ----------------------------
# Align shared features
# ----------------------------
common_genes <- intersect(rownames(W_single[["Transcriptomics"]]),
                          rownames(W_scot[["Transcriptomics"]]))

common_proteins <- intersect(rownames(W_single[["Proteomics"]]),
                             rownames(W_scot[["Proteomics"]]))

cat("Shared genes:", length(common_genes), "\n")
cat("Shared proteins:", length(common_proteins), "\n")

# ----------------------------
# Compute loading correlations
# ----------------------------
rna_corr  <- numeric(length(assignment))
prot_corr <- numeric(length(assignment))

for (i in seq_along(assignment)) {
  
  j <- assignment[i]  # matched SCOT+ factor
  
  # RNA
  w1_rna <- W_single[["Transcriptomics"]][common_genes, i]
  w2_rna <- W_scot[["Transcriptomics"]][common_genes, j]
  rna_corr[i] <- cor(abs(w1_rna), abs(w2_rna), use = "pairwise.complete.obs")
  
  # Protein
  w1_prot <- W_single[["Proteomics"]][common_proteins, i]
  w2_prot <- W_scot[["Proteomics"]][common_proteins, j]
  prot_corr[i] <- cor(abs(w1_prot), abs(w2_prot), use = "pairwise.complete.obs")
}

loading_df <- data.frame(
  Factor = paste0("F", seq_along(assignment)),
  RNA_loading_cor = rna_corr,
  Protein_loading_cor = prot_corr
)

print(loading_df)

# ----------------------------
# Heatmap of loading correlations
# ----------------------------
loading_mat <- cbind(Transcriptomics = rna_corr, Proteomics = prot_corr)
rownames(loading_mat) <- paste0("F", seq_along(assignment))

pheatmap(
  loading_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colors,
  breaks = breaks,
  main = "Concordance of MOFA+ factor loadings between reference and SCOT+ models)",
  display_numbers = TRUE,
  angle_col = 0
)

# ----------------------------
# Summary statistics
# ----------------------------
cat("Mean RNA loading correlation:", round(mean(rna_corr), 3), "\n")
cat("Mean protein loading correlation:", round(mean(prot_corr), 3), "\n")

cat("RNA factors > 0.5:", sum(rna_corr > 0.5), "\n")
cat("Protein factors > 0.5:", sum(prot_corr > 0.5), "\n")

