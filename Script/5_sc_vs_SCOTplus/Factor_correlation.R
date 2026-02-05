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
  main = "MOFA+ factor correlation: Single vs SCOT+",
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
