library(tidyverse)
library(readr)
library(org.Mm.eg.db)   # or org.Hs.eg.db if human
library(AnnotationDbi)
library(dplyr)
library(matrixStats)    # for rowSds

setwd("C:/Users/49152/Downloads/Multi-omics/")

# -------------------------------
# Step 1: Read raw MS intensities
# -------------------------------
raw <- read_tsv("Github_data/Proteomics/C10SVEC_singlecells_Protein_intensities.tsv",
                show_col_types = FALSE)

# Remove obvious contaminants
raw <- raw[!grepl("contam_sp", raw$PROTID), ]

# Extract sample columns
sample_cols <- setdiff(colnames(raw), "PROTID")

# Replace zeros with NA
raw[sample_cols] <- lapply(raw[sample_cols], function(x) ifelse(x == 0, NA, x))

# -------------------------------
# Step 2: Build initial expression matrix
# -------------------------------
expr_mat <- as.matrix(raw[, sample_cols])
rownames(expr_mat) <- raw$PROTID

# -------------------------------
# Step 3: Metadata
# -------------------------------
coldata <- data.frame(
  SampleID = colnames(expr_mat),
  Chip = substr(colnames(expr_mat), 1, 3),
  CellLine = ifelse(grepl("SVEC", colnames(expr_mat)), "SVEC",
                    ifelse(grepl("C10", colnames(expr_mat)), "C10", NA_character_)),
  stringsAsFactors = FALSE
)

# --------------------------------------------------------
# Step 4: Create pseudobulks (same logic as your pipeline)
# --------------------------------------------------------
set.seed(42)
cells_per_bulk <- 3  # adjust as needed

coldata <- coldata %>%
  group_by(CellLine, Chip) %>%
  group_modify(~{
    df <- .x
    group_keys <- .y
    n_cells <- nrow(df)
    n_full <- n_cells %/% cells_per_bulk
    n_leftover <- n_cells %% cells_per_bulk
    
    if (n_full > 0) {
      bulk_ids <- rep(1:n_full, each = cells_per_bulk)[1:n_cells]
    } else {
      bulk_ids <- rep(1, n_cells)
    }
    
    if (n_leftover > 0 && n_full > 0) {
      leftover_idx <- (n_full * cells_per_bulk + 1):n_cells
      redistribute <- rep(1:n_full, length.out = n_leftover)
      bulk_ids[leftover_idx] <- redistribute
    }
    
    df$BulkIndex <- bulk_ids
    df$PseudoBulkID <- paste(group_keys$CellLine, group_keys$Chip,
                             sprintf("bulk%02d", df$BulkIndex),
                             sep = "_")
    df
  }) %>%
  ungroup()

# -------------------------------------------------------------------------
# Step 4.5 (new): Aggregate RAW intensities by pseudobulk (sum, not mean)
# -------------------------------------------------------------------------
expr_mat_bulk <- expr_mat %>%
  as.data.frame() %>%
  rownames_to_column(var = "PROTID") %>%
  pivot_longer(-PROTID, names_to = "SampleID", values_to = "Intensity_raw") %>%
  left_join(coldata %>% dplyr::select(SampleID, PseudoBulkID), by = "SampleID") %>%
  group_by(PROTID, PseudoBulkID) %>%
  summarize(SumIntensity = sum(Intensity_raw, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = PseudoBulkID, values_from = SumIntensity) %>%
  column_to_rownames(var = "PROTID")

# Replace 0 (no detection) with NA
expr_mat_bulk[expr_mat_bulk == 0] <- NA_real_

# -------------------------------
# Step 4.6: Pseudobulk metadata
# -------------------------------
coldata_bulk <- coldata %>%
  group_by(PseudoBulkID) %>%
  summarize(
    CellLine = unique(CellLine)[1],
    Chip = unique(Chip)[1],
    nCells = n(),
    SingleCellIDs = paste(SampleID, collapse = ";"),
    .groups = "drop"
  ) %>%
  column_to_rownames(var = "PseudoBulkID")

# -------------------------------
# Step 5: PROTID → Gene symbol mapping
# -------------------------------
UniAc <- rownames(expr_mat_bulk)

geneSymbols <- select(org.Mm.eg.db, keys=UniAc, columns= c("SYMBOL","UNIPROT"),
                      keytype="UNIPROT", multiVals="first")

geneSymbols <- geneSymbols %>%
  mutate(Gene = SYMBOL) %>%
  dplyr::select(Gene,UNIPROT)

convert <- read_tsv("Github_data/convert_uniprot.tsv")

geneSymbols <- full_join(geneSymbols, convert, by = c("Gene","UNIPROT")) %>%
  filter(!is.na(Gene))

expr_mat_long <- expr_mat_bulk %>%
  as.data.frame() %>%
  rownames_to_column(var = "UNIPROT") %>%
  pivot_longer(-UNIPROT, names_to = "PseudoBulkID", values_to = "Intensity") %>%
  left_join(geneSymbols, by = "UNIPROT") %>%
  mutate(Gene = ifelse(is.na(Gene), UNIPROT, Gene)) %>%
  group_by(PseudoBulkID, Gene) %>%
  slice_max(order_by = Intensity, n = 1, with_ties = FALSE) %>%
  ungroup()

expr_mat_bulk <- expr_mat_long %>%
  dplyr::select(PseudoBulkID, Gene, Intensity) %>%
  pivot_wider(names_from = PseudoBulkID, values_from = Intensity) %>%
  column_to_rownames(var = "Gene") %>%
  as.matrix()

# -------------------------------
# Step 6: Log-transform + normalization (MOFA-friendly)
# -------------------------------
offset <- 1e-3
expr_log <- log2(expr_mat_bulk + offset)

# Median-center per pseudobulk (remove total intensity differences)
sample_medians <- apply(expr_log, 2, median, na.rm = TRUE)
expr_centered <- sweep(expr_log, 2, sample_medians, FUN = "-")

# Filter proteins detected in >=10% of pseudobulks
keep_genes <- rowMeans(!is.na(expr_centered)) >= 0.10
expr_centered <- expr_centered[keep_genes, , drop = FALSE]

# -------------------------------
# Step 7: Feature selection + scaling
# -------------------------------
feature_var <- rowVars(expr_centered, na.rm = TRUE)
feature_var <- feature_var[!is.na(feature_var) & feature_var > 0]
N_proteins <- min(1000, length(feature_var))
top_features <- names(sort(feature_var, decreasing = TRUE))[seq_len(N_proteins)]
prot_matrix <- expr_centered[top_features, , drop = FALSE]

# Z-score per protein (global center/scale), keep NAs
row_means <- rowMeans(prot_matrix, na.rm = TRUE)
row_sds <- rowSds(prot_matrix, na.rm = TRUE)
row_sds[row_sds == 0 | is.na(row_sds)] <- 1
prot_matrix_z <- sweep(prot_matrix, 1, row_means, FUN = "-")
prot_matrix_z <- sweep(prot_matrix_z, 1, row_sds, FUN = "/")

# -------------------------------
# Step 8: Export for MOFA+
# -------------------------------
write_csv(data.frame(Gene = rownames(prot_matrix_z), prot_matrix_z, check.names = FALSE),
          "MOFA/input/2_Pseudobulk_different_composition/proteomics_pdbk_different_for_MOFA.csv")

coldata_bulk_export <- coldata_bulk %>%
  rownames_to_column(var = "PseudoBulkID")

write_csv(coldata_bulk_export,
          "MOFA/input/2_Pseudobulk_different_composition/proteomics_pdbk_different_metadata.csv")

# -------------------------------
# Step 9: PCA diagnostic (impute only for plotting)
# -------------------------------
prot_imputed <- prot_matrix_z
for (i in 1:nrow(prot_imputed)) {
  row_vals <- prot_imputed[i, ]
  row_mean <- mean(row_vals, na.rm = TRUE)
  row_vals[is.na(row_vals)] <- row_mean
  prot_imputed[i, ] <- row_vals
}

pca <- prcomp(t(prot_imputed), center = FALSE, scale. = FALSE)
var_expl <- pca$sdev^2 / sum(pca$sdev^2)

plot(pca$x[,1], pca$x[,2],
     col = ifelse(coldata_bulk$CellLine == "C10", "red", "blue"),
     pch = 16,
     xlab = paste0("PC1 (", round(var_expl[1]*100,1), "%)"),
     ylab = paste0("PC2 (", round(var_expl[2]*100,1), "%)"),
     main = "Proteomics PCA (imputed): C10 (red) vs SVEC (blue)")

#### END PIPELINE ####
