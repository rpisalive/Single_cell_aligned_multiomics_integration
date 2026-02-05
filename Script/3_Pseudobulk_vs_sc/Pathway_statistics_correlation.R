library(dplyr)
library(ggplot2)
library(cowplot)
library(MOFAdata)
library(clusterProfiler)
library(org.Mm.eg.db)  # mouse
library(msigdbr)
library(MOFA2)

# ==========================
# Set file paths
# ==========================
setwd("C:/Users/49152/Downloads/Multi-omics/MOFA/")

# ==========================
# Load MOFA models
# ==========================
m_single <- load_model("output/single_cell_trained_model.hdf5")
m_pbulk  <- load_model("output/pdbk_different_trained_model.hdf5")

# ==========================
# Set parameters
# ==========================
view <- "Transcriptomics"
factors_to_analyze <- 1:7
sign <- "positive"
alpha <- 0.1

# Define suffix for removing from feature names
suffix <- if (view == "Transcriptomics") {
  "_RNA$"
} else if (view == "Proteomics") {
  "_PROT$"
}

# ==========================
# Load MSigDB C5 for mouse
# ==========================
data("MSigDB_v6.0_C5_mouse")

# ==========================
# Preprocess feature names
# ==========================
features_names(m_pbulk)[[view]] <- toupper(features_names(m_pbulk)[[view]])
features_names(m_pbulk)[[view]] <- sub(suffix, "", features_names(m_pbulk)[[view]])
head(features_names(m_pbulk)[[view]])

features_names(m_single)[[view]] <- toupper(features_names(m_single)[[view]])
features_names(m_single)[[view]] <- sub(suffix, "", features_names(m_single)[[view]])
head(features_names(m_single)[[view]])

# ==========================
# Run enrichment for single-cell model
# ==========================
enrichment_single <- run_enrichment(
  m_single,
  view = view,
  factors = factors_to_analyze,
  feature.sets = MSigDB_v6.0_C5_mouse,
  sign = sign,
  statistical.test = "parametric",
  alpha = alpha
)

# ==========================
# Run enrichment for pseudobulk model
# ==========================
enrichment_pbulk <- run_enrichment(
  m_pbulk,
  view = view,
  factors = factors_to_analyze,
  feature.sets = MSigDB_v6.0_C5_mouse,
  sign = sign,
  statistical.test = "parametric",
  alpha = alpha
)

plot_nes_correlation <- function(enrich_single, enrich_pbulk, fac = 1, title_suffix = "", alpha = 0.4){
  
  # Extract set statistics (NES-like)
  nes_single <- enrich_single$set.statistics[, fac]
  nes_pbulk  <- enrich_pbulk$set.statistics[, fac]
  
  # Keep only common pathways
  common_pathways <- intersect(names(nes_single), names(nes_pbulk))
  nes_single <- nes_single[common_pathways]
  nes_pbulk  <- nes_pbulk[common_pathways]
  
  # Data frame for plotting
  df <- data.frame(
    SingleCell = nes_single,
    Pseudobulk = nes_pbulk,
    Pathway = common_pathways
  )
  
  # Significant pathways in both
  sig_common <- intersect(
    enrich_single$sigPathways[[fac]],
    enrich_pbulk$sigPathways[[fac]]
  )
  df$Significant <- ifelse(df$Pathway %in% sig_common, "Yes", "No")
  
  # Plot
  p <- ggplot(df, aes(x = SingleCell, y = Pseudobulk, color = Significant)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    theme_minimal() +
    labs(
      x = "NES / Statistic (Single-Cell)",
      y = "NES / Statistic (Pseudobulk)",
      title = paste("NES correlation for Factor", fac, title_suffix)
    ) +
    annotate(
      "text",
      x = min(df$SingleCell, na.rm = TRUE),
      y = max(df$Pseudobulk, na.rm = TRUE),
      label = paste0("r = ", round(cor(df$SingleCell, df$Pseudobulk, use = "complete.obs"), 2),
                     ", alpha = ", alpha),
      hjust = 0, vjust = 1, size = 5
    ) +
    scale_color_manual(values = c("Yes" = "red", "No" = "grey"))
  
  return(p)
}


for(fac in factors_to_analyze){
  p <- plot_nes_correlation(enrichment_single, enrichment_pbulk,
                            fac = fac,
                            title_suffix = paste("(", view, " view)", sep=""),
                            alpha = alpha)
  ggsave(paste0("Model_correlation/NES_correlation_factor",fac,"_",view,".png"),
         p, width = 8, height = 6, dpi = 300)
}
