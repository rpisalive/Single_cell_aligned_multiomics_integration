library(dplyr)
library(ggplot2)
#library(cowplot)
library(MOFAdata)
#library(clusterProfiler)
#library(org.Mm.eg.db)  # mouse
#library(msigdbr)
library(MOFA2)

# ==========================
# Set file paths
# ==========================
setwd("C:/Users/49152/Downloads/Multi-omics/SCOT_plus/MOFA_output/")

# ==========================
# Load MOFA models
# ==========================
m_single <- load_model("C:/Users/49152/Downloads/Multi-omics/MOFA/output/single_cell_trained_model.hdf5")
m_scot  <- load_model("C:/Users/49152/Downloads/Multi-omics/SCOT_plus/MOFA_output/SCOTplus_trained_model.hdf5")

# ==========================
# Set parameters
# ==========================
view <- "Transcriptomics"
factors_to_analyze <- 1:3
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
features_names(m_scot)[[view]] <- toupper(features_names(m_scot)[[view]])
features_names(m_scot)[[view]] <- sub(suffix, "", features_names(m_scot)[[view]])
head(features_names(m_scot)[[view]])

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
# Run enrichment for SCOT+ model
# ==========================
enrichment_scot <- run_enrichment(
  m_scot,
  view = view,
  factors = factors_to_analyze,
  feature.sets = MSigDB_v6.0_C5_mouse,
  sign = sign,
  statistical.test = "parametric",
  alpha = alpha
)

plot_nes_correlation <- function(enrich_single, enrich_scot, fac = 1, title_suffix = "", alpha = alpha){
  
  # Extract set statistics (NES-like)
  nes_single <- enrich_single$set.statistics[, fac]
  nes_scot  <- enrich_scot$set.statistics[, fac]
  
  # Keep only common pathways
  common_pathways <- intersect(names(nes_single), names(nes_scot))
  nes_single <- nes_single[common_pathways]
  nes_scot  <- nes_scot[common_pathways]
  
  # Data frame for plotting
  df <- data.frame(
    SingleCell = nes_single,
    SCOT_plus = nes_scot,
    Pathway = common_pathways
  )
  
  # Significant pathways in both
  sig_common <- intersect(
    enrich_single$sigPathways[[fac]],
    enrich_scot$sigPathways[[fac]]
  )
  df$Significant <- ifelse(df$Pathway %in% sig_common, "Yes", "No")
  
  # Plot
  p <- ggplot(df, aes(x = SingleCell, y = SCOT_plus, color = Significant)) +
    geom_point(alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE, color = "black") +
    theme_minimal() +
    labs(
      x = "NES / Statistic (Single-Cell)",
      y = "NES / Statistic (SCOT+)",
      title = paste("NES correlation for Factor", fac, title_suffix)
    ) +
    annotate(
      "text",
      x = min(df$SingleCell, na.rm = TRUE),
      y = max(df$SCOT_plus, na.rm = TRUE),
      label = paste0("r = ", round(cor(df$SingleCell, df$SCOT_plus, use = "complete.obs"), 2),
                     ", alpha = ", alpha),
      hjust = 0, vjust = 1, size = 5
    ) +
    scale_color_manual(values = c("Yes" = "red", "No" = "grey"))
  
  return(p)
}


for(fac in factors_to_analyze){
  p <- plot_nes_correlation(enrichment_single, enrichment_scot,
                            fac = fac,
                            title_suffix = paste("(", view, " view)", sep=""),
                            alpha = alpha)
  ggsave(paste0("Model_correlation/NES_correlation_factor",fac,"_",view,".png"),
         p, width = 8, height = 6, dpi = 300)
}
