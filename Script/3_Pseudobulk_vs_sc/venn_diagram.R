library(dplyr)
library(VennDiagram)
library(MOFA2)

setwd("C:/Users/49152/Downloads/Multi-omics/MOFA/Model_correlation")

# Load models
m_single <- load_model("C:/Users/49152/Downloads/Multi-omics/MOFA/output/single_cell_trained_model.hdf5")
m_pbulk <- load_model("C:/Users/49152/Downloads/Multi-omics/MOFA/output/pdbk_different_trained_model.hdf5")

view <- "Proteomics"
fac <- 3

# Get loadings for factor 1 for transcriptomics
# You can change "RNA" to your view name if different
weights_single <- get_weights(m_single, views = view, factors = fac, as.data.frame = TRUE)
weights_pbulk  <- get_weights(m_pbulk,  views = view, factors = fac, as.data.frame = TRUE)

top_n <- 100

top_genes_single <- weights_single %>%
  arrange(desc(abs(value))) %>%
  dplyr::slice(1:top_n) %>%
  pull(feature)

top_genes_pbulk <- weights_pbulk %>%
  arrange(desc(abs(value))) %>%
  dplyr::slice(1:top_n) %>%
  pull(feature)

# Compare overlap
length(intersect(top_genes_single, top_genes_pbulk))
# optionally fraction of overlap
length(intersect(top_genes_single, top_genes_pbulk)) / top_n

# Set the title dynamically based on the view
venn_title <- paste0("Top 100 genes contributing to factor ", fac , " (", view, " view)")

venn.plot <- venn.diagram(
  x = list(SingleCell = top_genes_single, Pseudobulk = top_genes_pbulk),
  filename = paste0("f", fac , "_venn_", view, ".png"),
  height = 3000,
  width = 3000,
  resolution = 300,
  fill = c("red", "blue"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 1.5,
  main = venn_title,
  category.names = c("SingleCell", "Pseudobulk"),
  cat.pos = c(0, 180),
  cat.dist = c(0.03, 0.03)
)
