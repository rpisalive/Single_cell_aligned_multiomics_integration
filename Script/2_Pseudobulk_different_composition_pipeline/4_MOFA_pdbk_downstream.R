library(ggplot2)
library(dplyr)
library(MOFA2)

setwd("C:/Users/49152/Downloads/Multi-omics/MOFA/")
model <- load_model("output/pdbk_different_trained_model.hdf5")

plot_data_overview(model)
model@samples_metadata
model@samples_metadata <- model@samples_metadata %>%
  mutate(cell_line = sub("_.*", "", sample))
model@samples_metadata

# Total variance explained per view and group
model@cache$variance_explained$r2_total
model@cache$variance_explained$r2_per_factor
plot_variance_explained(model, x="view", y="factor", plot_total = FALSE)
plot_variance_explained(model, x="view", y="factor", plot_total = TRUE)

plot_factor(model, factor = 1:9, color_by = "cell_line")

p <- plot_factor(model, 
                 factors = c(1,2,3),
                 color_by = "cell_line",
                 dot_size = 3,        # change dot size
                 dodge = T,           # dodge points with different colors
                 legend = F,          # remove legend
                 add_violin = T,      # add violin plots,
                 violin_alpha = 0.25  # transparency of violin plots
)

# The output of plot_factor is a ggplot2 object that we can edit
p <- p + 
  scale_color_manual(values=c("C10"="black", "SVEC"="red")) +
  scale_fill_manual(values=c("C10"="black", "SVEC"="red"))

print(p)

plot_factors(model, 
             factors = 1:5,
             color_by = "cell_line"
)

plot_weights(model,
             view = "Proteomics",
             factor = 2,
             nfeatures = 20,     # Number of features to highlight
             scale = T,          # Scale weights from -1 to 1
             abs = F             # Take the absolute value?
)

plot_top_weights(model,
                 view = "Proteomics",
                 factor = 2,
                 nfeatures = 20
)

plot_data_heatmap(model,
                  view = "Proteomics",         # view of interest
                  factor = 1,             # factor of interest
                  features = 20,          # number of features to plot (they are selected by weight)
                  
                  # extra arguments that are passed to the `pheatmap` function
                  cluster_rows = TRUE, cluster_cols = TRUE,
                  show_rownames = TRUE, show_colnames = TRUE
)

plot_data_heatmap(model,
                  view = "Transcriptomics",         # view of interest
                  factor = 1,             # factor of interest
                  features = 15,          # number of features to plot (they are selected by weight)
                  
                  # extra arguments that are passed to the `pheatmap` function
                  cluster_rows = TRUE, cluster_cols = TRUE,
                  show_rownames = TRUE, show_colnames = TRUE
)

plot_data_scatter(model,
                  view = "Proteomics",         # view of interest
                  factor = 1,             # factor of interest
                  features = 5,           # number of features to plot (they are selected by weight)
                  add_lm = TRUE,          # add linear regression
                  color_by = "cell_line"
)

plot_variance_explained_per_feature(
  model,
  "Proteomics",
  features = 20)



set.seed(42)
model <- run_umap(model,n_neighbors = 25)
model <- run_umap(model,n_neighbors = 20, factors = c("Factor1","Factor2"))

plot_dimred(model,
            method = "UMAP",  # method can be either "TSNE" or "UMAP"
            color_by = "cell_line"
)

# "factors" is a list of matrices, one matrix per group with dimensions (nsamples, nfactors)
factors <- get_factors(model, factors = "all")
lapply(factors,dim)

# "weights" is a list of matrices, one matrix per view with dimensions (nfeatures, nfactors)
weights <- get_weights(model, views = "all", factors = "all")
lapply(weights,dim)

# "data" is a nested list of matrices, one matrix per view and group with dimensions (nfeatures, nsamples)
data <- get_data(model)
lapply(data, function(x) lapply(x, dim))[[1]]
