library(data.table)
library(purrr)
library(ggplot2)
library(cowplot)
library(MOFAdata)
library(dplyr)
library(MOFA2)

view <- "Proteomics"
fac <- 4
# Define suffix based on view
suffix <- if (view == "Transcriptomics") {
  "_RNA$"
} else if (view == "Proteomics") {
  "_PROT$"
}
sign <- "positive"
#alpha <- alpha
  
# C5: extracted from the Gene Ontology data.base
data("MSigDB_v6.0_C5_mouse")
head(rownames(MSigDB_v6.0_C5_mouse), n=5)
head(colnames(MSigDB_v6.0_C5_mouse), n=5)

setwd("C:/Users/49152/Downloads/Multi-omics/MOFA/")
model <- load_model("output/single_cell_trained_model.hdf5")
model@samples_metadata
model@samples_metadata <- model@samples_metadata %>%
  mutate(cell_line = sub("^[^_]+_([^_]+).*", "\\1", sample))
model@samples_metadata

features_names(model)[[view]] <- toupper(features_names(model)[[view]])
features_names(model)[[view]] <- sub(suffix, "", features_names(model)[[view]])
head(features_names(model)[[view]])

enrichment.parametric <- run_enrichment(model,
                                        view = view, factors = 1:6,
                                        feature.sets = MSigDB_v6.0_C5_mouse,
                                        sign = sign,
                                        statistical.test = "parametric")

names(enrichment.parametric)

enrichment.parametric$set.statistics[1:5,1]
enrichment.parametric$pval.adj[1:5,1]
plot_enrichment_heatmap(enrichment.parametric)
plot_enrichment(enrichment.parametric, 
                factor = fac, 
                max.pathways = 15)
plot_enrichment_detailed(enrichment.parametric, 
                         factor = fac, 
                         max.genes = 8, 
                         max.pathways = 5)

genes <- list("WNT5A","RNF2","MKNK1")

genes %>% map(~ plot_factors(model, 
                             factors = c(1,2), 
                             color_by = "cell_line", 
                             scale = T,
                             legend = T
)) %>% cowplot::plot_grid(plotlist=., nrow=1)

genes %>% map(~ plot_factors(model, 
                             factors = c(1,2), 
                             color_by = ., 
                             scale = T,
                             legend = T
)) %>% cowplot::plot_grid(plotlist=., nrow=1)

enrichment.parametric.adj <- run_enrichment(model,
                                            view = view, factors = 1:5,
                                            feature.sets = MSigDB_v6.0_C5_mouse,
                                            sign = sign,
                                            statistical.test = "cor.adj.parametric"
)

dt <- rbind(
  enrichment.parametric$pval[,1:3] %>% as.data.table %>% .[,c("test","pathway"):=list("parametric",1:.N)],
  enrichment.parametric.adj$pval[,1:3] %>% as.data.table %>% .[,c("test","pathway"):=list("parametric.adj",1:.N)]
) %>% melt(id.vars=c("test","pathway"), variable.name="factor")

ggplot(dt, aes(x=value, fill=test)) +
  facet_wrap(~factor, scales="free_y", nrow=1) +
  geom_histogram() +
  theme_bw() +
  theme(
    legend.position = "top",
    legend.title = element_blank()
  )
dt2 <- dt %>% dcast(factor+pathway~test)

ggplot(dt2, aes(x=parametric, y=parametric.adj)) +
  geom_point(size=0.5) +
  geom_abline(slope=1, intercept=0, color="orange") +
  facet_wrap(~factor, scales="free_y", nrow=1) +
  labs(x="Parametric p-value", y="Adjusted parametric p-value") +
  theme_bw() +
  theme(
    legend.position = "top"
  )
