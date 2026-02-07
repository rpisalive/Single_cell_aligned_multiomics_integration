# =========================================================
# LOAD LIBRARIES
# =========================================================
library(MOFA2)
library(dplyr)
library(MOFAdata)
library(ggplot2)
library(tidyr)
library(purrr)
library(pheatmap)
library(tibble)

# =========================================================
# LOAD MODELS
# =========================================================
m_single <- load_model("C:/Users/49152/Downloads/Multi-omics/MOFA/output/single_cell_trained_model.hdf5")
m_scot   <- load_model("C:/Users/49152/Downloads/Multi-omics/SCOT_plus/MOFA_output/SCOTplus_trained_model.hdf5")

# =========================================================
# LOAD MSigDB GENE SETS
# =========================================================
data("MSigDB_v6.0_C5_mouse")
gene_sets <- MSigDB_v6.0_C5_mouse

# =========================================================
# NORMALIZE FEATURE NAMES TO MATCH MSigDB
# =========================================================
normalize_features <- function(model){
  for(v in names(features_names(model))){
    feats <- features_names(model)[[v]]
    feats <- toupper(feats)                 # match MSigDB
    feats <- sub("_RNA$", "", feats)        # remove suffix
    feats <- sub("_PROT$", "", feats)
    features_names(model)[[v]] <- feats
  }
  return(model)
}

m_single <- normalize_features(m_single)
m_scot   <- normalize_features(m_scot)

# =========================================================
# PARAMETERS
# =========================================================
views <- c("Transcriptomics", "Proteomics")
signs <- c("positive", "negative")
n_features <- 100  # top features per factor
K_single <- get_dimensions(m_single)$K
K_scot   <- get_dimensions(m_scot)$K

# =========================================================
# HELPER FUNCTION: GET TOP FEATURES PER FACTOR
# =========================================================
top_features_factor <- function(model, view, factor, n_features = 100, sign = "positive"){
  w <- get_weights(model, views = view, factors = factor, as.data.frame = TRUE)
  if(nrow(w) == 0) return(character(0))
  
  w$value <- as.numeric(w$value)
  
  if(sign == "positive") w <- w[w$value > 0, ]
  if(sign == "negative") w <- w[w$value < 0, ]
  
  if(nrow(w) == 0) return(character(0))
  
  # Robust base R ranking
  ord <- order(abs(w$value), decreasing = TRUE)
  top_feats <- w[ord[seq_len(min(n_features, nrow(w)))], ]
  
  return(toupper(top_feats$feature))
}

# =========================================================
# FUNCTION TO RUN ENRICHMENT PER MODEL
# =========================================================
run_model_enrichment <- function(model, model_name){
  enrichment_list <- list()
  K <- get_dimensions(model)$K
  
  for(view in views){
    for(sign in signs){
      for(f in 1:K){
        features <- top_features_factor(model, view, f, n_features = n_features, sign = sign)
        if(length(features) == 0) next
        
        enrich_res <- tryCatch({
          run_enrichment(
            model,
            view = view,
            factors = f,
            feature.sets = gene_sets,
            sign = sign,
            statistical.test = "parametric"
          )
        }, error=function(e) NULL)
        
        if(!is.null(enrich_res)){
          enrichment_list[[paste(model_name, view, sign, f, sep="_")]] <- enrich_res
        }
      }
    }
  }
  
  return(enrichment_list)
}

# =========================================================
# RUN ENRICHMENT FOR BOTH MODELS
# =========================================================
enrich_single <- run_model_enrichment(m_single, "Single")
enrich_scot   <- run_model_enrichment(m_scot, "SCOTplus")

# =========================================================
# CHECK OUTPUT
# =========================================================
names(enrich_single)
names(enrich_scot)



# =========================================================
# HELPER FUNCTIONS
# =========================================================

# Get significant pathways from enrichment object
get_sig_pathways <- function(enrich_obj) {
  if(is.null(enrich_obj)) return(character(0))
  sig_list <- enrich_obj$sigPathways
  # Flatten to a single vector
  unique(unlist(sig_list))
}

# Compute Jaccard index between two sets
jaccard_index <- function(set1, set2) {
  if(length(set1) == 0 & length(set2) == 0) return(1)
  if(length(set1) == 0 | length(set2) == 0) return(0)
  length(intersect(set1, set2)) / length(union(set1, set2))
}

# =========================================================
# Flatten enrichment list
# =========================================================
flatten_enrichment <- function(enrich_list) {
  
  df_list <- lapply(names(enrich_list), function(nm) {
    obj <- enrich_list[[nm]]
    
    # Extract model, view, sign, factor from the list name
    parts <- strsplit(nm, "_")[[1]]
    model_name <- parts[1]
    view <- parts[2]
    sign <- parts[3]
    factor <- parts[4]
    
    paths <- get_sig_pathways(obj)
    if(length(paths) == 0) return(NULL)  # skip factors with no sig pathways
    
    data.frame(
      model  = model_name,
      view   = view,
      sign   = sign,
      factor = factor,
      pathway = paths,
      stringsAsFactors = FALSE
    )
  })
  
  df <- do.call(rbind, df_list)
  rownames(df) <- NULL
  return(df)
}


# =========================================================
# Flatten both enrichment lists
# =========================================================
df_single <- flatten_enrichment(enrich_single)
df_scot   <- flatten_enrichment(enrich_scot)

# =========================================================
# Combine and prepare for Jaccard calculation
# =========================================================
df_combined <- bind_rows(df_single, df_scot)

# Ensure factors are uniquely identified per model
df_combined <- df_combined %>%
  mutate(factor_id = paste(model, view, sign, factor, sep="_"))

# =========================================================
# Create pathway presence matrix
# =========================================================
pathway_mat <- df_combined %>%
  select(factor_id, pathway) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = pathway, values_from = value, values_fill = 0) %>%
  column_to_rownames("factor_id")

# =========================================================
# Compute Jaccard similarity between all factor pairs
# =========================================================
jaccard_similarity <- function(x, y) {
  intersection <- sum(x & y)
  union <- sum(x | y)
  if(union == 0) return(NA)
  intersection / union
}

factor_ids <- rownames(pathway_mat)
jaccard_mat <- matrix(0, nrow = length(factor_ids), ncol = length(factor_ids),
                      dimnames = list(factor_ids, factor_ids))

for(i in seq_along(factor_ids)) {
  for(j in seq_along(factor_ids)) {
    jaccard_mat[i,j] <- jaccard_similarity(pathway_mat[i,], pathway_mat[j,])
  }
}

# =========================================================
# Heatmap
# =========================================================

pheatmap(
  jaccard_mat,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Pathway overlap (Jaccard similarity) between Single vs SCOT+ factors",
  color = colorRampPalette(c("white", "red"))(50)
)

cross_mat <- jaccard_mat[grep("Single", rownames(jaccard_mat)),
                         grep("SCOTplus", colnames(jaccard_mat))]
pheatmap(cross_mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Cross-model pathway overlap (Jaccard)",
         color = colorRampPalette(c("white", "red"))(50))

binary_mat <- cross_mat
binary_mat[binary_mat < 0.1] <- 0
binary_mat[binary_mat >= 0.1] <- 1

pheatmap(binary_mat,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         main = "Conserved pathways (Jaccard > 0.1)",
         color = c("white", "red"))

df_jaccard <- as.data.frame(as.table(cross_mat))
df_jaccard %>%
  group_by(Var1, Var2) %>%
  summarise(jaccard = mean(Freq)) %>%
  ggplot(aes(x = Var1, y = Var2, fill = jaccard)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  theme(axis.text.x = element_text(angle=45, hjust=1))
