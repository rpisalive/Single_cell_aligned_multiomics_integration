library(dplyr)
library(ggplot2)
library(MOFA2)

setwd("C:/Users/49152/Downloads/Multi-omics/MOFA/")
model <- load_model("output/trained_model.hdf5")

# ---- Merge factors with your metadata ----
factors_obj <- get_factors(model)

if (is.list(factors_obj)) {
  factors_df <- as.data.frame(factors_obj[[1]])
} else {
  factors_df <- as.data.frame(factors_obj)
}
colnames(factors_df) <- paste0("Factor", seq_len(ncol(factors_df)))
factors_df$SampleID <- rownames(factors_df)

meta <- read.csv("input/transcriptomics_metadata.csv", stringsAsFactors = FALSE)
diag <- left_join(factors_df, meta, by="SampleID")

# ---- Correlation of factors with sizeFactor ----
cors <- sapply(colnames(factors_df)[1:(ncol(factors_df)-1)], function(f) {
  cor(diag[[f]], diag$sizeFactor, use="pairwise.complete.obs")
})
print(cors)

# ---- Plot Factor1 and Factor2 against sizeFactor ----
for(f in c("Factor1","Factor2","Factor3","Factor4","Factor5","Factor6","Factor7","Factor8","Factor9")) {
  if(f %in% colnames(diag)) {
    p <- ggplot(diag, aes(x=sizeFactor, y=.data[[f]], colour=CellLine)) +
      geom_point(size=3) + geom_smooth(method="lm", se=FALSE) +
      ggtitle(paste0(f, " vs sizeFactor"))
    print(p)
  }
}

# ---- Test association with chip/cell line ----
summary(lm(Factor1 ~ sizeFactor + Chip + CellLine, data=diag))
summary(lm(Factor2 ~ sizeFactor + Chip + CellLine, data=diag))
summary(lm(Factor3 ~ sizeFactor + Chip + CellLine, data=diag))
summary(lm(Factor4 ~ sizeFactor + Chip + CellLine, data=diag))
summary(lm(Factor5 ~ sizeFactor + Chip + CellLine, data=diag))
summary(lm(Factor6 ~ sizeFactor + Chip + CellLine, data=diag))
summary(lm(Factor7 ~ sizeFactor + Chip + CellLine, data=diag))
summary(lm(Factor8 ~ sizeFactor + Chip + CellLine, data=diag))
summary(lm(Factor9 ~ sizeFactor + Chip + CellLine, data=diag))
