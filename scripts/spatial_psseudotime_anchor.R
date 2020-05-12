# setRepositories(ind = 1:2)
# options(unzip = "internal")
# devtools::install_github("satijalab/seurat", ref = "spatial")
# library(Seurat, lib.loc = "/home/mm2937/miniconda3/r-seurat-3.1.1-r35h0357c0b_0/lib/R/library/")
library(Seurat)
packageVersion("Seurat")
library(Matrix); library(stringr)
library(readr); library(here)
library(fitdistrplus); library(dplyr)
library(SeuratData); library(ggplot2)
library(cowplot); library(reticulate)
library(pals); library(monocle)
setwd("/workdir/mm2937/chicken/") 

# load(here("robjs", "chicken_normalised_scanorama3.Robj"))
load("robjs/chicken_visium.4.prediction.1.Robj")

prediction.scores <- as.data.frame(t(GetAssayData(chicken_visium, assay = "predictions_cells")))

load("robjs/epicardial.Robj")
load("robjs/endocardial.Robj")
load("robjs/myocardial.Robj")

dim(prediction.scores)

for(i in 1:nrow(prediction.scores)){
  temp <- str_replace_all(str_replace_all(colnames(prediction.scores)[prediction.scores[i,1:22315] > 0.0], pattern = "V-", replacement = "V_") , pattern = "D4-", replacement = "D4_")
  temp <- temp[temp %in% colnames(epicardial)]
  prediction.scores$epicardial_PHATE1[i] <- mean(epicardial@reductions$scanorama_phate_gamma0@cell.embeddings[,"PHATE_1"][temp])
}
chicken_visium$epicardial_PHATE1 <- prediction.scores$epicardial_PHATE1

for(i in 1:nrow(prediction.scores)){
  temp <- str_replace_all(str_replace_all(colnames(prediction.scores)[prediction.scores[i,1:22315] > 0.0], pattern = "V-", replacement = "V_") , pattern = "D4-", replacement = "D4_")
  temp <- temp[temp %in% colnames(endocardial)]
  prediction.scores$endocardial_PHATE1[i] <- mean(endocardial@reductions$scanorama_phate_gamma0@cell.embeddings[,"PHATE_1"][temp])
}
chicken_visium$endocardial_PHATE1 <- prediction.scores$endocardial_PHATE1

for(i in 1:nrow(prediction.scores)){
  temp <- str_replace_all(str_replace_all(colnames(prediction.scores)[prediction.scores[i,1:22315] > 0.0], pattern = "V-", replacement = "V_") , pattern = "D4-", replacement = "D4_")
  temp <- temp[temp %in% colnames(myocardial)]
  prediction.scores$myocardial_PHATE1[i] <- mean(myocardial@reductions$scanorama_phate_gamma0@cell.embeddings[,"PHATE_1"][temp])
}
chicken_visium$myocardial_PHATE1 <- prediction.scores$myocardial_PHATE1

dim(prediction.scores)
prediction.scores$max <- NULL

load("robjs/epi_cds.Robj")
for(i in 1:nrow(prediction.scores)){
  temp <- str_replace_all(str_replace_all(colnames(prediction.scores)[prediction.scores[i,1:22315] > 0.0], pattern = "V-", replacement = "V_") , pattern = "D4-", replacement = "D4_")
  temp <- temp[temp %in% colnames(epicardial)]
  prediction.scores$epicardial_pseudotime[i] <- mean(pData(cds)[temp, "Pseudotime"])
}
chicken_visium$epicardial_pseudotime <- prediction.scores$epicardial_pseudotime

load("robjs/endo_cds.Robj")
for(i in 1:nrow(prediction.scores)){
  temp <- str_replace_all(str_replace_all(colnames(prediction.scores)[prediction.scores[i,1:22315] > 0.0], pattern = "V-", replacement = "V_") , pattern = "D4-", replacement = "D4_")
  temp <- temp[temp %in% colnames(endocardial)]
  prediction.scores$endocardial_pseudotime[i] <- mean(pData(cds)[temp, "Pseudotime"])
}
chicken_visium$endocardial_pseudotime <- prediction.scores$endocardial_pseudotime

load("robjs/myo_cds.Robj")
for(i in 1:nrow(prediction.scores)){
  temp <- str_replace_all(str_replace_all(colnames(prediction.scores)[prediction.scores[i,1:22315] > 0.0], pattern = "V-", replacement = "V_") , pattern = "D4-", replacement = "D4_")
  temp <- temp[temp %in% colnames(myocardial)]
  prediction.scores$myocardial_pseudotime[i] <- mean(pData(cds)[temp, "Pseudotime"])
}
chicken_visium$myocardial_pseudotime <- prediction.scores$myocardial_pseudotime

# save(chicken_visium, file="robjs/chicken_visium.4.prediction.1.Robj")
load("robjs/chicken_visium.4.prediction.1.Robj")

# These spatial PHATE and pseudotime can be used for visualisation on spatial maps using SpatialFeaturePlot