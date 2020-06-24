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

load(file = "cc.genes.rda") 
mito_genes = c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
samples = c("D4-A1-H1", "D4-A1-H2", "D4-A1-H3", "D4-A1-H4", "D4-A1-H5", "D7-B1-H1", "D7-B1-H2", "D7-B1-H3", "D7-B1-H4", "D10-C1-H1", "D10-C1-H2", "D14-D1")
samples = c("D4-A1", "D7-B1", "D10-C1", "D14-D1")
cc.genes <- NULL

prepare_datasets <- function(folder_name, sample_name, cc.genes, mito_genes){
  if(is.null(sample_name)){
    sample_name <- folder_name
  }
  seurat.object <- Load10X_Spatial(paste("spaceranger_count_output_6.0_round2/", folder_name, "/outs/", sep = ""));
  seurat.object$orig.ident <- sample_name
  seurat.object@project.name <- sample_name
  # Mitrocondria
  seurat.object <- PercentageFeatureSet(object = seurat.object, features = mito_genes, col.name = "percent.mito");
  return(seurat.object)
}

day4_visium <- prepare_datasets("D4-A1", "D4", cc.genes, mito_genes)
day7_visium <- prepare_datasets("D7-B1", "D7", cc.genes, mito_genes)
day10_visium <- prepare_datasets("D10-C1", "D10", cc.genes, mito_genes)
day14_visium <- prepare_datasets("D14-D1", "D14", cc.genes, mito_genes)
# save.image("robjs/all.visiums.4.RData")
# load("robjs/all.visiums.4.RData")

day4_h1 = prepare_datasets(samples[1], "D4-H1", cc.genes, mito_genes)
day4_h2 = prepare_datasets(samples[2], "D4-H2", cc.genes, mito_genes)
day4_h3 = prepare_datasets(samples[3], "D4-H3", cc.genes, mito_genes)
day4_h4 = prepare_datasets(samples[4], "D4-H4", cc.genes, mito_genes)
day4_h5 = prepare_datasets(samples[5], "D4-H5", cc.genes, mito_genes)
day7_h1 = prepare_datasets(samples[6], "D7-H1", cc.genes, mito_genes)
day7_h2 = prepare_datasets(samples[7], "D7-H2", cc.genes, mito_genes)
day7_h3 = prepare_datasets(samples[8], "D7-H3", cc.genes, mito_genes)
day7_h4 = prepare_datasets(samples[9], "D7-H4", cc.genes, mito_genes)
day10_h1 = prepare_datasets(samples[10], "D10-H1", cc.genes, mito_genes)
day10_h2 = prepare_datasets(samples[11], "D10-H2", cc.genes, mito_genes)
day14_h1 = prepare_datasets(samples[12], "D14-H1", cc.genes, mito_genes)
# save.image("robjs/all.visiums.12.RData")
# load("robjs/all.visiums.12.RData")

# Merge all visium seurat objects in one with the images
# chicken_visium <- merge(day4_h1, y = c(day4_h2, day4_h3, day4_h4, day4_h5, day7_h1, day7_h2, day7_h3, day7_h4, day10_h1, day10_h2, day14_h1), add.cell.ids = samples)
chicken_visium <- merge(day4_visium, y = c(day7_visium, day10_visium, day14_visium), add.cell.ids = c("D4-A1", "D7-B1", "D10-C1", "D14-D1"))

# Preprocessing of spatial RNAseq data
chicken_visium <- NormalizeData(chicken_visium) %>% FindVariableFeatures() %>% ScaleData()
dim(chicken_visium)
DefaultAssay(chicken_visium)

# Run PCA on spatial RNAseq data
chicken_visium <- RunPCA(object = chicken_visium)
ElbowPlot(chicken_visium)

# Run clustering and UMP reduction of spatial RNAseq data
chicken_visium <- FindNeighbors(object = chicken_visium, dims=1:20, force.recalc = TRUE)
chicken_visium <- FindClusters(object = chicken_visium, resolution=1.0)
chicken_visium <- RunUMAP(object = chicken_visium, dims = 1:20)
DimPlot(chicken_visium, reduction = "umap", group.by = "Spatial_snn_res.1", label = TRUE)
DimPlot(chicken_visium, reduction = "umap", group.by = "orig.ident")

save(chicken_visium, file = "robjs/chicken_visium.4.Robj")
load("robjs/chicken_visium.Robj")


################### This part includes SC Tranform based normalisation (optional) #######################

# Normalise using SC Transform method in Seurat v-3 
chicken_visium <- SCTransform(chicken_visium, assay = "Spatial")


# SCT creates a new SCT assay and run PCA on SCT values
DefaultAssay(chicken_visium)
chicken_visium <- RunPCA(object = chicken_visium, reduction.name = "SCT_pca")

# CLustering and UMAP dimesntion reduction using SCT values
chicken_visium <- FindNeighbors(object = chicken_visium, reduction = "SCT_pca", graph.name = "SCT_snn", dims=1:20, force.recalc = TRUE)
chicken_visium <- FindClusters(object = chicken_visium, graph.name = "SCT_snn", resolution=1.0)
table(chicken_visium$seurat_clusters)
chicken_visium <- RunUMAP(object = chicken_visium, reduction = "SCT_pca", reduction.name = "SCT_umap", dims = 1:20)
DimPlot(chicken_visium, reduction = "SCT_umap", group.by = "SCT_snn_res.1",label = TRUE)
DimPlot(chicken_visium, reduction = "SCT_umap", group.by = "orig.ident")


