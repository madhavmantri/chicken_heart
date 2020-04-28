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

load("robjs/chicken_visium.Robj")

dim(chicken_visium)
table(chicken_visium$orig.ident)

# Function to extract expression values for individual datasets
extractRNA_chicken <- function(seurat.object, sample_name, assay = "Spatial"){
  return(t(as.matrix(GetAssayData(seurat.object, assay = assay)))[colnames(seurat.object)[seurat.object$orig.ident == sample_name],])
}

# Prepare datasets and gene lists for scanorama integration
samples = unique(chicken_visium$orig.ident)
data = list(extractRNA_chicken(chicken_visium, samples[[1]]), extractRNA_chicken(chicken_visium, samples[[2]]), 
            extractRNA_chicken(chicken_visium, samples[[3]]), extractRNA_chicken(chicken_visium, samples[[4]]))
gene_list = list(rownames(chicken_visium), rownames(chicken_visium), rownames(chicken_visium), rownames(chicken_visium))


# Intregration starts here
scanorama = import('scanorama')
integrated.corrected.data = scanorama$correct(data, gene_list, return_dimred=TRUE, return_dense=TRUE, ds_names = samples, verbose = TRUE)

# save(integrated.corrected.data, file="robjs/corrected_Spatial_norm_scano.4.Robj")
load("robjs/corrected_Spatial_norm_scano.4.Robj")

corrected_scanorama <- t(do.call(rbind, integrated.corrected.data[[2]]))
colnames(corrected_scanorama) <- colnames(chicken_visium)
rownames(corrected_scanorama) <- integrated.corrected.data[[3]]
dim(corrected_scanorama)
corrected_scanorama_pca <- t(do.call(rbind, integrated.corrected.data[[1]]))
colnames(corrected_scanorama_pca) <- colnames(chicken_visium)
dim(corrected_scanorama_pca)

# Create assay from integrated values and save to seurat object
scanorama_assay <- CreateAssayObject(data = corrected_scanorama)
chicken_visium[["scanorama"]] <- scanorama_assay
DefaultAssay(chicken_visium) <- "scanorama"

# Preprocess scanorama values and perform PCA
chicken_visium <- FindVariableFeatures(chicken_visium, assay = "scanorama", selection.method = "vst")
chicken_visium <- ScaleData(chicken_visium)
chicken_visium <- RunPCA(object = chicken_visium, assay = "scanorama", reduction.name = "pca_scanorama")

# Clustering and UMAP reduction on scanorama values
chicken_visium <- FindNeighbors(object=chicken_visium,dims=1:20, reduction = "pca_scanorama", force.recalc = TRUE, graph.name = "scanorama_snn")
chicken_visium <- FindClusters(object=chicken_visium,graph.name = "scanorama_snn", resolution=1.0)
chicken_visium <- RunUMAP(object = chicken_visium, reduction = "pca_scanorama", dims = 1:20, reduction.name = "umap_scanorama")

# To visualise
DimPlot(chicken_visium, reduction = "umap_scanorama", group.by = "scanorama_snn_res.1", label = TRUE)
DimPlot(chicken_visium, reduction = "umap_scanorama", group.by = "orig.ident")

# save(chicken_visium, file="robjs/chicken_visium.4.Robj")
# load("robjs/chicken_visium.4.Robj")

