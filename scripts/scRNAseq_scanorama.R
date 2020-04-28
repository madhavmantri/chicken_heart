# Preprocessing
# Load the requisite packages and some additional helper functions.
library(Seurat, lib.loc = "/home/mm2937/miniconda3/pkgs/r-seurat-3.1.1-r35h0357c0b_0/lib/R/library/")
packageVersion("Seurat")
library(Matrix)
library(stringr)
library(readr)
library(here)
library(fitdistrplus)
library(dplyr)
library(monocle)
library(reticulate)
setwd("/workdir/mm2937/chicken/")

load(file = "cc.genes.rda")
samples = c("D4", "D7-LV", "D7-RV", "D10-LV", "D10-RV", "D14-LV", "D14-RV")

load("robjs/chicken_raw_SCT.Robj")

dim(chicken)
table(chicken$orig.ident)

# Function to extract expression values for individual datasets
extractRNA_chicken <- function(seurat.object, sample_name){
  return(t(as.matrix(GetAssayData(seurat.object)))[colnames(seurat.object)[seurat.object$orig.ident == sample_name],])
}

# Prepare datasets and gene lists for scanorama integration
data = list(extractRNA_chicken(chicken, samples[1]), extractRNA_chicken(chicken, samples[2]), extractRNA_chicken(chicken, samples[3]), extractRNA_chicken(chicken, samples[4]), extractRNA_chicken(chicken, samples[5]), extractRNA_chicken(chicken, samples[6]), extractRNA_chicken(chicken, samples[7]))
gene_list = list(rownames(chicken), rownames(chicken), rownames(chicken), rownames(chicken), rownames(chicken), rownames(chicken), rownames(chicken))


# Intregration starts here
scanorama = import('scanorama')
integrated.corrected.data = scanorama$correct(data, gene_list, return_dimred=TRUE, return_dense=TRUE, ds_names = samples, verbose = TRUE)

# save(integrated.corrected.data, file="robjs/corrected_SCT_counts_scano.Robj")
# load("robjs/corrected_norm_scano.Robj")

corrected_scanorama <- t(do.call(rbind, integrated.corrected.data[[2]]))
colnames(corrected_scanorama) <- colnames(chicken)
rownames(corrected_scanorama) <- integrated.corrected.data[[3]]
dim(corrected_scanorama)
corrected_scanorama_pca <- t(do.call(rbind, integrated.corrected.data[[1]]))
colnames(corrected_scanorama_pca) <- colnames(chicken)
dim(corrected_scanorama_pca)


# Create assay from integrated values and save to seurat object
scanorama_assay <- CreateAssayObject(data = corrected_scanorama)
chicken[["scanorama"]] <- scanorama_assay
DefaultAssay(chicken) <- "scanorama"


# Preprocess scanorama values and perform PCA
DefaultAssay(chicken)
chicken <- FindVariableFeatures(chicken, assay = "scanorama", selection.method = "vst")
chicken <- ScaleData(chicken)
chicken <- RunPCA(object = chicken, assay = "scanorama", reduction.name = "pca_scanorama")

# Clustering and UMAP reduction on scanorama values
chicken <- FindNeighbors(object=chicken, dims=1:20, reduction = "pca_scanorama", force.recalc = TRUE, graph.name = "scanorama_snn")
chicken <- FindClusters(object=chicken,graph.name = "scanorama_snn", resolution=0.5)
table(chicken$scanorama_snn_res.0.4)
chicken <- RunUMAP(object = chicken, reduction.use = "pca_scanorama", dims = 1:20, reduction.name = "umap_scanorama")

# To visualise
DimPlot(chicken, reduction = "umap_scanorama", label = TRUE, group.by = "scanorama_snn_res.0.4")
DimPlot(chicken, reduction = "umap", group.by = "orig.ident")

# save(chicken, file="robjs/chicken_normalised_scanorama2.Robj")
# load("robjs/chicken_normalised_scanorama2.Robj")
