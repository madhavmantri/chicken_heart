# Preprocessing
# Load the requisite packages and some additional helper functions.
library(Seurat)
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

load("robjs/chicken_raw.Robj")

dim(chicken)
table(chicken$orig.ident)
DefaultAssay(chicken)

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
chicken <- RunPCA(object = chicken, assay = "scanorama", verbose = F, reduction.name = "pca_scanorama")

# Clustering and UMAP reduction on scanorama values
chicken <- FindNeighbors(object=chicken, assay = "scanorama", reduction = "pca_scanorama", dims = 1:20, k.param = 30, force.recalc = TRUE)
chicken <- FindClusters(object=chicken, resolution=0.5)
table(chicken$scanorama_snn_res.0.5)
# This command should give seventeen clusters to replicate the analysis
# table(chicken$scanorama_snn_res.0.5)
# 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16 
# 3489 3487 2998 2272 1847 1844 1075 1053  790  741  686  615  471  352  265  264   66 

chicken <- RunUMAP(object = chicken, assay = "scanorama", reduction = "pca_scanorama", dims = 1:20, reduction.name = "umap_scanorama", metric = "correlation")

# To visualise
DimPlot(chicken, reduction = "umap_scanorama", label = TRUE, group.by = "scanorama_snn_res.0.5")
DimPlot(chicken, reduction = "umap_scanorama", group.by = "orig.ident")

# save(chicken, file="robjs/chicken_normalised_scanorama2.Robj")
load("robjs/chicken_normalised_scanorama2.Robj")


# Differnetial marker analysis
DefaultAssay(object = chicken) <- "RNA"
markers.all = FindAllMarkers(chicken, assay = "RNA", do.print = TRUE, logfc.threshold = 1.5, return.thresh = 0.1, min.pct = 0.5, only.pos = TRUE)
markers.top10 = markers.all %>% group_by(cluster) %>% top_n(10, avg_logFC)
markers.top20 = markers.all %>% group_by(cluster) %>% top_n(20, avg_logFC)


# Use differential marker analysis to label clusters with cell type names
Idents(chicken) <- chicken$scanorama_snn_res.0.5
DimPlot(object = chicken, reduction = "umap", label = TRUE)
chicken <- RenameIdents(chicken, `0` = "Fibroblast cells", `1` = "Cardiomyocyte cells", `2` = "Myocardial progenitor cells",
                                   `3` = "Endothelial cells", `4` = "Cardiomyocyte precursor cells", `5` = "Endocardial valve cells",
                                   `6` = "TMSB4X high cells",
                                   `7` = "Epicardial progenitor cells-1", `8` = "Erythrocytes", `9` = "Vascular endothelial cells", `10` = "Erythrocytes",
                                   `11` = "Smooth muscle cells", `12` = "Epicardial progenitor cells-2", `13` = "Mitochondria enriched cardiomyocytes", `14` = "Macrophages",
                                   `15` = "Erythrocytes", `16` = "Dendritic cells")
chicken$celltypes.0.5 <- Idents(chicken)

chicken.integrated <- chicken

# save(chicken.integrated, file = "robjs/chicken_normalised_scanorama3.Robj")
load("robjs/chicken_normalised_scanorama3.Robj")

#############################  This section includes PHATE reduction for the entire dataset (optional) ####################

library(phateR)
data <- t(GetAssayData(chicken.integrated, assay = "scanorama", slot = "data"))
tree_chicken <- phate(data, gamma = 0)
chicken.integrated@reductions[["scanorama_phate_gamma0"]] <- CreateDimReducObject(embeddings = tree_chicken$embedding, assay = "scanorama", key = "PHATE")
DimPlot(chicken.integrated, reduction = "scanorama_phate_gamma0", label = TRUE)
tree_chicken <- phate(data, gamma = 1)
chicken.integrated@reductions[["scanorama_phate_gamma1"]] <- CreateDimReducObject(embeddings = tree_chicken$embedding, assay = "scanorama", key = "PHATE")
DimPlot(chicken.integrated, reduction = "scanorama_phate_gamma1")


