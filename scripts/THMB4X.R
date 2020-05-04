# Preprocessing
# Load the requisite packages and some additional helper functions.
library(Seurat)
packageVersion("Seurat")
library(Matrix); library(stringr)
library(readr); library(here)
library(fitdistrplus); library(dplyr)
library("URD"); library(monocle)
setwd("/workdir/mm2937/chicken/")
library(future)
plan("multiprocess", workers = 64)
options(future.globals.maxSize = 3000 * 1024^2)

# load(here("robjs", "chicken_normalised_scanorama3.Robj"))

DimPlot(chicken.integrated, reduction = "umap", label = TRUE)
DimPlot(chicken.integrated, reduction = "umap", group.by = "orig.ident")
Idents(chicken.integrated) <- chicken.integrated$celltypes.0.5

cluster6_markers = FindMarkers(chicken.integrated, assay = "RNA", ident.1 = "TMSB4X high cells", logfc.threshold = 0.0, only.pos = TRUE)
# write.table(cluster6_markers, file = "Cluster6_all_markers.csv", sep = ",")
# cluster6_markers <- read.csv("Cluster6_all_markers.csv")

library(EnhancedVolcano)
EnhancedVolcano(cluster6_markers,
                lab = rownames(cluster6_markers),
                x = 'avg_logFC',
                y = 'p_val_adj', title = NULL, subtitle = NULL, caption = NULL, legendPosition = "none", pCutoff = 10^-15, transcriptLabSize = 2, 
                FCcutoff = 0.6, gridlines.major = FALSE, gridlines.minor = FALSE, colAlpha = 0.8) +
                theme_bw()+

  
# Creating a seuart object for TMSB4X high cells                  
subcluster <- subset(chicken.integrated, subset = celltypes.0.5 %in% c("TMSB4X high cells"))
DimPlot(object = subcluster, reduction = "umap", label = TRUE)

# Preprocess the expression vlaues for the subcluster 
DefaultAssay(subcluster) <- "RNA"
subcluster <- FindVariableFeatures(subcluster)
subcluster = ScaleData(subcluster)

# Run PCA on the TMSB4X cell sublcuster
subcluster <- RunPCA(object = subcluster)
ElbowPlot(object = subcluster, ndims = 50)
n.pcs <- 20;

# Run Clustering and UMAP reduction
subcluster <- FindNeighbors(object = subcluster, reduction = "pca", dims = 1:n.pcs, k = 30, force.recalc = TRUE, verbose = T)
subcluster <- FindClusters(object = subcluster, resolution = 0.15)
subcluster <- RunUMAP(object = subcluster, reduction.use = "pca", dims = 1:n.pcs)

# To visualise
DimPlot(subcluster, reduction = "umap", label = TRUE) +  theme(legend.position = "none")
DimPlot(subcluster, reduction = "umap", group.by = "orig.ident", label = TRUE)

# save(subcluster, file = "robjs/thsb4x_subcluster.Robj")
# load("robjs/thsb4x_subcluster.Robj")

# Differential marker analysis on subclusters in TMSB4X cells 
markers.all = FindAllMarkers(subcluster, assay = "RNA", do.print = TRUE, logfc.threshold = 0.5, return.thresh = 0.01, min.diff.pct = 0.5, only.pos = TRUE)
# markers.all <- subset(markers.all[!(rownames(markers.all) %in% grep("^ENSGAL", x = rownames(markers.all), value = TRUE)),])
markers.top10 = markers.all %>% group_by(cluster) %>% top_n(10, avg_logFC)
markers.top20 = markers.all %>% group_by(cluster) %>% top_n(20, avg_logFC)
write_csv(markers.all, "thymosin_subclusters_markers.csv")

subcluster <- RenameIdents(subcluster, "0" = "Vascular smooth muscle- like", "1" = "Vascular endothelial- like", "2" = "Cardiomyocytes_like")

DefaultAssay(subcluster) <- "RNA"
chicken.integrated$celltypes_manual <- as.character(chicken.integrated$celltypes.0.5)
DimPlot(subcluster)
subcluster <- RenameIdents(subcluster, "0" = "Vascular smooth muscle- like", "1" = "Vascular endothelial- like", "2" = "Cardiomyocytes_like")
subcluster$seurat_clusters <- Idents(subcluster)
subcluster$seurat_clusters <- as.character(paste("TMSB4X_high_", subcluster$seurat_clusters, sep = ""))
Idents(subcluster) <- subcluster$seurat_clusters
table(subcluster$seurat_clusters)
chicken.integrated$celltypes_manual[colnames(subcluster)] <- subcluster$seurat_clusters
