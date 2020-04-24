# Preprocessing
# Load the requisite packages and some additional helper functions.
library(Seurat, lib.loc = "/home/mm2937/miniconda3/pkgs/r-seurat-3.1.1-r35h0357c0b_0/lib/R/library")
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
mito_genes = c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
samples = c("D4", "D7-LV", "D7-RV", "D10-LV", "D10-RV", "D14-LV", "D14-RV")

cc.genes <- NULL

prepare_datasets <- function(sample_name, cc.genes, mito_genes){
  # Read data
  data <- Read10X(data.dir = paste("cellranger_count_output_6.0/", sample_name, "/outs/filtered_feature_bc_matrix/", sep = ""))
  seurat.object <- CreateSeuratObject(counts = data, min.cells = 1, min.features = 1, project = sample_name)
  # Mitrocondria
  seurat.object$percent.mito <- PercentageFeatureSet(object = seurat.object, features = mito_genes)
  
  return(seurat.object)
}

day4 = prepare_datasets(samples[1], cc.genes, mito_genes)
day7_lv = prepare_datasets(samples[2], cc.genes, mito_genes)
day7_rv = prepare_datasets(samples[3], cc.genes, mito_genes)
day10_lv = prepare_datasets(samples[4], cc.genes, mito_genes)
day10_rv = prepare_datasets(samples[5], cc.genes, mito_genes)
day14_lv = prepare_datasets(samples[6], cc.genes, mito_genes)
day14_rv = prepare_datasets(samples[7], cc.genes, mito_genes)
dim(day4)
dim(day7_lv)
dim(day7_rv)
dim(day10_lv)
dim(day10_rv)
dim(day14_lv)
dim(day14_rv)


# save.image("robjs/all.objs.RData")
load("robjs/all.objs.RData")


chicken = merge(day4, y = c(day7_lv, day7_rv, day10_lv, day10_rv, day14_lv, day14_rv), add.cell.ids = samples, project = "ChickenEmbryo")
# chicken = merge(day4, y = c(day7_lv, day7_rv), add.cell.ids = samples[1:3], project = "ChickenEmbryoEarly")
# chicken = merge(day10_lv, y = c(day10_rv, day14_lv, day14_rv), add.cell.ids = samples[4:7], project = "ChickenEmbryoLate")
# chicken = merge(day7_lv, y = c(day7_rv, day10_lv, day10_rv), add.cell.ids = samples[2:5], project = "ChickenEmbryoMiddle")
dim(chicken)
table(chicken$orig.ident)

VlnPlot(object = chicken, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), pt.size = 0.01, group.by = "orig.ident")
# Analysis
FeatureScatter(object = chicken, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", use.raw=T)
FeatureScatter(object = chicken, feature1 = "nFeature_RNA", feature2 = "percent.mito", do.return=TRUE) + geom_hline(yintercept = 20) + geom_vline(xintercept = 200)

# Filter out cells with few reads and few genes.
chicken <- subset(chicken, subset = nFeature_RNA >= 200 & percent.mito <= 30)
# Assign cell cycle score
chicken <- CellCycleScoring(object = chicken, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)

# Normalise and scale
chicken <- NormalizeData(object = chicken, scale.factor = 1e6)
chicken <- FindVariableFeatures(object = chicken, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
chicken <- ScaleData(object = chicken)
dim(GetAssayData(chicken, assay = "RNA", slot = "scale.data"))
chicken <- RunPCA(object = chicken, assay = "RNA")
FeaturePlot(chicken, reduction = "pca", c("nCount_RNA", "nFeature_RNA"))

library(sctransform)
chicken <- SCTransform(chicken, verbose = TRUE)
dim(GetAssayData(chicken, assay = "SCT", slot = "scale.data"))
DefaultAssay(chicken)
VlnPlot(object = chicken, features = c("nCount_RNA", "nFeature_RNA", "percent.mito", "nCount_SCT", "nFeature_SCT"), pt.size = 0.0, ncol = 3, group.by = "orig.ident")
chicken <- CellCycleScoring(object = chicken, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)

chicken <- RunPCA(object = chicken, assay = "SCT", reduction.name = "SCT_pca")
FeaturePlot(chicken, reduction = "SCT_pca", c("nCount_SCT", "nFeature_SCT"))

n.pcs = 20
# Clustering
chicken <- FindNeighbors(object = chicken, assay = "RNA", reduction = "pca", dims = 1:n.pcs, force.recalc = TRUE)
chicken <- FindClusters(object = chicken, resolution = 0.5)
table(Idents(chicken))

# chicken <- FindNeighbors(object = chicken, assay = "SCT", reduction = "SCT_pca", dims = 1:n.pcs, force.recalc = TRUE, graph.name = "SCT_graph")
chicken <- FindClusters(object = chicken, graph.name = "SCT_graph", resolution = 0.5)
table(Idents(chicken))

# To visualize
chicken <- RunUMAP(object = chicken, assay = "RNA", reduction = "pca", dims = 1:n.pcs)
DimPlot(chicken, reduction = "umap", label = TRUE, group.by = "RNA_snn_res.0.5")
DimPlot(chicken, reduction = "umap", group.by = "orig.ident")
  
chicken <- RunUMAP(object = chicken, assay = "SCT", reduction = "SCT_pca", dims = 1:n.pcs, reduction.name = "SCT_umap")
DimPlot(chicken, reduction = "SCT_umap", label = TRUE, group.by = "SCT_s")
DimPlot(chicken, reduction = "SCT_umap", group.by = "orig.ident")

# save(chicken, file="robjs/chicken_raw_SCT.Robj")
load(here("robjs", "chicken_raw.Robj"))

# Cluster17 - IFI6, GP9, BIN2, HPSE 
cluster17 = FindMarkers(chicken, ident.1 = 17, features = VariableFeatures(chicken), logfc.threshold = 3)

# Cluster15 - CSF1R, CCAH221, C1QC, LCP1, C1QB, GSTA3, CD83, LY86 (Macrophages)
cluster15 = FindMarkers(chicken, ident.1 = 15, features = VariableFeatures(chicken), logfc.threshold = 3)

chicken = SubsetData(chicken, ident.remove = c(8,11,15,16,17))
table(Idents(chicken))
# Rerun the pipeline starting from find variable genes

min(GetAssayData(chicken))

# save(chicken, file="robjs/chicken_raw_filtered.Robj")
load(here("robjs", "chicken_raw_filtered.Robj"))
FeaturePlot(chicken, c("TNNC1", "NKX2-5", "HBA1", "TMSB4X", "ACTB"))
