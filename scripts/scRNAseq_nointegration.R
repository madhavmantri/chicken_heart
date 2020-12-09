# Preprocessing
# Load the requisite packages and some additional helper functions.
library(Seurat)
packageVersion("Seurat");
library(Matrix); library(stringr); library(readr)
library(here); library(fitdistrplus)
library(dplyr); library(monocle); library(reticulate)
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
# load("robjs/all.objs.RData")

chicken = merge(day4, y = c(day7_lv, day7_rv, day10_lv, day10_rv, day14_lv, day14_rv), add.cell.ids = samples, project = "ChickenEmbryo")
dim(chicken)
table(chicken$orig.ident)

# Analysis starts here
VlnPlot(object = chicken, features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), pt.size = 0.01, group.by = "orig.ident")
FeatureScatter(object = chicken, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", use.raw=T)
FeatureScatter(object = chicken, feature1 = "nFeature_RNA", feature2 = "percent.mito", do.return=TRUE) + geom_hline(yintercept = 20) + geom_vline(xintercept = 200)

# Filter out cells with few reads and few genes.
chicken <- subset(chicken, subset = nFeature_RNA >= 200 & percent.mito <= 20)
dim(chicken)

# Assign cell cycle score
chicken <- CellCycleScoring(object = chicken, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)

# Normalise, scale, and run PCA
chicken <- NormalizeData(object = chicken, scale.factor = 1e6)
chicken <- FindVariableFeatures(object = chicken, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
chicken <- ScaleData(object = chicken)
dim(GetAssayData(chicken, assay = "RNA", slot = "scale.data"))
chicken <- RunPCA(object = chicken, assay = "RNA")
FeaturePlot(chicken, reduction = "pca", c("nCount_RNA", "nFeature_RNA"))
ElbowPlot(chicken, reduction = "pca")
n.pcs = 20

# CLustering and UMAP dimesntion reduction 
chicken <- FindNeighbors(object = chicken, assay = "RNA", reduction = "pca", dims = 1:n.pcs, force.recalc = TRUE)
chicken <- FindClusters(object = chicken, resolution = 0.5)
table(Idents(chicken))

# To visualize
chicken <- RunUMAP(object = chicken, assay = "RNA", reduction = "pca", dims = 1:n.pcs)
DimPlot(chicken, reduction = "umap", label = TRUE, group.by = "RNA_snn_res.0.5")
DimPlot(chicken, reduction = "umap", group.by = "orig.ident")

# save(chicken, file="robjs/chicken_raw.Robj")
# load(here("robjs", "chicken_raw.Robj"))


################### This part includes SC Tranform based normalisation (optional) #######################

# Normalise using SC Transform method in Seurat v-3 
library(sctransform)
chicken <- SCTransform(chicken, verbose = TRUE)

# SCT creates a new SCT assay and run PCA on SCT values
DefaultAssay(chicken)
chicken <- CellCycleScoring(object = chicken, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes)
chicken <- RunPCA(object = chicken, assay = "SCT", reduction.name = "SCT_pca")
FeaturePlot(chicken, reduction = "SCT_pca", c("nCount_SCT", "nFeature_SCT"))
ElbowPlot(chicken, reduction = "SCT_pca")
n.pcs = 20

# CLustering and UMAP dimesntion reduction using SCT values
chicken <- FindNeighbors(object = chicken, assay = "SCT", reduction = "SCT_pca", dims = 1:n.pcs, force.recalc = TRUE, graph.name = "SCT_graph")
chicken <- FindClusters(object = chicken, graph.name = "SCT_graph", resolution = 0.5)
table(Idents(chicken))

# To visualize
chicken <- RunUMAP(object = chicken, assay = "SCT", reduction = "SCT_pca", dims = 1:n.pcs, reduction.name = "SCT_umap")
DimPlot(chicken, reduction = "SCT_umap", label = TRUE, group.by = "SCT_graph_res.0.5")
DimPlot(chicken, reduction = "SCT_umap", group.by = "orig.ident")

# save(chicken, file="robjs/chicken_raw_SCT.Robj")
# load(here("robjs", "chicken_raw_SCT.Robj"))
