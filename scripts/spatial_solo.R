library(Seurat)
packageVersion("Seurat")
library(Matrix); library(stringr)
library(readr); library(here)
library(fitdistrplus); library(dplyr)
library(SeuratWrappers)
library(ggplot2)
library(cowplot)
library(reticulate)
library(pals)
# library(monocle)
setwd("/workdir/mm2937/chicken/")

load(file = "cc.genes.rda") 
mito_genes = c("ND1", "ND2", "COX1", "COII", "ATP8", "ATP6", "COX3", "ND3", "ND4L", "ND4", "ND5", "CYTB", "ND6")
samples = c("D4-A1-H1", "D4-A1-H2", "D4-A1-H3", "D4-A1-H4", "D4-A1-H5", "D7-B1-H1", "D7-B1-H2", "D7-B1-H3", "D7-B1-H4", "D10-C1-H1", "D10-C1-H2", "D14-D1")
samples = c("D4-A1", "D7-B1", "D10-C1", "D14-D1")
cc.genes <- NULL

# This script prepares and analyses visium datasets for individual stages and saves information from integration of scRNAseq and spatial datasets

load("robjs/all.visiums.4.RData")
load("robjs/chicken_visium.4.prediction.1.Robj")
load("robjs/chicken_normalised_scanorama3.Robj")

# # Testing clustering at high resolution on combined datasets to find analtomical regions
# DefaultAssay(chicken_visium) <- "scanorama"
# chicken_visium <- FindNeighbors(object = chicken_visium, dims=1:20, force.recalc = TRUE)
# chicken_visium <- FindClusters(object = chicken_visium, resolution=1.0)
# SpatialDimPlot(chicken_visium, group.by = "seurat_clusters", pt.size.factor = 1.0, crop = F, ncol = 2) + coord_cartesian()

##########################  This section preprocesses individual visium datasets and reduce to PCA dimensions  ##################################################
##########################  This sections also calculates spatially variable genes for visium datasets ##########################################################

day14_visium <- NormalizeData(day14_visium) %>% FindVariableFeatures() %>% ScaleData()  %>% RunPCA()
day14_visium <- FindSpatiallyVariableFeatures(day14_visium, assay = "Spatial", features = VariableFeatures(day14_visium)[1:2000], 
                                              selection.method = "markvariogram", verbose = T)
day10_visium <- NormalizeData(day10_visium) %>% FindVariableFeatures() %>% ScaleData()  %>% RunPCA()
day10_visium <- FindSpatiallyVariableFeatures(day10_visium, assay = "Spatial", features = VariableFeatures(day10_visium)[1:2000], 
                                              selection.method = "markvariogram", verbose = T)
day7_visium <- NormalizeData(day7_visium) %>% FindVariableFeatures() %>% ScaleData()  %>% RunPCA()
day7_visium <- FindSpatiallyVariableFeatures(day7_visium, assay = "Spatial", features = VariableFeatures(day7_visium)[1:2000], 
                                              selection.method = "markvariogram", verbose = T)
day4_visium <- NormalizeData(day4_visium) %>% FindVariableFeatures() %>% ScaleData()  %>% RunPCA()
day4_visium <- FindSpatiallyVariableFeatures(day4_visium, assay = "Spatial", features = VariableFeatures(day4_visium)[1:2000], 
                                             selection.method = "markvariogram", verbose = T)

# save(day4_visium, day7_visium, day10_visium, day14_visium, file = "all.visiums.4.solo.Robj")
load("robjs/all.visiums.4.solo.Robj")


################################### This section performs clustering on individual datasets and label anatomical regions ##############################################

day14_visium <- FindNeighbors(object = day14_visium, features = SpatiallyVariableFeatures(day4_visium), dims=1:20, force.recalc = TRUE)
day14_visium <- FindClusters(object = day14_visium, resolution=0.6)
SpatialDimPlot(day14_visium, crop = F, pt.size.factor = 1.0) + coord_cartesian()
# SpatialDimPlot(day14_visium, crop = F, pt.size.factor = 1.0) + coord_cartesian() + geom_vline(xintercept = 310) + geom_vline(xintercept = 400)
day14_visium <- RenameIdents(day14_visium, "0" = "Compact left ventricle and \ninter-ventricular septum", "2" = "Right ventricle", 
                             "4" = "Endothelium", "1" = "Atria", "5" = "Epicardium", "3" = "Valves")
day14_visium$region <- Idents(day14_visium)

day10_visium <- FindNeighbors(object = day10_visium, dims=1:20, force.recalc = TRUE)
day10_visium <- FindClusters(object = day10_visium, resolution=0.5)
SpatialDimPlot(day10_visium, crop = F, pt.size.factor = 1.0) + coord_cartesian()
SpatialDimPlot(day10_visium, crop = F, pt.size.factor = 1.0) + coord_cartesian() + geom_vline(xintercept = 200) + geom_vline(xintercept = 260) +
  geom_vline(xintercept = 325) + geom_vline(xintercept = 400) +  geom_vline(xintercept = 470)
day10_visium <- RenameIdents(day10_visium, "0" = "Compact left ventricle and \ninter-ventricular septum", "3" = "Right ventricle", 
                             "1" = "Atria", "5" = "Epicardium", "4" = "Valves",
                             "2" = 'Trabecular left ventricle and \nendocardium', "6" = "Outflow tract")
day10_visium$region <- Idents(day10_visium)

day7_visium <- FindNeighbors(object = day7_visium, dims=1:20, force.recalc = TRUE)
day7_visium <- FindClusters(object = day7_visium, resolution=0.6)
SpatialDimPlot(day7_visium, crop = F, pt.size.factor = 1.0) + coord_cartesian()
SpatialDimPlot(day7_visium, crop = F, pt.size.factor = 1.0) + coord_cartesian() + geom_vline(xintercept = 200) + geom_vline(xintercept = 260) +
  geom_vline(xintercept = 325) + geom_vline(xintercept = 400) +  geom_vline(xintercept = 470)
day7_visium <- RenameIdents(day7_visium, "0" = "Compact left ventricle and \ninter-ventricular septum", "1" = "Trabecular left ventricle and \nendocardium", "5" = "Right ventricle", "3" = "Endothelium", "2" = "Atria", "4" = "Epicardium", "6" = "Valves")
day7_visium$region <- Idents(day7_visium)

day4_visium <- FindNeighbors(object = day4_visium, dims=1:20, force.recalc = TRUE)
day4_visium <- FindClusters(object = day4_visium, resolution=0.5)
SpatialDimPlot(day4_visium, crop = F, pt.size.factor = 1.0) + coord_cartesian()
# SpatialDimPlot(day4_visium, crop = F, pt.size.factor = 1.0) + coord_cartesian() + geom_vline(xintercept = 200) + geom_vline(xintercept = 260) + 
#   geom_vline(xintercept = 325) + geom_vline(xintercept = 400) +  geom_vline(xintercept = 470) 
day4_visium <- RenameIdents(day4_visium, "0" = "Ventricle", "1" = "Epicardium- like", 
                            "2" = "Outflow tract", "4" = "Atria", "3" = "Valves")
day4_visium$region <- Idents(day4_visium)

# save(day4_visium, day7_visium, day10_visium, day14_visium, file = "all.visiums.4.solo1.Robj")
load("robjs/all.visiums.4.solo1.Robj")

chicken_visium$region <- NA
chicken_visium$region[paste("D4-A1", colnames(day4_visium), sep = "_")] <- as.character(day4_visium$region)
chicken_visium$region[paste("D7-B1", colnames(day7_visium), sep = "_")] <- as.character(day7_visium$region)
chicken_visium$region[paste("D10-C1", colnames(day10_visium), sep = "_")] <- as.character(day10_visium$region)
chicken_visium$region[paste("D14-D1", colnames(day14_visium), sep = "_")] <- as.character(day14_visium$region)

# save(chicken_visium, file = "robjs/chicken_visium.4.prediction.1.Robj")
load("robjs/chicken_visium.4.prediction.1.Robj")

######################################################### This region find differnetial genes expressed in anatomical regions ###########################################


# DefaultAssay(chicken_visium) <- "Spatial"
# chicken_visium_markers <- FindAllMarkers(chicken_visium, only.pos = T)
# chicken_visium_markers <- subset(chicken_visium_markers[!(rownames(chicken_visium_markers) %in% grep("^ENSGAL", x = rownames(chicken_visium_markers), value = TRUE)),])
# markers.top20 = chicken_visium_markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
# markers.top10 = chicken_visium_markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
# DoHeatmap(chicken_visium_markers, markers.top10$gene)


SpatialDimPlot(day4_visium, crop = F, pt.size.factor = 1.0) + coord_cartesian()
day4_spatial_markers <- FindAllMarkers(day4_visium, assay = "Spatial", only.pos = T)
day4_spatial_markers <- subset(day4_spatial_markers[!(rownames(day4_spatial_markers) %in% grep("^ENSGAL", x = rownames(day4_spatial_markers), value = TRUE)),])
markers.top20 = day4_spatial_markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
markers.top10 = day4_spatial_markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
markers.top5 = day4_spatial_markers %>% group_by(cluster) %>% top_n(5, avg_logFC)

SpatialDimPlot(day7_visium, crop = F, pt.size.factor = 1.0) + coord_cartesian()
day7_spatial_markers <- FindAllMarkers(day7_visium, assay = "Spatial",  only.pos = T, min.pct = 0.3, logfc.threshold = 0.5)
day7_spatial_markers <- subset(day7_spatial_markers[!(rownames(day7_spatial_markers) %in% grep("^ENSGAL", x = rownames(day7_spatial_markers), value = TRUE)),])
markers.top20 = day7_spatial_markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
markers.top10 = day7_spatial_markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
markers.top5 = day7_spatial_markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
DoHeatmap(day7_visium, markers.top10$gene)

SpatialDimPlot(day10_visium, crop = F, pt.size.factor = 1.0) + coord_cartesian()
day10_spatial_markers <- FindAllMarkers(day10_visium, assay = "Spatial", only.pos = T)
day10_spatial_markers <- subset(day10_spatial_markers[!(rownames(day10_spatial_markers) %in% grep("^ENSGAL", x = rownames(day10_spatial_markers), value = TRUE)),])
markers.top20 = day10_spatial_markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
markers.top10 = day10_spatial_markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
markers.top5 = day10_spatial_markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
DoHeatmap(day10_visium, markers.top10$gene)

SpatialDimPlot(day14_visium, crop = F, pt.size.factor = 1.0) + coord_cartesian()
day14_spatial_markers <- FindAllMarkers(day14_visium, assay = "Spatial", only.pos = T)
day14_spatial_markers <- subset(day14_spatial_markers[!(rownames(day14_spatial_markers) %in% grep("^ENSGAL", x = rownames(day14_spatial_markers), value = TRUE)),])
markers.top20 = day14_spatial_markers %>% group_by(cluster) %>% top_n(20, avg_logFC)
markers.top10 = day14_spatial_markers %>% group_by(cluster) %>% top_n(10, avg_logFC)
markers.top5 = day14_spatial_markers %>% group_by(cluster) %>% top_n(5, avg_logFC)
DoHeatmap(day10_visium, markers.top10$gene)

DotPlot(day14_visium, features = unique(markers.top5$gene), cols = c("gray", "brown"), scale.by = "size", dot.scale = 2.7)


####################################### This section analyses and visualises the top spatially variable/ restricted genes ###############################

top.features.day14 <- head(SpatiallyVariableFeatures(day14_visium, selection.method = "markvariogram"), 100)
top.features.day10 <- head(SpatiallyVariableFeatures(day10_visium, selection.method = "markvariogram"), 100)
top.features.day7 <- head(SpatiallyVariableFeatures(day7_visium, selection.method = "markvariogram"), 100)
top.features.day4 <- head(SpatiallyVariableFeatures(day4_visium, selection.method = "markvariogram"), 100)
all.genes = unique(c(rownames(day4_visium), rownames(day7_visium), rownames(day10_visium), rownames(day14_visium)))
# save(all.genes, top.features.day4, top.features.day7, top.features.day10, top.features.day10, top.features.day14, file = "robjs/top.features.100.Robj")
load("robjs/top.features.100.Robj")

SpatialFeaturePlot(day14_visium, features = c("PITX2"), crop = F, pt.size.factor = 1.0)+ coord_cartesian()
SpatialFeaturePlot(day10_visium, features = c("SPARC"), crop = F, pt.size.factor = 1.0)+ coord_cartesian()
SpatialFeaturePlot(day7_visium, features = c("BMP10"), crop = F, pt.size.factor = 1.0) + coord_cartesian()
SpatialFeaturePlot(day4_visium, features = c("CSRP2"), crop = F, pt.size.factor = 1.0)+ coord_cartesian()

common_genes <- intersect(top.features.day4, intersect(top.features.day7, intersect(top.features.day10, top.features.day14)))


################################### This section adds the cell type composition to individual visium objects ##############################################
################################### This section also checks for spatially variable restricted cell types ###############################################

temp <- GetAssayData(subset(chicken_visium, subset = orig.ident == "D14"), assay = "predictions")
temp <- temp[rownames(temp) != "max",]
rowSums(temp)
rownames(temp)
colnames(temp) <- str_split_fixed(colnames(temp), "_", 2)[,2]
day14_visium[["predictions"]] <- CreateAssayObject(data = temp)
sum(colnames(day14_visium) != colnames(temp))
DefaultAssay(day14_visium) <- "predictions"
day14_visium <- FindSpatiallyVariableFeatures(day14_visium, assay = "predictions", selection.method = "markvariogram", 
                                             features = rownames(day14_visium)[rownames(day14_visium) != "max"], r.metric = 5, slot = "data")

temp <- GetAssayData(subset(chicken_visium, subset = orig.ident == "D10"), assay = "predictions")
temp <- temp[rownames(temp) != "max",]
rowSums(temp)
rownames(temp)
colnames(temp) <- str_split_fixed(colnames(temp), "_", 2)[,2]
day10_visium[["predictions"]] <- CreateAssayObject(data = temp)
sum(colnames(day10_visium) != colnames(temp))
DefaultAssay(day10_visium) <- "predictions"
day10_visium <- FindSpatiallyVariableFeatures(day10_visium, assay = "predictions", selection.method = "markvariogram", 
                                              features = rownames(day10_visium)[rownames(day10_visium) != "max"], r.metric = 5, slot = "data")

temp <- GetAssayData(subset(chicken_visium, subset = orig.ident == "D7"), assay = "predictions")
temp <- temp[rownames(temp) != "max",]
rowSums(temp)
rownames(temp)
colnames(temp) <- str_split_fixed(colnames(temp), "_", 2)[,2]
day7_visium[["predictions"]] <- CreateAssayObject(data = temp)
sum(colnames(day7_visium) != colnames(temp))
DefaultAssay(day7_visium) <- "predictions"
day7_visium <- FindSpatiallyVariableFeatures(day7_visium, assay = "predictions", selection.method = "markvariogram", 
                                              features = rownames(day7_visium)[rownames(day7_visium) != "max"], r.metric = 5, slot = "data")

temp <- GetAssayData(subset(chicken_visium, subset = orig.ident == "D4"), assay = "predictions")
temp <- temp[rownames(temp) != "max",]
rowSums(temp)
rownames(temp)
colnames(temp) <- str_split_fixed(colnames(temp), "_", 2)[,2]
day4_visium[["predictions"]] <- CreateAssayObject(data = temp)
sum(colnames(day4_visium) != colnames(temp))
DefaultAssay(day4_visium) <- "predictions"
day4_visium <- FindSpatiallyVariableFeatures(day4_visium, assay = "predictions", selection.method = "markvariogram", 
                                              features = rownames(day4_visium)[rownames(day4_visium) != "max"], r.metric = 5, slot = "data")

# save(day4_visium, day7_visium, day10_visium, day14_visium, file = "all.visiums.4.solo2.Robj")
load("robjs/all.visiums.4.solo2.Robj")

################################### Calculate entropy in cell type prediction values ############################################################
data <- as.data.frame(GetAssayData(chicken_visium, assay = "predictions"))
data <- data[rownames(data) != "max",]
library(entropy)
temp <- lapply(data, FUN = entropy.empirical)
chicken_visium$num_celltypes_0.05 = colSums(data >= 0.05)
chicken_visium$emp_entropy = unlist(temp)


######################### This (optional) section creates a dimension reduciton for visium spatial maps (removes the images) ########################################

visium.embeddings <- as.matrix(day10_visium@images$slice1@coordinates[,c("col", "row")])
colnames(visium.embeddings) <- c("visium_1", "visium_2")
day10_visium[['visium']] <- CreateDimReducObject(
  embeddings = visium.embeddings,
  assay = "Spatial",
  key = "visium_"
)
SpatialDimPlot(day10_visium)
DimPlot(day10_visium, reduction = "visium")
