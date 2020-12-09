library(phateR)
library(Seurat, lib.loc = "/home/mm2937/miniconda3/pkgs/r-seurat-3.1.1-r35h0357c0b_0/lib/R/library")
packageVersion("Seurat")
library(Matrix); library(stringr)
library(readr); library(here)
library(fitdistrplus); library(dplyr)
library(monocle); library(reticulate)
setwd("/workdir/mm2937/chicken/")

samples = c("D4", "D7-LV", "D7-RV", "D10-LV", "D10-RV", "D14-LV", "D14-RV")

seurat3tomonocle <- function(otherCDS, assay, slot, lowerDetectionLimit = 0, import_all = FALSE) {
  if(class(otherCDS)[1] == 'Seurat') {
    requireNamespace("Seurat")
    data <- GetAssayData(otherCDS, assay = assay, slot = slot)
    data <- data[rowSums(as.matrix(data)) != 0,]
    pd <- new("AnnotatedDataFrame", data = otherCDS@meta.data)
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new("AnnotatedDataFrame", data = fData)
    if(all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    } else if(any(data < 0)){
      expressionFamily <- uninormal()
    } else {
      expressionFamily <- tobit()
    }
    valid_data <- data[, row.names(pd)]
    monocle_cds <- newCellDataSet(data,
                                  phenoData = pd, 
                                  featureData = fd,
                                  lowerDetectionLimit=lowerDetectionLimit,
                                  expressionFamily=expressionFamily)
    if(import_all) {
      if("Monocle" %in% names(otherCDS@misc)) {
        otherCDS@misc$Monocle@auxClusteringData$seurat <- NULL
        otherCDS@misc$Monocle@auxClusteringData$scran <- NULL
        monocle_cds <- otherCDS@misc$Monocle
        mist_list <- otherCDS
      } else {
        mist_list <- otherCDS
      }
    } else {
      mist_list <- list()
    } 
    # monocle_cds@auxClusteringData$seurat <- mist_list
  } 
  return(monocle_cds)
}


# Load full integrated dataset for lineage analysis
load("robjs/chicken_normalised_scanorama2.Robj")
DimPlot(object = chicken.integrated, reduction = "umap") 

# Subset the myocardial population out and preprocess the data
myocardial <- subset(chicken.integrated, subset = seurat_clusters %in% c(3))
DefaultAssay(myocardial) <- "RNA"
myocardial <- FindVariableFeatures(myocardial) %>% ScaleData()

# Perform PCA on myocardial cells
myocardial <- RunPCA(object = myocardial)
ElbowPlot(object = myocardial, ndims = 50)

# Perform clustering and UMAP reduction
myocardial <- FindNeighbors(object=myocardial, dims = 1:20, force.recalc = TRUE)
myocardial <- FindClusters(object=myocardial, resolution=0.30)
table(Idents(myocardial))
myocardial <- RunUMAP(object = myocardial, dims = 1:20)

# To visualise
DimPlot(myocardial, reduction = "umap", label = TRUE)
DimPlot(object = myocardial, reduction = "umap", group.by = "orig.ident") 
DimPlot(object = myocardial, reduction = "umap", group.by = "Phase")

# Integrate myocardial lineage again using scanorama
extractRNA_chicken <- function(seurat.object, assay = "RNA", slot = "data", sample_name = NULL){
  return(t(as.matrix(GetAssayData(seurat.object, assay = assay, slot = slot)))[colnames(seurat.object)[seurat.object$orig.ident == sample_name],])
}

chicken <- myocardial
samples <- unique(chicken$orig.ident)
data = list(extractRNA_chicken(chicken, sample_name = samples[1]), extractRNA_chicken(chicken, sample_name = samples[2]), 
            extractRNA_chicken(chicken, sample_name = samples[3]), extractRNA_chicken(chicken, sample_name = samples[4]),
            extractRNA_chicken(chicken, sample_name = samples[5]), extractRNA_chicken(chicken, sample_name = samples[6]), 
            extractRNA_chicken(chicken, sample_name = samples[7]))
gene_list = list(rownames(chicken), rownames(chicken), rownames(chicken), rownames(chicken), rownames(chicken), rownames(chicken), rownames(chicken))

# Scanorama integartion 
scanorama = import('scanorama')
integrated.corrected.data = scanorama$correct(data, gene_list, return_dimred=TRUE, return_dense=TRUE, ds_names = samples, verbose = TRUE)

corrected_scanorama <- t(do.call(rbind, integrated.corrected.data[[2]]))
colnames(corrected_scanorama) <- colnames(chicken)
rownames(corrected_scanorama) <- integrated.corrected.data[[3]]
dim(corrected_scanorama)

# Add scanorama values to original seurat object and set as default assay for subcluster analyis 
scanorama_assay <- CreateAssayObject(data = corrected_scanorama)
myocardial[["scanorama"]] <- scanorama_assay
DefaultAssay(myocardial) <- "scanorama"

# Find variable features and process data again
myocardial <- FindVariableFeatures(myocardial) %>% ScaleData() 

# Run PCA on scanorama values
myocardial <- RunPCA(object = myocardial, reduction.name = "pca_scanorama_data")
ElbowPlot(object = myocardial, ndims = 50)

# Perform clustering and UMAP reduction
myocardial <- FindNeighbors(object=myocardial, reduction = "pca_scanorama_data", dims = 1:20, force.recalc = TRUE, graph.name = "scanorama_snn_data")
myocardial <- FindClusters(object=myocardial, graph.name = "scanorama_snn_data", resolution=0.20)
table(Idents(myocardial))
myocardial <- RunUMAP(object = myocardial, reduction = "pca_scanorama_data", dims = 1:20, reduction.name = "umap_scanorama_data")

# To visualise
DimPlot(myocardial, reduction = "umap_scanorama_data", label = TRUE)
DimPlot(object = myocardial, reduction = "umap_scanorama_data", group.by = "orig.ident") 
DimPlot(object = myocardial, reduction = "umap_scanorama_data", group.by = "Phase")

# save(myocardial, file="robjs/myocardial.Robj")
# load("robjs/myocardial.Robj")

# Run 
data <- t(GetAssayData(myocardial, assay = "scanorama_data", slot = "data"))
tree_chicken <- phate(data, gamma = 0)
myocardial@reductions[["scanorama_phate_gamma0"]] <- CreateDimReducObject(embeddings = tree_chicken$embedding, assay = "scanorama", key = "PHATE")
DimPlot(myocardial, reduction = "scanorama_phate_gamma0", group.by = "orig.ident") 
DimPlot(myocardial, reduction = "scanorama_phate_gamma0") 

tree_chicken <- phate(data, gamma = 1)
myocardial@reductions[["scanorama_phate_gamma1"]] <- CreateDimReducObject(embeddings = tree_chicken$embedding, assay = "scanorama", key = "PHATE")
DimPlot(myocardial, reduction = "scanorama_phate_gamma1", group.by = "orig.ident") 

# save(myocardial, file="robjs/myocardial.Robj")
load("robjs/myocardial.Robj")

FeaturePlot(myocardial, reduction = "umap_scanorama_data", c("TNNC2"))
VlnPlot(myocardial, c("nFeature_RNA", "nCount_RNA"))

DefaultAssay(myocardial) <- "RNA"
table(myocardial$scanorama_snn_data_res.0.1)
Idents(myocardial) <- myocardial$scanorama_snn_data_res.0.1
myocardial <- RenameIdents(myocardial, "0" = "C1",  "2" = "C2", "1" = "C3")
DefaultAssay(myocardial) <- "RNA"
markers.myocardial.RNA <- FindAllMarkers(myocardial, assay = "RNA", logfc.threshold = 0.5, return.thresh = 0.01, only.pos = T, do.print = TRUE)
markers.myocardial.RNA <- subset(markers.myocardial.RNA[!(rownames(markers.myocardial.RNA) %in% grep("^ENSGAL", x = rownames(markers.myocardial.RNA), value = TRUE)),])
write_csv(markers.myocardial.RNA, path = "myocardial.markers.csv")
markers.top10 = markers.myocardial.RNA %>% group_by(cluster) %>% top_n(10, avg_logFC)

myocardial$subcluster <- factor(Idents(myocardial), levels = c("C1", "C2", "C3"))

pdf(file="MyoClusters.pdf",
    width=1.5, height=1.5, paper="special", bg="transparent",
    fonts="Helvetica", colormodel = "rgb", pointsize=5)
DimPlot(myocardial, reduction = "umap_scanorama_data", group.by = "subcluster") + scale_color_manual(values = as.vector(alphabet())[1:20]) + theme_nothing()
DimPlot(myocardial, reduction = "scanorama_phate_gamma0", group.by = "subcluster") + scale_color_manual(values = as.vector(alphabet())[1:20]) + theme_nothing()
dev.off()

pdf(file="MyoDay.pdf",
    width=1.5, height=1.5, paper="special", bg="transparent",
    fonts="Helvetica", colormodel = "rgb", pointsize=5)
DimPlot(myocardial, reduction = "umap_scanorama_data", group.by = "day") + theme_nothing()
DimPlot(myocardial, reduction = "scanorama_phate_gamma0", group.by = "day") + theme_nothing()
dev.off()

pdf(file="myo_EGFL7.pdf",
    width=1.5, height=1.5, paper="special", bg="transparent",
    fonts="Helvetica", colormodel = "rgb", pointsize=5)
FeaturePlot(myocardial, reduction = "umap_scanorama_data", features = c("EGFL7"), pt.size = 0.001) + theme_nothing()
dev.off()

pdf(file="Myo_nCount.pdf",
    width=1.3, height=1.3, paper="special", bg="transparent",
    fonts="Helvetica", colormodel = "rgb", pointsize=5)
VlnPlot(object = myocardial, features = c("nCount_RNA"), group.by = "subcluster", pt.size = 0) + scale_y_log10() + scale_fill_manual(values = as.vector(alphabet())[1:20]) +
  stat_summary(fun.y = median, geom='point', size = 5, colour = "black", shape = 95) +
  xlab("Sample")+
  ylab("Number of genes")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background= element_rect(fill = "transparent", color = NA),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y= element_blank(),
        axis.text = element_text(family = "Helvetica", size = 6, color = "black"),
        axis.text.x = element_text(colour = "black", size = 6, family = "Helvetica", angle = 0, color = as.vector(alphabet())[1:20]), 
        plot.title=element_blank())
dev.off()

# save(myocardial, file="robjs/myocardial.Robj")
load("robjs/myocardial.Robj")

## monocle analysis on reintegarted datasets
cds <- seurat3tomonocle(myocardial, assay = "RNA", slot = "counts")
cds <- estimateSizeFactors(cds)
cds <- detectGenes(cds)
DefaultAssay(myocardial)

seurat_variable_genes <- VariableFeatures(FindVariableFeatures(myocardial, assay = "RNA"), assay = "RNA")

DimPlot(myocardial, reduction = "umap_scanorama_data", group.by = "seurat_clusters")
day_genes <- differentialGeneTest(cds, fullModelFormulaStr = '~seurat_clusters', cores = 16)
# save(day_genes, file="robjs/myo_day_genes.Robj")
load("robjs/myo_day_genes.Robj")

day_genes1 <- row.names(day_genes)[order(day_genes$qval)][1:2000]

cds <- setOrderingFilter(cds, ordering_genes = day_genes1)
dim(cds)

# cds <- reduceDimension(cds, method = 'DDRTree', max_components = 2, norm_method = 'log', verbose = TRUE)
cds <- reduceDimension(cds,  method = 'DDRTree', max_components = 2, norm_method = "none", pseudo_expr = 0, verbose = TRUE)
cds <- orderCells(cds)

GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$orig.ident)[,"D4"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

cds <- orderCells(cds, root_state = GM_state(cds))


plot_cell_trajectory(cds)
plot_cell_trajectory(cds, color_by = "day")
plot_cell_trajectory(cds, color_by = "seurat_clusters") + facet_wrap(~seurat_clusters)
plot_cell_trajectory(cds, color_by = "Phase")
plot_cell_trajectory(cds, color_by = "Pseudotime")

# save(cds, file="robjs/myo_cds.Robj")
load("robjs/myo_cds.Robj")

pdf(file="MyoCDSClusters.pdf",
    width=1.5, height=1.5, paper="special", bg="transparent",
    fonts="Helvetica", colormodel = "rgb", pointsize=5)
plot_cell_trajectory(cds, color_by = "seurat_clusters", show_branch_points = FALSE, cell_size = 0.5) + theme_nothing() + scale_color_manual(values = as.vector(alphabet())[1:20]) 
dev.off()

pdf(file="MyoCDSDay.pdf",
    width=1.5, height=1.5, paper="special", bg="transparent",
    fonts="Helvetica", colormodel = "rgb", pointsize=5)
plot_cell_trajectory(cds, color_by = "day", show_branch_points = FALSE, cell_size = 0.5) + theme_nothing() 
dev.off()
