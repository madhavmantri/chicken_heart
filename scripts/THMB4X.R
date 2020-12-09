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
                theme_bw()
  
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
subcluster <- RenameIdents(subcluster, "0" = "Vascular smooth muscle- like", "1" = "Vascular endothelial- like", "2" = "Cardiomyocytes_like")
subcluster$seurat_clusters <- Idents(subcluster)

pdf(file="tmsb4x_exp_subclusters_day.pdf",
    width=1.0, height=0.8, paper="special", bg="white",
    fonts="Helvetica", colormodel = "rgb", pointsize=5)
library(pals)
library(ggplot2)
subcluster$day <- factor(subcluster$day, levels = c("D4", "D7", "D10", "D14"))
VlnPlot(object = subcluster, features = c("TMSB4X"), group.by = "day", pt.size = 0) + 
  stat_summary(fun.y = median, geom='point', size = 5, colour = "black", shape = 95) + 
  theme_bw() + labs(y = "Log-normalised TMSB4X expression") + 
  theme(plot.background=element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        legend.key.size = unit(1, "pt"),
        axis.title.x = element_blank() ,
        axis.title.y=element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 6, family = "Helvetica"),
        axis.ticks.length.x = unit(0, "pt"),
        plot.title=element_blank())
dev.off()

pdf(file="TMSB4X_Clusters.pdf",
    width=1.5, height=1.5, paper="special", bg="transparent",
    fonts="Helvetica", colormodel = "rgb", pointsize=5)
DimPlot(subcluster, reduction = "umap") + theme_nothing() # + scale_color_manual(values = as.vector(alphabet())[1:20])
dev.off()

pdf(file="TMSB4X_Day.pdf",
    width=1.5, height=1.5, paper="special", bg="transparent",
    fonts="Helvetica", colormodel = "rgb", pointsize=5)
DimPlot(subcluster, reduction = "umap", group.by = "day") + theme_nothing()
dev.off()

pdf(file="TMSB4X_PECAM1.pdf",
    width=1.5, height=1.5, paper="special", bg="transparent",
    fonts="Helvetica", colormodel = "rgb", pointsize=5)
FeaturePlot(subcluster, reduction = "umap", features = c("EGFL7"), pt.size = 0.001) + theme_nothing()
dev.off()

library(pals)
library(ggplot2)
temp <- SubsetData(chicken.integrated, ident.remove = c("Macrophages", "Dendritic cells", "Erythrocytes"))
temp1 <- cbind("TMSB4X" = as.matrix(GetAssayData(temp))["TMSB4X",], temp@meta.data)
cdata <- ddply(temp1, c("day", "celltypes.0.5"), summarise, 
               N    = length(TMSB4X),
               mean = mean(TMSB4X),
               sd   = sd(TMSB4X), 
               se = sd / sqrt(N)
)
pdf(file="TMSB4X_heatmap.pdf",
    width=4, height=1.5, paper="special", bg="white",
    fonts="Helvetica", colormodel = "rgb", pointsize=5)
ggplot(data = cdata, mapping = aes(x = celltypes.0.5, y = day, fill = mean)) + geom_bin2d() + theme_classic() + scale_fill_viridis(direction = -1) + 
  labs(x = "Cell types", y = "Stage", fill = "TMSB4X") + 
  theme(plot.background=element_blank(),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "right",
        legend.text = element_text(colour = "black", size = 6, family = "Helvetica"),
        axis.title.x =  element_text(colour = "black", size = 8, family = "Helvetica"),
        axis.title.y =  element_text(colour = "black", size = 8, family = "Helvetica"),
        legend.title =  element_text(colour = "black", size = 8, family = "Helvetica"),
        axis.text.x = element_text(colour = "black", size = 6, family = "Helvetica", angle = 30, vjust = 1, hjust = 1), # 
        axis.text.y = element_text(colour = "black", size = 6, family = "Helvetica"),
        plot.title=element_blank())
dev.off()


