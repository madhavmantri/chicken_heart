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

# load(here("robjs", "chicken_normalised_scanorama3.Robj"))
# load("robjs/chicken_visium.4.Robj")

DefaultAssay(chicken.integrated) <- "RNA"
DefaultAssay(chicken_visium) <- "Spatial"

# Find gene anchors between scRNAseq and spatila RNAseq datasets
anchors <- FindTransferAnchors(reference = chicken.integrated, query = chicken_visium, reduction = "cca")

# Transfer labels from scRNAseq to spatial using the anchor based approach
chicken.integrated$cellname <- colnames(chicken.integrated)
predictions.assay <- TransferData(anchorset = anchors, refdata = chicken.integrated$celltypes.0.5, prediction.assay = TRUE, 
                                  weight.reduction = chicken_visium[["pca"]])
# save(anchors, predictions.assay, file = "robjs/anchors_and_prediction_assay.RData")
# load("robjs/anchors_and_prediction_assay.RData")

# Adding cell type predictions to original seurat object
chicken_visium[["predictions"]] <- predictions.assay
dim(GetAssayData(chicken_visium, assay = "predictions"))

# Adding cell type predictions in meta data as well
chicken_visium <- AddMetaData(chicken_visium, metadata = as.data.frame(t(GetAssayData(chicken_visium, assay = "predictions"))))
head(chicken_visium@meta.data)

# Define cell type with maximum prediction score as spot type 
prediction.scores <- as.data.frame(t(GetAssayData(chicken_visium, assay = "predictions")))
prediction.scores$max <- NULL
sum(is.na(prediction.scores))
prediction.scores$celltype_prediction <- NA
dim(prediction.scores)
for(i in 1:nrow(prediction.scores)){
  prediction.scores$celltype_prediction[i] <- colnames(prediction.scores)[prediction.scores[i,1:15] == max(prediction.scores[i,1:15])]
}

table(prediction.scores$celltype_prediction)
chicken_visium$cellprediction <- prediction.scores$cellprediction_max

# save(chicken_visium, file="robjs/chicken_visium.4.prediction.1.Robj")
load("robjs/chicken_visium.4.prediction.1.Robj")

#####################  This sections calclulates the celltype pair neighborhood maps  ############################
# Example shown for D10 (Run this 4 times for individual stages
prediction.scores <- as.data.frame(t(GetAssayData(chicken_visium, assay = "predictions")))
prediction.scores$max <- NULL
dim(prediction.scores)
prediction.scores.1 <- prediction.scores[colnames(chicken_visium)[chicken_visium$orig.ident == "D10"],]
dim(prediction.scores.1)
interaction_matrix = matrix(0, ncol = length(unique(chicken_visium$celltype_prediction)), nrow = length(unique(chicken_visium$celltype_prediction)))
rownames(interaction_matrix) <- unique(chicken_visium$celltype_prediction)
colnames(interaction_matrix) <- unique(chicken_visium$celltype_prediction)
for(i in 1:nrow(prediction.scores.1)){
  temp <- colnames(sort(prediction.scores.1[i,prediction.scores.1[i,] > 0], decreasing = T))
  if(length(temp) == 2){
    interaction_matrix[temp[1], temp[2]] <- interaction_matrix[temp[1], temp[2]] + 1
  } else if(length(temp) == 3){
    interaction_matrix[temp[1], temp[2]] <- interaction_matrix[temp[1], temp[2]] + 1
    interaction_matrix[temp[2], temp[3]] <- interaction_matrix[temp[2], temp[3]] + 1
    interaction_matrix[temp[1], temp[3]] <- interaction_matrix[temp[1], temp[3]] + 1
  } else if(length(temp) >= 4){
    interaction_matrix[temp[1], temp[2]] <- interaction_matrix[temp[1], temp[2]] + 1
    interaction_matrix[temp[2], temp[3]] <- interaction_matrix[temp[2], temp[3]] + 1
    interaction_matrix[temp[3], temp[4]] <- interaction_matrix[temp[3], temp[4]] + 1
    interaction_matrix[temp[1], temp[3]] <- interaction_matrix[temp[1], temp[3]] + 1
    interaction_matrix[temp[1], temp[4]] <- interaction_matrix[temp[1], temp[4]] + 1
    interaction_matrix[temp[2], temp[4]] <- interaction_matrix[temp[2], temp[4]] + 1
  }
}

interaction_matrix <- interaction_matrix + t(interaction_matrix)
colnames(interaction_matrix)
temp <- colnames(interaction_matrix)[!colnames(interaction_matrix) %in% c("Erythrocytes", "Macrophages", "Mitochondria enriched cardiomyocytes")]
interaction_matrix <- interaction_matrix[temp, temp]

library(pals)
color_pelette <- rev(as.vector(kelly()[3:(2+length(levels(chicken.integrated$celltypes.0.5)))]))
names(color_pelette) <- levels(chicken.integrated$celltypes.0.5)

# write.csv(interaction_matrix, file = "interactions-D10.csv")
# interaction_matrix <- read.csv("interactions-D10.csv", row.names = 1)
interaction_matrix[lower.tri(interaction_matrix)] <- 0

library(circlize)
color_used <- color_pelette[colnames(interaction_matrix)]
row_col <- color_used
row_col[names(row_col) != "TMSB4X high cells"] <- "#cecece"

col <- matrix(rep(color_used, each = ncol(interaction_matrix), T), nrow = nrow(interaction_matrix), ncol = ncol(interaction_matrix))
rownames(col) <- rownames(interaction_matrix)
colnames(col) <- colnames(interaction_matrix)

chordDiagram(interaction_matrix, grid.col = color_used, col = col, annotationTrack = "grid")

col[rownames(col)!= "TMSB4X high cells", colnames(col) != "TMSB4X high cells"] <- "#cecece"
chordDiagram(interaction_matrix, grid.col = color_used, col = col, annotationTrack = "grid")


#####################  This section saves cell ids for visium samples  ############################

table(chicken_visium$orig.ident)
colnames(chicken_visium)
Images(chicken_visium)
sample_cell_id_map <- data.frame(sample = chicken_visium$orig.ident, cell_id = str_split_fixed(colnames(chicken_visium), "_", 2)[,2])
head(sample_cell_id_map)
# save(sample_cell_id_map, file="robjs/sample_cell_id_map.Robj")
load("robjs/sample_cell_id_map.Robj")


#####################  This section calculates the cell spot similarity map ############################
# Transfer cellnames from scRNAseq to spatial using the anchor based approach to get a cell spot similairy map
chicken.integrated$cellname <- colnames(chicken.integrated)
predictions.assay <- TransferData(anchorset = anchors, refdata = chicken.integrated$cellname, prediction.assay = TRUE, 
                                  weight.reduction = chicken_visium[["pca"]])


# Adding cellname predictions to original seurat object
chicken_visium[["predictions_cells"]] <- predictions.assay
dim(GetAssayData(chicken_visium, assay = "predictions_cells"))


#####################  This section uses the cell spot similarity map and defines spot type in 2 differet ways (optional: not used in manuscript) ############################

# Spot type defined by cell type with maxium 
prediction.scores <- as.data.frame(t(GetAssayData(chicken_visium, assay = "predictions_cells")))
# prediction.scores$max <- NULL
sum(is.na(prediction.scores))
prediction.scores$cellprediction_max <- NA
dim(prediction.scores)
for(i in 1:nrow(prediction.scores)){
  prediction.scores$cellprediction_max[i] <- colnames(prediction.scores)[prediction.scores[i,1:22315] == prediction.scores$max[i]]
}
prediction.scores$cellprediction_max
sum(is.na(prediction.scores$cellprediction_max))
temp <- str_replace_all(prediction.scores$cellprediction_max, pattern = "V-", replacement = "V_") 
temp <- str_replace_all(temp, pattern = "D4-", replacement = "D4_")
sum(!temp %in% rownames(chicken.integrated@meta.data))
prediction.scores$celltype_prediction_max <- chicken.integrated$celltypes.0.5[temp]
dim(chicken.integrated)

prediction.scores$celltype_prediction_mode <- NA
dim(prediction.scores)
for(i in 1:nrow(prediction.scores)){
  temp <- table(chicken.integrated$celltypes.0.5[prediction.scores[i,1:22315] > 0.0])
  prediction.scores$celltype_prediction_mode[i] <- names(temp)[temp == max(temp)]
}
table(prediction.scores$celltype_prediction_mode)

sum(prediction.scores$celltype_prediction_max == chicken_visium$celltype_prediction) 
sum(prediction.scores$celltype_prediction_mode == chicken_visium$celltype_prediction)

chicken_visium$celltype_prediction_max <- prediction.scores$celltype_prediction_max
chicken_visium$celltype_prediction_mode <- prediction.scores$celltype_prediction_mode

# save(chicken_visium, file="robjs/chicken_visium.4.prediction.1.Robj")
load("robjs/chicken_visium.4.prediction.1.Robj")



