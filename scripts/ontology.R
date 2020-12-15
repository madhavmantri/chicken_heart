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

setwd("/workdir/mm2937/chicken")
load(here("robjs", "chicken_normalised_scanorama3.Robj"))
load(here("robjs", "chicken_normalised_scanorama_markers_RNA.Robj"))
library(org.Gg.eg.db)
library(topGO)

cluster6_markers <- read.csv("csvs/Cluster6_all_markers.csv")
cluster6_markers <- cluster6_markers[cluster6_markers$avg_logFC > 0.6,]
geneList <- cluster6_markers$p_val_adj
names(geneList) <- rownames(cluster6_markers)

load("robjs/fibro_day_genes.Robj")
geneList <- day_genes$qval
names(geneList) <- rownames(day_genes)

names(geneList)
topDiffGenes <- function(allScore) {
    return(allScore < 10^-15)
}
sum(topDiffGenes(geneList))

GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = geneList,
              geneSel = topDiffGenes,
              annot = annFUN.org,
              mapping = "org.Gg.eg.db",
              ID = "symbol")


resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")
allRes <- GenTable(GOdata,  classicFisher = resultFisher, classicKS = resultKS, KS.elim = resultKS.elim,
                   orderBy = "classicFisher",
                   ranksOf = "classicFisher", topNodes = 50, numChar = 100)

# write.csv(allRes, file = "csvs/fibroblast-go.csv")
allRes <- read.csv("csvs/fibroblast-go.csv", row.names = 1)

library(ggplot2)
library(pals)

library(stringr)
allRes$newTerm = str_wrap(allRes$Term, width = 40)

pdf(file="fibroblast_GO.pdf",
    width=2.5, height=1.6, paper="special", bg="transparent",
    fonts="Helvetica", colormodel = "rgb", pointsize=5)
allRes$classicFisher <- as.numeric(allRes$classicFisher)
ggplot(data = allRes[1:10,]) + geom_bar(mapping = aes(x = reorder(newTerm, -classicFisher), y = -log10(classicFisher), fill = GO.ID), stat = "identity") +
  coord_flip() + theme_bw() + scale_fill_manual(values = as.vector(alphabet(20))) + 
  theme(panel.background = element_rect(fill = "transparent", color = NA),
        plot.background= element_rect(fill = "transparent", color = NA),
        plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y= element_blank(),
        axis.text = element_text(family = "Helvetica", size = 5, color = "black"),
        axis.text.x = element_text(family = "Helvetica", size = 5, color = "black"),
        plot.title=element_blank())
dev.off()  

pdf(file="SupFig3c.pdf",
    width=3.0, height=3.0, paper="special", bg="transparent",
    fonts="Helvetica", colormodel = "rgb", pointsize=5)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 6, useInfo = 'all')                                    
dev.off()


############################################### Optional script to get a mapp of GO terms for all chicken genes ###################################

library(stringr)
library(org.Gg.eg.db)
library(topGO)
selGenes <- sample(ls(org.Gg.eg.db), 50)
gene2GO <- lapply(mget(selGenes, envir = org.Gg.eg.db), names)

## Bimap interface:
x <- org.Gg.egGO
# Get the entrez gene identifiers that are mapped to a GO ID
mapped_genes <- mappedkeys(x)
# Convert to a list
gene2GO <- as.list(x[mapped_genes])
results <- data.frame("GOID" = NA, "Evidence" = NA, "Ontology" = NA)
results <- na.omit(results)
dim(results)
for(i in 1:length(gene2GO)){
  for (j in 1:length(gene2GO[[i]])){
    results[paste(names(gene2GO)[[i]], j, sep = "-"),] = gene2GO[[i]][[j]]
  }
}
results$Entrez = str_split_fixed(string = rownames(results), pattern = "-", n = 2)[,1]

## Bimap interface:
x <- org.Gg.egSYMBOL
# Get the gene symbol that are mapped to an entrez gene identifiers
mapped_genes <- mappedkeys(x)
# Convert to a list
xx <- as.list(x[mapped_genes])
results$GeneSymbol <- unlist(xx[results$Entrez])

## Bimap interface:
# Convert the object to a list
xx <- as.list(GOTERM)

temp <- results
temp$Term <- NA
for(i in 1:nrow(temp)){
  temp$Term[i] <- xx[[temp$GOID[i]]]@Term
}

write.csv(temp, file = "chicken_ontology.csv")

table(temp[temp$Ontology == "CC", "Term"])



a <- xx[temp$GOID]@`Term`
# For the reverse map:
# Convert to a list
GO2gene <- as.list(org.Gg.egGO2EG)
if(length(GO2gene) > 0){
  # Gets the entrez gene ids for the top 2nd and 3nd GO identifiers
  goids <- GO2gene[2:3]
  # Gets the entrez gene ids for the first element of goids
  goids[[1]]
  # Evidence code for the mappings
  names(goids[[1]])
}
# For org.Gg.egGO2ALLEGS
xx <- as.list(org.Gg.egGO2ALLEGS)
if(length(xx) > 0){
  # Gets the Entrez Gene identifiers for the top 2nd and 3nd GO identifiers
  goids <- xx[2:3]
  # Gets all the Entrez Gene identifiers for the first element of goids
  goids[[1]]
  # Evidence code for the mappings
  names(goids[[1]])
}

