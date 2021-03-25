**Title:** Spatiotemporal single-cell RNA sequencing of developing chicken hearts identifies interplay between cellular differentiation and morphogenesis.

**Author:** Madhav Mantri (mm2937@cornell.edu)

**Citation:** Mantri, M., Scuderi, G.J., Abedini-Nassab, R. et al. Spatiotemporal single-cell RNA sequencing of developing chicken hearts identifies interplay between cellular differentiation and morphogenesis. Nat Commun 12, 1771 (2021). https://doi.org/10.1038/s41467-021-21892-z

**DOI:** https://doi.org/10.1038/s41467-021-21892-z

**Contents:** This repository contains processed data, scripts and meta data required to analyse spatially resolved RNA-seq data with high-throughput time course scRNA-seq to study the spatiotemporal interactions and regulatory programs that drive embryonic development of the chicken hearts. The sequencing data needed for this analysis have been deposited in NCBI's Gene Expression Omnibus and are accessible through GEO Series accession number [GSE149457](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149457).

**Processed data:** The data/chicken_heart_spatial_RNAseq_processed directory contains processed data required to analyse the time course spatial RNAseq chicken heart datasets. The dircetory contains three files for each sptial RNA-seq sample.
1. <sample_name>_tissue_lowres_image.png: Low resolution H&E stained image file that can be loaded in the Seurat object.
2. <sample_name>_tissue_positions_list.csv: List of spot coordinates corvered with tissue section.
3. <sample_name>_scalefactors_json.json: Scaling factors for spots in the original tissue image to the low resolution image. 
Note: All three files for a spatial RNA-seq sample need to be loaded in the Seurat object for analysis.

**Meta data:** The data directory also contains meta data files with cell-type annotaitons for scRNAseq data and cell-type compositions for sptial RNAseq data for developing chicken hearts.

**Packages:** The analysis scripts written in R programming language use the following packages:
1. [Seurat-v3.2](https://satijalab.org/seurat/v3.1/spatial_vignette.html)
2. [monocle-v2](http://cole-trapnell-lab.github.io/monocle-release/docs/#installing-monocle)
3. [phateR](https://cran.r-project.org/web/packages/phateR/readme/README.html#installation)
4. [topGO](https://www.bioconductor.org/packages/release/bioc/html/topGO.html)

**R scripts:** The scripts directory contains the R scripts required for analysis of the transcriptomic datasets.
