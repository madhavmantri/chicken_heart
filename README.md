**Title:** Spatiotemporal single-cell RNA sequencing of chick hearts identifies novel stage-specific regulatory programs that govern cardiac development.

**Author:** Madhav Mantri (mm2937@cornell.edu)
Link to preprint: 

**Contents:** This repository contains processed data and scripts required to analyse spatially resolved RNA-seq with high-throughput time course scRNA-seq to study the spatiotemporal interactions and regulatory programs that drive fetal development of the chicken hearts. The sequencing data needed for this analysis have been deposited in NCBI's Gene Expression Omnibus and are accessible through GEO Series accession number [GSE149457](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149457).

**Processed data:** The data directory contains processed data required to analyse the time course spatial RNA-seq chicken heart datasets. The dircetory contains three files for each sptial RNA-seq sample.
1. <sample_name>_tissue_lowres_image.png: Low resolution H&E stained image file that can be loaded in the Seurat object.
2. <sample_name>_tissue_positions_list.csv: List of spot coordinates corvered with tissue section.
3. <sample_name>_scalefactors_json.json: Scaling factors for spots in the original tissue image to the low resolution image. 
Note: All three files for a spatial RNA-seq sample need to be loaded in the Seurat object for analysis.

**Packages:** The analysis scripts written in R programming language use the following packages:
1. [Seurat-v3.2](https://satijalab.org/seurat/v3.1/spatial_vignette.html)
2. [monocle-v2](http://cole-trapnell-lab.github.io/monocle-release/docs/#installing-monocle)
3. [phateR](https://cran.r-project.org/web/packages/phateR/readme/README.html#installation)
4. [topGO](https://www.bioconductor.org/packages/release/bioc/html/topGO.html)

**R scripts:** The scripts directory contains the R scripts required for analysis of the transcriptomic datasets.
