Cross-Tissue Gene Co-expression Analysis Using WGCNA
Overview
This repository contains an R pipeline to integrate transcriptomic data across multiple tissues (adipose, liver, aorta, heart, and bone) and associate gene modules with metabolic traits using Weighted Gene Co-expression Network Analysis (WGCNA).
Data Inputs
- Gene expression files from the HMDP mouse chow dataset:
- HMDP_chow_trx_<tissue>.txt
- Trait data: HMDP_chow_traits.txt
Key Analyses
- Data preprocessing and strain intersection across tissues
- Construction of a unified expression matrix with tissue-labeled gene symbols
- Trait matrix harmonization and outlier detection via sample clustering
- WGCNA-based module detection, soft-thresholding selection, and dendrogram visualization
- Module-trait association using bicor and trait heatmaps
- Graph-based visualization of module connectivity via qgraph
Outputs
- Module membership files:
- ME1 module members.txt
- M2 module members.txt
- Module-trait integration matrix
- Co-expression network visualized with tissue-informed labels
Dependencies
- WGCNA, bnstruct, qgraph, reshape2, dplyr
- R version ≥ 4.0 recommended
# HMDP-network-construction
