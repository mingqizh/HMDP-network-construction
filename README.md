# 🧬 Cross-Tissue Gene Co-expression Analysis Using WGCNA

## Overview
This repository contains an R pipeline to integrate transcriptomic data across multiple tissues (adipose, liver, aorta, heart, and bone) and associate gene modules with metabolic traits using **Weighted Gene Co-expression Network Analysis (WGCNA)**.

## 📂 Data Inputs
- Gene expression files from the HMDP chow dataset:
  - `HMDP_chow_trx_adipose.txt`
  - `HMDP_chow_trx_liver.txt`
  - `HMDP_chow_trx_aorta.txt`
  - `HMDP_chow_trx_heart.txt`
  - `HMDP_chow_trx_bone.txt`
- Trait file:
  - `HMDP_chow_traits.txt`

## 🧪 Key Analyses
- Tissue-wise strain matching and preprocessing
- Unified expression matrix construction with tissue-labeled gene symbols
- Trait matrix harmonization and outlier detection via sample clustering
- WGCNA-based module detection and soft-threshold selection
- Module–trait association analysis using bicorrelation
- Network visualization with `qgraph`

## 🗃️ Outputs
- Module membership files:
  - `ME1 module members.txt`
  - `M2 module members.txt`
- Module–trait integration matrix
- Co-expression network plot with gene–gene edges

## 🛠️ Dependencies
Make sure the following R packages are installed:
```r
install.packages(c("WGCNA", "reshape2", "dplyr", "qgraph"))
# For bnstruct (from GitHub):
# install.packages("devtools")
# devtools::install_github("sambofra/bnstruct")
