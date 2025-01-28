# Differential Gene Expression Analysis in Lung Squamous Cell Carcinoma (LUSC)

A comprehensive bioinformatics project to identify **differentially expressed genes (DEGs)** in LUSC using paired and independent hypothesis testing, fold change analysis, volcano plots, and gene set enrichment analysis (GSEA).

## 📌 Project Overview

This study analyzes paired gene expression (GE) data from TCGA to identify DEGs between cancerous and healthy lung tissues in LUSC. Three methods are employed:
1. **Hypothesis Testing** (paired & independent t-tests)
2. **Fold Change Analysis**
3. **Volcano Plot Integration** (combining p-values and fold change)
4. **Gene Set Enrichment Analysis (GSEA)** on DEGs from paired tests.

## 📂 Dataset
- **Cancer Data**: `lusc-rsem-fpkm-tcga-t_paired.txt` (FPKM values for tumor tissues)
- **Healthy Data**: `lusc-rsem-fpkm-tcga_paired.txt` (FPKM values for normal tissues)
- **Paired Design**: Samples are ordered identically in both files (same patients).

## ⚙️ Requirements
- **Python Packages**: `pandas`, `numpy`, `scipy`, `matplotlib`, `seaborn`
- **GSEA Software**: [Download here](https://www.gsea-msigdb.org/gsea/index.jsp)

## 📊 Methodology
### 1. Hypothesis Testing
- **Paired t-test**: Accounts for patient-matched tumor-normal pairs.
- **Independent t-test**: Treats samples as unpaired (for comparison).
- **Multiple Testing Correction**: Benjamini-Hochberg FDR adjustment.

### 2. Fold Change Analysis
- Log2 fold change threshold: **|log2FC| > 1** (2x up/downregulation).

### 3. Volcano Plot
- Integrates p-values (from paired tests) and fold change to visualize DEGs.

### 4. GSEA
- Input: Ranked gene list (using paired test results).
- Pathway enrichment analysis via GSEA software.

## 📈 Results
- **DEG Lists**: 
  - `paired_DEGs.csv`: DEGs from paired t-test (FDR < 0.05, |log2FC| > 1).
  - `independent_DEGs.csv`: DEGs from independent t-test.
- **Volcano Plot**: `volcano_plot.png` highlights significant DEGs.

  ![image](https://github.com/user-attachments/assets/a1253be3-40dd-47f4-93c6-2fc0c3d65641)
- **GSEA Output**:

  ![image](https://github.com/user-attachments/assets/c5ac4893-97f8-41c4-b412-603ccf98a42b)
### 🔬 Contributors  
- **[Hesham Tamer]**   
  [![GitHub](https://img.shields.io/badge/GitHub-heshamtamer-blue?logo=github)](https://github.com/heshamtamer)  
- **[Reem Adel]**   
  [![GitHub](https://img.shields.io/badge/GitHub-Reeem2001-blue?logo=github)](https://github.com/Reeem2001) 
- **[Mina Adel]**
- **[Mariem Magdy]** 
