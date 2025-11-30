# Single-Cell Transcriptomic Analysis of Immune Cell Composition in Systemic Lupus Erythematosus
Bioinformatics Final Project — GSE135779 Analysis

Jennifer Brow · BIOL-5340-001 · 1000853174

Overview

This project analyzes single-cell RNA-seq (scRNA-seq) PBMC data from pediatric SLE (cSLE) patients and healthy controls (cHD) using dataset GSE135779.
The goal is to characterize immune cell composition and cell-type–specific gene expression differences between disease and control samples.

A balanced ~5,000-cell subset was used for RAM-efficient analysis in Google Colab.

Workflow Summary
Data Processing

Load QC-filtered AnnData (adata_all_qc.h5ad)

Create balanced subset (~5k cells)

Compute QC metrics (genes, counts, mito %)

Normalize → log-transform

Select 2,000 highly variable genes

Dimensionality Reduction & Clustering

PCA → neighbor graph

UMAP visualization

Leiden clustering

Cell-Type Identification

Based on canonical markers:

Naive/CM T

Effector/Memory T

Classical & Non-classical Monocytes

B cells & Activated B cells

NK cells

Small populations: plasmablasts, DC-like, IFN-stimulated cells

Comparative Analyses

Composition analysis: relative cell-type abundance in cSLE vs cHD

Per–cell-type differential expression (DEG): Wilcoxon test (Scanpy)

Key Findings

cSLE and cHD contain the same major immune cell types

cSLE shows:

↑ Classical Monocytes

↓ Naive T cells

↓ NK cells

Strongest DEG seen in Classical Monocytes, with increased:

IFITM3, IFI6, TNFSF10, PLSCR1 (interferon/inflammatory genes)

Naive T cells show minimal transcriptional differences, suggesting abundance-driven changes

Overall: Pediatric lupus exhibits cell-type redistribution and modest monocyte interferon activation, which is best resolved using single-cell methods.

Tools Used

Python (Google Colab)

Scanpy / AnnData

NumPy, Pandas, SciPy

Matplotlib

UMAP, PCA

Leiden clustering

Wilcoxon DEG (rank_genes_groups)

Repository Structure
/notebooks/
    lupus_scRNAseq_pipeline.ipynb
    lupus_subset_analysis.ipynb

/data/
    adata_all_qc.h5ad
    adata_mixed_subset_5k.h5ad

/figures/
    umap_celltypes.png
    umap_diagnosis.png
    composition_barplot.png
    volcano_monocyte_naiveT.png

lupus_subset_pipeline.py
README.md

Dataset

GSE135779 — PBMC scRNA-seq from

33 pediatric cSLE patients

11 healthy controls

Technology: 10x Genomics (3’)

AI Usage Disclosure

ChatGPT assisted in code development, troubleshooting, figure interpretation, and text generation. All analysis was run and validated manually.
