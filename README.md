# Single-Cell Transcriptomic Analysis of Immune Cell Composition in Systemic Lupus Erythematosus
Bioinformatics Final Project — GSE135779 Analysis

Jennifer Brow + 1000853174 + BIOL-5340-001

Project Overview

This project analyzes publicly available single-cell RNA sequencing (scRNA-seq) data from the GEO dataset GSE135779, which includes PBMC samples from pediatric lupus patients (cSLE) and healthy controls (cHD). The goal is to identify differences in immune cell composition and gene expression between disease and healthy conditions using single-cell transcriptomics.

The analysis pipeline includes:

Data acquisition and formatting

QC filtering

Normalization and log-transformation

Highly variable gene selection

PCA, neighbor graph construction

UMAP visualization

Leiden clustering

Marker-based cell type annotation

Initial condition-specific comparisons (cSLE vs cHD)


All computations were performed in Google Colab using Python and Scanpy.


Repository Contents

/notebooks/
    lupus_scRNAseq_pipeline.ipynb     # Main Colab notebook (full workflow)
    lupus_subset_analysis.ipynb        # Subset-based testing notebook

/data/
    adata_test.h5ad                    # Small high-quality subset for fast testing
    adata_test_filtered_normalized.h5ad# Normalized subset (latest version)
    adata_all_qc.h5ad                  # Full merged + QC-filtered dataset

/figures/
    qc_violin_plots.png
    umap_clusters.png
    umap_conditions.png
    marker_dotplot.png
    ...
    

README.md                              # Project summary and repo guide


How to Use This Repository

Run the analysis in Google Colab


Open the main notebook:

https://colab.research.google.com/drive/1QcEpGQUus6s0nHrxckS6Xc8JQ29U91sV?usp=sharing


Mount Google Drive

Follow the pipeline from data loading → QC → clustering → annotation


Dependencies

Python 3.10+

Scanpy

anndata

numpy / pandas / scipy

matplotlib / seaborn

python-igraph (for Leiden)


All dependencies are installed automatically inside the Colab notebook.

Dataset Information

Source: GEO Accession GSE135779

Samples:

33 pediatric SLE (cSLE1–cSLE33)

11 healthy controls (cHD1–cHD11)

Technology: 10x Genomics scRNA-seq (PBMCs)

Raw files: barcodes, gene tables, and matrix files (.tsv.gz and .mtx.gz)


AI Usage Disclosure

ChatGPT was used as a coding and conceptual assistant for:

Designing the analysis workflow

Troubleshooting Scanpy operations and merging issues

Generating summary language for the report

Interpreting QC plots and UMAP clusters


All code was reviewed, executed, and validated manually.
