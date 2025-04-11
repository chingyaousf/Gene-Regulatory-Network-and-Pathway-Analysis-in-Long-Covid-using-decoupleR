# Transcription Factor Activity and Target Gene Annotation in Long Covid using decoupleR

## Overview

This repository contains the analysis of transcription factor (TF) activities across immune cell types in Long Covid using the `decoupleR` package. The study identifies key TFs driving immune dysregulation and links them to their target genes across Long Covid and Control samples.

## Objectives

-   Estimate transcription factor activity across 27 annotated immune cell types.
-   Compare TF activity between Long Covid and Control groups.
-   Identify the top 20 most variable TFs in each group.
-   Annotate top TFs with downstream target genes using the CollecTRI regulon.
-   Visualize TF activity and target gene profiles using heatmaps.

## Methodology

-   Single-cell RNA-seq data processed with Seurat.
-   TF activity inference using `run_ulm()` from `decoupleR` with CollecTRI regulons.
-   Group comparison across Long Covid vs. Control.
-   Visualization using heatmaps to highlight:
    -   Top 20 TFs globally.
    -   Top 20 TFs within each condition.
    -   Top 20 TFs with the largest differences between conditions.
-   Target gene annotation from CollecTRI regulons.
-   Export of annotated TF–target pairs to Excel for interpretation.

## Dependencies

-   R 4.3.0 or later\
-   Seurat (v4.3.0.1 or v5.1.0)\
-   decoupleR\
-   dplyr, tidyr, tibble\
-   pheatmap, RColorBrewer\
-   openxlsx

## Data

-   Processed Seurat object with 27 cell types:\
    `/home/cyang40/chingyao/long_covid_project/merged_samples_02/longcovid_post_labeltransfer/updated_seurat_object_with_CT_Sub_combined_02.rds`

-   TF activity results stored in:\
    `/home/cyang40/chingyao/long_covid_project/merged_samples_02/pathway_analysis_CT_Sub_combined_27celltypes/TF_analysis/`

## Workflow

### 1. TF Activity Inference

-   Subset cells by `CT_Sub_combined_02` and `disease`.
-   Infer TF activity using the Univariate Linear Model (ULM).
-   Repeat analysis for all 27 cell types across Long Covid and Control.

### 2. Merge & Summarize Activity

-   Combine all ULM result files into a single matrix.
-   Calculate the mean TF activity per cell type and condition.
-   Save cleaned matrix for visualization and analysis.

### 3. Visualization

-   Generate heatmaps for:
    -   Top 20 globally variable TFs.
    -   Top 20 most variable TFs in **Long Covid**.
    -   Top 20 most variable TFs in **Control**.
    -   Top 20 TFs with largest activity differences between groups.
-   Output figures include:
    -   `top20_TF_activities_heatmap_longcovid.png`
    -   `top20_TF_activities_heatmap_control.png`
    -   `top20_TF_activities_difference_heatmap.png`

### 4. Target Gene Annotation

-   Annotate top 20 TFs (global, per group, and difference-based) using CollecTRI regulons.
-   Filter TF–target pairs to match gene expression in Seurat object.
-   Save annotated targets to Excel files:
    -   `organized_top20_TF_targets_filtered.xlsx`
    -   `organized_top20_TF_targets_diff.xlsx`

### 5. Target Gene Activity

-   Join ULM scores with annotated target genes.
-   Calculate target gene activity across conditions.
-   Pivot to wide matrix format for visualization or export.

## Results

-   **Top connected TFs** showed high variability across immune cell types.
-   **Notable TFs in Long Covid**:
    -   CEBPB, STAT1, IRF1, EGR1, and SMAD3
-   **Key immune cell types** with distinct TF activity:
    -   CD8 T cells, Monocytes, and Plasma cells
-   **TFs with largest activity differences**:
    -   Highlighted transcriptional rewiring between Long Covid and Control.
-   **Target gene annotation** connects key TFs to inflammation, T cell exhaustion, and IFN signaling.
