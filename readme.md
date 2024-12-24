# Description

This repository contains a collection of bioinformatics scripts primarily used in my ongoing research. Most of these scripts are integral to my current work and have already been submitted to a journal for publication. Once the work is published, I will make the full set of scripts publicly available for the broader research community.
The scripts cover a variety of bioinformatics tasks, including but not limited to data analysis, statistical methods, omics data processing, pathway analysis, and network construction. These tools are designed to aid in genomic, transcriptomic, and proteomic data analysis.

**Note:** Please check back after the publication for full access to the repository. If you need any of the scripts urgently before the publication, please feel free to contact me directly, and I will personally send the requested scripts to you.
 # Email: aftab07alig@gmail.com
# ************************************************************************

# TCGA Data

- **[KM_survival.R]**: Survival analysis with TCGA data in R - Create Kaplan-Meier Curves.
- **[PAN-CANCER ANALYSIS TCGA.R]**: Pan-cancer analysis using TCGA data in R - Comprehensive analysis across multiple cancer types.
- **[tcga_data_Analysis.R]**: Analysis of TCGA data in R - General analysis workflow for TCGA datasets.
- **[tcga_data_Analysis---SURVIVAL.R]**: TCGA data analysis focused on survival - Specific analysis for survival data in TCGA.
- **[tcga_data_download.R]**: Download data from GDC Portal using TCGAbiolinks R Package.

# Immune_Infilteration

- **ImmuCellAI_mouse.R**: Implements ImmuCellAI (Immune Cell Abundance Identifier) for mouse datasets. Estimates the abundance of 24 immune cell types from transcriptomic data. Adapted for mouse gene expression datasets.
- **Immunedeconv.R**: Implements the immunedeconv R package for immune cell deconvolution. Supports multiple deconvolution methods like CIBERSORT, xCell, EPIC, and more, suitable for studying immune microenvironments.

# Data Formatting

- **Chemical name to standard names.R**: Converts chemical names to standardized forms (e.g., IUPAC or PubChem names).
- **Clustering of DEGs by K-mean.R**: Performs K-means clustering on differentially expressed genes (DEGs) to group genes with similar expression patterns.
- **Correct GMT files.R**: Cleans and corrects Gene Matrix Transposed (GMT) files for pathway analysis.
- **Count-matrix to FPKM matrix.R**: Converts raw count matrices (RNA-seq data) to FPKM values for normalization.
- **crop_all PNG files at once.R**: Batch processes PNG files by cropping them to remove unwanted margins.
- **csv_files to excel sheets.R**: Converts multiple CSV files into individual sheets within a single Excel file.
- **data tranform...multiple rows to single.R**: Combines multiple rows into a single row based on specific criteria.
- **date reshaping.R**: Formats date columns into a consistent structure (e.g., MM-DD-YYYY → YYYY-MM-DD).
- **dotplot by ordr.R**: Generates dot plots where data points are arranged in a specific order.
- **dotplot.R**: Creates dot plots showing gene/pathway significance.
- **download_proteins.R**: Automates downloading protein sequences from online resources like UniProt and NCBI.
- **Heatmap__simple.R**: Generates a basic heatmap from an input matrix.
- **human gene to mouse ortholog.R**: Maps human gene names to their corresponding mouse orthologs.
- **make logn data format.R**: Converts wide-format data to long-format for analysis.
- **Merge__ALL PDF.R**: Merges multiple PDF files into a single document.
- **attach_a col_similar.R**: Finds and attaches similar columns in datasets.
- **Merge__ALL txt or csv files.R**: Merges all TXT or CSV files in a directory into a single file.
- **Sorter__&_merger.R**: Sorts and merges data files based on specific conditions.
- **merge_all_text_file by file_name.R**: Merges text files by their filenames.
- **merging PDFs and grid image.R**: Merges PDFs and includes a grid image for visualization.
- **PCA.R**: Performs Principal Component Analysis (PCA) for dimensionality reduction.
- **Protien_names_2_Fastafiles.R**: Converts protein names or IDs to FASTA format.
- **Remove double qoutes.R**: Cleans text by removing double quotes from strings.
- **remove strings bw ().R**: Removes strings enclosed within parentheses from gene names or annotations.
- **row2Matrix_splitting.R**: Splits rows into matrix-like structures based on conditions.
- **data collapse and longer.R**: Collapses or expands data into longer formats for analysis.
- **GMT file maker.R**: Generates GMT files for pathway enrichment analyses.
- **human_GMT to Mouse GMT converterr.R**: Converts human GMT files into mouse-compatible GMT files using ortholog mapping.

# Data Visualization

- **BARPLOT WITH circle--gsea.R**: Generates a bar plot with circular markers for GSEA results.
- **Boxplot.R**: Creates boxplots for comparing data distributions across groups.
- **Clustering of DEGs by K-mean.R**: Performs K-means clustering on differentially expressed genes (DEGs).
- **dotplot.R**: Generates dot plots for pathway or enrichment results.
- **HC--clustering__with annotation.R**: Performs hierarchical clustering with annotations.
- **Heatmap split by groups.R**: Generates a heatmap with data split into predefined groups.
- **Heatmap__annotation_expression__ALLin1.R**: Creates an annotated heatmap integrating multiple elements.
- **Heatmap__simple if you have matrix.R**: A simplified heatmap script for quick visualization.
- **Heatmap__simple.R**: Another basic heatmap script for matrix-based data.
- **HEATMAP_4_Pathways& Ontologies.R**: Generates a heatmap for pathways and ontology data.
- **lollipop, bar plots.R**: Generates both lollipop and bar plots for data visualization.
- **PCA & Corr plot ##.R**: Performs PCA and generates correlation plots for sample clustering.
- **PCA & Corr plot.R**: Similar script with additional formatting options for PCA and correlation visualization.
- **Pheatmap.R**: Creates heatmaps using the pheatmap R package, widely used for gene expression visualization.
- **REACTOME & MSIGdb OVERrepresneted Analysis.R**: Performs overrepresentation analysis for Reactome and MSigDB pathways.
- **VOLCANO PLOT MAKER FROM LIST OF GENES.r**: Generates volcano plots for gene expression data, highlighting upregulated and downregulated genes.

# Microbiota Data Analysis

- **MicroBiomics_FULL PIPELINE #RANK-1.R**: Comprehensive pipeline for microbiome data analysis, including data import, normalization, filtering, and statistical analysis.
- **MicroBiomics_FULL PIPELINE-temp # Back-up.R**: Backup version of the main microbiome analysis pipeline.
- **Reltv abudnce plots # Violin, HM, bar.R**: Generates relative abundance plots using violin, heatmap, and bar plots.
- **LDA(LEFSE) # LDA_PLOTS.R**: Generates LDA plots for LEfSe (Linear Discriminant Analysis Effect Size) analysis.
- **Abundance_heatmap.R**: Creates heatmaps to visualize microbial abundance data across samples.
- **Abundance_heatmap--AGGGRETGATE SAMPLES.R**: Extended version of the heatmap script that aggregates samples based on experimental conditions.
- **lefse.pl**: Perl script for preprocessing microbiome data for LEfSe analysis.
- **lefse_format_input.py**: Python script to format input data for LEfSe.
- **lefse_run.py**: Executes the LEfSe pipeline to identify significantly enriched taxa.
- **raxmlHPC-PTHREADS-SSE3.exe**: RAxML executable for high-performance phylogenetic tree reconstruction.
- **raxml-ng**: Advanced version of RAxML for improved performance in phylogenetic analyses.
- **Optional script but useful.zip**: A collection of optional but useful scripts for microbiome analysis, including data preprocessing and visualization tools.

# RNA-Seq

- **DEGs Analysis using DESEQ2.R**: Performs Differential Expression Analysis using the DESeq2 package in R.
- **DESEQ2.r**: A foundational script for running DESeq2 analysis and generating DEGs.
- **DESEQ2------2.R**: Extended version of the DESeq2 analysis with additional filtering steps.
- **DESEQ2------3 (DEGS,VENN,HEATMAP).R**: Enhances DESeq2 analysis by adding Venn diagrams and heatmaps.
- **DSEQ2___new 16 SEP 2024.R**: A customized version of the DESeq2 analysis workflow.
## PPI Network
### PPI network__BioGRID+IntAcT.R
This script constructs a Protein-Protein Interaction (PPI) network by integrating data from the BioGRID and IntAct databases. It likely processes the raw interaction data to build a comprehensive PPI network for downstream analysis.

### INTACT database download.R
This script automates the process of downloading PPI data from the IntAct database. It may involve querying the database, retrieving relevant interaction datasets, and saving them locally for further use.

### attach_a col_similar.R
This script seems to deal with attaching or comparing column data (e.g., gene or pathway similarity) within a dataset. It might be used to identify similar attributes, such as co-expression or pathway overlaps.

### Community_enrichment.R
This script performs community detection within a network and analyzes the functional enrichment of these communities. Communities are likely clusters of highly interconnected nodes (genes/proteins), and enrichment could reveal shared biological pathways or functions.

### Complete____Rich Club script.R
This script analyzes the rich-club phenomenon in networks, which identifies a group of highly connected nodes (hubs) that form a tightly interconnected sub-network. It is used to study network topology and the importance of hubs in biological systems.

### DISGENET_in__R___An R package to explore the molecular underpinnings of human diseases.r
This script utilizes the DisGeNET R package to explore the molecular mechanisms underlying human diseases. It likely involves querying the DisGeNET database to identify disease-associated genes, pathways, and interactions.

### Hub_gene identification.R
This script identifies hub genes within a network based on centrality metrics (e.g., degree, betweenness, closeness). Hub genes are highly connected and often critical for maintaining network integrity and biological function.

### LEV-iGRAPH method.R
This script implements the LEV-iGraph method, which is likely a network analysis approach for identifying specific patterns, subgraphs, or key nodes within a large network. It could be used for graph partitioning or hierarchical analysis.

### make sub graph from main graph.R
This script generates subgraphs from a larger main graph (network). It likely focuses on extracting specific portions of the network, such as communities, pathways, or clusters, for targeted analysis.

## Pathway Enrichment Analysis
### Enrichmemnt_ANALYSIS.R
This script performs general pathway enrichment analysis using various methods (e.g., over-representation analysis, GSEA, or others). It processes gene lists and identifies enriched pathways, possibly using R packages like clusterProfiler or enrichR.

### Enrocher && EnricGO.R
This script combines enrichment analysis tools, specifically focusing on GO (Gene Ontology) enrichment and other pathway enrichment analyses. It likely uses packages like enrichGO for functional annotation.

### FISHER___ENRICHER.R
This script employs the Fisher's exact test for pathway enrichment analysis. It calculates over-representation of genes within specific pathways or functional categories.

### FULL GSEA_ENRCIHMNET ANALYSIS with gsea & Dot plot.R
This script performs a full GSEA (Gene Set Enrichment Analysis) workflow, including pathway ranking and visualization through dot plots. It provides a comprehensive analysis and visual output for enriched pathways.

### FULL GSEA_ENRCIHMNET ANALYSIS_with all plots.R
This script extends the GSEA analysis by generating multiple visualizations (e.g., enrichment plots, dot plots, bar charts, and heatmaps) for a complete overview of pathway enrichment results.

### FULL__ENRICHMENT.R
This is a comprehensive script that consolidates various pathway enrichment methods, likely including over-representation analysis (ORA), GSEA, and GSVA. It provides a full pathway enrichment workflow with relevant visualizations.

### Fisher'TEST by customised GMT___PATHWAYS.R
This script applies Fisher's exact test for enrichment analysis using a customized GMT file (a gene set database format). It enables testing for pathway enrichment tailored to specific user-defined gene sets.

### GSEA enrichment using GMT file.R
This script performs GSEA using a specific GMT file containing custom gene sets. It allows pathway enrichment analysis based on curated or user-defined gene collections.

### GSEA mouse__full MSIGDB.R
This script conducts GSEA on mouse gene expression data using the full MSigDB database. It likely involves converting mouse gene symbols to human orthologs where necessary.

### GSVA-using_countdata_2_HM.R
This script implements GSVA (Gene Set Variation Analysis) using count data as input and generates heatmaps to visualize pathway-level enrichment across samples.

### HEATMAP_4_Pathways& Ontologies.R
This script generates heatmaps to visualize pathway and ontology enrichment results. It displays the enrichment scores or p-values for multiple pathways across conditions.

### KEGG__PATH.R
This script focuses on pathway enrichment using the KEGG database to identify enriched biological pathways. It likely uses R packages like KEGGREST or clusterProfiler.

### Overrepresent Analysis 19 JUNE__GO.R
This script performs over-representation analysis specifically for Gene Ontology (GO) terms. It calculates enrichment for GO biological processes, cellular components, and molecular functions.

### PATHWAY___REACTOME by ORA.R
This script performs over-representation analysis (ORA) using the Reactome pathway database to identify enriched biological pathways.

### Pathway mapping by KEGG and Reactome.R
This script maps gene sets to pathways in both the KEGG and Reactome databases, providing dual analysis for pathway enrichment and functional annotation.

### REACTOME & MSIGdb OVERrepresneted Analysis #######perfect.R
This script conducts over-representation analysis using both Reactome pathways and the MSigDB gene set collections, ensuring robust and extensive pathway analysis.

### ReactomePathways__GSEA by GMT.R
This script performs GSEA for Reactome pathways using a GMT file as input, enabling pathway-level analysis with Reactome gene sets.

### Reactomemapping___mouse2Human_UD convertsion--perfect.R
This script maps mouse genes to human orthologs for Reactome pathway analysis, ensuring compatibility when analyzing pathways across species.

### SIMILARTYE_ENRICHER__Binary_cut Clustering.zip
This file likely contains scripts for similarity-based enrichment analysis and binary clustering of pathway or gene set data, identifying patterns of pathway activity across conditions.

### enrichPathawy_ABC.R
This script performs pathway enrichment analysis, likely in a customized manner, for a set of input gene lists or conditions.

### extract___ENTREZIDs_from MOSUE.r
This script extracts Entrez IDs from mouse gene data, converting gene symbols to standardized Entrez gene identifiers for compatibility with pathway enrichment tools.

### gprofiler2.R
This script uses the gprofiler2 package for gene set enrichment analysis. It queries databases like GO, KEGG, and Reactome for functional annotation and pathway enrichment.

### reactome mapping by GSEA (Set your parameters).R
This script performs Reactome pathway enrichment using GSEA, with options for the user to customize parameters (e.g., gene set size, p-value thresholds).

### simplifyEnrichment.R
This script simplifies and summarizes pathway enrichment results by grouping similar pathways, likely using clustering methods. It helps reduce redundancy and provides clearer interpretations of the results.

## MISC Scripts
### [dataManipulation.R]
Description: How to manipulate gene expression data from NCBI GEO in R using dplyr.

### [visualize_ggplot2.R]
Description: Visualize gene expression data in R using ggplot2.

### [dataformats.R]
Description: Reading single-cell data in R: H5AD, loom, MEX, AnnData formats.

### [singleCell_standard_workflow.R]
Description: How to analyze single-cell RNA-Seq data in R.

### [runDESeq2.R]
Description: DESeq2 workflow tutorial.

### [getData.R]
Description: DESeq2 workflow tutorial.

### [singleCell_integration.R]
Description: Integrate single-cell RNA-Seq datasets in R using Seurat (CCA).

### [singleCell_doublets.R]
Description: DoubletFinder: Detect doublets in single-cell RNA-Seq data in R.

### [microrray_RMA_normalize.R]
Description: How to read and normalize microarray data in R - RMA normalization.

### [singleCell_integrate_harmony.R]
Description: Integrate single-cell RNA-Seq data in R using Harmony.

### [singleCell_CI_markers.R]
Description: Find markers and cluster identification in single-cell RNA-Seq using Seurat.

### [singleCell_pseudoBulk.R]
Description: Pseudo-bulk analysis for single-cell RNA-Seq data.

### [convert_geneID_geneSymbols.R]
Description: 3 ways to convert Ensembl IDs to gene symbols.

### [RNASeqpipeline.sh]
Description: Setup RNA-Seq Pipeline from scratch: fastq (reads) to counts.

### [TI_monocle3.R]
Description: Single-cell Trajectory analysis using Monocle3 and Seurat.

### [WGCNA.R]
Description: Weighted Gene Co-expression Network Analysis (WGCNA) Step-by-step Tutorial.

### [variant_calling.sh]
Description: WGS Variant Calling: Variant calling with GATK - Part 1.

### [DESeq2ErrorFix.R]
Description: DESeq2 Error Fix: DESeqDataSetFromMatrix ncol(countData) == nrow(colData) is not TRUE.

### [variant_filtering_annotation.sh]
Description: WGS Variant Calling: Variant Filtering and Annotation - Part 2.

### [installPackages.R]
Description: How to install packages in R? What is CRAN? What is Bioconductor?

### [tcga_data_download.R]
Description: Download data from GDC Portal using TCGAbiolinks R Package.

### [KM_survival.R]
Description: Survival analysis with TCGA data in R - Create Kaplan-Meier Curves.

### [annotateSingleR.R]
Description: How to annotate single-cell RNA-Seq data in R using SingleR.

### [count_from_bamfile.R]
Description: Generate read count matrix from BAM file in R.

### [make_graph_from_vcf.R]
Description: Make a graph of genetic variants from VCF file in R.

### [tpm_to_counts.R]
Description: Convert TPM to counts in R.

### [singleCell_Clustering.R]
Description: How to perform clustering analysis for single-cell RNA-Seq data in R.
