# Clase_genOMICAFuncional - Ana Morilla a00840496


Hello this is my repository
### i find it important to note that on top of each code file theres a short note on what they are about, also for the pipeline, i ended using the WGCNA code thats called "WGCNAcode.R"

# Coffea arabica × Xylella fastidiosa — Extended Transcriptomic & Co-expression Analysis

This repository contains the downstream bioinformatics pipeline for analyzing the transcriptomic response of two *Coffea arabica* cultivars (Catuai and CR95) to *Xylella fastidiosa* infection. The workflow encompasses differential expression analysis (DESeq2), automated annotation rescue, functional enrichment (CAMERA and Hypergeometric tests), and weighted gene co-expression network analysis (WGCNA). Additionally, it provides specialized plotting scripts to visualize network hub genes and functional enrichment outputs.

## Repository Structure

```text
├── README.md
└── R_codes/
    ├── ANNOTATION SEARCH AUTOMATIZATION.R
    ├── Camera plots Up and DOWN regulation.R
    ├── Count and Log2fc graphs of cytoscape network genes.R
    ├── Modified Camera plots.R
    └── Most of the pipeline(PCA to WGCNA).R

```

## Requirements

* **R** (>= 4.2)
* **CRAN packages:** `ggplot2`, `pheatmap`, `VennDiagram`, `gridExtra`, `rentrez`, `patchwork`, `data.table`, `stringr`.
* **Bioconductor packages:** `DESeq2`, `edgeR`, `limma`, `GO.db`, `AnnotationDbi`, `WGCNA`, `pathview`.

## Input Data

To run this pipeline successfully, the working directory must contain the following input files (not included here due to size restrictions):

| File | Description | Used in |
| --- | --- | --- |
| `Matrix_hisat2.txt` | Raw gene count matrix (genes x samples) | Core pipeline, Gene plots |
| `metadata.txt` | Sample metadata (sample, cultivar, treatment) | Core pipeline, Gene plots |
| `fullAnnotation.txt` / `full_annotation_clean.csv` | Initial functional annotation matrix | Core pipeline, Gene plots |
| `DEG_summary_*.csv` | DESeq2 output files | Annotation search, KEGG mapping |
| `Cytoscape_node_attributes_*.txt` | WGCNA network outputs | Gene plots |

## Execution Order & Script Descriptions

The scripts are designed to be run in a logical progression, starting from the core pipeline down to specific visualizations and annotation rescues:

**1. `Most of the pipeline(PCA to WGCNA).R**`
This is the main workhorse script of the analysis. It loads the count matrices and metadata, filters low-expression genes, and runs a comprehensive factorial DESeq2 analysis to find differentially expressed genes (DEGs) across cultivars and treatments. It generates exploratory PCA, sample distance heatmaps, Volcano plots, and MA plots. Furthermore, it conducts CAMERA and Hypergeometric GO enrichment analysis, builds a WGCNA co-expression network, identifies trait-related modules, and exports the resulting network edges and node attributes for Cytoscape.

**2. `ANNOTATION SEARCH AUTOMATIZATION.R**`
A utility script designed to rescue missing gene annotations. It takes lists of unannotated genes (DEGs without functional descriptions) and uses the `rentrez` package to query the NCBI database automatically. It processes the gene IDs in batches to avoid server saturation and retrieves the official functional descriptions, saving the newly found annotations to text files.

**3. `Camera plots Up and DOWN regulation.R**`
A specialized visualization script for the functional enrichment results. It takes the output from the CAMERA analysis and generates custom bar plots for the top enriched GO terms. Crucially, it splits and color-codes the terms based on their direction of regulation (Upregulated vs. Downregulated) and includes a comparative plot highlighting GO terms that are shared between the Catuai and CR95 cultivars.

**4. `Modified Camera plots.R**`
An alternative iteration of the CAMERA visualization script. It produces clean bar plots for significant GO terms ranked by their False Discovery Rate (FDR). It also filters the top enriched terms to find intersections between the two cultivars, mapping the shared terms on a side-by-side grouped bar chart.

**5. `Count and Log2fc graphs of cytoscape network genes.R**`
A downstream plotting script intended to visualize specific hub genes identified in the Cytoscape network. It extracts specific `LOC` IDs from the WGCNA/Cytoscape node attributes and cross-references them with the full annotation file. Using `ggplot2` and `patchwork`, it generates a composite figure for each gene: the top panel shows raw read counts across all samples separated by cultivar and treatment, while the bottom panel displays the log2FoldChange values across all biological comparisons.

## Outputs

Each script outputs its results directly into the working directory or a designated output folder (e.g., `salido_complete/` and `plots_DEG/`). Outputs include high-resolution `.png` or `.pdf` plots (heatmaps, PCAs, Venn diagrams), tabular data (`.csv` / `.txt`) containing DESeq2 results, GO enrichment tables, and Cytoscape-ready network files.
