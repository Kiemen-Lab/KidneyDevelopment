# Figure 5 – Spatial Transcriptomics

This folder contains R scripts used for the spatial transcriptomics analyses presented in Figure 5, including clustering, marker identification, differential expression, pathway enrichment, and visualization.

## Requirements

* R ≥ 4.3
* Seurat
* SeuratObject
* hdf5r
* dplyr
* ggplot2
* clusterProfiler
* enrichplot
* fgsea
* msigdbr
* AnnotationDbi
* org.Hs.eg.db
* ggrepel
* patchwork
* pheatmap
* devEMF

## Early mouse

`early mouse/r_copy_macaqueanalysis.R`

Processes the MA40 spatial transcriptomics dataset. The script:

* Loads the spatial expression matrix and metanephros annotations.
* Creates and filters a Seurat object.
* Performs SCTransform normalization, PCA, UMAP, and clustering.
* Identifies cluster marker genes using the Wilcoxon test.
* Exports significant and top cluster markers.
* Saves the processed Seurat object.

The script requires the Space Ranger expression matrix and the corresponding metanephros annotation file.

## Late human

### `late human/r_Good.R`

Processes the human HK3 spatial transcriptomics dataset using imported UMAP and cluster annotations. The script normalizes and scales the spatial expression data, identifies marker genes for the imported clusters, and generates marker-gene heatmaps and related outputs.

Required inputs include the Space Ranger expression matrix, UMAP coordinates, and cluster annotation files.

### `late human/r_Plot_volcano_good.R`

Generates a volcano plot comparing differential expression between medullary and pelvic fibroblasts. Genes are classified according to log2 fold change and statistical significance, with the most significant genes labeled. The script also prepares developmental, kidney, epithelial, mesenchymal, and related gene sets for pathway visualization.

Input: differential-expression table containing gene names, log2 fold changes, and p-values.

### `late human/r_GSEA_good.R`

Performs gene set enrichment analysis of medullary versus pelvic fibroblast differential expression. Gene symbols are converted to Entrez IDs, ranked by log2 fold change, and tested against Reactome pathways using `clusterProfiler`. Significant pathways are summarized and visualized using bar plots and enrichment plots.

Input: differential-expression table containing gene symbols and log2 fold changes.

## Usage

1. Install the required R packages.
2. Download or generate the required spatial transcriptomics and differential-expression data.
3. Update the input and output paths in each script.
4. Run the scripts corresponding to the desired analysis.
5. The resulting processed objects, marker tables, pathway-enrichment results, and figures are saved to the specified output locations.

The input datasets are not included in this repository.
