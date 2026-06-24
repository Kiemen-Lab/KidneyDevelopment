# Load required libraries
library(Seurat)
library(SeuratObject)
library(hdf5r)
# Single-sample path
sample_path <- "\\Lucie Dequiedt\\Visium kidney\\SpaceRanger Outputs\\Human_MA40\\outs"
seg_path <- file.path(sample_path, "segmented_outputs")
# Output directory (early_mouse)
output_dir <- "\\Lucie Dequiedt\\Visium kidney\\Lucie R Analysis\\early_mouse\\7thMAYCHANGEMACAQUE\\mouse_rds"
cat("Loading counts from:", seg_path, "\n")
counts <- Read10X(file.path(seg_path, "filtered_feature_cell_matrix"))
cat("Creating Seurat object...\n")
adata <- CreateSeuratObject(
  counts = counts,
  project = "Human_MA40",
  assay = "Spatial"
)
cat("Total cells in full dataset:", ncol(adata), "\n")
cat("Loading annotations...\n")
annotations <- read.csv(file.path(seg_path, "metanephros_annotation.csv"), sep = ",")
metanephros_cells <- annotations$Barcode
matching_cells <- metanephros_cells[metanephros_cells %in% colnames(adata)]
cat("Total cells in annotation:", length(metanephros_cells), "\n")
cat("Matching cells found in data:", length(matching_cells), "\n")
adata <- subset(adata, cells = matching_cells)
adata$annotation <- "metanephros"
adata$sample_id <- "Human_MA40"
adata$orig.ident <- "Human_MA40"
cat("Final cell count after subsetting:", ncol(adata), "\n")
print(adata)
## QC filtering
par(mfrow = c(1, 1))
hist(adata$nCount_Spatial, xlim = c(0, 200), breaks = 5000,
     main = "Human_MA40", xlab = "UMI counts (nCount_Spatial)", col = "lightblue")
abline(v = 40, col = "red", lwd = 2)
adata <- subset(adata, subset = nCount_Spatial >= 40)
## Normalization / scaling
adata <- SCTransform(
  adata,
  assay = "Spatial",
  variable.features.n = 3000,
  verbose = TRUE
)
## PCA, UMAP, Clustering
adata <- RunPCA(adata, npcs = 50, verbose = TRUE)
ElbowPlot(adata, ndims = 50)
DimPlot(adata, reduction = "pca", group.by = "sample_id")
adata <- RunUMAP(adata, dims = 1:40, n.neighbors = 50, min.dist = 0.3, return.model = TRUE, verbose = FALSE)
DimPlot(adata, group.by = "sample_id")
adata <- FindNeighbors(adata, reduction = "pca", dims = 1:40)
adata <- FindClusters(adata, resolution = 0.5)
DimPlot(adata, group.by = "sample_id")
DimPlot(adata)
## Export clustering results
clusters_df <- data.frame(
  Barcode = colnames(adata),
  Cluster = as.character(Idents(adata)),
  stringsAsFactors = FALSE
)
write.table(clusters_df, file.path(output_dir, "integrated_clusters_for_Macaque_MA40.csv"),
            sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)
cat("Clusters exported:", nrow(clusters_df), "rows\n")
## Find markers
# Ensure SCT internal models reference the Spatial assay
if (!is.null(adata@assays$SCT@SCTModel.list$model1.1)) {
  adata@assays$SCT@SCTModel.list$model1.1@umi.assay <- "Spatial"
}
if (!is.null(adata@assays$SCT@SCTModel.list$model1)) {
  adata@assays$SCT@SCTModel.list$model1@umi.assay <- "Spatial"
}
adata <- PrepSCTFindMarkers(adata, assay = "SCT", verbose = TRUE)
all_markers <- FindAllMarkers(
  adata,
  assay = "SCT",
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.1,
  test.use = "wilcox"
)
significant_markers <- all_markers[all_markers$p_val_adj < 0.05, ]
head(significant_markers, 20)
write.csv(significant_markers, file.path(output_dir, "cluster_markers_significantMACAQUE.csv"), row.names = FALSE)
# Optional: Create TOP50 markers per cluster
library(dplyr)
TOP_markers <- significant_markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 50)
write.csv(TOP_markers, file.path(output_dir, "TOP50_cluster_markers_Macaque_MA40.csv"), row.names = FALSE)
## Save Seurat object
saveRDS(adata, file.path(output_dir, "Human_MA40_processed.rds"))
cat("Done. Outputs written to:", output_dir, "\n")