# ================================================================================
# SPATIAL TRANSCRIPTOMICS ANALYSIS - MOUSE E17 KIDNEY
# Combined Analysis Pipeline with Cluster Annotation and Heatmap Generation
# ================================================================================

# Load required libraries
library(Seurat)
library(SeuratObject)
library(hdf5r)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(patchwork)

# Set output directory
output_dir <- "\\Lucie Dequiedt\\Visium kidney\\Lucie R Analysis\\late_mouse"
setwd(output_dir)

# ================================================================================
# 1. LOAD AND MERGE SAMPLES
# ================================================================================

cat("\n================================================================================\n")
cat("LOADING AND PROCESSING SAMPLES\n")
cat("================================================================================\n")

# Define paths for both samples
sample_paths <- list(
  Mouse_E17_544M = "\\Lucie Dequiedt\\Visium kidney\\SpaceRanger Outputs\\Mouse_E17_544M\\outs",
  Mouse_E17_1044M = "\\Lucie Dequiedt\\Visium kidney\\SpaceRanger Outputs\\Mouse_E17_1044M_updated\\outs"
)

# Initialize list to store Seurat objects
seurat_list <- list()

# Loop through each sample
for (sample_name in names(sample_paths)) {
  
  cat("\n========================================\n")
  cat("Processing sample:", sample_name, "\n")
  cat("========================================\n")
  
  # Define paths
  sample_path <- sample_paths[[sample_name]]
  seg_path <- file.path(sample_path, "segmented_outputs")
  
  # Load the counts
  cat("Loading counts...\n")
  counts <- Read10X(file.path(seg_path, "filtered_feature_cell_matrix"))
  
  # Create Seurat object
  cat("Creating Seurat object...\n")
  adata <- CreateSeuratObject(
    counts = counts,
    project = sample_name,
    assay = "Spatial"
  )
  
  cat("Total cells in full dataset:", ncol(adata), "\n")
  
  # Load the annotation file
  cat("Loading annotations...\n")
  annotations <- read.csv(
    file.path(seg_path, "metanephros_annotation.csv"),
    sep = ","
  )
  
  # Get metanephros cell IDs
  metanephros_cells <- annotations$Barcode
  
  # Check matching
  matching_cells <- metanephros_cells[metanephros_cells %in% colnames(adata)]
  
  cat("Total cells in annotation:", length(metanephros_cells), "\n")
  cat("Matching cells found in data:", length(matching_cells), "\n")
  
  # Subset Seurat object to metanephros cells only
  adata <- subset(adata, cells = matching_cells)
  
  # Add metadata
  adata$annotation <- "metanephros"
  adata$sample_id <- sample_name
  adata$orig.ident <- sample_name
  
  cat("Final cell count after subsetting:", ncol(adata), "\n")
  
  # Store in list
  seurat_list[[sample_name]] <- adata
  
  print(adata)
}

# Print summary of both objects
cat("\n========================================\n")
cat("Summary of loaded samples:\n")
cat("========================================\n")
for (sample_name in names(seurat_list)) {
  cat(sample_name, ":", ncol(seurat_list[[sample_name]]), "cells\n")
}

# Merge the two Seurat objects
cat("\n========================================\n")
cat("Merging samples...\n")
cat("========================================\n")

adata_merged <- merge(
  x = seurat_list[["Mouse_E17_544M"]],
  y = seurat_list[["Mouse_E17_1044M"]],
  add.cell.ids = c("544M", "1044M"),
  project = "Mouse_E17_Combined"
)

# Verify the merge
cat("Total cells in merged object:", ncol(adata_merged), "\n")
cat("Samples in merged object:\n")
print(table(adata_merged$sample_id))

print(adata_merged)

# ================================================================================
# 2. QUALITY CONTROL FILTERING
# ================================================================================

cat("\n================================================================================\n")
cat("QUALITY CONTROL\n")
cat("================================================================================\n")

# Make histograms of transcript counts for each object
pdf("QC_histograms.pdf", width = 12, height = 6)
par(mfrow = c(1, 2))

# Plot for Mouse_E17_544M
hist(adata_merged$nCount_Spatial[adata_merged$sample_id == "Mouse_E17_544M"], 
     xlim = c(0, 200), 
     breaks = 5000,
     main = "Mouse_E17_544M",
     xlab = "UMI counts (nCount_Spatial)",
     col = "lightblue")
abline(v = 40, col = "red", lwd = 2)

# Plot for Mouse_E17_1044M
hist(adata_merged$nCount_Spatial[adata_merged$sample_id == "Mouse_E17_1044M"], 
     xlim = c(0, 200), 
     breaks = 5000,
     main = "Mouse_E17_1044M",
     xlab = "UMI counts (nCount_Spatial)",
     col = "lightgreen")
abline(v = 40, col = "red", lwd = 2)

par(mfrow = c(1, 1))
dev.off()

# Filter cells based on UMI counts
cat("\nFiltering cells with nCount_Spatial >= 40\n")
adata_merged <- subset(adata_merged, subset = nCount_Spatial >= 40)
cat("Cells after filtering:", ncol(adata_merged), "\n")

# ================================================================================
# 3. NORMALIZATION AND DIMENSIONALITY REDUCTION
# ================================================================================

cat("\n================================================================================\n")
cat("NORMALIZATION AND PCA\n")
cat("================================================================================\n")
# 
# # SCTransform - normalization and variance stabilization
# adata_merged <- SCTransform(
#   adata_merged, 
#   assay = "Spatial",
#   variable.features.n = 3000,
#   verbose = TRUE
# )

adata_merged <- NormalizeData(adata_merged,scale.factor=median(adata_merged$nCount_Spatial), verbose = FALSE)
DefaultAssay(adata_merged) <- "Spatial"
adata_merged <- FindVariableFeatures(adata_merged, verbose = FALSE)
adata_merged <- ScaleData(adata_merged, verbose = FALSE)

# Run PCA
adata_merged <- RunPCA(adata_merged, npcs = 50, verbose = TRUE)

# Save PCA plots
pdf("PCA_analysis.pdf", width = 12, height = 6)
print(ElbowPlot(adata_merged, ndims = 50))
print(DimPlot(adata_merged, reduction = "pca", group.by = "sample_id"))
dev.off()

# Check top PCs
cat("\nTop principal components:\n")
print(adata_merged[["pca"]], dims = 1:5, nfeatures = 5)

# ================================================================================
# 4. CLUSTERING
# ================================================================================

cat("\n================================================================================\n")
cat("CLUSTERING\n")
cat("================================================================================\n")

# Run UMAP
adata_merged <- RunUMAP(adata_merged, dims = 1:30, return.model = TRUE, verbose = FALSE)

# Find neighbors and clusters
adata_merged <- FindNeighbors(adata_merged, reduction = "pca", dims = 1:30)
adata_merged <- FindClusters(adata_merged, resolution = 0.5)

# Save UMAP plots
pdf("UMAP_initial_clustering.pdf", width = 12, height = 6)
p1 <- DimPlot(adata_merged, group.by = "sample_id") + ggtitle("By Sample")
p2 <- DimPlot(adata_merged, group.by = "seurat_clusters") + ggtitle("By Cluster")
print(p1 | p2)
dev.off()

# ================================================================================
# 5. EXPORT CLUSTERING RESULTS FOR LOUPE
# ================================================================================

cat("\n================================================================================\n")
cat("EXPORTING CLUSTERING RESULTS\n")
cat("================================================================================\n")

# Get all clustering results
clusters_all <- data.frame(
  cell_id = colnames(adata_merged),
  samples = adata_merged$sample_id,
  cluster_full = Idents(adata_merged),
  stringsAsFactors = FALSE
)

# Export for Mouse_E17_544M
clusters_544M <- clusters_all[(clusters_all$samples == "Mouse_E17_544M"), c("cell_id", "cluster_full")]
clusters_544M$cell_id <- sub("^544M_", "", clusters_544M$cell_id)
colnames(clusters_544M) <- c("Barcode", "Cluster")
write.table(clusters_544M, "integrated_clusters_for_Mouse_E17_544M.csv",
            sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)

# Export for Mouse_E17_1044M
clusters_1044M <- clusters_all[(clusters_all$samples == "Mouse_E17_1044M"), c("cell_id", "cluster_full")]
clusters_1044M$cell_id <- sub("^1044M_", "", clusters_1044M$cell_id)
colnames(clusters_1044M) <- c("Barcode", "Cluster")
write.table(clusters_1044M, "integrated_clusters_for_Mouse_E17_1044M.csv",
            sep = ",", col.names = TRUE, row.names = FALSE, quote = FALSE)

cat("\nClustering results exported:\n")
cat("Mouse_E17_544M:", nrow(clusters_544M), "cells\n")
cat("Mouse_E17_1044M:", nrow(clusters_1044M), "cells\n")
cat("Total clusters found:", length(unique(clusters_all$cluster_full)), "\n")

# ================================================================================
# 6. FIND MARKER GENES FOR IMPORTED CLUSTERS
# ================================================================================

cat("\n================================================================================\n")
cat("FINDING MARKER GENES FOR DEFINED CLUSTERS\n")
cat("================================================================================\n")

# Set default assay to Spatial (log-normalized data)
DefaultAssay(adata_merged) <- "Spatial"

# join cancer and normal matrices in sketch (50k cells)
adata_merged[["Spatial"]] <- JoinLayers(adata_merged[["Spatial"]])

# calculate marker genes for each cluster
all.markers <- FindAllMarkers(adata_merged, only.pos = TRUE)

cat("\nTotal markers found:", nrow(all.markers), "\n")

# Filter for significance
significant_markers <- all.markers[all.markers$p_val_adj < 0.05, ]
cat("Significant markers (p_adj < 0.05):", nrow(significant_markers), "\n\n")

# Markers per cluster
cat("Markers per cluster:\n")
print(table(significant_markers$cluster))

# View top markers
cat("\nTop 30 markers:\n")
top_markers <- significant_markers[order(significant_markers$p_val_adj), ]
print(head(top_markers[, c("gene", "cluster", "avg_log2FC", "p_val_adj", "pct.1", "pct.2")], 30))

cat("\n========================================\n")
cat("MARKER FINDING SUMMARY\n")
cat("========================================\n")
cat("Total markers found:", nrow(all.markers), "\n")
cat("Significant markers (p_val_adj < 0.05):", nrow(significant_markers), "\n")
cat("\nMarkers per cluster:\n")
print(table(significant_markers$cluster))

# Display top markers for each cluster
cat("\n========================================\n")
cat("TOP 5 MARKERS PER CLUSTER\n")
cat("========================================\n")

top_5_per_cluster <- significant_markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))


# Get top markers for heatmap
top_n_markers <- 5
top_markers <- significant_markers %>%
  group_by(cluster) %>%
  top_n(n = top_n_markers, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

genes_to_plot <- unique(top_markers$gene)
cat("\nTotal unique genes to plot:", length(genes_to_plot), "\n")

print(as.data.frame(top_5_per_cluster))

# Save all significant markers
write.csv(significant_markers, "cluster_markers_significant_E17.csv", row.names = FALSE)
cat("\nSaved all significant markers to: cluster_markers_significant_E17.csv\n")

# Save top 10 markers per cluster
top_10_per_cluster <- significant_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

write.csv(top_10_per_cluster, "cluster_markers_top10_E17.csv", row.names = FALSE)
cat("Saved top 10 markers per cluster to: cluster_markers_top10_E17.csv\n")


# ================================================================================
# 8. CLUSTER ANNOTATION
# ================================================================================

cat("\n================================================================================\n")
cat("ANNOTATING CLUSTERS\n")
cat("================================================================================\n")

# Check current clusters
cat("Current clusters:", levels(Idents(adata_merged)), "\n")

# Define cluster names based on marker genes
cluster_names <- c(
  "0" = "Medullary Fibroblasts",
  "1" = "Nephrogenic Zone Other",
  "2" = "Blastema",
  "3" = "Proximal Tubule",
  "4" = "Collecting Duct",
  "5" = "Developing Nephron",
  "6" = "Distal Tubule",
  "7" = "Thick Ascending Limb of Loop of Henle",
  "8" = "Podocytes",
  "9" = "Vascular Endothelial Cells",
  "10" = "Arteries",
  "11" = "Connecting Tubule",
  "12" = "Thin Descending Limb of Loop of Henle",
  "13" = "Urothelium"
)

# Apply cluster names
adata_merged <- RenameIdents(adata_merged, cluster_names)
adata_merged$cell_type <- Idents(adata_merged)

# Define RGB-based colors
rgb_to_hex <- function(r, g, b) {
  rgb(r, g, b, maxColorValue = 255)
}

# Define cluster colors in the DESIRED ORDER for the legend
cluster_colors <- c(
  "Urothelium" = rgb_to_hex(235, 134, 181),
  "Collecting Duct" = rgb_to_hex(235, 186, 134),
  "Distal Tubule" = rgb_to_hex(191, 188, 109),
  "Thick Ascending Limb of Loop of Henle" = rgb_to_hex(113, 191, 109),
  "Thin Descending Limb of Loop of Henle" = rgb_to_hex(184, 233, 134),
  "Proximal Tubule" = rgb_to_hex(135, 214, 193),
  "Podocytes" = rgb_to_hex(134, 166, 235),
  "Developing Nephron" = rgb_to_hex(100, 66, 168),
  "Connecting Tubule" = rgb_to_hex(244, 240, 144),
  "Blastema" = rgb_to_hex(144, 55, 148),
  "Nephrogenic Zone Other" = rgb_to_hex(79, 30, 82),
  "Medullary Fibroblasts" = rgb_to_hex(246, 232, 250),
  "Arteries" = rgb_to_hex(145, 29, 29),
  "Vascular Endothelial Cells" = rgb_to_hex(77, 77, 244)
)

cat("\nCluster colors (hex codes):\n")
print(cluster_colors)

# Set the factor levels of cell_type to match the order in cluster_colors
# This controls the legend order
adata_merged$cell_type <- factor(
  adata_merged$cell_type, 
  levels = names(cluster_colors)
)

# Also set Idents to the same order
Idents(adata_merged) <- adata_merged$cell_type

cat("\nCell type levels (legend order):\n")
print(levels(adata_merged$cell_type))

# Save annotated UMAP
pdf("UMAP_annotated_clusters.pdf", width = 16, height = 6)

p1 <- DimPlot(adata_merged, 
              reduction = "umap", 
              #label = TRUE, 
              label.size = 3,
              cols = cluster_colors,
              pt.size = 0.5) + 
  ggtitle("Annotated Cell Types") +
  guides(color = guide_legend(override.aes = list(size = 8), ncol = 1)) +
  theme(legend.key.size = unit(1, "lines"),
        legend.text = element_text(size = 10))

p2 <- DimPlot(adata_merged, 
              reduction = "umap", 
              group.by = "seurat_clusters",
              pt.size = 0.5) + 
  ggtitle("Unannotated Clusters")

print(p1 | p2)
dev.off()

cat("\nCell type counts:\n")
print(table(adata_merged$cell_type))

# ================================================================================
# 7. CREATE HEATMAPS
# ================================================================================

cat("\n================================================================================\n")
cat("CREATING HEATMAPS\n")
cat("================================================================================\n")

# DIAGNOSTIC: Check what cluster names we have
cat("\nDiagnostic - Cluster names in significant_markers:\n")
print(unique(significant_markers$cluster))

cat("\nDiagnostic - Cell type names in Seurat object:\n")
print(unique(adata_merged$cell_type))

# Get top markers maintaining cluster order
top_n_markers <- 5

# IMPORTANT: The cluster column needs to match cell_type names
# Check if significant_markers$cluster contains numeric IDs or cell type names
if (is.numeric(significant_markers$cluster) || all(grepl("^[0-9]+$", significant_markers$cluster))) {
  cat("\nCluster column contains numeric IDs - need to convert to cell type names\n")
  
  # Create a mapping from numeric cluster to cell type
  # Get the mapping from your Seurat object
  cluster_to_celltype <- data.frame(
    cluster = as.character(Idents(adata_merged)),
    cell_type = adata_merged$cell_type
  )
  cluster_to_celltype <- unique(cluster_to_celltype)
  
  cat("\nCluster to cell type mapping:\n")
  print(cluster_to_celltype)
  
  # If your original clusters were stored, use those
  # Otherwise we need to map based on the seurat_clusters metadata
  if ("seurat_clusters" %in% colnames(adata_merged@meta.data)) {
    cluster_mapping <- data.frame(
      seurat_cluster = adata_merged$seurat_clusters,
      cell_type = adata_merged$cell_type
    )
    cluster_mapping <- unique(cluster_mapping)
    cluster_mapping <- setNames(cluster_mapping$cell_type, cluster_mapping$seurat_cluster)
    
    cat("\nMapping from seurat_clusters to cell_type:\n")
    print(cluster_mapping)
    
    # Apply mapping to significant_markers
    significant_markers$cell_type <- cluster_mapping[as.character(significant_markers$cluster)]
  } else {
    stop("Cannot find mapping between numeric clusters and cell types")
  }
} else {
  # Cluster names are already cell types
  significant_markers$cell_type <- significant_markers$cluster
}

# Now use cell_type column for grouping
top_markers_ordered <- significant_markers %>%
  group_by(cell_type) %>%
  top_n(n = top_n_markers, wt = avg_log2FC) %>%
  arrange(cell_type, desc(avg_log2FC))

# Convert to data frame
top_markers_ordered <- as.data.frame(top_markers_ordered)

cat("\nTotal markers before filtering:", nrow(top_markers_ordered), "\n")
cat("Cell types with markers:\n")
print(table(top_markers_ordered$cell_type))

# ===========================================
# Heatmap: Average expression per cell type
# ===========================================
cat("\nGenerating average expression heatmap for annotated clusters...\n")

# First, define the desired order of cell types (matching cluster_colors)
desired_celltype_order <- names(cluster_colors)

# Calculate average expression per cell type
avg_exp_annotated <- AverageExpression(
  adata_merged,
  assays = "Spatial",
  slot = "data",
  group.by = "cell_type"
)

avg_exp_matrix_all <- avg_exp_annotated$Spatial

# Reorder columns to match desired order (only keep cell types that exist)
desired_col_order <- desired_celltype_order[desired_celltype_order %in% colnames(avg_exp_matrix_all)]
avg_exp_matrix_all <- avg_exp_matrix_all[, desired_col_order]

cat("\nCell types in data (ordered):\n")
print(colnames(avg_exp_matrix_all))

# Now build the gene list in the SAME order as cell types
genes_to_plot_ordered <- c()
gene_cluster_map_list <- list()

for (celltype in desired_col_order) {
  # Get markers for this cell type
  celltype_markers <- top_markers_ordered[top_markers_ordered$cell_type == celltype, ]
  
  if (nrow(celltype_markers) > 0) {
    # Get top genes for this cell type
    top_genes <- head(celltype_markers$gene, top_n_markers)
    
    # Only add genes that aren't already in the list (avoid duplicates)
    new_genes <- top_genes[!top_genes %in% genes_to_plot_ordered]
    
    if (length(new_genes) > 0) {
      genes_to_plot_ordered <- c(genes_to_plot_ordered, new_genes)
      
      # Record which cluster each gene belongs to
      for (g in new_genes) {
        gene_cluster_map_list[[g]] <- celltype
      }
    }
  }
}

cat("\nTotal unique genes to plot:", length(genes_to_plot_ordered), "\n")

# Create gene_cluster_map in the correct order
gene_cluster_map <- data.frame(
  gene = genes_to_plot_ordered,
  cluster = unlist(gene_cluster_map_list[genes_to_plot_ordered]),
  stringsAsFactors = FALSE
)

# Convert cluster to factor with correct levels
gene_cluster_map$cluster <- factor(gene_cluster_map$cluster, levels = desired_col_order)

cat("\nGene order verification (first 20):\n")
print(head(gene_cluster_map, 20))

# Subset the expression matrix to only include our selected genes
avg_exp_matrix_annotated <- avg_exp_matrix_all[genes_to_plot_ordered, ]

# Z-score scale
avg_exp_scaled_annotated <- t(scale(t(avg_exp_matrix_annotated)))

# Create gene annotation (which cluster each gene belongs to)
gene_annotation <- data.frame(
  Cluster = gene_cluster_map$cluster,
  row.names = gene_cluster_map$gene
)

# Create cell type annotation for columns
celltype_annotation <- data.frame(
  CellType = colnames(avg_exp_scaled_annotated),
  row.names = colnames(avg_exp_scaled_annotated)
)

# Define cluster colors for gene annotation
clusters_with_genes <- unique(as.character(gene_cluster_map$cluster))
gene_cluster_colors <- cluster_colors[clusters_with_genes]

annotation_colors_list <- list(
  CellType = cluster_colors[colnames(avg_exp_scaled_annotated)],
  Cluster = gene_cluster_colors
)

# Calculate dimensions
n_genes <- length(genes_to_plot_ordered)
n_clusters <- length(desired_col_order)
plot_height <- max(12, min(n_genes * 0.25, 30))
plot_width <- max(10, min(n_clusters * 1.5, 20))

# Find where cluster groups change for gaps (on ROWS)
cluster_changes <- which(diff(as.numeric(gene_cluster_map$cluster)) != 0)

cat("\nCluster order in heatmap:\n")
cat("Columns (cell types):", paste(colnames(avg_exp_scaled_annotated), collapse = ", "), "\n")
cat("Rows (gene clusters):", paste(levels(gene_cluster_map$cluster)[levels(gene_cluster_map$cluster) %in% gene_cluster_map$cluster], collapse = ", "), "\n")
cat("Gap positions (row indices):", paste(cluster_changes, collapse = ", "), "\n")

pdf("heatmap_annotated_average_expression.pdf", width = plot_width, height = plot_height)

pheatmap(
  avg_exp_scaled_annotated,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_row = gene_annotation,
  annotation_col = celltype_annotation,
  annotation_colors = annotation_colors_list,
  color = colorRampPalette(c("#1C5CFE", "#D858E0", "#FEB533"))(100),
  breaks = seq(-2, 2, length.out = 101),
  fontsize_row = 7,
  fontsize_col = 10,
  angle_col = 45,
  main = "Top 5 Marker Genes per Cluster - Mouse E17",
  border_color = "grey60",
  show_rownames = TRUE,
  show_colnames = TRUE,
  gaps_row = cluster_changes
)

dev.off()

cat("\nSaved: heatmap_annotated_average_expression.pdf\n")

umap_data <- data.frame(
    Barcode = colnames(adata_merged),
    x = Embeddings(adata_merged, reduction = "umap")[, 1],
    y = Embeddings(adata_merged, reduction = "umap")[, 2],
    Cluster = as.character(adata_merged$cell_type),
    stringsAsFactors = FALSE
 )
 
write.csv(
         umap_data,
         "UMAP_coordinates_E17.csv",
         row.names = FALSE
     )
# Save scaled average expression matrix used for the heatmap
   write.csv(
         avg_exp_scaled_annotated,
         "heatmap_scaled_average_expression_E17.csv",
         row.names = TRUE
     )
