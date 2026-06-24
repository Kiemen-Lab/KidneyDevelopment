# ================================================================================
# SPATIAL TRANSCRIPTOMICS ANALYSIS - HUMAN KIDNEY K3
# Analysis Pipeline Using Imported Cluster Annotations
# ================================================================================

# Load required libraries
library(Seurat)
library(SeuratObject)
library(hdf5r)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(patchwork)
library(devEMF)  # for EMF export
# Set output directory
output_dir <- "\\Lucie Dequiedt\\Visium kidney\\Lucie R Analysis\\late_human"
setwd(output_dir)

# ================================================================================
# 1. LOAD SAMPLE
# ================================================================================

cat("\n================================================================================\n")
cat("LOADING AND PROCESSING HUMAN KIDNEY K3 SAMPLE\n")
cat("================================================================================\n")

# Define path for Human K3 sample
sample_name <- "Human_HK3"
sample_path <- "\\Lucie Dequiedt\\Visium kidney\\SpaceRanger Outputs\\Human_HK3\\outs"

cat("\n========================================\n")
cat("Processing sample:", sample_name, "\n")
cat("========================================\n")

# Check if segmented_outputs exists, otherwise use main filtered_feature_bc_matrix
seg_path <- file.path(sample_path, "segmented_outputs")
if (dir.exists(seg_path)) {
  data_path <- file.path(seg_path, "filtered_feature_cell_matrix")
  cat("Using segmented outputs data\n")
} else {
  data_path <- file.path(sample_path, "filtered_feature_bc_matrix")
  cat("Using standard filtered feature matrix\n")
}

# Load the counts
cat("Loading counts from:", data_path, "\n")
counts <- Read10X(data_path)

# Create Seurat object
cat("Creating Seurat object...\n")
adata <- CreateSeuratObject(
  counts = counts,
  project = sample_name,
  assay = "Spatial"
)

cat("Total cells/spots in dataset:", ncol(adata), "\n")

# Add metadata
adata$sample_id <- sample_name
adata$orig.ident <- sample_name

cat("Final cell/spot count:", ncol(adata), "\n")
print(adata)

# ================================================================================
# 2. LOAD CLUSTER ANNOTATIONS
# ================================================================================

cat("\n================================================================================\n")
cat("LOADING CLUSTER ANNOTATIONS\n")
cat("================================================================================\n")

## Load the UMAP data
umap_data <- read.csv("\\Lucie Dequiedt\\Visium kidney\\SpaceRanger Outputs\\Human_HK3\\outs\\segmented_outputs\\Test_sameParam-UMAP-Projection1.csv")

## Load the clustering data
clusters <- read.csv("\\Lucie Dequiedt\\Visium kidney\\SpaceRanger Outputs\\Human_HK3\\outs\\segmented_outputs\\Test_sameParam1.csv")

# View the first few rows
cat("\nFirst few rows of cluster data:\n")
print(head(clusters))

# Check column names to verify
cat("\nColumn names in clusters file:\n")
print(colnames(clusters))

# Merge UMAP data with cluster assignments
umap_clustered <- merge(umap_data, clusters, by = "Barcode")

# View the merged data
cat("\nFirst few rows of merged data:\n")
print(head(umap_clustered))

# Convert cluster column to factor
umap_clustered$Test_sameParam <- as.factor(umap_clustered$Test_sameParam)

# Check how many cells per cluster
cat("\nCells per cluster:\n")
print(table(umap_clustered$Test_sameParam))

# Define custom color palette for each cluster
cluster_colors <- c("Urothelium" = rgb(234, 134, 181, maxColorValue = 255),
                    "Collecting Duct" = rgb(235, 186, 134, maxColorValue = 255),
                    "Loop of Henle" = rgb(113, 191, 109, maxColorValue = 255),
                    "Proximal Tubule" = rgb(135, 214, 193, maxColorValue = 255),
                    "Podocytes" = rgb(134, 166, 235, maxColorValue = 255),
                    "Developing nephron" = rgb(100, 66, 168, maxColorValue = 255),
                    "Blastema" = rgb(144,55,148, maxColorValue = 255),
                    "Arteries" = rgb(145, 29, 29, maxColorValue = 255),
                    "Medullary Fibroblast" = rgb(169, 136, 209, maxColorValue = 255),
                    "Pelvic Fibroblasts" = rgb(246,232,250, maxColorValue = 255),
                    "Bladder Smooth Muscle" = rgb(148, 120, 156, maxColorValue = 255),
                    "Sympathetic Nerves" = rgb(64, 3, 3, maxColorValue = 255),
                    "Neuroendocrine cells" = rgb(74, 74, 74, maxColorValue = 255),
                    "UA" = rgb(10, 10, 10, maxColorValue = 255))

# Print the color palette to verify
cat("\nCustom color palette:\n")
print(cluster_colors)

# ================================================================================
# 3. CREATE UMAP VISUALIZATION WITH IMPORTED CLUSTERS
# ================================================================================

cat("\n================================================================================\n")
cat("CREATING UMAP VISUALIZATION\n")
cat("================================================================================\n")

# Set factor levels to match the order in cluster_colors
umap_clustered$Test_sameParam <- factor(
  umap_clustered$Test_sameParam,
  levels = names(cluster_colors)  # This forces the legend order
)

# Create the plot
p1 <- ggplot(umap_clustered, aes(x = X.Coordinate, y = Y.Coordinate, color = Test_sameParam)) +
  geom_point(size = 0.5, alpha = 1) +
  scale_color_manual(
    values = cluster_colors,
    breaks = names(cluster_colors)  # Explicitly set legend order
  ) +
  guides(color = guide_legend(override.aes = list(size = 8))) +  # Make legend dots bigger
  theme_minimal() +
  labs(title = "UMAP Projection - Human_HK3 Clusters",
       x = "UMAP 1",
       y = "UMAP 2",
       color = "Cell Type") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 11, face = "bold"),
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank()   # Remove minor grid lines
  )

# Display the plot
print(p1)
# Save main bar plot as EMF
emf("\\Lucie Dequiedt\\Visium kidney\\Lucie R Analysis\\Final plots\\UMAP_human.emf",
    width = 597/96,
    height = 431/96)
print(p1)
dev.off()
# ================================================================================
# 4. ADD IMPORTED CLUSTERS TO SEURAT OBJECT
# ================================================================================

cat("\n================================================================================\n")
cat("ADDING CLUSTER ANNOTATIONS TO SEURAT OBJECT\n")
cat("================================================================================\n")

# Create a named vector for easy matching
cluster_vector <- setNames(umap_clustered$Test_sameParam, umap_clustered$Barcode)

# Check matching between Seurat object and cluster annotations
matching_cells <- intersect(colnames(adata), names(cluster_vector))

cat("\nCells in Seurat object:", ncol(adata), "\n")
cat("Cells in cluster file:", length(cluster_vector), "\n")
cat("Matching cells:", length(matching_cells), "\n")

# Subset Seurat object to only cells with cluster annotations
adata <- subset(adata, cells = matching_cells)

cat("Cells after subsetting to annotated cells:", ncol(adata), "\n")

# Add cluster annotations to Seurat object
adata$imported_clusters <- cluster_vector[colnames(adata)]

# Set imported clusters as active identity
Idents(adata) <- "imported_clusters"

# Verify cluster assignments
cat("\nCells per cluster in Seurat object:\n")
print(table(adata$imported_clusters))

# ================================================================================
# 5. NORMALIZATION FOR DOWNSTREAM ANALYSIS
# ================================================================================

cat("\n================================================================================\n")
cat("NORMALIZATION AND PREPARATION FOR MARKER FINDING\n")
cat("================================================================================\n")

# # Log-normalization workflow
# cat("Starting log-normalization...\n")
adata <- NormalizeData(adata,
                       assay = "Spatial",
                       scale.factor=median(adata$nCount_Spatial), verbose = FALSE)
cat("✓ Normalization complete\n")

# Scale data
cat("Scaling data...\n")
adata <- ScaleData(adata,
                   assay = "Spatial",
                   verbose = TRUE)
cat("✓ Data scaling complete\n")
# 
# cat("\nNormalization workflow complete. Ready for downstream analysis.\n")

# ================================================================================
# 6. FIND MARKER GENES FOR IMPORTED CLUSTERS
# ================================================================================

cat("\n================================================================================\n")
cat("FINDING MARKER GENES FOR IMPORTED CLUSTERS\n")
cat("================================================================================\n")

# Set default assay to Spatial (log-normalized data)
DefaultAssay(adata) <- "Spatial"

# join cancer and normal matrices in sketch (50k cells)
adata[["Spatial"]] <- JoinLayers(adata[["Spatial"]])

# calculate marker genes for each cluster
all.markers <- FindAllMarkers(adata, only.pos = TRUE)

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
write.csv(significant_markers, "cluster_markers_significant_HK3_1.csv", row.names = FALSE)
cat("\nSaved all significant markers to: cluster_markers_significant_HK3.csv\n")

# Save top 10 markers per cluster
top_10_per_cluster <- significant_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

write.csv(top_10_per_cluster, "cluster_markers_top10_HK3_1.csv", row.names = FALSE)
cat("Saved top 10 markers per cluster to: cluster_markers_top10_HK3.csv\n")

# ================================================================================
# 7. CREATE HEATMAPS
# ================================================================================
# Get top markers maintaining cluster order
top_n_markers <- 5
top_markers_ordered <- significant_markers %>%
  group_by(cluster) %>%
  top_n(n = top_n_markers, wt = avg_log2FC) %>%
  arrange(cluster, desc(avg_log2FC))

# Convert to data frame to avoid Rle issues
top_markers_ordered <- as.data.frame(top_markers_ordered)

# Create ordered gene list (DON'T use unique() - it destroys order!)
genes_to_plot <- top_markers_ordered$gene

# Remove duplicates while maintaining order
genes_to_plot <- genes_to_plot[!duplicated(genes_to_plot)]

cat("\nTotal unique genes to plot:", length(genes_to_plot), "\n")

# Create a gene-to-cluster mapping for annotation (using base R)
gene_cluster_map <- do.call(rbind, lapply(genes_to_plot, function(g) {
  # Find first occurrence of this gene
  idx <- which(top_markers_ordered$gene == g)[1]
  data.frame(
    gene = g,
    cluster = top_markers_ordered$cluster[idx],
    stringsAsFactors = FALSE
  )
}))

cat("Gene order verification:\n")
print(head(gene_cluster_map, 20))

# ===========================================
# Heatmap: Average expression per cell type
# ===========================================
cat("\nGenerating average expression heatmap for annotated clusters...\n")

# Calculate average expression per cell type
avg_exp_annotated <- AverageExpression(
  adata,
  features = genes_to_plot,
  assays = "Spatial",
  slot = "data",
  group.by = "imported_clusters"
)

avg_exp_matrix_annotated <- avg_exp_annotated$Spatial

# Ensure genes are in correct order
avg_exp_matrix_annotated <- avg_exp_matrix_annotated[genes_to_plot, ]

# REORDER COLUMNS to match cluster_colors order
# Get the desired column order (only clusters present in the data)
desired_col_order <- names(cluster_colors)[names(cluster_colors) %in% colnames(avg_exp_matrix_annotated)]

# Reorder the matrix columns
avg_exp_matrix_annotated <- avg_exp_matrix_annotated[, desired_col_order]

# REORDER ROWS (genes) to match cluster_colors order
# First, update gene_cluster_map to use the ordered cluster levels
gene_cluster_map$cluster <- factor(gene_cluster_map$cluster, levels = names(cluster_colors))

# Sort genes by cluster order, then by original position within each cluster
gene_cluster_map <- gene_cluster_map[order(gene_cluster_map$cluster), ]

# Reorder genes_to_plot based on the sorted gene_cluster_map
genes_to_plot_ordered <- gene_cluster_map$gene

# Reorder the matrix rows
avg_exp_matrix_annotated <- avg_exp_matrix_annotated[genes_to_plot_ordered, ]

# Z-score scale
avg_exp_scaled_annotated <- t(scale(t(avg_exp_matrix_annotated)))

# Create gene annotation (which cluster each gene belongs to) - now ordered
gene_annotation <- data.frame(
  Cluster = gene_cluster_map$cluster,
  row.names = gene_cluster_map$gene
)

# Create cell type annotation for columns (now in correct order)
celltype_annotation <- data.frame(
  CellType = colnames(avg_exp_scaled_annotated),
  row.names = colnames(avg_exp_scaled_annotated)
)

# Define cluster colors for gene annotation
# Only include clusters that have genes
clusters_with_genes <- unique(as.character(gene_cluster_map$cluster))
gene_cluster_colors <- cluster_colors[clusters_with_genes]

annotation_colors_list <- list(
  CellType = cluster_colors[colnames(avg_exp_scaled_annotated)],
  Cluster = gene_cluster_colors
)

# Calculate dimensions
n_genes <- length(genes_to_plot_ordered)
n_clusters <- length(unique(colnames(avg_exp_scaled_annotated)))
plot_height <- max(12, min(n_genes * 0.25, 30))
plot_width <- max(10, min(n_clusters * 1.5, 20))

# Find where cluster groups change for gaps (on ROWS)
cluster_changes <- which(diff(as.numeric(gene_cluster_map$cluster)) != 0)

cat("Cluster order in heatmap:\n")
cat("Columns:", paste(colnames(avg_exp_scaled_annotated), collapse = ", "), "\n")
cat("Gene clusters (rows):", paste(unique(as.character(gene_cluster_map$cluster)), collapse = ", "), "\n")
cat("Gap positions:", paste(cluster_changes, collapse = ", "), "\n")

pdf("heatmap_annotated_average_expression.pdf", width = plot_width, height = plot_height)

pheatmap(
  avg_exp_scaled_annotated,
  cluster_rows = FALSE,          # Don't cluster - keep gene order
  cluster_cols = FALSE,          # Don't cluster - keep cluster order
  annotation_row = gene_annotation,     # Show which cluster each gene belongs to
  annotation_col = celltype_annotation, # Show cell types
  annotation_colors = annotation_colors_list,
  color = colorRampPalette(c("#1C5CFE", "#D858E0", "#FEB533"))(100),
  breaks = seq(-2, 2, length.out = 101),
  fontsize_row = 7,
  fontsize_col = 10,
  angle_col = 45,
  main = "Top 5 Marker Genes per Cluster - Human HK3",
  border_color = "grey60",
  show_rownames = TRUE,
  show_colnames = TRUE,
  gaps_row = cluster_changes  # Add gaps between cluster groups
)

dev.off()

cat("Saved: heatmap_annotated_average_expression.pdf\n")
