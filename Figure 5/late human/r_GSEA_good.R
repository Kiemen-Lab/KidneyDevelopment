# Load required libraries
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(msigdbr)
library(dplyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(devEMF)  # for EMF export

# Read your data
df <- read.csv("\\Lucie Dequiedt\\Visium kidney\\SpaceRanger Outputs\\Human_HK3\\outs\\segmented_outputs\\DE_Fib2.csv")

cat("Number of genes:", nrow(df), "\n")

# Convert gene symbols to ENTREZ IDs
entrez_data <- AnnotationDbi::select(
  org.Hs.eg.db, 
  keys = df$FeatureName,
  columns = c("SYMBOL", "ENTREZID"),
  keytype = "SYMBOL"
)

# Merge with your differential expression data
anno_result <- entrez_data %>%
  filter(!is.na(ENTREZID)) %>%
  inner_join(
    df %>% dplyr::select(FeatureName, Medullary.Fibroblast.Log2.Fold.Change),
    by = c("SYMBOL" = "FeatureName"),
    relationship = "many-to-many"
  )

# Handle duplicate ENTREZ IDs by averaging Log2FC
anno_result_unique <- anno_result %>%
  group_by(ENTREZID) %>%
  summarize(Log2FC = mean(Medullary.Fibroblast.Log2.Fold.Change, na.rm = TRUE))

cat("Number of unique genes with ENTREZ IDs:", nrow(anno_result_unique), "\n")

# Create ranked gene list
geneList <- with(anno_result_unique, setNames(Log2FC, ENTREZID))
geneList <- sort(geneList, decreasing = TRUE)

cat("Gene list range:", min(geneList), "to", max(geneList), "\n")

# ============================================
# Get ALL Reactome pathways (NO FILTERING)
# ============================================
cat("\n=== ALL REACTOME PATHWAYS ===\n")

# Get ALL Reactome pathways - NO FILTER
reactome_all <- msigdbr(
  species = "Homo sapiens", 
  category = "C2",
  subcollection = "CP:REACTOME"
) %>%
  dplyr::select(gs_name, entrez_gene)

cat("Total number of Reactome pathways:", length(unique(reactome_all$gs_name)), "\n")

# Run GSEA with ALL Reactome pathways
cat("\nRunning GSEA on all Reactome pathways...\n")
set.seed(42)
gsea_reactome <- GSEA(
  geneList, 
  TERM2GENE = reactome_all,
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  minGSSize = 10,
  maxGSSize = 500
)

cat("GSEA complete!\n")
cat("Total pathways tested:", nrow(gsea_reactome@result), "\n")
cat("Significant pathways (p.adj < 0.05):", 
    sum(gsea_reactome@result$p.adjust < 0.05), "\n")

# Convert ENTREZ IDs back to gene symbols
if(nrow(gsea_reactome@result) > 0) {
  gsea_reactome <- setReadable(gsea_reactome, 'org.Hs.eg.db', 'ENTREZID')
  
  # View top significant results
  cat("\nTop 20 significant pathways:\n")
  sig_results <- gsea_reactome@result[gsea_reactome@result$p.adjust < 0.05, ]
  sig_results <- sig_results[order(abs(sig_results$NES), decreasing = TRUE), ]
  print(head(sig_results[, c("Description", "setSize", "NES", "pvalue", "p.adjust")], 20))
  
  # Prepare data for plotting
  sorted_reactome <- gsea_reactome@result[order(gsea_reactome@result$NES, decreasing = FALSE), ]
  sorted_reactome$color <- ifelse(
    sorted_reactome$NES < 0, 
    "Enriched in Pelvic", 
    "Enriched in Medullary"
  )
  
  # Simplify pathway names
  sorted_reactome$Description_short <- gsub("REACTOME_", "", sorted_reactome$Description)
  sorted_reactome$Description_short <- gsub("_", " ", sorted_reactome$Description_short)
  sorted_reactome$Description_short <- tolower(sorted_reactome$Description_short)
  sorted_reactome$Description_short <- tools::toTitleCase(sorted_reactome$Description_short)
  
  # Create bar plot with top 15 pathways from each direction
  plot_data <- sorted_reactome %>%
    filter(p.adjust < 0.05) %>%
    dplyr::group_by(color) %>%
    dplyr::arrange(desc(abs(NES))) %>%
    slice_head(n = 15) %>%
    ungroup()
  
  p_reactome <- plot_data %>%
    ggplot(aes(x = NES, y = reorder(Description_short, NES), fill = color)) +
    geom_bar(stat = "identity") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    labs(
      y = NULL,
      x = "Normalized Enrichment Score (NES)",
      fill = NULL,
      title = "Reactome Pathways - Medullary vs Pelvic Fibroblasts"
    ) +
    theme_classic(base_size = 11) +
    scale_fill_manual(values = c(
      "Enriched in Pelvic" = "#cebed1", 
      "Enriched in Medullary" = "#A988D1"
    )) +
    theme(
      legend.position = "top",
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 9),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      panel.grid.major.y = element_blank()
    )
  
  print(p_reactome)
  
  # Save main bar plot as EMF
  emf("\\\\10.99.134.183\\kiemen-lab-data\\Lucie Dequiedt\\Visium kidney\\Lucie R Analysis\\Final plots\\GSEA_Reactome_barplot.emf",
      width = 14,
      height = 12)
  print(p_reactome)
  dev.off()
  
  # Create enrichment plots for top 3 pathways in each direction
  cat("\n=== Creating enrichment plots for top pathways ===\n")
  
  # Top 3 pathways enriched in Medullary
  top_medullary <- sorted_reactome %>%
    filter(NES > 0, p.adjust < 0.05) %>%
    arrange(desc(NES)) %>%
    slice_head(n = 3)
  
  if(nrow(top_medullary) > 0) {
    for(i in 1:min(3, nrow(top_medullary))) {
      idx <- which(gsea_reactome@result$ID == top_medullary$ID[i])
      p_med <- enrichplot::gseaplot2(
        gsea_reactome, 
        geneSetID = idx,
        title = top_medullary$Description_short[i],
        color = "#A988D1",
        base_size = 11,
        pvalue_table = TRUE
      )
      print(p_med)
      
      # Save as EMF
      emf(paste0("\\Lucie Dequiedt\\Visium kidney\\Lucie R Analysis\\Final plots\\GSEA_Medullary_top", i, "_enrichment.emf"),
          width = 10,
          height = 6)
      print(p_med)
      dev.off()
    }
  }
  
  # Top 3 pathways enriched in Pelvic
  top_pelvic <- sorted_reactome %>%
    filter(NES < 0, p.adjust < 0.05) %>%
    arrange(NES) %>%
    slice_head(n = 3)
  
  if(nrow(top_pelvic) > 0) {
    for(i in 1:min(3, nrow(top_pelvic))) {
      idx <- which(gsea_reactome@result$ID == top_pelvic$ID[i])
      p_pel <- enrichplot::gseaplot2(
        gsea_reactome, 
        geneSetID = idx,
        title = top_pelvic$Description_short[i],
        color = "#cebed1",
        base_size = 11,
        pvalue_table = TRUE
      )
      print(p_pel)
      
      # Save as EMF
      emf(paste0("\\Lucie Dequiedt\\Visium kidney\\Lucie R Analysis\\Final plots\\GSEA_Pelvic_top", i, "_enrichment.emf"),
          width = 10,
          height = 6)
      print(p_pel)
      dev.off()
    }
  }
  
  # Save all results
  write.csv(
    gsea_reactome@result, 
    "\\Lucie Dequiedt\\Visium kidney\\SpaceRanger Outputs\\Human_HK3\\outs\\segmented_outputs\\GSEA_All_Reactome_results.csv",
    row.names = FALSE
  )
  
  # Save only significant results
  write.csv(
    gsea_reactome@result[gsea_reactome@result$p.adjust < 0.05, ], 
    "\\Lucie Dequiedt\\Visium kidney\\SpaceRanger Outputs\\Human_HK3\\outs\\segmented_outputs\\GSEA_All_Reactome_significant.csv",
    row.names = FALSE
  )
  
  cat("\nReactome pathways analysis saved successfully!\n")
  
  # Summary statistics by category
  cat("\n=== SUMMARY ===\n")
  cat("Total Reactome pathways tested:", nrow(gsea_reactome@result), "\n")
  cat("Significant pathways (p.adj < 0.05):", sum(gsea_reactome@result$p.adjust < 0.05), "\n")
  cat("Enriched in Medullary:", sum(gsea_reactome@result$NES > 0 & gsea_reactome@result$p.adjust < 0.05), "\n")
  cat("Enriched in Pelvic:", sum(gsea_reactome@result$NES < 0 & gsea_reactome@result$p.adjust < 0.05), "\n")
  
  # Show pathway categories
  cat("\n=== Top pathway categories ===\n")
  sig_pathways <- gsea_reactome@result[gsea_reactome@result$p.adjust < 0.05, ]
  
  # Extract general categories
  sig_pathways$category <- gsub("REACTOME_", "", sig_pathways$Description)
  sig_pathways$category <- sapply(strsplit(sig_pathways$category, "_"), function(x) paste(x[1:min(3, length(x))], collapse = " "))
  
  cat("\nPathway themes in Medullary:\n")
  medullary_themes <- sig_pathways[sig_pathways$NES > 0, ]
  if(nrow(medullary_themes) > 0) {
    print(head(medullary_themes[order(medullary_themes$NES, decreasing = TRUE), 
                                c("Description", "NES", "p.adjust")], 10))
  }
  
  cat("\nPathway themes in Pelvic:\n")
  pelvic_themes <- sig_pathways[sig_pathways$NES < 0, ]
  if(nrow(pelvic_themes) > 0) {
    print(head(pelvic_themes[order(pelvic_themes$NES), 
                             c("Description", "NES", "p.adjust")], 10))
  }
  
} else {
  cat("No pathways found in GSEA results.\n")
}

cat("\n=== GSEA ANALYSIS COMPLETE ===\n")