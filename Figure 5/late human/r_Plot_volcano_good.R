# Load libraries
library(fgsea)
library(msigdbr)
library(ggplot2)
library(dplyr)
library(ggrepel)  # for gene labels
library(patchwork)  # for combining plots
library(devEMF)  # for EMF export

# Read data
df <- read.csv("\\Lucie Dequiedt\\Visium kidney\\SpaceRanger Outputs\\Human_HK3\\outs\\segmented_outputs\\DE_Fib2.csv")

# Get gene sets (developmental + canonical pathways)
gene_sets <- msigdbr(species = "Homo sapiens", 
                     category = "C2", subcategory = "CP:REACTOME") %>%
  bind_rows(msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP")) %>%
  filter(grepl("develop|morphogen|nephro|kidney|epithelial|mesenchym", 
               gs_name, ignore.case = TRUE)) %>%
  split(x = .$gene_symbol, f = .$gs_name)


# 1) VOLCANO PLOT WITH LABELS
df_volcano <- df %>%
  mutate(
    log2FC = Medullary.Fibroblast.Log2.Fold.Change,
    pval = Medullary.Fibroblast.P.Value,
    neg_log10p = -log10(pval),
    significant = ifelse(abs(log2FC) > 1 & pval < 0.05, 
                         ifelse(log2FC > 0, "Medullary", "Pelvic"), "NS")
  ) %>%
  filter(!is.na(log2FC) & !is.na(pval))

# Get top 5 genes for each population
top_medullary <- df_volcano %>% 
  filter(significant == "Medullary") %>% 
  arrange(pval) %>% 
  head(5)

top_pelvic <- df_volcano %>% 
  filter(significant == "Pelvic") %>% 
  arrange(pval) %>% 
  head(5)

genes_to_label <- bind_rows(top_medullary, top_pelvic)

# Create the plot
p <- ggplot(df_volcano, aes(x = log2FC, y = neg_log10p, color = significant)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text_repel(data = genes_to_label, 
                  aes(label = FeatureName),
                  size = 3.5,
                  max.overlaps = 20,
                  box.padding = 0.5,
                  point.padding = 0.3,
                  family = "Arial",
                  fontface = "bold",
                  show.legend = FALSE) +
  scale_color_manual(values = c("Medullary" = "#A988D1", "Pelvic" = "#cebed1", "NS" = "grey70"),
                     guide = guide_legend(override.aes = list(size = 5))) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed", color = "grey50", linewidth = 0.5) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.5) +
  labs(x = expression(bold("Log"[2]*" Fold Change (Medullary)")), 
       y = expression(bold("-Log"[10]*"(P-value)")), 
       color = "Enriched in") +
  theme_minimal(base_family = "Arial", base_size = 12) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold", size = 11),
    legend.text = element_text(size = 10),
    axis.title = element_text(face = "bold", size = 11),
    axis.text = element_text(color = "grey20", size = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

# Save as EMF in the specified folder
emf("\\Lucie Dequiedt\\Visium kidney\\Lucie R Analysis\\Final plots\\volcano_plot.emf", 
    width = 504/96, 
    height = 348/96)
print(p)
dev.off()

# Display the plot
print(p)