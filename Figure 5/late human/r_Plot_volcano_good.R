
# Load libraries
library(ggplot2)
library(dplyr)
library(ggrepel)
library(devEMF)

# Read data
df <- read.csv(
  "\\Lucie Dequiedt\\Visium kidney\\SpaceRanger Outputs\\Human_HK3\\outs\\segmented_outputs\\DE_Fib2.csv"
)

# Prepare data and apply Benjamini-Hochberg correction
df_volcano <- df %>%
  mutate(
    log2FC = Medullary.Fibroblast.Log2.Fold.Change,
    pval = Medullary.Fibroblast.P.Value
  ) %>%
  filter(!is.na(log2FC), !is.na(pval), is.finite(pval),
         pval >= 0, pval <= 1) %>%
  mutate(
    padj = p.adjust(pval, method = "BH"),
    neg_log10padj = -log10(pmax(padj, .Machine$double.xmin)),
    significant = case_when(
      padj < 0.05 & log2FC > 1  ~ "Medullary",
      padj < 0.05 & log2FC < -1 ~ "Pelvic",
      TRUE ~ "NS"
    )
  )

# Select the five most significant genes in each direction
top_medullary <- df_volcano %>%
  filter(significant == "Medullary") %>%
  arrange(padj) %>%
  head(5)

top_pelvic <- df_volcano %>%
  filter(significant == "Pelvic") %>%
  arrange(padj) %>%
  head(5)

genes_to_label <- bind_rows(top_medullary, top_pelvic)

# Create volcano plot
p <- ggplot(
  df_volcano,
  aes(x = log2FC, y = neg_log10padj, color = significant)
) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text_repel(
    data = genes_to_label,
    aes(label = FeatureName),
    size = 3.5,
    max.overlaps = 20,
    box.padding = 0.5,
    point.padding = 0.3,
    family = "Arial",
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "Medullary" = "#A988D1",
      "Pelvic" = "#cebed1",
      "NS" = "grey70"
    ),
    guide = guide_legend(override.aes = list(size = 5))
  ) +
  geom_vline(
    xintercept = c(-1, 1),
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.5
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey50",
    linewidth = 0.5
  ) +
  labs(
    x = expression(bold("Log"[2] * " Fold Change (Medullary)")),
    y = expression(bold("-Log"[10] * "(Adjusted P-value)")),
    color = "Enriched in"
  ) +
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

# Save plot
emf(
  "\\\\Lucie Dequiedt\\Visium kidney\\Lucie R Analysis\\Final plots\\volcano_plot_BH.emf",
  width = 504/96,
  height = 348/96
)
print(p)
dev.off()

# Display plot
print(p)

# Export results, including adjusted P-values
write.csv(
  df_volcano,
  "\\\\\\Lucie Dequiedt\\Visium kidney\\Lucie R Analysis\\Final plots\\DE_Fib2_BH_corrected.csv",
  row.names = FALSE
)

# Report number of significant genes
cat("Significant medullary genes:",
    sum(df_volcano$significant == "Medullary"), "\n")
cat("Significant pelvic genes:",
    sum(df_volcano$significant == "Pelvic"), "\n")
