# =====================================================================
# Script 07 - Visualización: Volcano plot
# Proyecto:   Análisis RNA-seq Obeso2 vs Normopeso (Los Simpson)
# Criterio:   2.1 - Visualización de resultados y etiquetado (15%)
# =====================================================================

library(ggplot2)
library(ggrepel)   # para etiquetas que no se solapen

library(DESeq2)

# --- 1. Cargar los resultados de DESeq2 ---
res <- readRDS("03_resultados/res_object.rds")
res_df <- as.data.frame(res)

# Añadir columna 'gene' desde los rownames
res_df$gene <- rownames(res_df)

# Verificación
cat("Columnas:", colnames(res_df), "\n")
cat("Dimensiones:", dim(res_df), "\n")
cat("Primeros genes:", head(res_df$gene), "\n\n")

# --- 2. Definir categorías para colorear ---
res_df$categoria <- "No significativo"
# Umbral de cambio: |log2FC| > 0.58 equivale a un cambio de 1.5x
# Adaptado al tamaño de efecto de los datos simulados (~1.75x)
LFC_THRESHOLD <- 0.58

res_df$categoria[res_df$padj < 0.05 & res_df$log2FoldChange >  LFC_THRESHOLD] <- "UP en Obeso2"
res_df$categoria[res_df$padj < 0.05 & res_df$log2FoldChange < -LFC_THRESHOLD] <- "DOWN en Obeso2"

res_df$categoria <- factor(res_df$categoria,
                           levels = c("UP en Obeso2", "DOWN en Obeso2", "No significativo"))

cat("Recuento por categoría:\n")
print(table(res_df$categoria))

# --- 3. Identificar genes a etiquetar (los DEGs) ---
res_df$etiqueta <- ifelse(res_df$categoria != "No significativo",
                          res_df$gene, "")

cat("\nGenes a etiquetar:\n")
print(res_df$gene[res_df$etiqueta != ""])

# --- 4. El gráfico ---
volcano <- ggplot(res_df, aes(x = log2FoldChange,
                              y = -log10(padj),
                              color = categoria)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "gray50", linewidth = 0.5) +
  geom_vline(xintercept = c(-LFC_THRESHOLD, LFC_THRESHOLD), 
             linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_text_repel(aes(label = etiqueta),
                  size = 4,
                  fontface = "bold",
                  box.padding = 0.5,
                  point.padding = 0.3,
                  max.overlaps = Inf,
                  show.legend = FALSE) +
  scale_color_manual(values = c(
    "UP en Obeso2"      = "#E63946",
    "DOWN en Obeso2"    = "#2A9D8F",
    "No significativo"  = "gray70"
  )) +
  labs(
    title    = "Volcano plot: Obeso2 vs Normopeso",
    subtitle = "Genes diferencialmente expresados (padj < 0.05, |log2FC| > 0.58)",
    x        = expression(log[2]~"Fold Change"),
    y        = expression(-log[10]~"(p ajustado)"),
    color    = "Categoría"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title       = element_text(face = "bold", size = 16),
    plot.subtitle    = element_text(color = "gray40"),
    legend.position  = "right",
    panel.grid.minor = element_blank()
  )

print(volcano)

# --- 5. Guardar ---
dir.create("03_resultados/figuras", showWarnings = FALSE)
ggsave("03_resultados/figuras/volcano_plot.png",
       plot = volcano, width = 9, height = 7, dpi = 300)

cat("\n✔ Volcano plot guardado en: 03_resultados/figuras/volcano_plot.png\n")