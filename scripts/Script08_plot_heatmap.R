# =====================================================================
# Script 08 - Heatmap centrado en Obeso2 vs Normopeso (6 muestras)
# Proyecto:   Análisis RNA-seq Obeso2 vs Normopeso (Los Simpson)
# Criterio:   2.1 - Visualización de resultados (15%)
# =====================================================================

library(DESeq2)
library(pheatmap)
library(RColorBrewer)

# --- 1. Cargar datos ---
dds <- readRDS("03_resultados/dds_object.rds")
res <- readRDS("03_resultados/res_object.rds")

# --- 2. Filtrar a solo Obeso2 + Normopeso ---
muestras_seleccionadas <- colData(dds)$grupo %in% c("Normopeso", "Obeso2")
dds_filtrado <- dds[, muestras_seleccionadas]

cat("Muestras seleccionadas:\n")
print(data.frame(
  muestra = colnames(dds_filtrado),
  grupo   = colData(dds_filtrado)$grupo
))

# --- 3. Transformación VST ---
vsd <- varianceStabilizingTransformation(dds_filtrado, 
                                         blind = TRUE,
                                         fitType = "mean")
expr_matrix <- assay(vsd)

# --- 4. Seleccionar los DEGs ---
LFC_THRESHOLD <- 0.58

degs <- rownames(res)[which(
  !is.na(res$padj) &
    res$padj < 0.05 &
    abs(res$log2FoldChange) > LFC_THRESHOLD
)]
cat("\nDEGs seleccionados:", length(degs), "\n")
print(degs)

expr_degs <- expr_matrix[degs, ]

# --- 5. Z-score por gen ---
expr_z <- t(scale(t(expr_degs)))

# --- 6. Anotación de columnas ---
# Reordenamos los niveles del factor para que en la leyenda
# aparezca Normopeso primero y Obeso2 después
grupo_filtrado <- droplevels(colData(dds_filtrado)$grupo)

annotation_col <- data.frame(
  Grupo = grupo_filtrado,
  row.names = colnames(expr_z)
)

annotation_colors <- list(
  Grupo = c(
    Normopeso = "#90BE6D",   # verde
    Obeso2    = "#F94144"    # rojo
  )
)

# --- 7. Heatmap ---
pheatmap(
  expr_z,
  annotation_col    = annotation_col,
  annotation_colors = annotation_colors,
  
  color = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100),
  
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  
  fontsize_row    = 12,
  fontsize_col    = 11,
  cellwidth       = 40,
  cellheight      = 25,
  border_color    = "white",
  
  main = "Expresión diferencial: Obeso2 vs Normopeso\n(Z-score por gen, transformación VST)",
  
  filename = "03_resultados/figuras/heatmap_DEGs_obeso2_vs_normopeso.png",
  width    = 8,
  height   = 6
)

cat("\n✔ Heatmap guardado en: 03_resultados/figuras/heatmap_DEGs_obeso2_vs_normopeso.png\n")