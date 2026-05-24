# =====================================================================
# Script 09 - Tabla de expresión por persona y gen
# Proyecto:   Análisis RNA-seq Obeso2 vs Normopeso (Los Simpson)
# Criterio:   2.1 - Visualización de resultados (15%)
# =====================================================================
install.packages("openxlsx")
library(DESeq2)
library(openxlsx)  

# --- 1. Cargar datos ---
dds <- readRDS("03_resultados/dds_object.rds")
res <- readRDS("03_resultados/res_object.rds")

# --- 2. Filtrar a Obeso2 + Normopeso ---
muestras_seleccionadas <- colData(dds)$grupo %in% c("Normopeso", "Obeso2")
dds_filtrado <- dds[, muestras_seleccionadas]

# --- 3. Seleccionar los DEGs ---
LFC_THRESHOLD <- 0.58

degs <- rownames(res)[which(
  !is.na(res$padj) &
    res$padj < 0.05 &
    abs(res$log2FoldChange) > LFC_THRESHOLD
)]
cat("DEGs seleccionados:", length(degs), "\n")
print(degs)

# --- 4. Obtener cuentas normalizadas ---
# DESeq2 ajusta las cuentas por los size factors (aunque en tu caso
# todos son 1, esto es lo metodológicamente correcto).
norm_counts <- counts(dds_filtrado, normalized = TRUE)

# Filtrar a solo los DEGs y redondear a 1 decimal
expr_table <- round(norm_counts[degs, ], 1)

# Reordenar las columnas: primero Normopeso, luego Obeso2
orden_columnas <- c(
  "BartSimpson", "LisaSimpson", "MaggieSimpson",     # Normopeso
  "MargeSimpson", "PattyBouvier", "SelmaBouvier"     # Obeso2
)
expr_table <- expr_table[, orden_columnas]

# Mostrar
cat("\n=== TABLA A: Cuentas normalizadas de DEGs ===\n")
print(expr_table)

# --- 5. Añadir info de DESeq2 a la tabla ---
res_df <- as.data.frame(res)
res_df <- res_df[degs, ]

# Construir la tabla final con todas las columnas relevantes
tabla_completa <- data.frame(
  Gen          = rownames(expr_table),
  
  # Cuentas en cada muestra
  Bart         = expr_table[, "BartSimpson"],
  Lisa         = expr_table[, "LisaSimpson"],
  Maggie       = expr_table[, "MaggieSimpson"],
  Marge        = expr_table[, "MargeSimpson"],
  Patty        = expr_table[, "PattyBouvier"],
  Selma        = expr_table[, "SelmaBouvier"],
  
  # Medias por grupo
  Media_Normopeso = round(rowMeans(expr_table[, 1:3]), 1),
  Media_Obeso2    = round(rowMeans(expr_table[, 4:6]), 1),
  
  # Estadísticas DESeq2
  log2FC       = round(res_df$log2FoldChange, 2),
  padj         = formatC(res_df$padj, format = "e", digits = 2),
  
  # Dirección (interpretación)
  Direccion    = ifelse(res_df$log2FoldChange > 0,
                        "UP en Obeso2", "DOWN en Obeso2"),
  
  stringsAsFactors = FALSE
)

# Ordenar: primero los UP, luego los DOWN, por significancia
tabla_completa <- tabla_completa[order(tabla_completa$Direccion,
                                       -abs(tabla_completa$log2FC)), ]

cat("\n=== TABLA B: Resumen completo de DEGs ===\n")
print(tabla_completa)

# --- 6. Guardar como CSV ---
dir.create("03_resultados/tablas", showWarnings = FALSE)

write.csv(tabla_completa,
          "03_resultados/tablas/tabla_DEGs_completa.csv",
          row.names = FALSE)

cat("\n✔ Tabla guardada en: 03_resultados/tablas/tabla_DEGs_completa.csv\n")

# --- 7. (Opcional) Guardar también como Excel con formato ---
wb <- createWorkbook()
addWorksheet(wb, "DEGs")
writeData(wb, "DEGs", tabla_completa)

# Formato bonito: cabecera en negrita, fondo gris
header_style <- createStyle(textDecoration = "bold",
                            fgFill = "#4472C4",
                            fontColour = "white",
                            halign = "center")
addStyle(wb, "DEGs", header_style, rows = 1, cols = 1:ncol(tabla_completa))

# Ancho de columnas automático
setColWidths(wb, "DEGs", cols = 1:ncol(tabla_completa), widths = "auto")

# Filas de UP en rojo claro, DOWN en azul claro
up_rows   <- which(tabla_completa$Direccion == "UP en Obeso2") + 1
down_rows <- which(tabla_completa$Direccion == "DOWN en Obeso2") + 1

up_style   <- createStyle(fgFill = "#FFE5E5")
down_style <- createStyle(fgFill = "#E5F0FF")

addStyle(wb, "DEGs", up_style,   rows = up_rows,   cols = 1:ncol(tabla_completa), gridExpand = TRUE)
addStyle(wb, "DEGs", down_style, rows = down_rows, cols = 1:ncol(tabla_completa), gridExpand = TRUE)

saveWorkbook(wb, "03_resultados/tablas/tabla_DEGs_completa.xlsx", overwrite = TRUE)

cat("✔ Tabla en Excel: 03_resultados/tablas/tabla_DEGs_completa.xlsx\n")