# =====================================================================
# Script 06 - Análisis de expresión diferencial con DESeq2
# Proyecto:   Análisis RNA-seq Obeso2 vs Normopeso (Los Simpson)
#
# Objetivo:   Identificar genes diferencialmente expresados (DEGs)
#             entre el grupo Obeso2 y el grupo Normopeso.
#
# Input:      03_resultados/featurecounts_resultado.rds
# Output:     03_resultados/deseq2_resultados.csv
#             03_resultados/dds_object.rds
#
# Criterio:   1.4 - Análisis de expresión diferencial (15%)
# =====================================================================
#BiocManager::install("DESeq2")
library(DESeq2)

# --- 1. Cargar la matriz de cuentas ---
fc <- readRDS("03_resultados/featurecounts_resultado.rds")
cuentas <- fc$counts

# Limpiar nombres de columnas
colnames(cuentas) <- sub("\\.bam$", "", colnames(cuentas))

# DESeq2 necesita ENTEROS (no decimales). Como usamos fraction=TRUE
# en featureCounts, las cuentas son decimales. Redondeamos.
cuentas <- round(cuentas)

cat("Dimensiones de la matriz:", dim(cuentas), "\n")
head(cuentas)

# --- 2. Construir la tabla de metadatos (coldata) ---

coldata <- data.frame(
  row.names = colnames(cuentas),
  grupo     = c(
    "Obeso1",     # AbrahamSimpson
    "Normopeso",  # BartSimpson
    "Obeso1",     # HomerSimpson
    "Normopeso",  # LisaSimpson
    "Normopeso",  # MaggieSimpson
    "Obeso2",     # MargeSimpson
    "Obeso2",     # PattyBouvier
    "Obeso2"      # SelmaBouvier
  )
)

# Convertir a factor con orden de niveles explícito

coldata$grupo <- factor(coldata$grupo,
                        levels = c("Normopeso", "Obeso1", "Obeso2"))

print(coldata)

# --- 3. Verificación crítica ---

if (!all(rownames(coldata) == colnames(cuentas))) {
  stop("Error: las muestras de coldata no coinciden con las de cuentas")
}
cat("✔ Las muestras coinciden entre coldata y cuentas\n")


# --- 4. Crear el objeto DESeqDataSet ---
# Es la estructura central que DESeq2 usa para todo el análisis.
# Combina la matriz de cuentas + los metadatos + la fórmula del diseño.
#
# La fórmula "~ grupo" significa: "modela la expresión en función del grupo"
# Es lo que le dice a DESeq2 qué variable usar para comparar.
dds <- DESeqDataSetFromMatrix(
  countData = cuentas,
  colData   = coldata,
  design    = ~ grupo
)

cat("Objeto DESeqDataSet creado:\n")
print(dds)

# --- 5. Se podria filtrar genes con muy pocas cuentas ---
# Genes con casi 0 lecturas en todas las muestras no aportan información
# y solo aumentan el número de tests. Aquí mantenemos solo genes con al
# menos 10 cuentas totales.
mantener <- rowSums(counts(dds)) >= 10
cat("Genes antes del filtrado:", nrow(dds), "\n")
dds <- dds[mantener, ]
cat("Genes después del filtrado:", nrow(dds), "\n")

# --- 6. Ejecutar DESeq2  ---
# Necesario porque con pocos genes (37) la estimación automática
# de la curva de dispersión falla. Hacemos los 3 pasos por separado.

# Paso A: estimar los size factors (normalización entre muestras)
dds <- estimateSizeFactors(dds)
cat("✔ Size factors estimados:\n")
print(sizeFactors(dds))

# Paso B: estimar la dispersión gen a gen (sin ajustar curva)
dds <- estimateDispersionsGeneEst(dds)
# Usamos las estimaciones gen a gen como valores finales
dispersions(dds) <- mcols(dds)$dispGeneEst
cat("\n✔ Dispersiones estimadas gen a gen.\n")

# Paso C: test estadístico (Wald test, el default de DESeq2)
dds <- nbinomWaldTest(dds)
cat("\n✔ Test de Wald completado.\n")


# --- 7. Extraer los resultados de Obeso2 vs Normopeso ---
# La función results() extrae log2FC, p-value, padj, etc.
# El argumento contrast = c("variable", "nivel_A", "nivel_B") significa:
#   "compara nivel_A frente a nivel_B"
# Es decir: log2FC POSITIVO = más expresión en Obeso2 que en Normopeso
#           log2FC NEGATIVO = menos expresión en Obeso2 que en Normopeso

res <- results(dds,
               contrast = c("grupo", "Obeso2", "Normopeso"),
               alpha    = 0.05)  # umbral de significancia para padj

# Resumen automático que da DESeq2
summary(res)


# Tabla ordenada por padj
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[, c("gene", "baseMean", "log2FoldChange",
                     "lfcSE", "stat", "pvalue", "padj")]
res_df <- res_df[order(res_df$padj), ]
print(res_df)


# Regenerar res_df y degs desde dds y res que sí están en memoria
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[, c("gene", "baseMean", "log2FoldChange",
                     "lfcSE", "stat", "pvalue", "padj")]
res_df <- res_df[order(res_df$padj), ]

degs <- res_df[!is.na(res_df$padj) &
                 res_df$padj < 0.05 &
                 abs(res_df$log2FoldChange) > 1, ]

cat("DEGs encontrados:", nrow(degs), "\n")

# --- 10. Guardar los resultados ---
write.csv(res_df, "03_resultados/deseq2_resultados.csv", row.names = FALSE)
write.csv(degs, "03_resultados/deseq2_DEGs.csv", row.names = FALSE)
saveRDS(dds, "03_resultados/dds_object.rds")
saveRDS(res, "03_resultados/res_object.rds")

cat("\n✔ Resultados guardados en 03_resultados/\n")