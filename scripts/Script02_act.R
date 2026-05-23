# =====================================================================
# Script 02 - Construcción de la referencia para alineamiento
# Proyecto:   Análisis RNA-seq Obeso2 vs Normopeso (Los Simpson)
# Autor:      Analia Pastrana
#
# Objetivo:   Cargar el mapeo transcrito→gen y combinar todas las
#             secuencias rna.fna en un único archivo FASTA que servirá
#             de referencia para el alineamiento.
#
# Input:      00_data/referencia/tx2gene.tsv
#             00_data/referencia/carpetas_genes/*/rna.fna
# Output:     02_alineamiento/referencia_combinada.fa
#             02_alineamiento/tx2gene.rds  (el TSV ya procesado)
#
# Criterio:   1.2 - Alineamiento de secuencias (10%)
# =====================================================================

# --- BLOQUE 1: cargar el TSV de transcritos ---

# Definimos la ruta al archivo. Si tu archivo se llama diferente,
# cambia el nombre aquí (por ejemplo "transcritos.tsv").
tx2gene_path <- "00_data/referencia/Transcrito_a_Gen.tsv"

# Verificación rápida: el archivo debe existir
if (!file.exists(tx2gene_path)) {
  stop("No se encuentra el TSV en: ", tx2gene_path,
       "\nRevisa el nombre y la ubicación del archivo.")
}

# Leemos el TSV.
# - sep = "\t":          separador es tabulación (TSV)
# - header = FALSE:      tu archivo NO tiene fila de cabecera (lo vimos en la captura)
# - col.names = c(...):  le ponemos nombres a las columnas nosotros
# - stringsAsFactors = FALSE: dejar las columnas como texto, no como "factor"
tx2gene <- read.table(
  tx2gene_path,
  sep              = "\t",
  header           = FALSE,
  col.names        = c("transcript_id", "gene_symbol"),
  stringsAsFactors = FALSE
)

head(tx2gene)


# --- BLOQUE 2: Explorar la tabla tx2gene ---

# ¿Cuántas filas tiene? 
cat("Total de transcritos:", nrow(tx2gene), "\n")

# ¿Cuántos genes únicos?
cat("Total de genes únicos:", length(unique(tx2gene$gene_symbol)), "\n")

# Lista de todos los genes únicos (los símbolos)
genes_unicos <- unique(tx2gene$gene_symbol)
print(genes_unicos)

# Cuántos transcritos tiene cada gen, ordenados de más a menos
transcritos_por_gen <- sort(table(tx2gene$gene_symbol), decreasing = TRUE)
print(transcritos_por_gen)

# Estadísticos generales
cat("\nTranscritos por gen:\n")
cat("  Mínimo:  ", min(transcritos_por_gen), "\n")
cat("  Máximo:  ", max(transcritos_por_gen), "\n")
cat("  Mediana: ", median(transcritos_por_gen), "\n")
cat("  Media:   ", round(mean(transcritos_por_gen), 1), "\n")