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


# --- BLOQUE 3: Revisar la carpeta Genes ---

# Ruta a la carpeta principal donde están las carpetas de cada gen
carpetas_genes_dir <- "00_data/referencia/Genes"

# Verificación: la carpeta debe existir
if (!dir.exists(carpetas_genes_dir)) {
  stop("No se encuentra la carpeta: ", carpetas_genes_dir,
       "\nRevisa la ruta.")
}

# Listar todas las subcarpetas
carpetas_existentes <- list.dirs(
  carpetas_genes_dir,
  recursive  = FALSE,
  full.names = FALSE
)

cat("Subcarpetas encontradas:", length(carpetas_existentes), "\n\n")
print(carpetas_existentes)


genes_en_carpetas <- sub("_datasets$", "", carpetas_existentes)

# Comparación con la lista de genes del TSV
cat("\n--- Comparación ---\n")
cat("Genes en TSV:        ", length(genes_unicos), "\n")
cat("Carpetas encontradas:", length(genes_en_carpetas), "\n")

# Genes del TSV que NO tienen carpeta 
faltantes_carpeta <- setdiff(genes_unicos, genes_en_carpetas)
if (length(faltantes_carpeta) > 0) {
  cat("\n⚠ Genes SIN carpeta:\n")
  print(faltantes_carpeta)
} else {
  cat("\n✔ Todos los genes del TSV tienen su carpeta.\n")
}

# Carpetas que NO están en el TSV 
sobrantes_carpeta <- setdiff(genes_en_carpetas, genes_unicos)
if (length(sobrantes_carpeta) > 0) {
  cat("\n⚠ Carpetas que NO aparecen en el TSV:\n")
  print(sobrantes_carpeta)
}


# --- BLOQUE 4: Combinar todos los rna.fna en un solo archivo ---

# Vamos a usar Biostrings para leer y escribir FASTA correctamente
library(Biostrings)

# Carpeta de salida (debe existir; la creé arriba si no)
dir.create("02_alineamiento", showWarnings = FALSE)

# Ruta del archivo combinado de salida
output_fasta <- "02_alineamiento/referencia_combinada.fa"

# Listamos las rutas completas a todos los rna.fna
# (uno por cada carpeta de gen)
rna_fna_paths <- file.path(
  carpetas_genes_dir,
  carpetas_existentes,
  "rna.fna"
)

# Verificamos que todos existen antes de empezar
existen <- file.exists(rna_fna_paths)
if (!all(existen)) {
  cat("⚠ rna.fna faltantes en:\n")
  print(rna_fna_paths[!existen])
  stop("Hay archivos rna.fna que no se encuentran. Revisa antes de seguir.")
}
cat("✔ Los 37 archivos rna.fna existen.\n\n")

# Leemos las secuencias de cada archivo y las acumulamos
# DNAStringSet es el tipo de dato de Biostrings para secuencias de ADN.
todas_secuencias <- DNAStringSet()  # contenedor vacío

for (i in seq_along(rna_fna_paths)) {
  ruta <- rna_fna_paths[i]
  gen  <- carpetas_existentes[i]
  
  # Leer el FASTA de ese gen
  secs <- readDNAStringSet(ruta)
  
  # Concatenar al contenedor general
  todas_secuencias <- c(todas_secuencias, secs)
  
  cat(sprintf("  [%2d/%d] %-25s → %3d transcritos\n",
              i, length(rna_fna_paths), gen, length(secs)))
}

cat("\nTotal de transcritos combinados:", length(todas_secuencias), "\n")

# Guardar todo en un único archivo FASTA
writeXStringSet(todas_secuencias, filepath = output_fasta)

cat("✔ Referencia combinada guardada en:", output_fasta, "\n")
cat("  Tamaño:", round(file.info(output_fasta)$size / (1024^2), 2), "MB\n")