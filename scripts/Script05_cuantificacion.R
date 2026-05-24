r# =====================================================================
# Script 05 - Cuantificación con featureCounts
# Proyecto:   Análisis RNA-seq Obeso2 vs Normopeso (Los Simpson)
#
# Objetivo:   Contar lecturas por gen en cada una de las 8 muestras,
#             generando una matriz gen × muestra para DESeq2.
#
# Input:      02_alineamiento/bam/*.bam (8 archivos)
#             00_data/referencia/tx2gene.tsv  (mapeo transcrito→gen)
#             02_alineamiento/referencia_combinada.fa  (para longitudes)
# Output:     03_resultados/cuentas_genes.csv
#             03_resultados/cuentas_genes.rds
#
# Criterio:   1.3 - Cuantificación y procesamiento (10%)
# =====================================================================

library(Rsubread)
library(Biostrings)

# --- 1. Cargar el TSV transcrito→gen ---
tx2gene <- read.table(
  "00_data/referencia/Transcrito_a_Gen.tsv",
  sep = "\t", header = FALSE,
  col.names = c("transcript_id", "gene_symbol"),
  stringsAsFactors = FALSE
)
cat("Transcritos en TSV:", nrow(tx2gene), "\n")

# --- 2. Cargar las secuencias de referencia para conocer la longitud
#       de cada transcrito  ---
ref <- readDNAStringSet("02_alineamiento/referencia_combinada.fa")

# Limpiar los nombres del FASTA (quitar la descripción larga,
# quedarnos solo con el ID del transcrito)
ref_ids <- sub(" .*$", "", names(ref))

# Crear tabla con ID y longitud de cada transcrito
ref_info <- data.frame(
  transcript_id = ref_ids,
  length        = width(ref),
  stringsAsFactors = FALSE
)
cat("Transcritos en referencia:", nrow(ref_info), "\n")

# --- 3. Construir la anotación SAF ---
# Combinamos tx2gene con la longitud de cada transcrito
saf <- merge(tx2gene, ref_info, by = "transcript_id")

# Reorganizar columnas en el orden que featureCounts espera:
# GeneID, Chr, Start, End, Strand
saf <- data.frame(
  GeneID = saf$gene_symbol,    # ← agrupará por gen automáticamente
  Chr    = saf$transcript_id,  # ← cada transcrito es como un "cromosoma"
  Start  = 1,                  # empieza en la base 1
  End    = saf$length,         # termina al final del transcrito
  Strand = "+",                # transcritos siempre se leen forward
  stringsAsFactors = FALSE
)

cat("\n--- Resumen de la anotación SAF ---\n")
cat("Total de filas (transcritos):", nrow(saf), "\n")
cat("Genes únicos:                ", length(unique(saf$GeneID)), "\n")
head(saf)

# --- 4. Ejecutar featureCounts sobre los 8 BAM construidos---
bam_files <- list.files("02_alineamiento/bam",
                        pattern = "\\.bam$",
                        full.names = TRUE)
cat("\nArchivos BAM encontrados:", length(bam_files), "\n")
print(basename(bam_files))

fc <- featureCounts(
  files                = bam_files,
  annot.ext            = saf,
  isGTFAnnotationFile  = FALSE,    # nuestra anotación es SAF, no GTF
  isPairedEnd          = FALSE,    # single-end 
  nthreads             = 2,
  
  # ¿Qué hacer con lecturas multi-mapping?
  # Como vimos, 85% de las lecturas son multi-mapping (transcritos de un
  # mismo gen comparten exones). Las contamos pero las "reparten":
  countMultiMappingReads = TRUE,
  fraction               = TRUE   # 1 lectura en 3 sitios 
)

# --- 5. La matriz de cuentas ---
cuentas <- fc$counts
cat("\n--- Matriz de cuentas ---\n")
cat("Dimensiones (genes × muestras):", dim(cuentas)[1], "×", dim(cuentas)[2], "\n")

# Limpiar nombres de columnas 
colnames(cuentas) <- sub("\\.bam$", "", colnames(cuentas))

# Mostrar las primeras filas
head(cuentas)

# --- 6. Guardar la matriz ---
dir.create("03_resultados", showWarnings = FALSE)

# Guardar como CSV (para inspección humana)
write.csv(cuentas, "03_resultados/cuentas_genes.csv")

# Guardar como RDS (para cargarlo rápido en R en pasos posteriores)
saveRDS(fc, "03_resultados/featurecounts_resultado.rds")

cat("\n✔ Matriz guardada en: 03_resultados/cuentas_genes.csv\n")
cat("✔ Resultado completo (fc) en: 03_resultados/featurecounts_resultado.rds\n")

# --- 7. Resumen final ---
cat("\n--- Estadísticas del conteo ---\n")
print(fc$stat)



# Cargar la matriz
fc <- readRDS("03_resultados/featurecounts_resultado.rds")
cuentas <- fc$counts
colnames(cuentas) <- sub("\\.bam$", "", colnames(cuentas))

# Ver toda la matriz (son solo 37 filas, cabe entera)
print(round(cuentas, 1))

# Dimensiones
dim(cuentas)