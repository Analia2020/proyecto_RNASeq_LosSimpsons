# =====================================================================
# Script 04 - Alineamiento single-end (solo R1)
# Proyecto:   Análisis RNA-seq Obeso2 vs Normopeso (Los Simpson)
#
# Justificación: Se detectó que los archivos R2 del grupo Obeso2
# (Marge, Patty, Selma) no eran parejas reales de R1 (los IDs base
# divergían a partir de la lectura 289). Para mantener homogeneidad
# metodológica entre las 8 muestras, se procesan todas como single-end
# usando exclusivamente los archivos R1.
#
# Input:      00_data/fastq/*_R1.fastq.gz (8 archivos)
#             02_alineamiento/indice/ref_obesidad.*  (índice ya construido)
# Output:     02_alineamiento/bam/*.bam (8 archivos)
#
# Criterio:   1.2 - Alineamiento de secuencias (10%)
# =====================================================================

library(Rsubread)

# --- Rutas ---
indice_prefijo <- "02_alineamiento/indice/ref_obesidad"

bam_dir <- "02_alineamiento/bam"
dir.create(bam_dir, showWarnings = FALSE, recursive = TRUE)

# --- Listar solo los R1 ---
r1_files <- list.files("00_data/fastq",
                       pattern = "_R1\\.fastq\\.gz$",
                       full.names = TRUE)

sample_names <- sub("_R1\\.fastq\\.gz$", "", basename(r1_files))

cat("Muestras a alinear:", length(sample_names), "\n")
print(sample_names)

# --- Rutas de salida ---
bam_files <- file.path(bam_dir, paste0(sample_names, ".bam"))

# --- Lanzar el alineamiento ---
# Nota la diferencia con la versión paired-end:
#   - readfile1 = r1_files (los R1)
#   - readfile2 = NULL     (no hay segundo archivo: single-end)
align(
  index             = indice_prefijo,
  readfile1         = r1_files,
  readfile2         = NULL,           # single-end
  output_file       = bam_files,
  type              = "rna",
  nthreads          = 2,
  unique            = FALSE,
  nBestLocations    = 1
)

cat("\n✔ Alineamiento terminado.\n")
cat("Archivos BAM generados:\n")
print(list.files(bam_dir))